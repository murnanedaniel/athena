/*
  Copyright (C) 2002-2021 CERN for the benefit of the ATLAS collaboration
*/

/** @file AthenaRootSharedWriterSvc.cxx
 *  @brief This file contains the implementation for the AthenaRootSharedWriterSvc class.
 *  @author Peter van Gemmeren <gemmeren@anl.gov>
 **/

#include "AthenaRootSharedWriterSvc.h"

#include "TBranch.h"
#include "TClass.h"
#include "TFile.h"
#include "TFileMerger.h"
#include "TKey.h"
#include "TLeaf.h"
#include "TMemFile.h"
#include "TMessage.h"
#include "TMonitor.h"
#include "TServerSocket.h"
#include "TSocket.h"
#include "TString.h"
#include "TTree.h"

#include <set>
#include <map>

/// Definiton of a branch descriptor from RootTreeContainer
struct BranchDesc {
public:
   TClass* clazz;
   using dummy_ptr_t = std::unique_ptr<void, std::function<void(void*)> >;
   std::unique_ptr<void, std::function<void(void*)> > dummyptr;
   void* dummy = 0;

   BranchDesc(TClass* cl) : clazz(cl) {}

   void*     dummyAddr()
   {
      if (clazz) {
         void(TClass::*dxtor)(void*, Bool_t) = &TClass::Destructor;
         std::function<void(void*)> del = std::bind(dxtor, clazz, std::placeholders::_1, false);
         dummyptr = std::unique_ptr<void, std::function<void(void*)> >(clazz->New(), std::move(del));
         dummy = dummyptr.get();
         return &dummy;
      }
      return nullptr;
   }
};

/* Code from ROOT tutorials/net/parallelMergeServer.C, reduced to handle TTrees only */

struct ParallelFileMerger : public TObject
{
   TString       fFilename;
   TFileMerger   fMerger;

   ParallelFileMerger(const char *filename, int compress = ROOT::RCompressionSetting::EDefaults::kUseCompiledDefault) : fFilename(filename), fMerger(kFALSE, kTRUE)
   {
      fMerger.OutputFile(filename, "RECREATE", compress);
   }

   ~ParallelFileMerger()
   {
   }

   ULong_t Hash() const
   {
      return fFilename.Hash();
   }

   const char* GetName() const
   {
      return fFilename;
   }

// Add missing branches to client tree and BackFill before merging
   bool syncBranches(TTree* fromTree, TTree* toTree)
   {
      bool updated = false;
      const TObjArray* fromBranches = fromTree->GetListOfBranches();
      const TObjArray* toBranches = toTree->GetListOfBranches();
      int nBranches = fromBranches->GetEntriesFast();
      for (int k = 0; k < nBranches; ++k) {
         TBranch* branch = static_cast<TBranch*>(fromBranches->UncheckedAt(k));
         if (toBranches->FindObject(branch->GetName()) == nullptr) {
            TBranch* newBranch = nullptr;
            TClass* cl = TClass::GetClass(branch->GetClassName());
            BranchDesc desc(cl);
            void* empty = desc.dummyAddr();
            char buff[32];
            if (strlen(branch->GetClassName()) > 0) {
               newBranch = toTree->Branch(branch->GetName(), branch->GetClassName(), nullptr, branch->GetBasketSize(), branch->GetSplitLevel());
               newBranch->SetAddress(empty);
            } else {
               TObjArray* outLeaves = branch->GetListOfLeaves();
               TLeaf* leaf = static_cast<TLeaf*>(outLeaves->UncheckedAt(0));
               std::string type = leaf->GetTypeName();
               std::string attr = leaf->GetName();
               if (type == "Int_t") type = attr + "/I";
               else if (type == "Short_t") type = attr + "/S";
               else if (type == "Long_t") type = attr + "/L";
               else if (type == "UInt_t") type = attr + "/i";
               else if (type == "UShort_t") type = attr + "/s";
               else if (type == "UShort_t") type = attr + "/s";
               else if (type == "Float_t") type = attr + "/F";
               else if (type == "Double_t") type = attr + "/D";
               else if (type == "Char_t") type = attr + "/B";
               else if (type == "UChar_t") type = attr + "/b";
               else if (type == "Bool_t") type = attr + "/O";
               newBranch = toTree->Branch(branch->GetName(), buff, type.c_str(), 2048);
            }
            int nEntries = toTree->GetEntries();
            for (int m = 0; m < nEntries; ++m) {
               newBranch->BackFill();
            }
            updated = true;
         }
      }
      return updated;
   }

   Bool_t MergeTrees(TFile *input)
   {
      fMerger.AddFile(input);
      TTree* outCollTree = static_cast<TTree*>(fMerger.GetOutputFile()->Get("CollectionTree"));
      TTree* inCollTree = static_cast<TTree*>(input->Get("CollectionTree"));
      if (inCollTree != nullptr && outCollTree != nullptr) {
         if (syncBranches(outCollTree, inCollTree)) {
            input->Write();
         }
         syncBranches(inCollTree, outCollTree);
      }
      Bool_t result = fMerger.PartialMerge(TFileMerger::kIncremental | TFileMerger::kResetable | TFileMerger::kKeepCompression);
      TIter nextKey(input->GetListOfKeys());
      while (TKey* key = static_cast<TKey*>(nextKey())) {
         TClass* cl = TClass::GetClass(key->GetClassName());
         if (0 != cl->GetResetAfterMerge()) {
            key->Delete();
            input->GetListOfKeys()->Remove(key);
            delete key;
         }
      }
      return result;
   }
};

//___________________________________________________________________________
AthenaRootSharedWriterSvc::AthenaRootSharedWriterSvc(const std::string& name, ISvcLocator* pSvcLocator)
  : AthService(name, pSvcLocator)
  , m_rootServerSocket(nullptr), m_rootMonitor(nullptr), m_rootMergers(), m_rootClientIndex(0), m_rootClientCount(0) {
}
//___________________________________________________________________________
StatusCode AthenaRootSharedWriterSvc::initialize() {
   ATH_MSG_INFO("in initialize()");

   // Initialize IConversionSvc
   ATH_CHECK(m_cnvSvc.retrieve());
   IProperty* propertyServer = dynamic_cast<IProperty*>(m_cnvSvc.get());
   if (propertyServer == nullptr) {
      ATH_MSG_ERROR("Unable to cast conversion service to IProperty");
      return StatusCode::FAILURE;
   } else {
      std::string propertyName = "ParallelCompression";
      bool parallelCompression(false);
      BooleanProperty parallelCompressionProp(propertyName, parallelCompression);
      if (propertyServer->getProperty(&parallelCompressionProp).isFailure()) {
         ATH_MSG_INFO("Conversion service does not have ParallelCompression property");
      } else if (parallelCompressionProp.value()) {
         int streamPort = 0;
         propertyName = "StreamPortString";
         std::string streamPortString("");
         StringProperty streamPortStringProp(propertyName, streamPortString);
         if (propertyServer->getProperty(&streamPortStringProp).isFailure()) {
            ATH_MSG_INFO("Conversion service does not have StreamPortString property, using default: " << streamPort);
         } else {
            streamPort = atoi(streamPortStringProp.value().substr(streamPortStringProp.value().find(':') + 1).c_str());
         }
         m_rootServerSocket = new TServerSocket(streamPort, (streamPort == 0 ? false : true), 100);
         if (m_rootServerSocket == nullptr || !m_rootServerSocket->IsValid()) {
            ATH_MSG_FATAL("Could not create ROOT TServerSocket: " << streamPort);
            return StatusCode::FAILURE;
         }
         streamPort = m_rootServerSocket->GetLocalPort();
         const std::string newStreamPortString{streamPortStringProp.value().substr(0,streamPortStringProp.value().find(':')+1) + std::to_string(streamPort)};
         if(propertyServer->setProperty(propertyName,newStreamPortString).isFailure()) {
            ATH_MSG_FATAL("Could not set Conversion Service property " << propertyName << " from " << streamPortString << " to " << newStreamPortString);
            return StatusCode::FAILURE;
         }
         m_rootMonitor = new TMonitor;
         m_rootMonitor->Add(m_rootServerSocket);
         ATH_MSG_DEBUG("Successfully created ROOT TServerSocket and added it to TMonitor: ready to accept connections, " << streamPort);
      }
   }
   return StatusCode::SUCCESS;
}
//___________________________________________________________________________
StatusCode AthenaRootSharedWriterSvc::share(int/* numClients*/, bool motherClient) {
   ATH_MSG_DEBUG("Start commitOutput loop");
   StatusCode sc = m_cnvSvc->commitOutput("", false);

   // Allow ROOT clients to start up (by setting active clients)
   // and wait to stop the ROOT server until all clients are done and metadata is written (commitOutput fail).
   bool anyActiveClients = (m_rootServerSocket != nullptr);
   while (sc.isSuccess() || sc.isRecoverable() || anyActiveClients) {
      if (sc.isSuccess()) {
         ATH_MSG_VERBOSE("Success in commitOutput loop");
      } else if (m_rootMonitor != nullptr) {
         TSocket* socket = m_rootMonitor->Select(1);
         if (socket != nullptr && socket != (TSocket*)-1) {
            ATH_MSG_DEBUG("ROOT Monitor got: " << socket);
            if (socket->IsA() == TServerSocket::Class()) {
               TSocket* client = ((TServerSocket*)socket)->Accept();
               client->Send(m_rootClientIndex, 0);
               client->Send(1, 1);
               ++m_rootClientIndex;
               ++m_rootClientCount;
               m_rootMonitor->Add(client);
               ATH_MSG_INFO("ROOT Monitor add client: " << m_rootClientIndex << ", " << client);
            } else {
               TMessage* message = nullptr;
               Int_t result = socket->Recv(message);
               if (result < 0) {
                  ATH_MSG_ERROR("ROOT Monitor got an error while receiving the message from the socket: " << result);
                  return StatusCode::FAILURE;
               }
               if (message == nullptr) {
                  ATH_MSG_WARNING("ROOT Monitor got no message from socket: " << socket);
               } else if (message->What() == kMESS_STRING) {
                  char str[64];
                  message->ReadString(str, 64);
                  ATH_MSG_INFO("ROOT Monitor client: " << socket << ", " << str);
                  m_rootMonitor->Remove(socket);
                  ATH_MSG_DEBUG("ROOT Monitor client: " << socket << ", " << socket->GetBytesRecv() << ", " << socket->GetBytesSent());
                  socket->Close();
                  --m_rootClientCount;
                  if (m_rootMonitor->GetActive() == 0 || m_rootClientCount == 0) {
                     if (!motherClient) {
                        anyActiveClients = false;
                        ATH_MSG_INFO("ROOT Monitor: No more active clients...");
                     } else {
                        motherClient = false;
                        ATH_MSG_INFO("ROOT Monitor: Mother process is done...");
                        if (!m_cnvSvc->commitCatalog().isSuccess()) {
                           ATH_MSG_FATAL("Failed to commit file catalog.");
                           return StatusCode::FAILURE;
                        }
                     }
                  }
               } else if (message->What() == kMESS_ANY) {
                  long long length;
                  TString filename;
                  int clientId;
                  message->ReadInt(clientId);
                  message->ReadTString(filename);
                  message->ReadLong64(length);
                  ATH_MSG_DEBUG("ROOT Monitor client: " << socket << ", " << clientId << ": " << filename << ", " << length);
                  std::unique_ptr<TMemFile> transient(new TMemFile(filename, message->Buffer() + message->Length(), length, "UPDATE"));
                  message->SetBufferOffset(message->Length() + length);
                  ParallelFileMerger* info = static_cast<ParallelFileMerger*>(m_rootMergers.FindObject(filename));
                  if (!info) {
                     info = new ParallelFileMerger(filename, transient->GetCompressionSettings());
                     m_rootMergers.Add(info);
                     ATH_MSG_INFO("ROOT Monitor ParallelFileMerger: " << info << ", for: " << filename);
                  }
                  info->MergeTrees(transient.get());
               }
               delete message; message = nullptr;
            }
         }
      } else if (m_rootMonitor == nullptr) {
         usleep(100);
      }
      // Once commitOutput failed all legacy clients are finished (writing metadata), do not call again.
      if (sc.isSuccess() || sc.isRecoverable()) {
         sc = m_cnvSvc->commitOutput("", false);
         if (sc.isFailure() && !sc.isRecoverable()) {
            ATH_MSG_INFO("commitOutput failed, metadata done.");
            if (anyActiveClients && m_rootClientCount == 0) {
              ATH_MSG_INFO("ROOT Monitor: No clients, terminating the loop...");
              anyActiveClients = false;
            }
         }
      }
   }
   ATH_MSG_INFO("End commitOutput loop");
   return StatusCode::SUCCESS;
}
//___________________________________________________________________________
StatusCode AthenaRootSharedWriterSvc::stop() {
   m_rootMergers.Delete();
   return StatusCode::SUCCESS;
}
//___________________________________________________________________________
StatusCode AthenaRootSharedWriterSvc::finalize() {
   ATH_MSG_INFO("in finalize()");
   delete m_rootMonitor; m_rootMonitor = nullptr;
   delete m_rootServerSocket; m_rootServerSocket = nullptr;
   return StatusCode::SUCCESS;
}
//___________________________________________________________________________
StatusCode AthenaRootSharedWriterSvc::queryInterface(const InterfaceID& riid, void** ppvInterface) {
   if ( IAthenaSharedWriterSvc::interfaceID().versionMatch(riid) ) {
      *ppvInterface = (IAthenaSharedWriterSvc*)this;
   } else {
      // Interface is not directly available: try out a base class
      return(AthService::queryInterface(riid, ppvInterface));
   }
   addRef();
   return(StatusCode::SUCCESS);
}
