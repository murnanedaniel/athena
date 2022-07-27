/*
   Copyright (C) 2002-2021 CERN for the benefit of the ATLAS collaboration
 */

// System include(s):
#include <memory>
// ROOT include(s):
#include <TFile.h>
// Infrastructure include(s):
#ifdef XAOD_STANDALONE
#include "xAODRootAccess/Init.h"
#include "xAODRootAccess/TEvent.h"
#include "xAODRootAccess/TStore.h"
#endif // XAOD_STANDALONE

// EDM include(s):
#include "ElectronEfficiencyCorrection/AsgElectronEfficiencyCorrectionTool.h"
#include "Messaging.h"
#include "SFHelpers.h"
#include "xAODEgamma/Electron.h"
#include "xAODEgamma/ElectronContainer.h"

// To disable sending data
#include "xAODRootAccess/tools/TFileAccessTracer.h"
int
main(int argc, char* argv[])
{

  xAOD::TFileAccessTracer::enableDataSubmission(false);
  // The application's name:
  const char* APP_NAME = argv[0];

  MSG::Level mylevel = MSG::INFO;
  MSGHELPERS::getMsgStream().msg().setLevel(mylevel);
  MSGHELPERS::getMsgStream().msg().setName(APP_NAME);

  // Check if we received a file name:
  if (argc < 2) {
    MSG_ERROR(APP_NAME << "No file name received!");
    MSG_ERROR(
      APP_NAME
      << "  Usage: <<APP_NAME <<  [xAOD file name] [Num of events to use]");
    return 1;
  }

  // Initialise the application:
  CHECK(xAOD::Init(APP_NAME));

  // Open the input file:
  const TString fileName = argv[1];
  MSG_INFO("Opening file: " << fileName.Data());
  std::unique_ptr<TFile> ifile(TFile::Open(fileName, "READ"));
  CHECK(ifile.get());

  // Create a TEvent object:
  // xAOD::TEvent event( xAOD::TEvent::kBranchAccess );
  xAOD::TEvent event;

  // Then the tools
  AsgElectronEfficiencyCorrectionTool ElEffCorrectionTool(
    "ElEffCorrectionTool");
  CHECK(ElEffCorrectionTool.setProperty("IdKey", "Medium"));
  CHECK(ElEffCorrectionTool.setProperty("ForceDataType", 1));
  CHECK(ElEffCorrectionTool.setProperty("OutputLevel", mylevel));
  CHECK(ElEffCorrectionTool.setProperty("CorrelationModel", "FULL"));
  CHECK(ElEffCorrectionTool.setProperty("UseRandomRunNumber", false));
  CHECK(ElEffCorrectionTool.initialize());

  // Then open the file(s)
  CHECK(event.readFrom(ifile.get()));
  MSG_INFO("Number of available events to read in:  "
           << static_cast<long long int>(event.getEntries()));

  // Decide how many events to run over:
  long long int entries = event.getEntries();
  if (argc > 2) {
    const long long int e = atoll(argv[2]);
    if (e < entries) {
      entries = e;
    }
  }
  MSG_INFO("Number actual events to read in:  " << entries);

  // Loop over the events:
  for (long long int entry = 0; entry < entries; ++entry) {
    event.getEntry(entry);
    MSG_INFO(" \n ==> Event " << entry);

    const xAOD::ElectronContainer* electrons = 0;
    CHECK(event.retrieve(electrons, "Electrons"));

    for (const xAOD::Electron* el : *electrons) {
      const xAOD::CaloCluster* cluster = el->caloCluster();
      if (!cluster) {
	MSG_ERROR("ERROR no cluster associated to the Electron \n");
	return CP::CorrectionCode::Error;
      }
      if (el->pt() < 7000)
        continue; // skip electrons outside of recommendations
      if (std::abs(cluster->etaBE(2))>=2.47) {
	continue;
      }

      int index = ElEffCorrectionTool.systUncorrVariationIndex(*el);
      /*
       * Set up the systematic variations
       */
      bool isToys = false;
      double nominalSF{};
      double totalNeg{};
      double totalPos{};
      CHECK(
        SFHelpers::result(
          ElEffCorrectionTool, *el, nominalSF, totalPos, totalNeg, isToys) ==
        0);

      MSG_INFO("===> electron : Pt = "
               << el->pt() << " : eta = " << el->eta()
               << " : Bin index = " << index << " : SF = " << nominalSF << " + "
               << totalPos << " - " << totalNeg << " <===");
    }
  }

  MSG_INFO("===> DONE <===\n");
  return 0;
}
