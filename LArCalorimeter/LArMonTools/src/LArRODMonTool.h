/*
  Copyright (C) 2002-2022 CERN for the benefit of the ATLAS collaboration
*/

//Dear emacs, this is -*-c++-*-

/**
 * @class LArRODMonTool
 * @author Mathieu Aurousseau <mathieu.aurousseau@lapp.in2p3.fr>
 *
 */

#ifndef LARMONTOOLS_LARRODMONTOOL_H
#define LARMONTOOLS_LARRODMONTOOL_H

#include "AthenaMonitoring/ManagedMonitorToolBase.h"
#include "LArElecCalib/ILArPedestal.h"  //in template
#include "LArElecCalib/ILArOFC.h"  //in template
#include "LArElecCalib/ILArShape.h"  //in template
#include "LArElecCalib/ILArHVScaleCorr.h"  //in template
#include "LArIdentifier/LArOnlineID.h"  //inline usage
#include "LArRawConditions/LArADC2MeV.h"  //in template
#include "StoreGate/ReadCondHandleKey.h"
#include "StoreGate/ReadHandleKey.h"
#include "CaloDetDescr/CaloDetDescrManager.h"  //in template

//Events infos:

#include "LArRawEvent/LArFebHeaderContainer.h"  //in template
//Events infos:
#include "xAODEventInfo/EventInfo.h"  //in template

#include "LArRecConditions/LArBadChannelMask.h"  //member
#include "LArRecConditions/LArBadChannelCont.h"  //in template

#include <iosfwd>
#include <cstdio>  //FILE* typedef
#include <ios>  //for std::hex
#include <fstream> //for std::ofstream
#include <set>
#include <vector>

class LArRawChannel;
class HWIdentifier;
class LArOnOffIdMapping;
class LArRawChannelContainer;
class LArDigitContainer;

class LArEM_ID;
class TH1F;
class TH2I;
class TH2F;
class TProfile2D;
class LArDigit;

//using namespace std;

class LArRODMonTool: public ManagedMonitorToolBase
{
 public:
  LArRODMonTool(const std::string& type, 
		const std::string& name,
		const IInterface* parent);

  /** @brief Default destructor */
  virtual ~LArRODMonTool();

  virtual StatusCode initialize() override;

  virtual StatusCode finalize() override;

  // Book general histograms
  virtual StatusCode bookHistograms() override;

  // Called each event
  virtual StatusCode fillHistograms() override;

  //StatusCode procHistograms(bool isEndOfEventsBlock, bool isEndOfLumiBlock, bool isEndOfRun);

  //bool Checksum_Check(std::vector<int>& myFEBvector);

  //bool RodStatus_Check(std::vector<int>& myStatusvector);

  /** @brief For a given cell, dump info that can be used by DspTest, into two txt files and one dat file */
  void DumpCellEvent(int count,                            // counter
                     ILArOFC::OFCRef_t ofc,                // vector containing OFC's
                     ILArOFC::OFCRef_t ofcb,               // vector containing OFCb's
		     ILArShape::ShapeRef_t ofch,           // vector containing OFCh's
		     ILArShape::ShapeRef_t ofcd,           // vector containing OFCd's
		     float Escale,                         // ADC2MeV factor
		     float Ramp0,                          // ramp intercept
		     float pedestal,                       // pedestal
		     const std::vector<short>* digits,     // vector containing samples (digits)
		     int E_off,                            // Energy calculated offline (i.e. fromDigits)
		     int E_on,                             // Energy calculated online (i.e. fromBytestream)
		     FILE* ofcfile,                        // output file containing ofc's
		     std::ofstream& digitsfile,                 // output file containing digits
		     std::ofstream& energyfile,                 // output file containing energies
		     const HWIdentifier chid,                   // FEB HW ID
		     int event);                           // Event #
  
  
  /** @brief Dump a cell's information and calculated energies into a txt file */
  void DumpCellInfo(HWIdentifier chid,                    // Channel HW ID
		    const char* partc,
		    int gain,
		    int event, 
		    int E_off,                             // Energy calculated offline (i.e. fromDigits)  
		    int E_on,                              // Energy calculated online (i.e. fromBytestream)
		    float T_off,                             // Time calculated offline (i.e. fromDigits)  
		    float T_on,                              // Time calculated online (i.e. fromBytestream)
		    float Q_off,                             // Quality Factor calculated offline (i.e. fromDigits)  
		    float Q_on,                              // Quality Factor calculated online (i.e. fromBytestream)
		    FILE* dumpfile);
  
private:
  const LArOnlineID* m_LArOnlineIDHelper;
  unsigned m_counter;
  unsigned m_eventsCounter;
  unsigned m_ndump;
  enum PARTITION {EMBC=0,EMBA,EMECC,EMECA,HECC,HECA,FCALC,FCALA,N_PARTITIONS};

  


  PARTITION getPartition(const HWIdentifier chid) const;

  struct HistGroup{
  public:
    HistGroup(): 
      m_hDE(NULL),m_hDT(NULL), m_hDQ(NULL),		
      m_hOut_E_FT_vs_SLOT(NULL), m_hOut_T_FT_vs_SLOT(NULL),m_hOut_Q_FT_vs_SLOT(NULL),	
      m_hOut_E_FT_vs_SLOT_shadow(NULL), m_hOut_T_FT_vs_SLOT_shadow(NULL),m_hOut_Q_FT_vs_SLOT_shadow(NULL),	
      m_hEon_VS_Eoff(NULL),m_hTon_VS_Toff(NULL),	
      m_hQon_VS_Qoff(NULL),m_Sweet(NULL),
      m_hDE_ranges(0),
      m_monGroup(0)
    {};

    TH1F* m_hDE;		
    TH1F* m_hDT;		
    TH1F* m_hDQ;		
    TH2F* m_hOut_E_FT_vs_SLOT;	
    TH2F* m_hOut_T_FT_vs_SLOT;	
    TH2F* m_hOut_Q_FT_vs_SLOT;	
    TH2F** m_hOut_E_FT_vs_SLOT_shadow;	
    TH2F** m_hOut_T_FT_vs_SLOT_shadow;	
    TH2F** m_hOut_Q_FT_vs_SLOT_shadow;	
    TH2F* m_hEon_VS_Eoff;	
    TH2F* m_hTon_VS_Toff;	
    TH2F* m_hQon_VS_Qoff;	
    TH1F* m_Sweet;
    TH2F* m_hDE_ranges;
    std::string m_partName;
    MonGroup* m_monGroup;
  };

  std::vector<HistGroup> m_histos;

  std::vector<unsigned> m_errsPerFEB;

  // initialize monitoring bookkeeping info
  //StatusCode initMonInfo();

  bool FebStatus_Check();

  StatusCode compareChannels(const CaloDetDescrManager* ddman,
                             const LArOnOffIdMapping* cabling,
                             const ILArOFC* ofcs,
                             const ILArShape* shapes,
                             const ILArHVScaleCorr* hvScaleCorrs,
                             const ILArPedestal* pedestals,
                             const LArADC2MeV* adc2mev,
                             const HWIdentifier chid, 
			     const LArRawChannel& rcDig, 
			     const LArRawChannel& rcBS, 
			     const LArDigit* dig=nullptr

			     );

  class ERRCOUNTER {
  public:
    ERRCOUNTER() {
      clear();
    };
    void clear();
    unsigned errors_E[3];
    unsigned errors_T[3];
    unsigned errors_Q[3];
  };
  
  std::vector<ERRCOUNTER> m_errcounters;

  unsigned m_count_gain[3]{};

  TH2F* m_hSummaryErrors_Energy;
  TH2F* m_hSummaryErrors_Time;
  TH2F* m_hSummaryErrors_Quality;

  TH1F* m_hE_all;
  TH1F* m_hT_all;
  TH1F* m_hQ_all;

  TH2F* m_hE_ranges_all;

  TH2I* m_hEErrors_LB_part;
  TH2I* m_hTErrors_LB_part;
  TH2I* m_hQErrors_LB_part;

  TH2I* m_hEErrors_LB_stream;
  TH2I* m_hTErrors_LB_stream;
  TH2I* m_hQErrors_LB_stream;

 
  SG::ReadHandleKey<xAOD::EventInfo> m_eventInfoKey{this,"EventInfoKey","EventInfo","SG Key of EventInfo object"};
  SG::ReadHandleKey<LArFebHeaderContainer> m_febContKey{this, "LArFebHeaderContainerKey", "LArFebHeader"};
  // Keys for LArRawChannels containers
  SG::ReadHandleKey<LArRawChannelContainer> m_channelKey_fromBytestream{this,"LArRawChannelKey_fromBytestream","LArRawChannels","SG key of LArRawChannels produced by teh DSP"};
  SG::ReadHandleKey<LArRawChannelContainer> m_channelKey_fromDigits{this,"LArRawChannelKey_fromDigits","LArRawChannels_FromDigits","SG key of LArRawChannels produced offline"};

  SG::ReadHandleKey<LArDigitContainer> m_digitContainerKey{this,"LArDigitContainerKey","FREE","SG key of LArDigitContainer read from Bytestream"};


  SG::ReadCondHandleKey<ILArOFC>         m_keyOFC{this,"KeyOFC","LArOFC","SG key of LArOFC CDO"};
  SG::ReadCondHandleKey<ILArShape>       m_keyShape{this,"KeyShape","LArShape","SG key of LArShape CDO"};
  SG::ReadCondHandleKey<ILArHVScaleCorr> m_keyHVScaleCorr{this,"KeyHVScaleCorr","LArHVScaleCorr","SG key of LArHVScaleCorr CDO"};
  SG::ReadCondHandleKey<ILArPedestal>    m_keyPedestal{this,"LArPedestalKey","LArPedestal","SG key of LArPedestal CDO"};

  SG::ReadCondHandleKey<LArADC2MeV> m_adc2mevKey{this,"LArADC2MeVKey","LArADC2MeV","SG Key of the LArADC2MeV CDO"};


  LArBadChannelMask m_bcMask;
  SG::ReadCondHandleKey<LArBadChannelCont> m_bcContKey {this, "BadChanKey", "LArBadChannel", "SG key for LArBadChan object"};
  Gaudi::Property<std::vector<std::string> > m_problemsToMask{this,"ProblemsToMask",{}, "Bad-Channel categories to mask"}; 
  

  SG::ReadCondHandleKey<LArOnOffIdMapping> m_cablingKey{this,"CablingKey","LArOnOffIdMap","SG Key of LArOnOffIdMapping CDO"};
  SG::ReadCondHandleKey<CaloDetDescrManager> m_caloMgrKey {this, "CaloDetDescrManager", "CaloDetDescrManager", "SG Key for CaloDetDescrManager in the Condition Store" };


  // Output files names
  std::string m_DigitsFileName;
  std::string m_EnergyFileName;
  std::string m_AiFileName;
  std::string m_DumpCellsFileName;

  // Ofstream / files for test outputs
  std::ofstream m_fdig, m_fen; // ,ai
  FILE* m_fdump;
  FILE* m_fai;

  void closeDumpfiles();

  bool m_doDspTestDump;
  bool m_doCellsDump;
  bool m_doCheckSum;
  bool m_doRodStatus;
  bool m_printEnergyErrors;
  bool m_removeNoiseBursts;
  bool m_useEvtCounter;

  int m_range_E_0;
  int m_range_E_1;
  int m_range_E_2;
  int m_range_E_3;

  int m_range_T_0 = 0;
  int m_range_T_1 = 0;
  int m_range_T_2 = 0;
  int m_range_T_3 = 0;

  int m_range_Q_0 = 0;
  int m_range_Q_1 = 0;
  int m_range_Q_2 = 0;
  int m_range_Q_3 = 0;

  int m_precision_E_0;
  int m_precision_E_1;
  int m_precision_E_2;
  int m_precision_E_3;
  int m_precision_E_max;

  int m_precision_T_0;
  int m_precision_T_1;
  int m_precision_T_2;
  int m_precision_T_3;
  int m_precision_T_max;

  int m_precision_Q_0;
  int m_precision_Q_1;
  int m_precision_Q_2;
  int m_precision_Q_3;
  int m_precision_Q_max;

  int m_history_size;
  int m_history_granularity;
  int m_hsize;
  //bool *m_hdone;
  std::vector<bool> m_hdone;

  bool m_skipKnownProblematicChannels;
  bool m_skipNullPed;
  bool m_skipNullQT;

  float m_unit_online;
  float m_unit_offline;
  float m_timeOffset;
  short m_adc_th;
  float m_peakTime_cut;
  float m_nb_lb;

  float m_BC; // value of 1 bunch-crossing = 25ns

  bool m_IsOnline;
  int m_last_lb;
  int m_curr_lb;

  std::vector<std::string> m_streams;

  std::set<HWIdentifier> m_ignoreFEBs;

  bool m_dumpDigits;

};


inline 
LArRODMonTool::PARTITION LArRODMonTool::getPartition(const HWIdentifier chid) const{
  const int side=m_LArOnlineIDHelper->pos_neg(chid);
  if (m_LArOnlineIDHelper->isEMBchannel(chid)) {
    if (side==0)
      return EMBC;
    else
      return EMBA;
  }
  if (m_LArOnlineIDHelper->isEMECchannel(chid)) {
     if (side==0)
       return EMECC;
     else
       return EMECA;
  }
  if (m_LArOnlineIDHelper->isHECchannel(chid)){
    if (side==0)
      return HECC;
    else
      return HECA;
  }
  if (m_LArOnlineIDHelper->isFCALchannel(chid)) {
    if (side==0)
      return FCALC;
    else
      return FCALA;
  }

  msg(MSG::FATAL) << "Channel 0x "<< std::hex << chid.get_identifier32().get_compact() << std::dec << " neither EMB nor EMEC nor HEC nor FCAL???" << std::endl;
  return N_PARTITIONS;
}


#endif
