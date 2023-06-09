#!/bin/sh

get_files -jo TrigInDetValidation_RTT_topOptions_MuonSlice.py
athena.py  -c 'XMLDataSet="TrigInDetValidation_mu_bphys";EventMax=-1;runMergedChain=True' TrigInDetValidation_RTT_topOptions_MuonSlice.py 


get_files -data TIDAdata11-rtt.dat
get_files -data Test_bin.dat
TIDArdict.exe TIDAdata11-rtt.dat -f data-muon-bphys-merge.root -p 13 -b Test_bin.dat

get_files -data data-mu_bphys-reference.root
TIDAcomparitor.exe data-muon-bphys-merge.root data-mu_bphys-reference.root HLT_mu6_idperf_InDetTrigTrackingxAODCnv_Muon_FTF HLT_mu6_L2Star_idperf_TrigL2SiTrackFinder_Muon_1 HLT_mu6_L2Star_idperf_TrigL2SiTrackFinder_Muon_2 -d HLTL2-plots

get_files -data data-mu_bphys-reference.root
TIDAcomparitor.exe data-muon-bphys-merge.root data-mu_bphys-reference.root HLT_mu6_idperf_InDetTrigTrackingxAODCnv_Muon_FTF HLT_mu6_idperf_InDetTrigTrackingxAODCnv_Muon_EFID HLT_mu6_L2Star_idperf_InDetTrigTrackingxAODCnv_Muon_EFID -d HLTEF-plots

get_files -data expert-monitoring-mu_bphys-ref.root
TIDAcpucost.exe expert-monitoring.root expert-monitoring-mu_bphys-ref.root --auto -o times

RunTrigCostD3PD.exe -f trig_cost.root --outputTagFromAthena --monitorAllChainSeqAlg --monitorROI --linkOutputDir

TIDAcpucost.exe costMon/TrigCostRoot_Results.root costMon/TrigCostRoot_Results.root -o cost-perCall --auto -d "/Algorithm" -p "_Time_perCall"

TIDAcpucost.exe costMon/TrigCostRoot_Results.root costMon/TrigCostRoot_Results.root -o cost-perEvent --auto -d "/Algorithm" -p "_Time_perEvent"

TIDAcpucost.exe costMon/TrigCostRoot_Results.root costMon/TrigCostRoot_Results.root -o cost-perCall-chain --auto -d "/Chain_Algorithm" -p "_Time_perCall"

TIDAcpucost.exe costMon/TrigCostRoot_Results.root costMon/TrigCostRoot_Results.root -o cost-perEvent-chain --auto -d "/Chain_Algorithm" -p "_Time_perEvent" check_log.pl --conf checklogTriggerTest.conf %JOBLOG


