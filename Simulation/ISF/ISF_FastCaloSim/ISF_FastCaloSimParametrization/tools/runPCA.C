/*
  Copyright (C) 2002-2017 CERN for the benefit of the ATLAS collaboration
*/

{

 //!!!remember to source the hacked tmva!!!
 //export LD_LIBRARY_PATH=/afs/cern.ch/user/s/schaarsc/public/fastcalo/trash/tests/mywork/TMVA-v4.2.0/lib:${LD_LIBRARY_PATH}
 //to do: migrate to root6, where this hack is implemented

 gSystem->AddIncludePath(" -I.. ");
 gROOT->LoadMacro("../Root/IntArray.cxx+");
 gROOT->LoadMacro("../Root/TreeReader.cxx+");
 gROOT->LoadMacro("../Root/firstPCA.cxx+");
 gROOT->LoadMacro("../Root/TFCS1DFunction.cxx+");
 gROOT->LoadMacro("../Root/TFCS1DFunctionRegression.cxx+");
 gROOT->LoadMacro("../Root/TFCS1DFunctionRegressionTF.cxx+");
 gROOT->LoadMacro("../Root/TFCS1DFunctionHistogram.cxx+");
 gROOT->LoadMacro("../Root/TFCSFunction.cxx+");
 gROOT->LoadMacro("../Root/secondPCA.cxx+");
 
 string label;
 vector<string> input;
 
 //label="pions";
 //input.push_back("root://eosatlas//eos/atlas/user/z/zhubacek/FastCaloSim/NTUP_090315/ISF_HitAnalysis_evgen_calo__211_E50000_50000_eta20_25_Evts0-5500_vz_0_origin_calo.standard.matched.pool.root");
  
 label="s2864";
 input.push_back("root://eosatlas//eos/atlas/user/s/schaarsc/FCS/user.fladias.428137.FastCalo_pid11_E65536_etam35_35_zv_m100.e4001_s2864_r7736.w0_162706_matched_output.root/*.root");
 
 //label="s2865";
 //input.push_back("root://eosatlas//eos/atlas/user/s/schaarsc/FCS/user.fladias.428137.FastCalo_pid11_E65536_etam35_35_zv_m100.e4001_s2865_r7736.w0_162706_matched_output.root/*.root");
 
 cout<<"*** Preparing to run on "<<label<<" ***"<<endl;
 TChain* mychain= new TChain("FCS_ParametrizationInput");
 for(int i=0;i<input.size();i++)
 {
 	cout<<"input: "<<input[i]<<endl;
  mychain->Add(input[i].c_str());
 }
 cout<<"TChain entries: "<<mychain->GetEntries()<<endl;
 
 system(("mkdir /afs/cern.ch/user/s/schaarsc/public/fastcalo/pca/tools/output/"+label).c_str());
 string pca1_outfilename="/afs/cern.ch/user/s/schaarsc/public/fastcalo/pca/tools/output/"+label+"/firstPCA.root";
 string pca2_outfilename="/afs/cern.ch/user/s/schaarsc/public/fastcalo/pca/tools/output/"+label+"/secondPCA.root";
 
 firstPCA *myfirstPCA=new firstPCA(mychain,pca1_outfilename);
 myfirstPCA->set_cumulativehistobins(5000);
 myfirstPCA->set_edepositcut(0.001);
 myfirstPCA->set_etacut(0.2,0.25);
 myfirstPCA->apply_etacut(0); //this flag is for the old files, which are already sliced in eta
 myfirstPCA->set_pcabinning(5,1);
 myfirstPCA->run();
 
 secondPCA* mysecondPCA=new secondPCA(pca1_outfilename,pca2_outfilename);
 mysecondPCA->set_PCAbin(0); //all bins
 mysecondPCA->set_storeDetails(0);
 mysecondPCA->set_cumulativehistobins(5000);
 mysecondPCA->set_cut_maxdeviation_regression(5);
 mysecondPCA->set_cut_maxdeviation_smartrebin(5);
 mysecondPCA->set_Ntoys(1000);
 mysecondPCA->set_neurons_iteration(2,10);
 mysecondPCA->set_skip_regression(1);
 mysecondPCA->run();
 
 //add validation class
 
 
}

