/*
  Copyright (C) 2002-2021 CERN for the benefit of the ATLAS collaboration
*/

#define CaloD3PDClass_cxx
#include "LArSamplesMon/CaloD3PDClass.h"
#include <TH2.h>
#include <TStyle.h>
#include <TCanvas.h>

void CaloD3PDClass::Loop()
{
//   In a ROOT session, you can do:
//      Root > .L CaloD3PDClass.C
//      Root > CaloD3PDClass t
//      Root > t.GetEntry(12); // Fill t data members with entry number 12
//      Root > t.Show();       // Show values of entry 12
//      Root > t.Show(16);     // Read and show values of entry 16
//      Root > t.Loop();       // Loop on all entries
//

//     This is the loop skeleton where:
//    jentry is the global entry number in the chain
//    ientry is the entry number in the current Tree
//  Note that the argument to GetEntry must be:
//    jentry for TChain::GetEntry
//    ientry for TTree::GetEntry and TBranch::GetEntry
//
//       To read only selected branches, Insert statements like:
// METHOD1:
//    fChain->SetBranchStatus("*",0);  // disable all branches
//    fChain->SetBranchStatus("branchname",1);  // activate branchname
// METHOD2: replace line
//    fChain->GetEntry(jentry);       //read all branches
//by  b_branchname->GetEntry(ientry); //read only this branch
   if (fChain == nullptr) return;

   Long64_t nentries = fChain->GetEntriesFast();

   for (Long64_t jentry=0; jentry<nentries;jentry++) {
      Long64_t ientry = LoadTree(jentry);
      if (ientry < 0) break;
      fChain->GetEntry(jentry);
      // if (Cut(ientry) < 0) continue;
   }
}
