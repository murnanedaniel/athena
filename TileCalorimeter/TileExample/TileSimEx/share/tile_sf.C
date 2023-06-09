//-----------------------------------------------------------------------
//     
//    root macro for RTT tests - sampling fraction value
//    
//
//-----------------------------------------------------------------------

#include <stdlib.h>
#include <vector>

#include "TSystem.h"
#include "TFile.h"
#include "TTree.h"
#include "TChain.h"
#include "TH1.h"
#include "TH2.h"
#include "TProfile.h"
#include "TF1.h"
#include "TROOT.h"
#include "TStyle.h"
#include "TString.h"
#include "TCanvas.h"
#include "TVector3.h"
#include "TPaveStats.h"


void Make_plots(TString infile_hit, TString infile_truth, TString outfile, bool dofit) {

  TTree *t = new TTree();
  t->AddFriend("h32",infile_hit);
  t->AddFriend("truth",infile_truth);
  //t->AddFriend("truth=tile_ntuple",infile_truth);

   //function - dependence of the total energy in scintillators on the impact point
  Float_t thefun_p0   = 3000.;
  Float_t thefun_p1   = 50.;
  Float_t thefun_p2   = 18.;
  Float_t thefun_p3   = 0.;
    
  TF1 * thefun = new TF1("thefun","[0]+ [1]*cos((x-[3])*(2*TMath::Pi()/[2]))");
   
    //************************************************************************
    // Histo definitions
    //************************************************************************
   
    TH1F* h_totalE = new TH1F("h_totalE","Total hit energy in scintillators",600,0.,6000.);
    TH1F* h_ptgen = new TH1F("h_ptgen","Pt of the generated electron",100,0.,120000.);
    TProfile2D* hp_EvsYvsZ = new TProfile2D("hp_EvsYvsZ","Total hit energy vs Y & Z",50,-25.,25.,50,-25.,25.,0.,100000.);   
    TProfile* hp_totalEvsZ = new TProfile("hp_totalEvsZ","Total hit energy vs Z",40,-20.,20.,0.,100000.);   
    
    TObjArray HistArray(0);
 
    TCanvas *c_g4=new TCanvas("c_g4","Tile TB simulation plots");
    c_g4->Divide(2,2);

    c_g4->cd(1);
    t->Draw("truth.mcpart_pt>>h_ptgen","","histo");
    h_ptgen->SetXTitle("Energy (MeV)");
    HistArray.Add(h_ptgen);
    
    c_g4->cd(2);
    t->Draw("h32.totalE>>h_totalE","","histo");
    thefun_p0 = h_totalE->GetMean();
    thefun_p1 = h_totalE->GetRMS();
    
    //h_totalE->Fit("gaus","","");
    h_totalE->SetXTitle("Energy (MeV)");
    HistArray.Add(h_totalE);   
   
    c_g4->cd(3);
    t->Draw("(h32.totalE):(truth.mcvtx_y):(truth.mcvtx_z)>>hp_EvsYvsZ","","prof");
    hp_EvsYvsZ->SetXTitle("ZVGen (mm)");
    hp_EvsYvsZ->SetYTitle("YVGen (mm)");
    hp_EvsYvsZ->Draw("COLZ");
    //hp_EvsYvsZ->GetZaxis()->SetRangeUser(1000,5000);
    HistArray.Add(hp_EvsYvsZ);   

    c_g4->cd(4);
    t->Draw("(h32.totalE):(truth.mcvtx_z)>>hp_totalEvsZ","","prof");
    hp_totalEvsZ->SetXTitle("ZVGen (mm)");
    hp_totalEvsZ->SetYTitle("Energy (MeV)");
    hp_totalEvsZ->GetYaxis()->SetRangeUser(1000,5000);

    thefun->SetParameter(0,thefun_p0);
    thefun->SetParameter(1,thefun_p1);
    thefun->SetParameter(2,thefun_p2);
    thefun->SetParameter(3,thefun_p3);

    thefun->SetParLimits(0,0.,9000.);
    thefun->SetParLimits(1,0.,2000.);
    thefun->SetParLimits(2,1.,20.);
    thefun->SetParLimits(3,0.,20.);
    
    if (dofit) {
        hp_totalEvsZ->Fit("thefun","","",-15.,15.);
    }
    
    cout << "file name " << infile_hit << endl;
    cout << "beamE/total_E(hit) - no z correction " <<   (h_ptgen->GetMean())/(h_totalE->GetMean()) << endl;
    if (dofit) {
        cout << "total_E(hits) with z correction " << thefun->GetParameter(0) << "+/-" << thefun->GetParError(0) << " dErr % " << thefun->GetParError(0)/thefun->GetParameter(0)*100 << endl;
        cout << "1/SF value: " << (h_ptgen->GetMean())/(thefun->GetParameter(0)) << "+/-" << (h_ptgen->GetMean()/thefun->GetParameter(0))*(thefun->GetParError(0)/thefun->GetParameter(0)) << " dErr % " <<  (thefun->GetParError(0)/thefun->GetParameter(0))*100 << endl;   
    }
    
    HistArray.Add(hp_totalEvsZ);   

    c_g4->Update();
    c_g4->cd();

    if (outfile.Length()>0) {
      c_g4->SaveAs(outfile+".png");
      TFile fout(outfile+".root","recreate");
      HistArray.Write(0);
      fout.Close();
    }

}

int tile_sf(TString suff="", bool dofit=true ) {

  gROOT->SetStyle("Plain");
  gStyle->SetPalette(51,0);    
  gStyle->SetPadBottomMargin(0.13);
  gStyle->SetPadTopMargin   (0.20);
  gStyle->SetStatFontSize  (0.05);
  gStyle->SetStatStyle     (0);
  gStyle->SetOptStat(111111);  
  gStyle->SetOptFit(1111);  
  gStyle->SetTitleH(0.1);
  gStyle->SetCanvasDefW(900);
  gStyle->SetCanvasDefH(600);

  TString directory = "./";
  TString file_hit = directory+"tiletb"+suff+".ntup.root";
  TString file_truth = directory+"tiletb"+suff+".d3pd.root";

  cout << file_hit << endl;
  cout << file_truth << endl;

  directory = "./";
  TString output_filename = directory+"hist"+suff;
  cout << output_filename << endl;
  
  Make_plots(file_hit, file_truth, output_filename, dofit);

  return 1;  

}



