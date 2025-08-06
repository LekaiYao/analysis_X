#define analysis_PSI2S_cxx
#include "analysis_PSI2S.h"
#include <TH2.h>
#include <TStyle.h>
#include <TCanvas.h>

void analysis_PSI2S::Loop()
{
//   In a ROOT session, you can do:
//      root> .L analysis_PSI2S.C
//      root> analysis_PSI2S t
//      root> t.GetEntry(12); // Fill t data members with entry number 12
//      root> t.Show();       // Show values of entry 12
//      root> t.Show(16);     // Read and show values of entry 16
//      root> t.Loop();       // Loop on all entries
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
   if (fChain == 0) return;
   TFile* outputFile = new TFile("MC_PSI2S.root", "RECREATE");
   TTree* tree = new TTree("tree", "tree");

   //define the variables
   Int_t           B_type;
   Float_t         B_mass;
   Float_t         B_pt;
   Float_t         B_y;
   Float_t         B_vtxX;
   Float_t         B_vtxY;
   Float_t         B_d0;
   Float_t         B_d0Err;
   Float_t         B_dxyz;
   Float_t         B_dxyzErr;
   Float_t         B_chi2ndf;
   Float_t         B_chi2cl;
   Float_t         B_dtheta;
   Float_t         B_cos_dtheta;
   Float_t         B_alpha;
   Float_t         B_svpvDistance;
   Float_t         B_svpvDisErr;
   Float_t         B_svpvDistance_2D;
   Float_t         B_svpvDisErr_2D;
   Float_t         B_norm_svpvDistance;
   Float_t         B_norm_svpvDistance_2D;
   Float_t         B_Qvalue;
   Float_t         B_Qvalueuj;
   Float_t         B_Qvaluemumu;
   Float_t         B_trk1Pt;
   Float_t         B_trk2Pt;
   Float_t         B_trkPtimb;
   Float_t         B_trk1Eta;
   Float_t         B_trk2Eta;
   Float_t         B_trk1PtErr;
   Float_t         B_trk2PtErr;
   Float_t         B_trk1Dz;
   Float_t         B_trk2Dz;
   Float_t         B_trk1DzError;
   Float_t         B_trk2DzError;
   Float_t         B_trk1Dxy;
   Float_t         B_trk2Dxy;
   Float_t         B_norm_trk1Dxy;
   Float_t         B_norm_trk2Dxy;
   Float_t         B_trk1DxyError;
   Float_t         B_trk2DxyError;
   Float_t         B_trk1nPixelLayer;
   Float_t         B_trk2nPixelLayer;
   Float_t         B_trk1nStripLayer;
   Float_t         B_trk2nStripLayer;
   Float_t         B_trk1Chi2ndf;
   Float_t         B_trk2Chi2ndf;
   Bool_t          B_trk1highPurity;
   Bool_t          B_trk2highPurity;
   Float_t         B_trk1dR;
   Float_t         B_trk2dR;
   Float_t         B_tktkmass;
   Float_t         B_tktkmassKK;
   Float_t         B_tktkvProb;
   Float_t         B_tktkpt;
   Float_t         B_doubletmass;
   Float_t         B_mu1pt;
   Float_t         B_mu2pt;
   Bool_t          B_mu1highPurity;
   Bool_t          B_mu2highPurity;
   Bool_t          B_mu1isAcc;
   Bool_t          B_mu2isAcc;
   Bool_t          B_mu1SoftMuID;
   Bool_t          B_mu2SoftMuID;
   Bool_t          B_mu1HybridSoftMuID;
   Bool_t          B_mu2HybridSoftMuID;
   Bool_t          B_mu1isTriggered;
   Bool_t          B_mu2isTriggered;
   Float_t         B_mumumass;
   Float_t         B_ujmass;
   Float_t         B_ujvProb;
   Float_t         B_gen;
   Int_t           B_genIndex;

   //define the branches
   tree->Branch("B_type", &B_type, "B_type/I");
   tree->Branch("B_mass", &B_mass, "B_mass/F");
   tree->Branch("B_pt", &B_pt, "B_pt/F");
   tree->Branch("B_y", &B_y, "B_y/F");
   tree->Branch("B_vtxX", &B_vtxX, "B_vtxX/F");
   tree->Branch("B_vtxY", &B_vtxY, "B_vtxY/F");
   tree->Branch("B_d0", &B_d0, "B_d0/F");
   tree->Branch("B_d0Err", &B_d0Err, "B_d0Err/F");
   tree->Branch("B_dxyz", &B_dxyz, "B_dxyz/F");
   tree->Branch("B_dxyzErr", &B_dxyzErr, "B_dxyzErr/F");
   tree->Branch("B_chi2ndf", &B_chi2ndf, "B_chi2ndf/F");
   tree->Branch("B_chi2cl", &B_chi2cl, "B_chi2cl/F");
   tree->Branch("B_dtheta", &B_dtheta, "B_dtheta/F");
   tree->Branch("B_cos_dtheta", &B_cos_dtheta, "B_cos_dtheta/F");
   tree->Branch("B_alpha", &B_alpha, "B_alpha/F");
   tree->Branch("B_svpvDistance", &B_svpvDistance, "B_svpvDistance/F");
   tree->Branch("B_svpvDisErr", &B_svpvDisErr, "B_svpvDisErr/F");
   tree->Branch("B_svpvDistance_2D", &B_svpvDistance_2D, "B_svpvDistance_2D/F");
   tree->Branch("B_svpvDisErr_2D", &B_svpvDisErr_2D, "B_svpvDisErr_2D/F");
   tree->Branch("B_norm_svpvDistance", &B_norm_svpvDistance, "B_norm_svpvDistance/F");
   tree->Branch("B_norm_svpvDistance_2D", &B_norm_svpvDistance_2D, "B_norm_svpvDistance_2D/F");
   tree->Branch("B_Qvalue", &B_Qvalue, "B_Qvalue/F");
   tree->Branch("B_Qvalueuj", &B_Qvalueuj, "B_Qvalueuj/F");
   tree->Branch("B_Qvaluemumu", &B_Qvaluemumu, "B_Qvaluemumu/F");
   tree->Branch("B_trk1Pt", &B_trk1Pt, "B_trk1Pt/F");
   tree->Branch("B_trk2Pt", &B_trk2Pt, "B_trk2Pt/F");
   tree->Branch("B_trkPtimb", &B_trkPtimb, "B_trkPtimb/F");
   tree->Branch("B_trk1Eta", &B_trk1Eta, "B_trk1Eta/F");
   tree->Branch("B_trk2Eta", &B_trk2Eta, "B_trk2Eta/F");
   tree->Branch("B_trk1PtErr", &B_trk1PtErr, "B_trk1PtErr/F");
   tree->Branch("B_trk2PtErr", &B_trk2PtErr, "B_trk2PtErr/F");
   tree->Branch("B_trk1Dz", &B_trk1Dz, "B_trk1Dz/F");
   tree->Branch("B_trk2Dz", &B_trk2Dz, "B_trk2Dz/F");
   tree->Branch("B_trk1DzError", &B_trk1DzError, "B_trk1DzError/F");
   tree->Branch("B_trk2DzError", &B_trk2DzError, "B_trk2DzError/F");
   tree->Branch("B_trk1Dxy", &B_trk1Dxy, "B_trk1Dxy/F");
   tree->Branch("B_trk2Dxy", &B_trk2Dxy, "B_trk2Dxy/F");
   tree->Branch("B_norm_trk1Dxy", &B_norm_trk1Dxy, "B_norm_trk1Dxy/F");
   tree->Branch("B_norm_trk2Dxy", &B_norm_trk2Dxy, "B_norm_trk2Dxy/F");
   tree->Branch("B_trk1DxyError", &B_trk1DxyError, "B_trk1DxyError/F");
   tree->Branch("B_trk2DxyError", &B_trk2DxyError, "B_trk2DxyError/F");
   tree->Branch("B_trk1nPixelLayer", &B_trk1nPixelLayer, "B_trk1nPixelLayer/F");
   tree->Branch("B_trk2nPixelLayer", &B_trk2nPixelLayer, "B_trk2nPixelLayer/F");
   tree->Branch("B_trk1nStripLayer", &B_trk1nStripLayer, "B_trk1nStripLayer/F");
   tree->Branch("B_trk2nStripLayer", &B_trk2nStripLayer, "B_trk2nStripLayer/F");
   tree->Branch("B_trk1Chi2ndf", &B_trk1Chi2ndf, "B_trk1Chi2ndf/F");
   tree->Branch("B_trk2Chi2ndf", &B_trk2Chi2ndf, "B_trk2Chi2ndf/F");
   tree->Branch("B_trk1highPurity", &B_trk1highPurity, "B_trk1highPurity/O");
   tree->Branch("B_trk2highPurity", &B_trk2highPurity, "B_trk2highPurity/O");
   tree->Branch("B_trk1dR", &B_trk1dR, "B_trk1dR/F");
   tree->Branch("B_trk2dR", &B_trk2dR, "B_trk2dR/F");
   tree->Branch("B_tktkmass", &B_tktkmass, "B_tktkmass/F");
   tree->Branch("B_tktkmassKK", &B_tktkmassKK, "B_tktkmassKK/F");
   tree->Branch("B_tktkvProb", &B_tktkvProb, "B_tktkvProb/F");
   tree->Branch("B_tktkpt", &B_tktkpt, "B_tktkpt/F");
   tree->Branch("B_doubletmass", &B_doubletmass, "B_doubletmass/F");
   tree->Branch("B_mu1pt", &B_mu1pt, "B_mu1pt/F");
   tree->Branch("B_mu2pt", &B_mu2pt, "B_mu2pt/F");
   tree->Branch("B_mu1highPurity", &B_mu1highPurity, "B_mu1highPurity/O");
   tree->Branch("B_mu2highPurity", &B_mu2highPurity, "B_mu2highPurity/O");
   tree->Branch("B_mu1isAcc", &B_mu1isAcc, "B_mu1isAcc/O");
   tree->Branch("B_mu2isAcc", &B_mu2isAcc, "B_mu2isAcc/O");
   tree->Branch("B_mu1SoftMuID", &B_mu1SoftMuID, "B_mu1SoftMuID/O");
   tree->Branch("B_mu2SoftMuID", &B_mu2SoftMuID, "B_mu2SoftMuID/O");
   tree->Branch("B_mu1HybridSoftMuID", &B_mu1HybridSoftMuID, "B_mu1HybridSoftMuID/O");
   tree->Branch("B_mu2HybridSoftMuID", &B_mu2HybridSoftMuID, "B_mu2HybridSoftMuID/O");
   tree->Branch("B_mu1isTriggered", &B_mu1isTriggered, "B_mu1isTriggered/O");
   tree->Branch("B_mu2isTriggered", &B_mu2isTriggered, "B_mu2isTriggered/O");
   tree->Branch("B_mumumass", &B_mumumass, "B_mumumass/F");
   tree->Branch("B_ujmass", &B_ujmass, "B_ujmass/F");
   tree->Branch("B_ujvProb", &B_ujvProb, "B_ujvProb/F");
   tree->Branch("B_gen", &B_gen, "B_gen/F");
   tree->Branch("B_genIndex", &B_genIndex, "B_genIndex/I");

   Long64_t nentries = fChain->GetEntriesFast();

   Long64_t nbytes = 0, nb = 0;
   for (Long64_t jentry=0; jentry<nentries;jentry++) {
      Long64_t ientry = LoadTree(jentry);
      if (ientry < 0) break;
      nb = fChain->GetEntry(jentry);   nbytes += nb;
      /////////////////////////////////////////////////////////////////////
      if(Bsize == 0) continue;
      Long64_t i_B = 0;
      for(i_B=0;i_B<Bsize;i_B++){
         //////////////////// set all the cuts (below)//////////////////
         //isMCsignal
         if(!((Bgen[i_B] == 23333) || (Bgen[i_B] == 24333) 
         || (Bgen[i_B] == 23433) || (Bgen[i_B] == 24433))) continue;
      
         //FIDreg, _RAW
         if(!((Bpt[i_B] >= 2 )&& (By[i_B]  <= 2.4))) continue;

         //ACCcuts_ppRef, _ACC
         if(!((Btrk1Eta[i_B] < 2.4 && Btrk2Eta[i_B] < 2.4)
         && (Btrk1Pt[i_B]  > 0.5 && Btrk2Pt[i_B]  > 0.5)
         && (Bmu1isAcc[i_B] == 1 && Bmu2isAcc[i_B] == 1))) continue;

         //SELcuts_ppRef, _SEL
         if(!((Btrk1PtErr[i_B] / Btrk1Pt[i_B]) < 0.1 && (Btrk2PtErr[i_B] / Btrk2Pt[i_B]) < 0.1
         && (Btrk1nPixelLayer[i_B] + Btrk1nStripLayer[i_B]) > 10
         && (Btrk2nPixelLayer[i_B] + Btrk2nStripLayer[i_B]) > 10
         && (Btrk1Chi2ndf[i_B]/(Btrk1nPixelLayer[i_B] + Btrk1nStripLayer[i_B])) < 0.18
         && (Btrk2Chi2ndf[i_B]/(Btrk2nPixelLayer[i_B] + Btrk2nStripLayer[i_B])) < 0.18
         && Btrk1highPurity[i_B] && Btrk2highPurity[i_B]
         && (Bmu1SoftMuID[i_B] == 1 && Bmu2SoftMuID[i_B] == 1)
         && abs(Bujmass[i_B] - 3.096916) < 0.15 
         && BujvProb[i_B] > 0.01)) continue;

         //TRGmatching, _TRG
         if(!(Bmu1isTriggered[i_B] == 1 && Bmu2isTriggered[i_B] == 1)) continue;

         //////////////////// set all the cuts (above)//////////////////

         //////////////////// fill all the variables //////////////////
         B_type = Btype[i_B];
         B_mass = Bmass[i_B];
         B_pt = Bpt[i_B];
         B_y = By[i_B];
         B_vtxX = BvtxX[i_B];
         B_vtxY = BvtxY[i_B];
         B_d0 = Bd0[i_B];
         B_d0Err = Bd0Err[i_B];
         B_dxyz = Bdxyz[i_B];
         B_dxyzErr = BdxyzErr[i_B];
         B_chi2ndf = Bchi2ndf[i_B];
         B_chi2cl = Bchi2cl[i_B];
         B_dtheta = Bdtheta[i_B];
         B_cos_dtheta = Bcos_dtheta[i_B];
         B_alpha = Balpha[i_B];
         B_svpvDistance = BsvpvDistance[i_B];
         B_svpvDisErr = BsvpvDisErr[i_B];
         B_svpvDistance_2D = BsvpvDistance_2D[i_B];
         B_svpvDisErr_2D = BsvpvDisErr_2D[i_B];
         B_norm_svpvDistance = Bnorm_svpvDistance[i_B];
         B_norm_svpvDistance_2D = Bnorm_svpvDistance_2D[i_B];
         B_Qvalue = BQvalue[i_B];
         B_Qvalueuj = BQvalueuj[i_B];
         B_Qvaluemumu = BQvaluemumu[i_B];
         B_trk1Pt = Btrk1Pt[i_B];
         B_trk2Pt = Btrk2Pt[i_B];
         B_trkPtimb = BtrkPtimb[i_B];
         B_trk1Eta = Btrk1Eta[i_B];
         B_trk2Eta = Btrk2Eta[i_B];
         B_trk1PtErr = Btrk1PtErr[i_B];
         B_trk2PtErr = Btrk2PtErr[i_B];
         B_trk1Dz = Btrk1Dz[i_B];
         B_trk2Dz = Btrk2Dz[i_B];
         B_trk1DzError = Btrk1DzError[i_B];
         B_trk2DzError = Btrk2DzError[i_B];
         B_trk1Dxy = Btrk1Dxy[i_B];
         B_trk2Dxy = Btrk2Dxy[i_B];
         B_norm_trk1Dxy = Bnorm_trk1Dxy[i_B];
         B_norm_trk2Dxy = Bnorm_trk2Dxy[i_B];
         B_trk1DxyError = Btrk1DxyError[i_B];
         B_trk2DxyError = Btrk2DxyError[i_B];
         B_trk1nPixelLayer = Btrk1nPixelLayer[i_B];
         B_trk2nPixelLayer = Btrk2nPixelLayer[i_B];
         B_trk1nStripLayer = Btrk1nStripLayer[i_B];
         B_trk2nStripLayer = Btrk2nStripLayer[i_B];
         B_trk1Chi2ndf = Btrk1Chi2ndf[i_B];
         B_trk2Chi2ndf = Btrk2Chi2ndf[i_B];
         B_trk1highPurity = Btrk1highPurity[i_B];
         B_trk2highPurity = Btrk2highPurity[i_B];
         B_trk1dR = Btrk1dR[i_B];
         B_trk2dR = Btrk2dR[i_B];
         B_tktkmass = Btktkmass[i_B];
         B_tktkmassKK = BtktkmassKK[i_B];
         B_tktkvProb = BtktkvProb[i_B];
         B_tktkpt = Btktkpt[i_B];
         B_doubletmass = Bdoubletmass[i_B];
         B_mu1pt = Bmu1pt[i_B];
         B_mu2pt = Bmu2pt[i_B];
         B_mu1highPurity = Bmu1highPurity[i_B];
         B_mu2highPurity = Bmu2highPurity[i_B];
         B_mu1isAcc = Bmu1isAcc[i_B];
         B_mu2isAcc = Bmu2isAcc[i_B];
         B_mu1SoftMuID = Bmu1SoftMuID[i_B];
         B_mu2SoftMuID = Bmu2SoftMuID[i_B];
         B_mu1HybridSoftMuID = Bmu1HybridSoftMuID[i_B];
         B_mu2HybridSoftMuID = Bmu2HybridSoftMuID[i_B];
         B_mu1isTriggered = Bmu1isTriggered[i_B];
         B_mu2isTriggered = Bmu2isTriggered[i_B];
         B_mumumass = Bmumumass[i_B];
         B_ujmass = Bujmass[i_B];
         B_ujvProb = BujvProb[i_B];
         B_gen = Bgen[i_B];
         B_genIndex = BgenIndex[i_B];

         tree->Fill();
      }
   }
   tree->Write();
   outputFile->Close();
   std::cout << "OK" << std::endl;
}
