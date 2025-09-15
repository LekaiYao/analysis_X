//////////////////////////////////////////////////////////
// This class has been automatically generated on
// Thu Jul 31 16:37:39 2025 by ROOT version 6.30/04
// from TTree ntmix/
// found on file: /lstore/cms/hlegoinha/X3872/MC_DATA/DATA_ppRef_X3872.root
//////////////////////////////////////////////////////////

#ifndef analysis_DATA_h
#define analysis_DATA_h

#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>

// Header file for the classes stored in the TTree if any.

class analysis_DATA {
public :
   TTree          *fChain;   //!pointer to the analyzed TTree or TChain
   Int_t           fCurrent; //!current Tree number in a TChain

// Fixed size dimensions of array or collections stored in the TTree if any.

   // Declaration of leaf types
   Int_t           RunNo;
   Int_t           CentBin;
   Int_t           EvtNo;
   Int_t           LumiNo;
   Float_t         PVx;
   Float_t         PVy;
   Float_t         PVz;
   Float_t         PVxE;
   Float_t         PVyE;
   Float_t         PVzE;
   Float_t         PVnchi2;
   Int_t           nChargedTracks;
   Int_t           nSelectedChargedTracks;
   Int_t           Bsize;
   Int_t           Bindex[252];   //[Bsize]
   Int_t           Btype[252];   //[Bsize]
   Float_t         Bmass[252];   //[Bsize]
   Float_t         Bpt[252];   //[Bsize]
   Float_t         By[252];   //[Bsize]
   Float_t         BvtxX[252];   //[Bsize]
   Float_t         BvtxY[252];   //[Bsize]
   Float_t         Bd0[252];   //[Bsize]
   Float_t         Bd0Err[252];   //[Bsize]
   Float_t         Bdxyz[252];   //[Bsize]
   Float_t         BdxyzErr[252];   //[Bsize]
   Float_t         Bchi2ndf[252];   //[Bsize]
   Float_t         Bchi2cl[252];   //[Bsize]
   Float_t         Bdtheta[252];   //[Bsize]
   Float_t         Bcos_dtheta[252];   //[Bsize]
   Float_t         Balpha[252];   //[Bsize]
   Float_t         BsvpvDistance[252];   //[Bsize]
   Float_t         BsvpvDisErr[252];   //[Bsize]
   Float_t         BsvpvDistance_2D[252];   //[Bsize]
   Float_t         BsvpvDisErr_2D[252];   //[Bsize]
   Float_t         Bnorm_svpvDistance[252];   //[Bsize]
   Float_t         Bnorm_svpvDistance_2D[252];   //[Bsize]
   Float_t         BQvalue[252];   //[Bsize]
   Float_t         BQvalueuj[252];   //[Bsize]
   Float_t         BQvaluemumu[252];   //[Bsize]
   Float_t         Btrk1Pt[252];   //[Bsize]
   Float_t         Btrk2Pt[252];   //[Bsize]
   Float_t         BtrkPtimb[252];   //[Bsize]
   Float_t         Btrk1Eta[252];   //[Bsize]
   Float_t         Btrk2Eta[252];   //[Bsize]
   Float_t         Btrk1PtErr[252];   //[Bsize]
   Float_t         Btrk2PtErr[252];   //[Bsize]
   Float_t         Btrk1Dz[252];   //[Bsize]
   Float_t         Btrk2Dz[252];   //[Bsize]
   Float_t         Btrk1DzError[252];   //[Bsize]
   Float_t         Btrk2DzError[252];   //[Bsize]
   Float_t         Btrk1Dxy[252];   //[Bsize]
   Float_t         Btrk2Dxy[252];   //[Bsize]
   Float_t         Bnorm_trk1Dxy[252];   //[Bsize]
   Float_t         Bnorm_trk2Dxy[252];   //[Bsize]
   Float_t         Btrk1DxyError[252];   //[Bsize]
   Float_t         Btrk2DxyError[252];   //[Bsize]
   Float_t         Btrk1nPixelLayer[252];   //[Bsize]
   Float_t         Btrk2nPixelLayer[252];   //[Bsize]
   Float_t         Btrk1nStripLayer[252];   //[Bsize]
   Float_t         Btrk2nStripLayer[252];   //[Bsize]
   Float_t         Btrk1Chi2ndf[252];   //[Bsize]
   Float_t         Btrk2Chi2ndf[252];   //[Bsize]
   Bool_t          Btrk1highPurity[252];   //[Bsize]
   Bool_t          Btrk2highPurity[252];   //[Bsize]
   Float_t         Btrk1dR[252];   //[Bsize]
   Float_t         Btrk2dR[252];   //[Bsize]
   Float_t         Btktkmass[252];   //[Bsize]
   Float_t         BtktkmassKK[252];   //[Bsize]
   Float_t         BtktkvProb[252];   //[Bsize]
   Float_t         Btktkpt[252];   //[Bsize]
   Float_t         Bdoubletmass[252];   //[Bsize]
   Float_t         Bmu1pt[252];   //[Bsize]
   Float_t         Bmu2pt[252];   //[Bsize]
   Bool_t          Bmu1highPurity[252];   //[Bsize]
   Bool_t          Bmu2highPurity[252];   //[Bsize]
   Bool_t          Bmu1isAcc[252];   //[Bsize]
   Bool_t          Bmu2isAcc[252];   //[Bsize]
   Bool_t          Bmu1SoftMuID[252];   //[Bsize]
   Bool_t          Bmu2SoftMuID[252];   //[Bsize]
   Bool_t          Bmu1HybridSoftMuID[252];   //[Bsize]
   Bool_t          Bmu2HybridSoftMuID[252];   //[Bsize]
   Bool_t          Bmu1isTriggered[252];   //[Bsize]
   Bool_t          Bmu2isTriggered[252];   //[Bsize]
   Float_t         Bmumumass[252];   //[Bsize]
   Float_t         Bujmass[252];   //[Bsize]
   Float_t         BujvProb[252];   //[Bsize]
   Float_t         Bgen[252];   //[Bsize]
   Int_t           BgenIndex[252];   //[Bsize]

   // List of branches
   TBranch        *b_RunNo;   //!
   TBranch        *b_CentBin;   //!
   TBranch        *b_EvtNo;   //!
   TBranch        *b_LumiNo;   //!
   TBranch        *b_PVx;   //!
   TBranch        *b_PVy;   //!
   TBranch        *b_PVz;   //!
   TBranch        *b_PVxE;   //!
   TBranch        *b_PVyE;   //!
   TBranch        *b_PVzE;   //!
   TBranch        *b_PVnchi2;   //!
   TBranch        *b_nChargedTracks;   //!
   TBranch        *b_nSelectedChargedTracks;   //!
   TBranch        *b_Bsize;   //!
   TBranch        *b_Bindex;   //!
   TBranch        *b_Btype;   //!
   TBranch        *b_Bmass;   //!
   TBranch        *b_Bpt;   //!
   TBranch        *b_By;   //!
   TBranch        *b_BvtxX;   //!
   TBranch        *b_BvtxY;   //!
   TBranch        *b_Bd0;   //!
   TBranch        *b_Bd0Err;   //!
   TBranch        *b_Bdxyz;   //!
   TBranch        *b_BdxyzErr;   //!
   TBranch        *b_Bchi2ndf;   //!
   TBranch        *b_Bchi2cl;   //!
   TBranch        *b_Bdtheta;   //!
   TBranch        *b_Bcos_dtheta;   //!
   TBranch        *b_Balpha;   //!
   TBranch        *b_BsvpvDistance;   //!
   TBranch        *b_BsvpvDisErr;   //!
   TBranch        *b_BsvpvDistance_2D;   //!
   TBranch        *b_BsvpvDisErr_2D;   //!
   TBranch        *b_Bnorm_svpvDistance;   //!
   TBranch        *b_Bnorm_svpvDistance_2D;   //!
   TBranch        *b_BQvalue;   //!
   TBranch        *b_BQvalueuj;   //!
   TBranch        *b_BQvaluemumu;   //!
   TBranch        *b_Btrk1Pt;   //!
   TBranch        *b_Btrk2Pt;   //!
   TBranch        *b_BtrkPtimb;   //!
   TBranch        *b_Btrk1Eta;   //!
   TBranch        *b_Btrk2Eta;   //!
   TBranch        *b_Btrk1PtErr;   //!
   TBranch        *b_Btrk2PtErr;   //!
   TBranch        *b_Btrk1Dz;   //!
   TBranch        *b_Btrk2Dz;   //!
   TBranch        *b_Btrk1DzError;   //!
   TBranch        *b_Btrk2DzError;   //!
   TBranch        *b_Btrk1Dxy;   //!
   TBranch        *b_Btrk2Dxy;   //!
   TBranch        *b_Bnorm_trk1Dxy;   //!
   TBranch        *b_Bnorm_trk2Dxy;   //!
   TBranch        *b_Btrk1DxyError;   //!
   TBranch        *b_Btrk2DxyError;   //!
   TBranch        *b_Btrk1nPixelLayer;   //!
   TBranch        *b_Btrk2nPixelLayer;   //!
   TBranch        *b_Btrk1nStripLayer;   //!
   TBranch        *b_Btrk2nStripLayer;   //!
   TBranch        *b_Btrk1Chi2ndf;   //!
   TBranch        *b_Btrk2Chi2ndf;   //!
   TBranch        *b_Btrk1highPurity;   //!
   TBranch        *b_Btrk2highPurity;   //!
   TBranch        *b_Btrk1dR;   //!
   TBranch        *b_Btrk2dR;   //!
   TBranch        *b_Btktkmass;   //!
   TBranch        *b_BtktkmassKK;   //!
   TBranch        *b_BtktkvProb;   //!
   TBranch        *b_Btktkpt;   //!
   TBranch        *b_Bdoubletmass;   //!
   TBranch        *b_Bmu1pt;   //!
   TBranch        *b_Bmu2pt;   //!
   TBranch        *b_Bmu1highPurity;   //!
   TBranch        *b_Bmu2highPurity;   //!
   TBranch        *b_Bmu1isAcc;   //!
   TBranch        *b_Bmu2isAcc;   //!
   TBranch        *b_Bmu1SoftMuID;   //!
   TBranch        *b_Bmu2SoftMuID;   //!
   TBranch        *b_Bmu1HybridSoftMuID;   //!
   TBranch        *b_Bmu2HybridSoftMuID;   //!
   TBranch        *b_Bmu1isTriggered;   //!
   TBranch        *b_Bmu2isTriggered;   //!
   TBranch        *b_Bmumumass;   //!
   TBranch        *b_Bujmass;   //!
   TBranch        *b_BujvProb;   //!
   TBranch        *b_Bgen;   //!
   TBranch        *b_BgenIndex;   //!

   analysis_DATA(TTree *tree=0);
   virtual ~analysis_DATA();
   virtual Int_t    Cut(Long64_t entry);
   virtual Int_t    GetEntry(Long64_t entry);
   virtual Long64_t LoadTree(Long64_t entry);
   virtual void     Init(TTree *tree);
   virtual void     Loop();
   virtual Bool_t   Notify();
   virtual void     Show(Long64_t entry = -1);
};

#endif

#ifdef analysis_DATA_cxx
analysis_DATA::analysis_DATA(TTree *tree) : fChain(0) 
{
// if parameter tree is not specified (or zero), connect the file
// used to generate this class and read the Tree.
   if (tree == 0) {
      TFile *f = (TFile*)gROOT->GetListOfFiles()->FindObject("/lstore/cms/hlegoinha/X3872/MC_DATA/DATA_ppRef_X3872.root");
      if (!f || !f->IsOpen()) {
         f = new TFile("/lstore/cms/hlegoinha/X3872/MC_DATA/DATA_ppRef_X3872.root");
      }
      TDirectory * dir = (TDirectory*)f->Get("/lstore/cms/hlegoinha/X3872/MC_DATA/DATA_ppRef_X3872.root:/Bfinder");
      dir->GetObject("ntmix",tree);

   }
   Init(tree);
}

analysis_DATA::~analysis_DATA()
{
   if (!fChain) return;
   delete fChain->GetCurrentFile();
}

Int_t analysis_DATA::GetEntry(Long64_t entry)
{
// Read contents of entry.
   if (!fChain) return 0;
   return fChain->GetEntry(entry);
}
Long64_t analysis_DATA::LoadTree(Long64_t entry)
{
// Set the environment to read one entry
   if (!fChain) return -5;
   Long64_t centry = fChain->LoadTree(entry);
   if (centry < 0) return centry;
   if (fChain->GetTreeNumber() != fCurrent) {
      fCurrent = fChain->GetTreeNumber();
      Notify();
   }
   return centry;
}

void analysis_DATA::Init(TTree *tree)
{
   // The Init() function is called when the selector needs to initialize
   // a new tree or chain. Typically here the branch addresses and branch
   // pointers of the tree will be set.
   // It is normally not necessary to make changes to the generated
   // code, but the routine can be extended by the user if needed.
   // Init() will be called many times when running on PROOF
   // (once per file to be processed).

   // Set branch addresses and branch pointers
   if (!tree) return;
   fChain = tree;
   fCurrent = -1;
   fChain->SetMakeClass(1);

   fChain->SetBranchAddress("RunNo", &RunNo, &b_RunNo);
   fChain->SetBranchAddress("CentBin", &CentBin, &b_CentBin);
   fChain->SetBranchAddress("EvtNo", &EvtNo, &b_EvtNo);
   fChain->SetBranchAddress("LumiNo", &LumiNo, &b_LumiNo);
   fChain->SetBranchAddress("PVx", &PVx, &b_PVx);
   fChain->SetBranchAddress("PVy", &PVy, &b_PVy);
   fChain->SetBranchAddress("PVz", &PVz, &b_PVz);
   fChain->SetBranchAddress("PVxE", &PVxE, &b_PVxE);
   fChain->SetBranchAddress("PVyE", &PVyE, &b_PVyE);
   fChain->SetBranchAddress("PVzE", &PVzE, &b_PVzE);
   fChain->SetBranchAddress("PVnchi2", &PVnchi2, &b_PVnchi2);
   fChain->SetBranchAddress("nChargedTracks", &nChargedTracks, &b_nChargedTracks);
   fChain->SetBranchAddress("nSelectedChargedTracks", &nSelectedChargedTracks, &b_nSelectedChargedTracks);
   fChain->SetBranchAddress("Bsize", &Bsize, &b_Bsize);
   fChain->SetBranchAddress("Bindex", Bindex, &b_Bindex);
   fChain->SetBranchAddress("Btype", Btype, &b_Btype);
   fChain->SetBranchAddress("Bmass", Bmass, &b_Bmass);
   fChain->SetBranchAddress("Bpt", Bpt, &b_Bpt);
   fChain->SetBranchAddress("By", By, &b_By);
   fChain->SetBranchAddress("BvtxX", BvtxX, &b_BvtxX);
   fChain->SetBranchAddress("BvtxY", BvtxY, &b_BvtxY);
   fChain->SetBranchAddress("Bd0", Bd0, &b_Bd0);
   fChain->SetBranchAddress("Bd0Err", Bd0Err, &b_Bd0Err);
   fChain->SetBranchAddress("Bdxyz", Bdxyz, &b_Bdxyz);
   fChain->SetBranchAddress("BdxyzErr", BdxyzErr, &b_BdxyzErr);
   fChain->SetBranchAddress("Bchi2ndf", Bchi2ndf, &b_Bchi2ndf);
   fChain->SetBranchAddress("Bchi2cl", Bchi2cl, &b_Bchi2cl);
   fChain->SetBranchAddress("Bdtheta", Bdtheta, &b_Bdtheta);
   fChain->SetBranchAddress("Bcos_dtheta", Bcos_dtheta, &b_Bcos_dtheta);
   fChain->SetBranchAddress("Balpha", Balpha, &b_Balpha);
   fChain->SetBranchAddress("BsvpvDistance", BsvpvDistance, &b_BsvpvDistance);
   fChain->SetBranchAddress("BsvpvDisErr", BsvpvDisErr, &b_BsvpvDisErr);
   fChain->SetBranchAddress("BsvpvDistance_2D", BsvpvDistance_2D, &b_BsvpvDistance_2D);
   fChain->SetBranchAddress("BsvpvDisErr_2D", BsvpvDisErr_2D, &b_BsvpvDisErr_2D);
   fChain->SetBranchAddress("Bnorm_svpvDistance", Bnorm_svpvDistance, &b_Bnorm_svpvDistance);
   fChain->SetBranchAddress("Bnorm_svpvDistance_2D", Bnorm_svpvDistance_2D, &b_Bnorm_svpvDistance_2D);
   fChain->SetBranchAddress("BQvalue", BQvalue, &b_BQvalue);
   fChain->SetBranchAddress("BQvalueuj", BQvalueuj, &b_BQvalueuj);
   fChain->SetBranchAddress("BQvaluemumu", BQvaluemumu, &b_BQvaluemumu);
   fChain->SetBranchAddress("Btrk1Pt", Btrk1Pt, &b_Btrk1Pt);
   fChain->SetBranchAddress("Btrk2Pt", Btrk2Pt, &b_Btrk2Pt);
   fChain->SetBranchAddress("BtrkPtimb", BtrkPtimb, &b_BtrkPtimb);
   fChain->SetBranchAddress("Btrk1Eta", Btrk1Eta, &b_Btrk1Eta);
   fChain->SetBranchAddress("Btrk2Eta", Btrk2Eta, &b_Btrk2Eta);
   fChain->SetBranchAddress("Btrk1PtErr", Btrk1PtErr, &b_Btrk1PtErr);
   fChain->SetBranchAddress("Btrk2PtErr", Btrk2PtErr, &b_Btrk2PtErr);
   fChain->SetBranchAddress("Btrk1Dz", Btrk1Dz, &b_Btrk1Dz);
   fChain->SetBranchAddress("Btrk2Dz", Btrk2Dz, &b_Btrk2Dz);
   fChain->SetBranchAddress("Btrk1DzError", Btrk1DzError, &b_Btrk1DzError);
   fChain->SetBranchAddress("Btrk2DzError", Btrk2DzError, &b_Btrk2DzError);
   fChain->SetBranchAddress("Btrk1Dxy", Btrk1Dxy, &b_Btrk1Dxy);
   fChain->SetBranchAddress("Btrk2Dxy", Btrk2Dxy, &b_Btrk2Dxy);
   fChain->SetBranchAddress("Bnorm_trk1Dxy", Bnorm_trk1Dxy, &b_Bnorm_trk1Dxy);
   fChain->SetBranchAddress("Bnorm_trk2Dxy", Bnorm_trk2Dxy, &b_Bnorm_trk2Dxy);
   fChain->SetBranchAddress("Btrk1DxyError", Btrk1DxyError, &b_Btrk1DxyError);
   fChain->SetBranchAddress("Btrk2DxyError", Btrk2DxyError, &b_Btrk2DxyError);
   fChain->SetBranchAddress("Btrk1nPixelLayer", Btrk1nPixelLayer, &b_Btrk1nPixelLayer);
   fChain->SetBranchAddress("Btrk2nPixelLayer", Btrk2nPixelLayer, &b_Btrk2nPixelLayer);
   fChain->SetBranchAddress("Btrk1nStripLayer", Btrk1nStripLayer, &b_Btrk1nStripLayer);
   fChain->SetBranchAddress("Btrk2nStripLayer", Btrk2nStripLayer, &b_Btrk2nStripLayer);
   fChain->SetBranchAddress("Btrk1Chi2ndf", Btrk1Chi2ndf, &b_Btrk1Chi2ndf);
   fChain->SetBranchAddress("Btrk2Chi2ndf", Btrk2Chi2ndf, &b_Btrk2Chi2ndf);
   fChain->SetBranchAddress("Btrk1highPurity", Btrk1highPurity, &b_Btrk1highPurity);
   fChain->SetBranchAddress("Btrk2highPurity", Btrk2highPurity, &b_Btrk2highPurity);
   fChain->SetBranchAddress("Btrk1dR", Btrk1dR, &b_Btrk1dR);
   fChain->SetBranchAddress("Btrk2dR", Btrk2dR, &b_Btrk2dR);
   fChain->SetBranchAddress("Btktkmass", Btktkmass, &b_Btktkmass);
   fChain->SetBranchAddress("BtktkmassKK", BtktkmassKK, &b_BtktkmassKK);
   fChain->SetBranchAddress("BtktkvProb", BtktkvProb, &b_BtktkvProb);
   fChain->SetBranchAddress("Btktkpt", Btktkpt, &b_Btktkpt);
   fChain->SetBranchAddress("Bdoubletmass", Bdoubletmass, &b_Bdoubletmass);
   fChain->SetBranchAddress("Bmu1pt", Bmu1pt, &b_Bmu1pt);
   fChain->SetBranchAddress("Bmu2pt", Bmu2pt, &b_Bmu2pt);
   fChain->SetBranchAddress("Bmu1highPurity", Bmu1highPurity, &b_Bmu1highPurity);
   fChain->SetBranchAddress("Bmu2highPurity", Bmu2highPurity, &b_Bmu2highPurity);
   fChain->SetBranchAddress("Bmu1isAcc", Bmu1isAcc, &b_Bmu1isAcc);
   fChain->SetBranchAddress("Bmu2isAcc", Bmu2isAcc, &b_Bmu2isAcc);
   fChain->SetBranchAddress("Bmu1SoftMuID", Bmu1SoftMuID, &b_Bmu1SoftMuID);
   fChain->SetBranchAddress("Bmu2SoftMuID", Bmu2SoftMuID, &b_Bmu2SoftMuID);
   fChain->SetBranchAddress("Bmu1HybridSoftMuID", Bmu1HybridSoftMuID, &b_Bmu1HybridSoftMuID);
   fChain->SetBranchAddress("Bmu2HybridSoftMuID", Bmu2HybridSoftMuID, &b_Bmu2HybridSoftMuID);
   fChain->SetBranchAddress("Bmu1isTriggered", Bmu1isTriggered, &b_Bmu1isTriggered);
   fChain->SetBranchAddress("Bmu2isTriggered", Bmu2isTriggered, &b_Bmu2isTriggered);
   fChain->SetBranchAddress("Bmumumass", Bmumumass, &b_Bmumumass);
   fChain->SetBranchAddress("Bujmass", Bujmass, &b_Bujmass);
   fChain->SetBranchAddress("BujvProb", BujvProb, &b_BujvProb);
   fChain->SetBranchAddress("Bgen", Bgen, &b_Bgen);
   fChain->SetBranchAddress("BgenIndex", BgenIndex, &b_BgenIndex);
   Notify();
}

Bool_t analysis_DATA::Notify()
{
   // The Notify() function is called when a new file is opened. This
   // can be either for a new TTree in a TChain or when when a new TTree
   // is started when using PROOF. It is normally not necessary to make changes
   // to the generated code, but the routine can be extended by the
   // user if needed. The return value is currently not used.

   return kTRUE;
}

void analysis_DATA::Show(Long64_t entry)
{
// Print contents of entry.
// If entry is not specified, print current entry
   if (!fChain) return;
   fChain->Show(entry);
}
Int_t analysis_DATA::Cut(Long64_t entry)
{
// This function may be called from Loop.
// returns  1 if entry is accepted.
// returns -1 otherwise.
   return 1;
}
#endif // #ifdef analysis_DATA_cxx
