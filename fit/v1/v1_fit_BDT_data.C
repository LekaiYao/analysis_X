#include "TFile.h"
#include "TTree.h"
#include "TCanvas.h"
#include "TLegend.h"
#include "TH2F.h"
#include "RooAbsPdf.h"
#include "RooExtendPdf.h"
#include "RooGaussian.h"
#include "RooExponential.h"
#include "RooRealVar.h"
#include "RooPlot.h"
#include "RooDataSet.h"
#include "RooAddPdf.h"
#include "RooArgList.h"
#include "RooFitResult.h"
#include "RooArgSet.h"
#include "RooFit.h"   

using namespace RooFit;

void fit_BDT_data() {
    //TFile *f = TFile::Open("../selection/root_files/DATA_XPSI_nonan.root");
	TFile *f = TFile::Open("../TMVA/BDT_root/test3/DATA_XPSI_BDT_test3.root");
    TTree *tree = (TTree*)f->Get("tree");
    
    float lowBmass=3.5,highBmass=4.1;
    Int_t N=100000;//N is the max of events
	Int_t BinNum=60;
    
    //define variables
    RooRealVar B_mass("B_mass", "B_mass", lowBmass, highBmass);
	RooRealVar BDT_score("BDT_score", "BDT_score", -1.0, 1.0);
    
    //define the RooArgSet of the Dataset
	RooArgSet vars;
		vars.add(B_mass);
		vars.add(BDT_score);
    RooDataSet data("data", "data", vars,Import(*tree));
	RooDataSet* newdata = (RooDataSet*)data.reduce("(BDT_score>0.211)");
	
	//signal of Psi(2S)
	RooRealVar mean_psi2s("mean_psi2s","mean_psi2s",3.69,3.68,3.70);
	RooRealVar sigma1_psi2s("sigma1_psi2s","sigma1_psi2s",0.006,1e-04,0.01);
	RooGaussian Gaus1_psi2s("Gaus1_psi2s","Gaus1_psi2s",B_mass,mean_psi2s,sigma1_psi2s);
	RooRealVar sigma2_psi2s("sigma2_psi2s","sigma2_psi2s",0.003,5e-05,0.005);
	RooGaussian Gaus2_psi2s("Gaus2_psi2s","Gaus2_psi2s",B_mass,mean_psi2s,sigma2_psi2s);
	RooRealVar frac_psi2s("frac_psi2s","frac_psi2s",0.5,0.2,0.8);
	RooAddPdf sig_psi2s("sig_psi2s","sig_psi2s",RooArgList(Gaus1_psi2s,Gaus2_psi2s),frac_psi2s);

	//signal of X(3872)
	RooRealVar mean_X("mean_X","mean_X",3.87,3.86,3.88);
	//RooRealVar sigma1_X("sigma1_X","sigma1_X",0.006,1e-04,0.01);
	RooGaussian Gaus1_X("Gaus1_X","Gaus1_X",B_mass,mean_X,sigma1_psi2s);
	//RooRealVar sigma2_X("sigma2_X","sigma2_X",0.003,5e-05,0.005);
	RooGaussian Gaus2_X("Gaus2_X","Gaus2_X",B_mass,mean_X,sigma2_psi2s);
	//RooRealVar frac_X("frac_X","frac_X",0.5,0.2,0.8);
	RooAddPdf sig_X("sig_X","sig_X",RooArgList(Gaus1_X,Gaus2_X),frac_psi2s);

	// Chebychev polynomial background
	RooRealVar c1("c1", "coefficient #1", 0.1, -1.0, 1.0);
	RooRealVar c2("c2", "coefficient #2", -0.1, -1.0, 1.0);
	// You can add more coefficients (c3, c4, ...) for higher-order polynomials

	RooChebychev bkg("bkg", "Chebychev background", B_mass, RooArgList(c1, c2));

	//normalization to get n_events
	RooRealVar n_sig_psi2s = RooRealVar("n_sig_psi2s","n_sig_psi2s",1000,0,N);
	RooRealVar n_sig_X = RooRealVar("n_sig_X","n_sig_X",100,0,N);
	RooRealVar n_bkg = RooRealVar("n_bkg","n_bkg",1000,0,N);

	RooExtendPdf sige_psi2s("sige_psi2s","sige_psi2s",sig_psi2s,n_sig_psi2s);
	RooExtendPdf sige_X("sige_X","sige_X",sig_X,n_sig_X);
	RooExtendPdf bkge("bkge","bkge",bkg,n_bkg);
	
	//add them together
	RooAddPdf tot("total","total",RooArgList(sige_psi2s,sige_X,bkge),RooArgList(n_sig_psi2s,n_sig_X,n_bkg));

	//fit
	//RooFitResult* result=tot.fitTo(data,Save(),PrintLevel(1),Verbose(true),Minimizer("Minuit", "minimize"));
	//RooFitResult* result=Gaus_ext.fitTo(data,Save(),Minimizer("Minuit2", "Migrad"));
	//RooFitResult* result=Gaus_ext.fitTo(*newdata,Save(),Minimizer("Minuit2", "Migrad"));
	RooFitResult* result=tot.fitTo(*newdata,Save(),Minimizer("Minuit2", "Migrad"));

    //draw
    /*
    TCanvas *c1=new TCanvas("c1","c1",10,10,800,600);
    c1->cd();
	TH1F *hist_Bmass=new TH1F("hist_Bmass","hist_Bmass",BinNum,lowBmass,highBmass);
	hist_Bmass=data->createHistogram(Bmass);
	hist_Bmass->Draw();
	c1->SaveAs("MCplot.pdf");
    */

	//draw the fitX/Y and the leg
	TCanvas *c_Bmass=new TCanvas("c_Bmass","c_Bmass",10,10,800,600);
	RooPlot *frame_Bmass=B_mass.frame(RooFit::Title("Fit #{Psi}(2S) & X(3872)"),Bins(BinNum));
	//data.plotOn(frame_Bmass,DataError(RooAbsData::SumW2));
	newdata->plotOn(frame_Bmass,DataError(RooAbsData::SumW2));
	tot.plotOn(frame_Bmass,Name("all"));
	tot.plotOn(frame_Bmass,Components(sig_psi2s),LineColor(1),LineStyle(1),Name("sig_psi2s"));
	tot.plotOn(frame_Bmass,Components(sig_X),LineColor(1),LineStyle(1),Name("sig_X"));
	tot.plotOn(frame_Bmass,Components(bkg),LineColor(3),LineStyle(1),Name("bkg"));
	
	TLegend leg1(0.7,0.7,0.9,0.9);
	leg1.AddEntry(frame_Bmass->findObject("all"),"data","L");
	leg1.AddEntry(frame_Bmass->findObject("sig_psi2s"),"signal #{Psi}(2S)","L");
	leg1.AddEntry(frame_Bmass->findObject("sig_X"),"signal X(3872)","L");
	leg1.AddEntry(frame_Bmass->findObject("bkg"),"background","L");
	frame_Bmass->Draw();
	leg1.DrawClone();
	//c_Bmass->SaveAs("./output_fit/1DfitBMass_MC_v8.pdf");
	c_Bmass->SaveAs("./pdf_BDT_data/1DfitPsiX_BDTdata_v2.pdf");

    f->Close();
}