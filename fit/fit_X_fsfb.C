#include "TFile.h"
#include "TTree.h"
#include "TCanvas.h"
#include "TLegend.h"
#include "TH2F.h"
#include "RooAbsPdf.h"
#include "RooExtendPdf.h"
#include "RooGaussian.h"
#include "RooExponential.h"
#include "RooChebychev.h"
#include "RooRealVar.h"
#include "RooPlot.h"
#include "RooDataSet.h"
#include "RooAddPdf.h"
#include "RooArgList.h"
#include "RooFitResult.h"
#include "RooArgSet.h"
#include "RooFit.h"    

using namespace RooFit;

void fit_X_fsfb() {
    TFile *f = TFile::Open("/user/l/lekai/work/ppRef/analysis_X/TMVA/DATA_XPSI_BDT_trainX4.root");
    TTree *tree = (TTree*)f->Get("tree");
    
    float lowBmass=3.75,highBmass=4.0,lowBDTscore=-1.0,highBDTscore=1.0;
    Int_t N=1000000;//N is the max of events
	Int_t BinNum=100;
    
    //define variables
    RooRealVar B_mass("B_mass", "B_mass", lowBmass, highBmass);
	RooRealVar BDT_score("BDT_score", "BDT_score", lowBDTscore, highBDTscore);
    
    //define the RooArgSet of the Dataset
	RooArgSet vars;
		vars.add(B_mass);
		vars.add(BDT_score);
    RooDataSet data("data", "data", vars,Import(*tree));
	RooDataSet* newdata = (RooDataSet*)data.reduce("(BDT_score>0.0)");

	//signal of X(3872)
	RooRealVar mean("mean","mean",3.87,3.83,3.91);
	RooRealVar sigma1("sigma1","sigma1",0.008,1e-05,0.02);
	RooGaussian Gaus1("Gaus1","Gaus1",B_mass,mean,sigma1);
	RooRealVar sigma2("sigma2","sigma2",0.004,1e-05,0.01);
	RooGaussian Gaus2("Gaus2","Gaus2",B_mass,mean,sigma2);
	RooRealVar frac("frac","frac",0.5,0.2,0.8);
	RooAddPdf sig("sig","sig",RooArgList(Gaus1,Gaus2),frac);

	//background 
	//RooRealVar lambda("lambda","lambda",1,1e-03,10);
	//RooExponential bkg("bkg","bkg",B_mass,lambda,true);
	// Chebychev polynomial background
	RooRealVar c1("c1", "coefficient #1", 0.1, -1.0, 1.0);
	RooRealVar c2("c2", "coefficient #2", -0.1, -1.0, 1.0);
	// You can add more coefficients (c3, c4, ...) for higher-order polynomials

	RooChebychev bkg("bkg", "Chebychev background", B_mass, RooArgList(c1, c2));


	//normalization to get n_events
	RooRealVar n_sig = RooRealVar("n_sig","n_sig",1000,0,N);
	RooRealVar n_bkg = RooRealVar("n_bkg","n_bkg",1000,0,N);

	RooExtendPdf sige("sige","sige",sig,n_sig);
	RooExtendPdf bkge("bkge","bkge",bkg,n_bkg);
	
	//add them together
	RooAddPdf tot("total","total",RooArgList(sige,bkge),RooArgList(n_sig,n_bkg));

	//fit
	//RooFitResult* result=tot.fitTo(data,Save(),Minimizer("Minuit2", "Migrad"));
	//RooFitResult* result=Gaus_ext.fitTo(data,Save(),Minimizer("Minuit2", "Migrad"));
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
	RooPlot *frame_Bmass=B_mass.frame(RooFit::Title("Fit PSI2S"),Bins(BinNum));
	//data.plotOn(frame_Bmass,DataError(RooAbsData::SumW2));
	newdata->plotOn(frame_Bmass,DataError(RooAbsData::SumW2));
	tot.plotOn(frame_Bmass,Name("all"));
	tot.plotOn(frame_Bmass,Components(sig),LineColor(1),LineStyle(1),Name("sig"));
	tot.plotOn(frame_Bmass,Components(bkg),LineColor(3),LineStyle(1),Name("bkg"));
	
	TLegend leg1(0.7,0.7,0.9,0.9);
	leg1.AddEntry(frame_Bmass->findObject("all"),"data","L");
	leg1.AddEntry(frame_Bmass->findObject("sig"),"signal","L");
	leg1.AddEntry(frame_Bmass->findObject("bkg"),"background","L");
	frame_Bmass->Draw();
	leg1.DrawClone();

	c_Bmass->SaveAs("./pdf_X_fsfb/Xmass_fsfb_trainX4_v1.pdf");

    f->Close();
}