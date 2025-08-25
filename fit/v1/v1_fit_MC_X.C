#include "TFile.h"
#include "TTree.h"
#include "TCanvas.h"
#include "TLegend.h"
#include "TH2F.h"
#include "RooAbsPdf.h"
#include "RooExtendPdf.h"
#include "RooGaussian.h"
#include "RooRealVar.h"
#include "RooPlot.h"
#include "RooDataSet.h"
#include "RooAddPdf.h"
#include "RooArgList.h"
#include "RooFitResult.h"
#include "RooArgSet.h"
#include "RooFit.h"   

using namespace RooFit;

void fit_MC_X() {
    TFile *f = TFile::Open("/user/u/u25lekai/work/ppRef/analysis_X/selection/root_files/MC_X3872.root");
    TTree *tree = (TTree*)f->Get("tree");
    
    float lowBmass=3.66,highBmass=3.72;
    Int_t N=10000;//N is the max of events
	Int_t BinNum=60;
    
    //define variables
    RooRealVar B_mass("B_mass", "B_mass", lowBmass, highBmass);
    
    //define the RooArgSet of the Dataset
	RooArgSet vars;
		vars.add(B_mass);
    RooDataSet data("data", "data", vars,Import(*tree));
	//RooDataSet* newdata = (RooDataSet*)data.reduce("(B_mass > 5.2) && (B_mass < 5.4)");
	
	//signal of B+
	RooRealVar mean("mean","mean",3.69,3.68,3.70);
	RooRealVar sigma1("sigma1","sigma1",0.005,1e-05,0.01);
	RooGaussian Gaus1("Gaus1","Gaus1",B_mass,mean,sigma1);
	RooRealVar sigma2("sigma2","sigma2",0.0003,1e-05,0.006);
	RooGaussian Gaus2("Gaus2","Gaus2",B_mass,mean,sigma2);
	RooRealVar frac("frac","frac",0.5,0.2,0.8);
	RooAddPdf Gaus_B("Gaus_B","Gaus_B",RooArgList(Gaus1,Gaus2),frac);

	
	//normalization to get n_events
	RooRealVar n_Gaus1=RooRealVar("n_Gaus1","n_Gaus1",300,0,N);
	RooRealVar n_Gaus2=RooRealVar("n_Gaus2","n_Gaus2",200,0,N);

	RooExtendPdf Gaus1e("Gaus1e","Gaus1e",Gaus1,n_Gaus1);
	RooExtendPdf Gaus2e("Gaus2e","Gaus2e",Gaus2,n_Gaus2);

	//add them together
	RooAddPdf tot("total","total",RooArgList(Gaus1e,Gaus2e),RooArgList(n_Gaus1,n_Gaus2));

	//fit
	RooFitResult* result=tot.fitTo(data,Save(),Minimizer("Minuit2", "Migrad"));
	//RooFitResult* result=Gaus_ext.fitTo(data,Save(),Minimizer("Minuit2", "Migrad"));
	//RooFitResult* result=Gaus_ext.fitTo(*newdata,Save(),Minimizer("Minuit2", "Migrad"));

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
	RooPlot *frame_Bmass=B_mass.frame(RooFit::Title("Fit B+"),Bins(BinNum));
	data.plotOn(frame_Bmass,DataError(RooAbsData::SumW2));
	//newdata->plotOn(frame_Bmass,DataError(RooAbsData::SumW2));
	tot.plotOn(frame_Bmass,Name("all"));
	tot.plotOn(frame_Bmass,Components(Gaus1),LineColor(1),LineStyle(1),Name("Gaus1"));
	tot.plotOn(frame_Bmass,Components(Gaus2),LineColor(3),LineStyle(1),Name("Gaus2"));
	
	TLegend leg1(0.7,0.7,0.9,0.9);
	leg1.AddEntry(frame_Bmass->findObject("all"),"Gaus","L");
	leg1.AddEntry(frame_Bmass->findObject("Gaus1"),"Gaus1","L");
	leg1.AddEntry(frame_Bmass->findObject("Gaus2"),"Gaus2","L");
	frame_Bmass->Draw();
	leg1.DrawClone();
	//c_Bmass->SaveAs("./output_fit/1DfitBMass_MC_v8.pdf");
	c_Bmass->SaveAs("./pdf_MC_PSI/1Dfitpsi2Smass_MC_v1.pdf");

    f->Close();
}