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

void fit_data() {
    TFile *f = TFile::Open("../selection/root_files/DATA_XPSI_cut0.root");
    TTree *tree = (TTree*)f->Get("tree");
    
    float lowBmass=3.6,highBmass=4.0;
	//float lowDis=0,highDis=4;//1.92
	//float lowDis_2D=0,highDis_2D=0.92;//1.92
	//float lowBalpha=0.0,highBalpha=3.2;
	//float lowBchi2cl=0.0,highBchi2cl=1.0;
	float lowBQvalueuj=-0.6,highBQvalueuj=0.095;
	float lowBtrk1dR=0.0,highBtrk1dR=0.639;
	float lowBtrk2dR=0.0,highBtrk2dR=0.64;
    Int_t N=2000000;//N is the max of events
	Int_t BinNum=80;
    
    //define variables
    RooRealVar B_mass("B_mass", "B_mass", lowBmass, highBmass);
	//RooRealVar B_norm_svpvDistance("B_norm_svpvDistance","B_norm_svpvDistance",lowDis,highDis);
	//RooRealVar B_norm_svpvDistance_2D("B_norm_svpvDistance_2D","B_norm_svpvDistance_2D",lowDis_2D,highDis_2D);
	//RooRealVar B_alpha("B_alpha", "B_alpha", lowBalpha, highBalpha);
	//RooRealVar B_chi2cl("B_chi2cl","B_chi2cl",lowBchi2cl,highBchi2cl);
	RooRealVar B_Qvalueuj("B_Qvalueuj","B_Qvalueuj",lowBQvalueuj,highBQvalueuj);
	RooRealVar B_trk1dR("B_trk1dR","B_trk1dR",lowBtrk1dR,highBtrk1dR);
	RooRealVar B_trk2dR("B_trk2dR","B_trk2dR",lowBtrk2dR,highBtrk2dR);
    
    //define the RooArgSet of the Dataset
	RooArgSet vars;
		vars.add(B_mass);
		//vars.add(B_norm_svpvDistance);
		//vars.add(B_norm_svpvDistance_2D);
		//vars.add(B_alpha);
		//vars.add(B_chi2cl);
		vars.add(B_Qvalueuj);
		vars.add(B_trk1dR);
		vars.add(B_trk2dR);
    RooDataSet data("data", "data", vars,Import(*tree));
	//RooDataSet* newdata = (RooDataSet*)data.reduce("(B_norm_svpvDistance<1.92)");
	
	/*
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
	RooRealVar sigma1_X("sigma1_X","sigma1_X",0.006,1e-04,0.01);
	RooGaussian Gaus1_X("Gaus1_X","Gaus1_X",B_mass,mean_X,sigma1_X);
	RooRealVar sigma2_X("sigma2_X","sigma2_X",0.003,5e-05,0.005);
	RooGaussian Gaus2_X("Gaus2_X","Gaus2_X",B_mass,mean_X,sigma2_X);
	RooRealVar frac_X("frac_X","frac_X",0.5,0.2,0.8);
	RooAddPdf sig_X("sig_X","sig_X",RooArgList(Gaus1_X,Gaus2_X),frac_X);
	*/
	//signal of Psi(2S)
	RooRealVar mean_psi2s("mean_psi2s","mean_psi2s",3.69,3.68,3.70);
	RooRealVar sigma1("sigma1","sigma1",0.009,1e-03,0.015);
	RooGaussian Gaus1_psi2s("Gaus1_psi2s","Gaus1_psi2s",B_mass,mean_psi2s,sigma1);
	RooRealVar sigma2("sigma2","sigma2",0.004,1e-04,0.006);
	RooGaussian Gaus2_psi2s("Gaus2_psi2s","Gaus2_psi2s",B_mass,mean_psi2s,sigma2);
	RooRealVar frac("frac","frac",0.5,0.2,0.8);
	RooAddPdf sig_psi2s("sig_psi2s","sig_psi2s",RooArgList(Gaus1_psi2s,Gaus2_psi2s),frac);

	//signal of X(3872)
	RooRealVar mean_X("mean_X","mean_X",3.87,3.86,3.88);
	RooGaussian Gaus1_X("Gaus1_X","Gaus1_X",B_mass,mean_X,sigma1);
	RooGaussian Gaus2_X("Gaus2_X","Gaus2_X",B_mass,mean_X,sigma2);
	RooAddPdf sig_X("sig_X","sig_X",RooArgList(Gaus1_X,Gaus2_X),frac);


	// Chebychev polynomial background
	RooRealVar c1("c1", "coefficient #1", 0.07, -0.5, 0.5);
	RooRealVar c2("c2", "coefficient #2", -0.1, -0.3, 0.3);
	// You can add more coefficients (c3, c4, ...) for higher-order polynomials

	RooChebychev bkg("bkg", "Chebychev background", B_mass, RooArgList(c1, c2));

	//normalization to get n_events
	RooRealVar n_sig_psi2s = RooRealVar("n_sig_psi2s","n_sig_psi2s",30000,0,N);
	RooRealVar n_sig_X = RooRealVar("n_sig_X","n_sig_X",3000,0,N);
	RooRealVar n_bkg = RooRealVar("n_bkg","n_bkg",40000,0,N);

	RooExtendPdf sige_psi2s("sige_psi2s","sige_psi2s",sig_psi2s,n_sig_psi2s);
	RooExtendPdf sige_X("sige_X","sige_X",sig_X,n_sig_X);
	RooExtendPdf bkge("bkge","bkge",bkg,n_bkg);
	
	//add them together
	RooAddPdf tot("total","total",RooArgList(sige_psi2s,sige_X,bkge),RooArgList(n_sig_psi2s,n_sig_X,n_bkg));

	//fit
	//RooFitResult* result=tot.fitTo(data,Save(),Minimizer("Minuit2", "Migrad"),Eps(0.1));
	//RooFitResult* result=tot.fitTo(*newdata,Save(),Minimizer("Minuit2", "Migrad"));
	
	// run Minuit2 with a looser EDM tolerance, e.g. 0.1
	RooAbsReal* nll = tot.createNLL(data, Extended(true));
	RooMinimizer m(*nll);

	m.setEps(1000.0); 
	m.setPrintLevel(1); 
	gErrorIgnoreLevel = kInfo;

	m.minimize("Minuit2","Migrad");
	m.hesse();
	RooFitResult* result = m.save();
	result->Print("v"); 

	// ============================================================
	// Drawing section with pull distribution and chi2/ndf
	// (uses: B_mass, BinNum, lowBmass, highBmass, newdata, tot,
	//        sig_psi2s, sig_X, bkg, and RooFitResult* result)
	// ============================================================

	TCanvas *c_Bmass = new TCanvas("c_Bmass","c_Bmass",10,10,800,700);

	// 1) Main frame WITHOUT title (avoid default "A RooPlot of ...")
	RooPlot *frame_Bmass = B_mass.frame(Bins(BinNum));
	frame_Bmass->SetTitle(""); // remove default title

	//define colors
	int brightAzure = TColor::GetColor(51, 153, 255);
	int green  = TColor::GetColor(0,   128,   0); 
	int pink   = TColor::GetColor(255, 105, 180); 
	int purple = TColor::GetColor(128,   0, 128); 
	int yellow = TColor::GetColor(255, 255,   0); 
	int orange = TColor::GetColor(255, 165,   0); 


	// Plot data and model (name them so pull can find them)
	data.plotOn(frame_Bmass, DataError(RooAbsData::SumW2), Name("data"));
	tot.plotOn(frame_Bmass, LineColor(kBlue), LineStyle(1), Name("all"));
	tot.plotOn(frame_Bmass, Components(sig_psi2s), LineColor(orange), LineStyle(3), Name("sig_psi2s"));
	tot.plotOn(frame_Bmass, Components(sig_X),     LineColor(pink),  LineStyle(2), Name("sig_X"));
	tot.plotOn(frame_Bmass, Components(bkg),       LineColor(brightAzure),LineStyle(1), Name("bkg"));

	// Override Y-axis label with fixed bin width (avoid 0.000999999 formatting)
	frame_Bmass->GetYaxis()->SetTitle("Events / (0.005)");

	// 2) Pull frame WITHOUT title
	RooHist* hpull = frame_Bmass->pullHist("data", "all"); // (data - pdf)/sigma per bin
	RooPlot* frame_pull = B_mass.frame();
	frame_pull->SetTitle(""); // remove default title
	frame_pull->addPlotable(hpull, "P");
	frame_pull->SetYTitle("Pull");
	frame_pull->SetMinimum(-5.0);
	frame_pull->SetMaximum(+5.0);

	// Axis formatting for the pull panel
	frame_pull->GetYaxis()->SetTitleSize(0.10);
	frame_pull->GetYaxis()->SetTitleOffset(0.45);
	frame_pull->GetYaxis()->SetLabelSize(0.09);
	frame_pull->GetXaxis()->SetTitleSize(0.12);
	frame_pull->GetXaxis()->SetLabelSize(0.10);

	// 3) Split canvas into two pads: top 80%, bottom 20%
	c_Bmass->Divide(1,2,0,0);
	TPad* p_top = (TPad*)c_Bmass->cd(1);
	p_top->SetPad(0.0, 0.20, 1.0, 1.0);   // top 80%
	p_top->SetTopMargin(0.10);            // move the plot slightly down
	p_top->SetBottomMargin(0.02);
	p_top->SetLeftMargin(0.12);
	p_top->SetRightMargin(0.04);

	// ---------------- Legend and TPaveText CREATED BEFORE DRAW ----------------
	// Compute data-space coordinates for legend/pave (inside plot box)
	double x0 = lowBmass + 0.60*(highBmass - lowBmass);
	double x1 = lowBmass + 0.72*(highBmass - lowBmass);
	double x2 = lowBmass + 1.00*(highBmass - lowBmass);
	double y0 = 0.80*frame_Bmass->GetMaximum();
	double y1 = 0.70*frame_Bmass->GetMaximum();
	double y2 = 1.00*frame_Bmass->GetMaximum();

	// Legend: compact, inside plot box
	TLegend leg1(x0, y0, x1, y2, "", "br");
	leg1.SetTextSize(0.018);
	leg1.SetBorderSize(1);
	leg1.SetFillStyle(0);
	leg1.AddEntry(frame_Bmass->findObject("data"),       "data",               "PE");
	leg1.AddEntry(frame_Bmass->findObject("all"),        "total fit",          "L");
	leg1.AddEntry(frame_Bmass->findObject("sig_psi2s"),  "signal #Psi(2S)",    "L");
	leg1.AddEntry(frame_Bmass->findObject("sig_X"),      "signal X(3872)",     "L");
	leg1.AddEntry(frame_Bmass->findObject("bkg"),        "background",         "L");

	// PaveText with fit parameters (attach to frame BEFORE drawing)
	TPaveText *pave = new TPaveText(x1, y1, x2, y2, "BR");
	pave->SetFillColor(0);
	pave->SetFillStyle(0);
	pave->SetBorderSize(1);
	pave->SetTextSize(0.02);
	pave->SetTextAlign(12); // left alignment

	if (result) {
		RooArgList finalPars = result->floatParsFinal();
		for (int i = 0; i < finalPars.getSize(); i++) {
			RooRealVar* var = (RooRealVar*)&finalPars[i];
			TString line = Form("%-10s = %.4g #pm %.2g",
								var->GetName(), var->getVal(), var->getError());
			pave->AddText(line);
		}
		// chi2/ndf (uses number of floating parameters)
		const int nFloat = finalPars.getSize();
		const double chi2ndf = frame_Bmass->chiSquare("all", "data", nFloat);
		pave->AddText(Form("#chi^{2}/ndf = %.2f", chi2ndf));
	} else {
		pave->AddText("Fit result not available");
	}

	// IMPORTANT: attach pave BEFORE drawing the frame so it appears
	frame_Bmass->addObject(pave);

	// Now draw the top frame and legend
	frame_Bmass->GetXaxis()->SetLabelSize(0);
	frame_Bmass->GetXaxis()->SetTitleSize(0);
	frame_Bmass->Draw();
	leg1.DrawClone();

	// ---------------- Draw the pull panel on the bottom pad ----------------
	TPad* p_bot = (TPad*)c_Bmass->cd(2);
	p_bot->SetPad(0.0, 0.0, 1.0, 0.20);   // bottom 20%
	p_bot->SetTopMargin(0.02);
	p_bot->SetBottomMargin(0.35);
	p_bot->SetLeftMargin(0.12);
	p_bot->SetRightMargin(0.04);

	frame_pull->Draw();

	// Add a horizontal reference line at y=0 on the pull panel
	const double xmin = B_mass.getMin();
	const double xmax = B_mass.getMax();
	TLine* line0 = new TLine(xmin, 0.0, xmax, 0.0);
	line0->SetLineStyle(2);
	line0->Draw("same");

	// Save with a filename that indicates pull is included
	c_Bmass->SaveAs("./pdf_data/PsiX_data_OPT_v1.pdf");

    f->Close();
}