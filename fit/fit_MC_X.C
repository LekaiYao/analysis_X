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
    TFile *f = TFile::Open("/user/l/lekai/work/ppRef/analysis_X/selection/root_files/MC_X3872.root");
    TTree *tree = (TTree*)f->Get("tree");
    
    float lowBmass=3.84,highBmass=3.90;
    Int_t N=10000;//N is the max of events
	Int_t BinNum=60;
    
    //define variables
    RooRealVar B_mass("B_mass", "B_mass", lowBmass, highBmass);
    
    //define the RooArgSet of the Dataset
	RooArgSet vars;
		vars.add(B_mass);
    RooDataSet data("data", "data", vars,Import(*tree));
	
	//signal
	RooRealVar mean("mean","mean",3.87,3.86,3.88);
	RooRealVar sigma1("sigma1","sigma1",0.005,1e-05,0.02);
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

    // ============================================================
    // Drawing section with pull distribution and chi2/ndf
    // ============================================================

    TCanvas *c_Bmass = new TCanvas("c_Bmass","c_Bmass",10,10,800,700);

    // 1) Main frame WITHOUT title (avoid default "A RooPlot of ...")
    RooPlot *frame_Bmass = B_mass.frame(Bins(BinNum));
    frame_Bmass->SetTitle(""); // remove default title
    data.plotOn(frame_Bmass, DataError(RooAbsData::SumW2), Name("data"));
    tot.plotOn(frame_Bmass, Name("all"));
    tot.plotOn(frame_Bmass, Components(Gaus1), LineColor(1), LineStyle(1), Name("Gaus1"));
    tot.plotOn(frame_Bmass, Components(Gaus2), LineColor(3), LineStyle(1), Name("Gaus2"));

    // Override Y-axis label with fixed bin width (avoid 0.000999999 float formatting)
    frame_Bmass->GetYaxis()->SetTitle("Events / (0.001)");//change

    // 2) Pull frame WITHOUT title
    RooHist* hpull = frame_Bmass->pullHist("data", "all"); // (data - pdf) / sigma per bin
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
    frame_pull->GetXaxis()->SetTitle("m_{J/#psi#pi^{+}#pi^{-}}[GeV/c^{2}]");
    frame_pull->GetXaxis()->SetTitleSize(0.12);
    frame_pull->GetXaxis()->SetLabelSize(0.10);

    // 3) Split canvas into two pads: top 80%, bottom 20%
    c_Bmass->Divide(1,2,0,0);
    TPad* p_top = (TPad*)c_Bmass->cd(1);
    p_top->SetPad(0.0, 0.20, 1.0, 1.0);   // top 80%
    p_top->SetTopMargin(0.10);            // <-- increase top margin to move plot down
    p_top->SetBottomMargin(0.02);
    p_top->SetLeftMargin(0.12);
    p_top->SetRightMargin(0.04);

    // ---------------- Legend and TPaveText CREATED BEFORE DRAW ----------------
    // Compute data-space coordinates for legend/pave
    double x0 = lowBmass + 0.6*(highBmass-lowBmass);//change
    double x1 = lowBmass + 0.72*(highBmass-lowBmass);
    double x2 = lowBmass + 1.0*(highBmass-lowBmass);
    double y0 = 0.8*frame_Bmass->GetMaximum();
    double y1 = 0.8*frame_Bmass->GetMaximum();
    double y2 = 1.0*frame_Bmass->GetMaximum();

    // Legend: smaller, placed inside plot box
    TLegend leg1(x0, y0, x1, y2, "", "br");
    leg1.SetTextSize(0.018);
    leg1.AddEntry(frame_Bmass->findObject("all"),   "Gaus",  "L");
    leg1.AddEntry(frame_Bmass->findObject("Gaus1"), "Gaus1", "L");
    leg1.AddEntry(frame_Bmass->findObject("Gaus2"), "Gaus2", "L");

    // PaveText with fit parameters (attach to frame BEFORE drawing)
    TPaveText *pave = new TPaveText(x1, y1, x2, y2, "BR");
    pave->SetFillColor(0);
    pave->SetBorderSize(1);
    pave->SetTextSize(0.02);
    pave->SetTextAlign(12); // left alignment

    RooArgList finalPars = result->floatParsFinal();
    for (int i = 0; i < finalPars.getSize(); i++) {
        RooRealVar* var = (RooRealVar*)&finalPars[i];
        TString line = Form("%-6s = %.2e #pm %.2e",
                            var->GetName(), var->getVal(), var->getError());
        pave->AddText(line);
    }
    // chi2/ndf
    int nFloat = result ? result->floatParsFinal().getSize() : 0;
    double chi2ndf = frame_Bmass->chiSquare("all", "data", nFloat);
    pave->AddText(Form("#chi^{2}/ndf = %.2f", chi2ndf));

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
    double xmin = B_mass.getMin();
    double xmax = B_mass.getMax();
    TLine* line0 = new TLine(xmin, 0.0, xmax, 0.0);
    line0->SetLineStyle(2);
    line0->Draw("same");

    // Save with a new filename to indicate the pull panel is included
    c_Bmass->SaveAs("./pdf_MC_X/Xmass_MC_v1.pdf");

    f->Close();
}