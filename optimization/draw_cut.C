#include "TFile.h"
#include "TTree.h"
#include "TCanvas.h"
#include "TH1F.h"
#include "TStyle.h"
#include "TPaveStats.h"
#include "TError.h"
#include "TSystem.h"

void draw_cut() {
    // ---------------- Global style (apply before any drawing) ----------------
    gStyle->SetOptTitle(0);        // remove default histogram/pad title
    gStyle->SetOptStat(1110);      // show Entries/Mean/RMS, but hide the "name" line
    gStyle->SetStatW(0.24);        // width of stats box (NDC)
    gStyle->SetStatH(0.16);        // height of stats box (NDC)
    gStyle->SetStatFontSize(0.03);
    // Note: we will place stats box manually relative to the frame (plot area), so SetStatX/Y is not critical

    // ---------------- Open file and get tree ----------------
    //TFile *f = TFile::Open("/user/l/lekai/work/ppRef/analysis_X/selection/root_files/DATA_XPSI.root");
    TFile *f = TFile::Open("/user/l/lekai/work/ppRef/analysis_X/TMVA/DATA_XPSI_BDT_trainX4.root");
    if (!f || f->IsZombie()) { ::Error("draw_cut","Cannot open input ROOT file."); return; }
    TTree *tree = (TTree*)f->Get("tree");
    if (!tree) { ::Error("draw_cut","Cannot find TTree 'tree' in the file."); return; }

    // ---------------- User cut and histogram ----------------
    //TString cut = "B_chi2cl > 0.003 && B_Qvalueuj < 0.2 && B_Qvalueuj < 0.095 && B_trk1dR<0.639 && B_trk2dR<0.64";
    TString cut = "BDT_score > -0.000999987";
    const int nBins = 160;
    const double xMin = 3.6, xMax = 4.0;
    TH1F *h_Bmass = new TH1F("h_Bmass", "", nBins, xMin, xMax); // empty title to avoid top title

    // ---------------- Canvas before drawing (so stats box belongs to this pad) ----------------
    TCanvas *c1 = new TCanvas("c1","B_mass distribution",800,600);
    c1->SetLeftMargin(0.10);   // keep your margins; adjust if y-title needs more room
    c1->SetRightMargin(0.06);
    c1->SetTopMargin(0.06);
    c1->SetBottomMargin(0.10);

    // ---------------- Fill histogram (Tree::Draw projection) ----------------
    // Use "goff" to avoid implicit drawing by TTree::Draw
    tree->Draw("B_mass >> h_Bmass", cut, "goff");

    // ---------------- Axis titles and text sizes ----------------
    h_Bmass->GetXaxis()->SetTitle("m_{J/#Psi#pi+#pi-} [GeV/c^{2}]");
    h_Bmass->GetYaxis()->SetTitle("Entries/(0.0025)");

    h_Bmass->GetXaxis()->SetLabelSize(0.03);
    h_Bmass->GetXaxis()->SetTitleSize(0.035);
    h_Bmass->GetXaxis()->SetTitleOffset(1.10);

    h_Bmass->GetYaxis()->SetLabelSize(0.027);
    h_Bmass->GetYaxis()->SetTitleSize(0.03);
    h_Bmass->GetYaxis()->SetTitleOffset(1.30);

    int brightAzure = TColor::GetColor(51, 153, 255);
    h_Bmass->SetLineColor(brightAzure);
    h_Bmass->SetLineWidth(2);
    h_Bmass->SetFillColor(brightAzure);
    h_Bmass->SetFillStyle(3345);

    // ---------------- Y-axis range control ----------------
    // Option A: manual range (set use_manual=true and set yMin/yMax).
    bool   use_manual = false;
    double yMin = 14000.0;  // set e.g. 0.0
    double yMax = 59000.0;  // set e.g. 1200.0

    // Option B: automatic range with headroom (default if use_manual==false).
    double headroom = 1.2; // 20% room above the tallest bin

    if (use_manual) {
        h_Bmass->SetMinimum(yMin);
        h_Bmass->SetMaximum(yMax);
        // Alternatively: h_Bmass->GetYaxis()->SetRangeUser(yMin, yMax);
    } else {
        double ymax = h_Bmass->GetBinContent(h_Bmass->GetMaximumBin());
        if (ymax <= 0) ymax = 1.0;
        h_Bmass->SetMaximum(ymax * headroom);
    }

    // ---------------- Draw histogram ----------------
    h_Bmass->SetLineWidth(2);
    h_Bmass->Draw("HIST");

    // First update so that TPaveStats is created
    c1->Update();

    // ---------------- Precisely position the stats box at the frame's top-right ----------------
    // Align to the plotting rectangle (frame) instead of pad edge: use pad margins
    if (TPaveStats *st = (TPaveStats*)h_Bmass->FindObject("stats")) {
        const double xr = 1.0 - c1->GetRightMargin(); // frame right edge in pad-NDC
        const double yt = 1.0 - c1->GetTopMargin();   // frame top edge in pad-NDC

        const double w = gStyle->GetStatW();          // stats box width (NDC)
        const double h = gStyle->GetStatH();          // stats box height (NDC)

        st->SetX2NDC(xr);                     // anchor outer top-right to frame corner
        st->SetY2NDC(yt);
        st->SetX1NDC(st->GetX2NDC() - w);
        st->SetY1NDC(st->GetY2NDC() - h);
        st->Draw();
    }

    c1->Update();

    // ---------------- Ensure output directory exists and save ----------------
    gSystem->mkdir("pdf_cut", kTRUE); // create directory if missing
    c1->SaveAs("pdf_cut/mass_BDTcuts.pdf");
}
