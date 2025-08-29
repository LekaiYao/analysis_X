#include "TFile.h"
#include "TTree.h"
#include "TCanvas.h"
#include "TH1F.h"
#include "TLegend.h"
#include "TStyle.h"

void plot_var_zoom() {
    // --------- User-defined variable name ---------
    TString varName = "B_chi2cl";  // Variable to inspect

    // --------- Automatically determined variable range ---------
    Float_t var_min = 1e10;
    Float_t var_max = -1e10;

    // Open ROOT files
    TFile *f_data = TFile::Open("/user/u/u25lekai/work/ppRef/analysis_X/selection/root_files/sideband.root");
    TFile *f_mc   = TFile::Open("/user/u/u25lekai/work/ppRef/analysis_X/selection/root_files/MC_PSI2S.root");
    TFile *f_mc2  = TFile::Open("/user/u/u25lekai/work/ppRef/analysis_X/selection/root_files/MC_X3872.root");

    TTree *tree_data = (TTree*)f_data->Get("tree");
    TTree *tree_mc   = (TTree*)f_mc->Get("tree");
    TTree *tree_mc2  = (TTree*)f_mc2->Get("tree");

    Float_t var_value;
    tree_data->SetBranchAddress(varName, &var_value);

    // Determine min/max from data
    Long64_t nentries = tree_data->GetEntries();
    for (Long64_t i = 0; i < nentries; ++i) {
        tree_data->GetEntry(i);
        if (!std::isnan(var_value) && !std::isinf(var_value)) {
            if (var_value < var_min) var_min = var_value;
            if (var_value > var_max) var_max = var_value;
        }
    }

    // Also determine min/max from MC
    tree_mc->SetBranchAddress(varName, &var_value);
    Long64_t nentries_mc = tree_mc->GetEntries();
    for (Long64_t i = 0; i < nentries_mc; ++i) {
        tree_mc->GetEntry(i);
        if (!std::isnan(var_value) && !std::isinf(var_value)) {
            if (var_value < var_min) var_min = var_value;
            if (var_value > var_max) var_max = var_value;
        }
    }

    // Also determine min/max from MC2
    tree_mc2->SetBranchAddress(varName, &var_value);
    Long64_t nentries_mc2 = tree_mc2->GetEntries();
    for (Long64_t i = 0; i < nentries_mc2; ++i) {
        tree_mc2->GetEntry(i);
        if (!std::isnan(var_value) && !std::isinf(var_value)) {
            if (var_value < var_min) var_min = var_value;
            if (var_value > var_max) var_max = var_value;
        }
    }

    int nBins = 100;
    TH1F *h_data = new TH1F("h_data", varName, nBins, var_min, var_max);
    TH1F *h_mc   = new TH1F("h_mc",   varName, nBins, var_min, var_max);
    TH1F *h_mc2  = new TH1F("h_mc2",  varName, nBins, var_min, var_max);

    // Fill histograms
    tree_data->Draw(varName + ">>h_data", "", "goff");
    tree_mc->Draw(varName + ">>h_mc", "", "goff");
    tree_mc2->Draw(varName + ">>h_mc2", "", "goff");

    // Normalize histograms
    if (h_data->Integral() > 0) h_data->Scale(1.0 / h_data->Integral());
    if (h_mc->Integral() > 0)   h_mc->Scale(1.0 / h_mc->Integral());
    if (h_mc2->Integral() > 0)  h_mc2->Scale(1.0 / h_mc2->Integral());

    // Set styles
    int brightAzure = TColor::GetColor(51, 153, 255);
    h_data->SetLineColor(brightAzure);
    h_data->SetLineWidth(2);
    h_data->SetFillColor(brightAzure);
    h_data->SetFillStyle(3345);  // diagonal hatch

    h_mc->SetLineColor(kOrange + 7);
    h_mc->SetLineWidth(2);
    h_mc->SetFillColorAlpha(kOrange + 7, 0.4); // semi-transparent orange
    h_mc->SetFillStyle(1001);

    int lightBrightYellow = TColor::GetColor(255, 255, 0);  
    h_mc2->SetLineColor(lightBrightYellow);
    h_mc2->SetLineWidth(2);
    h_mc2->SetFillColorAlpha(lightBrightYellow, 0.2);// semi-transparent yellow
    h_mc2->SetFillStyle(1001);

    float max_val = std::max({h_data->GetMaximum(), h_mc->GetMaximum(), h_mc2->GetMaximum()});
    h_data->SetMaximum(max_val * 1.2);

    // Draw to canvas
    TCanvas *c1 = new TCanvas("c1", "Comparison", 800, 600);
    gStyle->SetOptStat(0);
    h_data->SetTitle("Comparison of " + varName + "; " + varName + "; Normalized Entries");
    h_data->Draw("hist");
    h_mc->Draw("hist same");
    h_mc2->Draw("hist same");

    // Add legend
    TLegend *leg = new TLegend(0.65, 0.70, 0.88, 0.88);
    leg->AddEntry(h_data, "sideband", "lf");
    leg->AddEntry(h_mc,   "#Psi(2S)",   "lf");
    leg->AddEntry(h_mc2,  "X(3872)",  "lf");
    leg->Draw();

    // ---------- Inset zoom of [0, 0.05] with dedicated binning ----------
    c1->cd();

    // Create a transparent inset pad in NDC of the main canvas
    TPad* inset = new TPad("inset","inset", 0.2, 0.2, 0.65, 0.65);
    inset->SetBorderMode(1);
    inset->SetTopMargin(0.08);
    inset->SetBottomMargin(0.25);
    inset->SetLeftMargin(0.18);
    inset->SetRightMargin(0.05);
    inset->Draw();
    inset->cd();

    // Define binning for zoom region
    const int    zoom_bins = 50;   // you can change to 40, 50, etc.
    const double zoom_low  = 0.0;
    const double zoom_high = 0.05;

    // New histograms for zoomed region
    TH1F* h_data_zoom = new TH1F("h_data_zoom", "", zoom_bins, zoom_low, zoom_high);
    TH1F* h_mc_zoom   = new TH1F("h_mc_zoom",   "", zoom_bins, zoom_low, zoom_high);
    TH1F* h_mc2_zoom  = new TH1F("h_mc2_zoom",  "", zoom_bins, zoom_low, zoom_high);

    // Fill zoom histograms
    tree_data->Draw(varName + ">>h_data_zoom", "", "goff");
    tree_mc  ->Draw(varName + ">>h_mc_zoom",   "", "goff");
    tree_mc2 ->Draw(varName + ">>h_mc2_zoom",  "", "goff");

    // Normalize
    if (h_data_zoom->Integral() > 0) h_data_zoom->Scale(1.0 / h_data_zoom->Integral());
    if (h_mc_zoom->Integral()   > 0) h_mc_zoom->Scale(1.0 / h_mc_zoom->Integral());
    if (h_mc2_zoom->Integral()  > 0) h_mc2_zoom->Scale(1.0 / h_mc2_zoom->Integral());

    // Inherit styles from main histograms
    h_data_zoom->SetLineColor(h_data->GetLineColor());
    h_data_zoom->SetLineWidth(2);
    h_data_zoom->SetFillColor(h_data->GetFillColor());
    h_data_zoom->SetFillStyle(h_data->GetFillStyle());

    h_mc_zoom->SetLineColor(h_mc->GetLineColor());
    h_mc_zoom->SetLineWidth(2);
    h_mc_zoom->SetFillColorAlpha(h_mc->GetFillColor(), 0.4);
    h_mc_zoom->SetFillStyle(h_mc->GetFillStyle());

    h_mc2_zoom->SetLineColor(h_mc2->GetLineColor());
    h_mc2_zoom->SetLineWidth(2);
    h_mc2_zoom->SetFillColorAlpha(h_mc2->GetFillColor(), 0.2);
    h_mc2_zoom->SetFillStyle(h_mc2->GetFillStyle());

    // Adjust axis
    h_data_zoom->SetTitle("");
    h_data_zoom->GetXaxis()->SetTitle(varName);
    h_data_zoom->GetYaxis()->SetTitle("Normalized Entries");
    h_data_zoom->GetXaxis()->SetLabelSize(0.04);
    h_data_zoom->GetYaxis()->SetLabelSize(0.04);
    h_data_zoom->GetXaxis()->SetTitleSize(0.06);
    h_data_zoom->GetYaxis()->SetTitleSize(0.06);
    h_data_zoom->GetYaxis()->SetTitleOffset(0.9);

    // Auto-scale maximum
    double ymax = std::max({h_data_zoom->GetMaximum(),
                            h_mc_zoom->GetMaximum(),
                            h_mc2_zoom->GetMaximum()});
    h_data_zoom->SetMaximum(ymax * 1.2);

    // Draw zoomed histograms
    h_data_zoom->Draw("hist");
    h_mc_zoom  ->Draw("hist same");
    h_mc2_zoom ->Draw("hist same");

    // ---------- Vertical cut line at B_chi2cl = 0.003 ----------
    double cut_x = 0.003;
    TLine *cutLine = new TLine(cut_x, 0.0, cut_x, h_data_zoom->GetMaximum());
    cutLine->SetLineColor(TColor::GetColor(34,139,34));
    cutLine->SetLineStyle(2); // dashed line
    cutLine->SetLineWidth(2);
    cutLine->Draw("same");


    c1->SaveAs("pre-cut/"+varName + "_compare_zoom.pdf");
}
