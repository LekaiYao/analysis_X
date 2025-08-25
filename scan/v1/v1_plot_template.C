#include "TFile.h"
#include "TTree.h"
#include "TCanvas.h"
#include "TH1F.h"
#include "TLegend.h"
#include "TStyle.h"

void plot_@VARNAME@() {
    // --------- User-defined variable name ---------
    TString varName = "@VARNAME@";  // Variable to inspect

    // --------- Automatically determined variable range ---------
    Float_t var_min = 1e10;
    Float_t var_max = -1e10;

    // Open ROOT files
    TFile *f_data = TFile::Open("/user/u/u25lekai/work/ppRef/analysis_X/selection/root_files/DATA_XPSI.root");
    TFile *f_mc   = TFile::Open("/user/u/u25lekai/work/ppRef/analysis_X/selection/root_files/MC_PSI2S.root");

    TTree *tree_data = (TTree*)f_data->Get("tree");
    TTree *tree_mc   = (TTree*)f_mc->Get("tree");

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

    // --------- Optional: manually set variable range ---------
    // var_min = 5.0;
    // var_max = 5.6;

    int nBins = 100;
    TH1F *h_data = new TH1F("h_data", varName, nBins, var_min, var_max);
    TH1F *h_mc   = new TH1F("h_mc",   varName, nBins, var_min, var_max);

    // Fill histograms
    tree_data->Draw(varName + ">>h_data", "", "goff");
    tree_mc->Draw(varName + ">>h_mc", "", "goff");

    // Normalize histograms
    if (h_data->Integral() > 0) h_data->Scale(1.0 / h_data->Integral());
    if (h_mc->Integral() > 0)   h_mc->Scale(1.0 / h_mc->Integral());

    // Set styles
    h_data->SetLineColor(kBlue);
    h_data->SetLineWidth(2);
    h_mc->SetLineColor(kOrange + 7);
    h_mc->SetLineWidth(2);

    float max_val = std::max(h_data->GetMaximum(), h_mc->GetMaximum());
    h_data->SetMaximum(max_val * 1.1);  // ensure all content fits in y-axis

    // Draw to canvas
    TCanvas *c1 = new TCanvas("c1", "Comparison", 800, 600);
    h_data->SetTitle("Comparison of " + varName + "; " + varName + "; Normalized Entries");
    h_data->Draw("hist");
    h_mc->Draw("hist same");

    // Add legend
    TLegend *leg = new TLegend(0.65, 0.75, 0.88, 0.88);
    leg->AddEntry(h_data, "Data", "l");
    leg->AddEntry(h_mc,   "MC",   "l");
    leg->Draw();

    c1->SaveAs("output_scan/"+varName + "_compare.pdf");
}
