#include "TFile.h"
#include "TTree.h"
#include "TCanvas.h"
#include "TH1F.h"
#include "TLegend.h"
#include "TStyle.h"

void plot_var() {
    // --------- User-defined variable name ---------
    TString varName = "BDT_score";  // Variable to inspect

    // --------- Automatically determined variable range ---------
    Float_t var_min = 1e10;
    Float_t var_max = -1e10;

    // Open ROOT files
    //TFile *f_data = TFile::Open("/user/l/lekai/work/ppRef/analysis_X/selection/root_files/sideband.root");
    //TFile *f_mc   = TFile::Open("/user/l/lekai/work/ppRef/analysis_X/selection/root_files/MC_PSI2S.root");
    //TFile *f_mc2  = TFile::Open("/user/l/lekai/work/ppRef/analysis_X/selection/root_files/MC_X3872.root");
    TFile *f_data = TFile::Open("/user/l/lekai/work/ppRef/analysis_X/TMVA/sideband_X3872_BDT_trainX4.root");
    TFile *f_mc   = TFile::Open("/user/l/lekai/work/ppRef/analysis_X/TMVA/MC_PSI2S_BDT_trainX4.root");
    TFile *f_mc2  = TFile::Open("/user/l/lekai/work/ppRef/analysis_X/TMVA/MC_X3872_BDT_trainX4.root");

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
    //Define the region by hand
    var_min = -0.25;
    var_max = 0.15;

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

    // ---------- Vertical cut line at x = 0.2 on main plot by hand----------
    double cut_x = 0.2;
    TLine *cutLine_main = new TLine(cut_x, 0.0, cut_x, h_data->GetMaximum());
    cutLine_main->SetLineColor(TColor::GetColor(34,139,34)); 
    cutLine_main->SetLineStyle(2); // dashed
    cutLine_main->SetLineWidth(2);
    cutLine_main->Draw("same");

    // Add legend
    TLegend *leg = new TLegend(0.65, 0.70, 0.88, 0.88);
    leg->AddEntry(h_data, "sideband", "lf");
    leg->AddEntry(h_mc,   "#Psi(2S)",   "lf");
    leg->AddEntry(h_mc2,  "X(3872)",  "lf");
    leg->Draw();

    c1->SaveAs("pre-cut/"+varName + "_compare.pdf");
}
