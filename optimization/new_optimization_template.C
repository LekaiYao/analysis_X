#include <TFile.h>
#include <TTree.h>
#include <TCanvas.h>
#include <TH1F.h>
#include <TBox.h>
#include <TLegend.h>
#include <iostream>
#include <fstream>
#include <cmath>
#include <TStyle.h>
#include <TGraph.h>
#include <TLine.h>
#include <TMarker.h>

void optimization_@VAR@() {
    //TFile *file_data = TFile::Open("/user/u/u25lekai/work/ppRef/analysis_X/TMVA/sideband_PSI_BDT_v3test1.root");
    TFile *file_data = TFile::Open("/user/u/u25lekai/work/ppRef/analysis_X/selection/root_files/sideband.root");
    if (!file_data || file_data->IsZombie()) {
        std::cerr << "Failed to open data ROOT file." << std::endl;
        return;
    }
    TTree *tree_data = (TTree*)file_data->Get("tree");
    if (!tree_data) {
        std::cerr << "Failed to get 'tree' from data file." << std::endl;
        return;
    }
    int yellow = TColor::GetColor(204, 204,   0); 
	int orange = TColor::GetColor(204, 102,   0); 

    // --- MC for X(3872) (original) ---
    TFile *file_MC = TFile::Open("/user/u/u25lekai/work/ppRef/analysis_X/selection/root_files/MC_X3872.root");
    if (!file_MC || file_MC->IsZombie()) {
        std::cerr << "Failed to open MC_X3872 ROOT file." << std::endl;
        return;
    }
    TTree *tree_MC = (TTree*)file_MC->Get("tree");
    if (!tree_MC) {
        std::cerr << "Failed to get 'tree' from MC_X3872 file." << std::endl;
        return;
    }

    // --- MC for PSI(2S) (added) ---
    TFile *file_MC_psi = TFile::Open("/user/u/u25lekai/work/ppRef/analysis_X/selection/root_files/MC_PSI2S.root");
    if (!file_MC_psi || file_MC_psi->IsZombie()) {
        std::cerr << "Failed to open MC_PSI2S ROOT file." << std::endl;
        return;
    }
    TTree *tree_MC_psi = (TTree*)file_MC_psi->Get("tree");
    if (!tree_MC_psi) {
        std::cerr << "Failed to get 'tree' from MC_PSI2S file." << std::endl;
        return;
    }

    // Canvas for plotting
    TCanvas *c1 = new TCanvas("c1", "Entry Scan", 900, 700);
    c1->SetGrid();

    float min = @MIN@;
    float max = @MAX@;
    float step = @STEP@;
    int n_points = int((max - min) / step) + 1;
    int nbin = int((max - min) / step) + 1;

    TString var_name = "@VAR@";

    // l for "x < cut" (less-than), g for "x > cut" (greater-than)
    Double_t n_data_l = 0, n_MC_l = 0, FOM_l = 0;
    Double_t n_data_g = 0, n_MC_g = 0, FOM_g = 0;

    // for PSI(2S) overlay (added)
    Double_t n_MCpsi_l = 0, n_MCpsi_g = 0, FOMpsi_l = 0, FOMpsi_g = 0;

    TString flag = "";
    // scale factors (background and signal)
    Double_t fb = 0.329285, fs = 0.45116248;

    Double_t FOM_max_l = 0, cut_best_l = min;
    Double_t FOM_max_g = 0, cut_best_g = min;

    // track max for PSI(2S) (added)
    Double_t FOMpsi_max_l = 0, FOMpsi_max_g = 0;

    // Allocate graphs
    TGraph *graph_l = new TGraph(); // X(3872) < cut
    TGraph *graph_g = new TGraph(); // X(3872) > cut

    // Added: graphs for PSI(2S)
    TGraph *graph_l_psi = new TGraph(); // PSI(2S) < cut
    TGraph *graph_g_psi = new TGraph(); // PSI(2S) > cut

    // Histograms for cumulative integrals
    TH1F *hist_data     = new TH1F("hist_data",     "hist_data",     nbin, min, max);
    TH1F *hist_MC       = new TH1F("hist_MC",       "hist_MC",       nbin, min, max);
    TH1F *hist_MC_psi   = new TH1F("hist_MC_psi",   "hist_MC_psi",   nbin, min, max); // added

    // Fill histograms from trees
    TString cutStr = "B_chi2cl > 0.003 && B_Qvalueuj < 0.2";
    tree_data->Draw(Form("%s >> hist_data", var_name.Data()), cutStr);
    tree_MC  ->Draw(Form("%s >> hist_MC",   var_name.Data()), cutStr);
    tree_MC_psi->Draw(Form("%s >> hist_MC_psi", var_name.Data()), cutStr); // added

    // Loop over cut positions and compute FOM for both directions
    int ipt = 0;
    for (int i = 1; i < n_points; ++i) {
        float cut_val = min + i * step;

        // ------- LEFT: x < cut -------
        n_data_l   = hist_data->Integral(0, i);
        n_MC_l     = hist_MC->Integral(0, i);
        n_MCpsi_l  = hist_MC_psi->Integral(0, i);         // added

        FOM_l      = (n_data_l==0.0 && n_MC_l==0.0) ? 0.0
                     : (n_MC_l * fs) / std::sqrt(n_MC_l * fs + n_data_l * fb);
        FOMpsi_l   = (n_data_l==0.0 && n_MCpsi_l==0.0) ? 0.0     // uses same data sideband
                     : (n_MCpsi_l * fs) / std::sqrt(n_MCpsi_l * fs + n_data_l * fb); // added

        if (FOM_l > FOM_max_l) { FOM_max_l = FOM_l; cut_best_l = cut_val; }
        if (FOMpsi_l > FOMpsi_max_l) FOMpsi_max_l = FOMpsi_l;   // added

        graph_l->SetPoint(ipt, cut_val, FOM_l);
        graph_l_psi->SetPoint(ipt, cut_val, FOMpsi_l);          // added

        // ------- RIGHT: x > cut -------
        n_data_g   = hist_data->Integral(i, nbin + 1);
        n_MC_g     = hist_MC->Integral(i, nbin + 1);
        n_MCpsi_g  = hist_MC_psi->Integral(i, nbin + 1);        // added

        FOM_g      = (n_data_g==0.0 && n_MC_g==0.0) ? 0.0
                     : (n_MC_g * fs) / std::sqrt(n_MC_g * fs + n_data_g * fb);
        FOMpsi_g   = (n_data_g==0.0 && n_MCpsi_g==0.0) ? 0.0    // added
                     : (n_MCpsi_g * fs) / std::sqrt(n_MCpsi_g * fs + n_data_g * fb);

        if (FOM_g > FOM_max_g) { FOM_max_g = FOM_g; cut_best_g = cut_val; }
        if (FOMpsi_g > FOMpsi_max_g) FOMpsi_max_g = FOMpsi_g;   // added

        graph_g->SetPoint(ipt, cut_val, FOM_g);
        graph_g_psi->SetPoint(ipt, cut_val, FOMpsi_g);          // added

        ++ipt;
    }

    // Common axis/line/marker styling helpers (extended with color)
    auto style_graph = [&](TGraph* g, Color_t lc=kBlue+1, Style_t ms=20){
        g->SetLineColor(lc);
        g->SetMarkerColor(lc);
        g->SetLineWidth(2);
        g->SetMarkerStyle(ms);
        g->SetMarkerSize(0.9);
        g->SetTitle("");
        g->GetXaxis()->SetTitle(Form("%s", var_name.Data()));
        g->GetYaxis()->SetTitle("FOM");
        g->GetXaxis()->SetTitleSize(0.04);
        g->GetXaxis()->SetTitleOffset(1.05);
        g->GetYaxis()->SetTitleSize(0.04);
        g->GetYaxis()->SetTitleOffset(1.05);
    };

    // Y-range based on BOTH X(3872) and PSI(2S) (so overlay fits)
    const double yMax_combined = 1.3 * std::max(std::max(FOM_max_g, FOM_max_l),
                                                std::max(FOMpsi_max_g, FOMpsi_max_l));

    // Decide which direction is better (based on X(3872) ONLY) and draw
    if (FOM_max_g > FOM_max_l) {
        std::cout << "The best cut for " << var_name << " is greater than " << cut_best_g << std::endl;
        std::cout << "FOM is " << FOM_max_g << std::endl;
        flag = "g";

        style_graph(graph_g, yellow, 20);            // X(3872)
        graph_g->GetYaxis()->SetRangeUser(0, yMax_combined);
        graph_g->Draw("AL");

        // overlay PSI(2S)
        style_graph(graph_g_psi, orange, 24);       // PSI(2S)
        graph_g_psi->Draw("L SAME");

        // vertical line and best marker for X(3872)
        TLine *bestLine = new TLine(cut_best_g, 0, cut_best_g, FOM_max_g * 1.05);
        bestLine->SetLineColor(kOrange+7);
        bestLine->SetLineStyle(2);
        bestLine->SetLineWidth(2);
        bestLine->Draw("SAME");

        TMarker *bestMark = new TMarker(cut_best_g, FOM_max_g, 29);
        bestMark->SetMarkerColor(kMagenta+1);
        bestMark->SetMarkerSize(1.2);
        bestMark->Draw("SAME");

        // legend
        auto leg = new TLegend(0.58, 0.75, 0.90, 0.90); // NDC
        leg->SetTextSize(0.032);
        leg->SetTextFont(62);

        leg->AddEntry(graph_g,     "X(3872)",  "l");
        leg->AddEntry(graph_g_psi, "(#psi(2S))", "l");
        leg->AddEntry((TObject*)0, Form("Best cut : x > %.3f", cut_best_g), "h");
        leg->AddEntry((TObject*)0, Form("Max FOM = %.3f", FOM_max_g), "h");

        leg->Draw();

    } else {
        std::cout << "The best cut for " << var_name << " is less than " << cut_best_l << std::endl;
        std::cout << "FOM is " << FOM_max_l << std::endl;
        flag = "l";

        style_graph(graph_l, yellow, 20);            // X(3872)
        graph_l->GetYaxis()->SetRangeUser(0, yMax_combined);
        graph_l->Draw("AL");

        // overlay PSI(2S)
        style_graph(graph_l_psi, orange, 24);       // PSI(2S)
        graph_l_psi->Draw("L SAME");

        // vertical line and best marker for X(3872)
        TLine *bestLine = new TLine(cut_best_l, 0, cut_best_l, FOM_max_l * 1.05);
        bestLine->SetLineColor(kOrange+7);
        bestLine->SetLineStyle(2);
        bestLine->SetLineWidth(2);
        bestLine->Draw("SAME");

        TMarker *bestMark = new TMarker(cut_best_l, FOM_max_l, 29);
        bestMark->SetMarkerColor(kMagenta+1);
        bestMark->SetMarkerSize(1.2);
        bestMark->Draw("SAME");

        // legend
        auto leg = new TLegend(0.58, 0.75, 0.89, 0.90); // NDC
        leg->SetTextSize(0.032);
        leg->SetTextFont(62);
        leg->AddEntry(graph_l,     "X(3872)", "l");
        leg->AddEntry(graph_l_psi, "(#psi(2S))", "l");
        leg->AddEntry((TObject*)0, Form("Best cut : x < %.3f", cut_best_l), "");
        leg->AddEntry((TObject*)0, Form("Max FOM = %.3f", FOM_max_l), "");
        leg->Draw();
    }

    // Save
    gPad->RedrawAxis();
    c1->SaveAs(Form("./new_FOM/%s_%d_v1.pdf", var_name.Data(), nbin));

    // Summary (ONLY for X(3872), unchanged logic)
    std::ofstream fout("new_FOM/FOM_summary_v1.txt", std::ios::app);
    if (fout.is_open()) {
        if (flag == "g"){
            fout << "The best cut for " << var_name.Data() << " is greater than " << cut_best_g
                 << " FOM_max=" << FOM_max_g << std::endl;
        } else {
            fout << "The best cut for " << var_name.Data() << " is less than " << cut_best_l
                 << " FOM_max=" << FOM_max_l << std::endl;
        }
        fout << "The range of " << var_name.Data() << " is [ " << min << " , " << max << " ]" << std::endl;
        fout.close();
    } else {
        std::cerr << "Error opening FOM_summary_v1.txt for writing!" << std::endl;
    }

    file_data->Close();
    file_MC->Close();
    file_MC_psi->Close(); // added
}
