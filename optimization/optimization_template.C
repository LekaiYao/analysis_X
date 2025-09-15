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
    TFile *file_data = TFile::Open("/user/l/lekai/work/ppRef/analysis_X/TMVA/sideband_X3872_BDT_trainX4.root");
    //TFile *file_data = TFile::Open("/user/l/lekai/work/ppRef/analysis_X/selection/root_files/sideband.root");
    if (!file_data || file_data->IsZombie()) {
        std::cerr << "Failed to open data ROOT file." << std::endl;
        return;
    }
    TTree *tree_data = (TTree*)file_data->Get("tree");
    if (!tree_data) {
        std::cerr << "Failed to get 'tree' from data file." << std::endl;
        return;
    }

    TFile *file_MC = TFile::Open("/user/l/lekai/work/ppRef/analysis_X/TMVA/MC_X3872_BDT_trainX4.root");
    //TFile *file_MC = TFile::Open("/user/l/lekai/work/ppRef/analysis_X/selection/root_files/MC_PSI2S.root");
    //TFile *file_MC = TFile::Open("/user/l/lekai/work/ppRef/analysis_X/selection/root_files/MC_X3872.root");
    if (!file_MC || file_MC->IsZombie()) {
        std::cerr << "Failed to open MC ROOT file." << std::endl;
        return;
    }
    TTree *tree_MC = (TTree*)file_MC->Get("tree");
    if (!tree_MC) {
        std::cerr << "Failed to get 'tree' from MC file." << std::endl;
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
    TString flag = "";
    // scale factors (background and signal)
    Double_t fb = 0.329285, fs = 0.45116248;

    Double_t FOM_max_l = 0, cut_best_l = min;
    Double_t FOM_max_g = 0, cut_best_g = min;

    // Allocate graphs
    TGraph *graph_l = new TGraph(); // will set size after filling
    TGraph *graph_g = new TGraph();

    // Histograms for cumulative integrals
    TH1F *hist_data = new TH1F("hist_data", "hist_data", nbin, min, max);
    TH1F *hist_MC   = new TH1F("hist_MC",   "hist_MC",   nbin, min, max);

    // Fill histograms from trees, pre-cuts
    //tree_data->Draw(Form("%s >> hist_data", var_name.Data()),"B_chi2cl > 0.003 && B_Qvalueuj < 0.2");
    //tree_MC->Draw(Form("%s >> hist_MC", var_name.Data()),"B_chi2cl > 0.003 && B_Qvalueuj < 0.2");
    tree_data->Draw(Form("%s >> hist_data", var_name.Data()));
    tree_MC->Draw(Form("%s >> hist_MC", var_name.Data()));

    // Loop over cut positions and compute FOM for both directions
    // Note: we set points starting from index 0; total points used = n_points-1
    int ipt = 0;
    for (int i = 1; i < n_points; ++i) {
        float cut_val = min + i * step;

        // LEFT: x < cut (use integral from underflow up to bin i)
        n_data_l = hist_data->Integral(0, i);
        n_MC_l   = hist_MC->Integral(0, i);
        FOM_l    = (n_data_l == 0.0 && n_MC_l == 0.0) ? 0.0
                    : (n_MC_l * fs) / std::sqrt(n_MC_l * fs + n_data_l * fb);
        if (FOM_l > FOM_max_l) {
            FOM_max_l  = FOM_l;
            cut_best_l = cut_val;
        }
        graph_l->SetPoint(ipt, cut_val, FOM_l);

        // RIGHT: x > cut (use integral from bin i to overflow)
        n_data_g = hist_data->Integral(i, nbin + 1); // include overflow
        n_MC_g   = hist_MC->Integral(i, nbin + 1);
        FOM_g    = (n_data_g == 0.0 && n_MC_g == 0.0) ? 0.0
                    : (n_MC_g * fs) / std::sqrt(n_MC_g * fs + n_data_g * fb);
        if (FOM_g > FOM_max_g) {
            FOM_max_g  = FOM_g;
            cut_best_g = cut_val;
        }
        graph_g->SetPoint(ipt, cut_val, FOM_g);

        ++ipt;
    }

    // Common axis/line/marker styling helpers
    auto style_graph = [&](TGraph* g){
        g->SetLineColor(kBlue+1);
        g->SetLineWidth(2);
        g->SetMarkerStyle(20);
        g->SetMarkerSize(0.9);
        g->SetTitle("");
        g->GetXaxis()->SetTitle(Form("Cut Value (%s)", var_name.Data()));
        g->GetYaxis()->SetTitle("FOM (S / #sqrt{S+B})");
        g->GetXaxis()->SetTitleSize(0.045);
        g->GetXaxis()->SetTitleOffset(1.1);
        g->GetYaxis()->SetTitleSize(0.045);
        g->GetYaxis()->SetTitleOffset(1.2);
    };

    // Y-range based on the better of the two
    const double yMax = 1.3 * std::max(FOM_max_g, FOM_max_l);

    // Decide which direction is better and draw accordingly on c1
    if (FOM_max_g > FOM_max_l) {
        std::cout << "The best cut for " << var_name << " is greater than " << cut_best_g << std::endl;
        std::cout << "FOM is " << FOM_max_g << std::endl;
        flag = "g";

        style_graph(graph_g);
        graph_g->GetYaxis()->SetRangeUser(0, yMax);
        graph_g->Draw("AL");

        // vertical line and best marker
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
        auto leg = new TLegend(0.60, 0.78, 0.88, 0.90); // (x1,y1,x2,y2) NDC
        leg->SetTextSize(0.035);
        leg->SetTextFont(62);
        leg->AddEntry(graph_g, "FOM vs cut", "l");
        leg->AddEntry((TObject*)0, Form("Best cut: x > %.4f", cut_best_g), "");
        leg->AddEntry((TObject*)0, Form("Max FOM = %.4f", FOM_max_g), "");
        leg->Draw();

    } else {
        std::cout << "The best cut for " << var_name << " is less than " << cut_best_l << std::endl;
        std::cout << "FOM is " << FOM_max_l << std::endl;
        flag = "l";

        style_graph(graph_l);
        graph_l->GetYaxis()->SetRangeUser(0, yMax);
        graph_l->Draw("AL");

        // vertical line and best marker
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
        auto leg = new TLegend(0.60, 0.78, 0.88, 0.90); // (x1,y1,x2,y2) NDC
        leg->SetTextSize(0.035);
        leg->SetTextFont(62);
        leg->AddEntry(graph_l, "FOM vs cut", "l");
        leg->AddEntry((TObject*)0, Form("Best cut: x < %.4f", cut_best_l), "");
        leg->AddEntry((TObject*)0, Form("Max FOM = %.4f", FOM_max_l), "");
        leg->Draw();
    }

    // Keep original saving style and name
    gPad->RedrawAxis();
    c1->SaveAs(Form("./FOM_X/%s_%d_v1.pdf", var_name.Data(), nbin));

    // Keep original summary file output
    std::ofstream fout("FOM_X/FOM_summary_v1.txt", std::ios::app);
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
}
