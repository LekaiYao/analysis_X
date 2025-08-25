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


void optimization_@VAR@() {
    //TFile *file_data = TFile::Open("/user/u/u25lekai/work/ppRef/analysis_X/TMVA/sideband_PSI_BDT_v3test1.root");
    TFile *file_data = TFile::Open("/user/u/u25lekai/work/ppRef/analysis_X/selection/root_files/sideband.root");
    TTree *tree_data = (TTree*)file_data->Get("tree");

    //TFile *file_MC = TFile::Open("/user/u/u25lekai/work/ppRef/analysis_X/TMVA/MC_PSI2S_BDT_v3test1.root");
    TFile *file_MC = TFile::Open("/user/u/u25lekai/work/ppRef/analysis_X/selection/root_files/MC_PSI2S.root");
    //TFile *file_MC = TFile::Open("/user/u/u25lekai/work/ppRef/analysis_X/selection/root_files/MC_X3872.root");
    TTree *tree_MC = (TTree*)file_MC->Get("tree");

    TCanvas *c1 = new TCanvas("c1", "Entry Scan", 800, 600);

    float min = @MIN@;
    float max = @MAX@;
    float step = @STEP@;
    int n_points = int((max - min) / step) + 1;
    int nbin = int((max - min) / step) + 1;

    TString var_name = "@VAR@";
    Double_t n_data_l, n_MC_l, FOM_l, n_data_g, n_MC_g, FOM_g;//l for less than, g for greater than
    TString flag = "";
    Double_t fb = 0.559271, fs = 13.657809;
    Double_t FOM_max_l = 0, cut_best_l = min, FOM_max_g = 0, cut_best_g = min;

    TGraph *graph_l = new TGraph(n_points);
    TGraph *graph_g = new TGraph(n_points);
    TH1F *hist_data = new TH1F("hist_data", "hist_data", nbin, min, max);
    TH1F *hist_MC = new TH1F("hist_MC", "hist_MC", nbin, min, max);

    tree_data->Draw(Form("%s >> hist_data", var_name.Data()));
    tree_MC->Draw(Form("%s >> hist_MC", var_name.Data()));

    for (int i = 1; i < n_points; ++i) {
        float cut_val = min + i * step;
        n_data_l = hist_data->Integral(0, i);
        n_MC_l = hist_MC->Integral(0, i);
        FOM_l = (n_data_l == 0 && n_MC_l == 0) ? 0 : n_MC_l * fs / sqrt(n_MC_l * fs + n_data_l * fb);
        if (FOM_l > FOM_max_l) {
            FOM_max_l = FOM_l;
            cut_best_l = cut_val;
        }
        graph_l->SetPoint(i-1, cut_val, FOM_l);

        n_data_g = hist_data->Integral(i-1, nbin);
        n_MC_g = hist_MC->Integral(i-1, nbin);
        FOM_g = (n_data_g == 0 && n_MC_g == 0) ? 0 : n_MC_g * fs / sqrt(n_MC_g * fs + n_data_g * fb);

        if (FOM_g > FOM_max_g) {
            FOM_max_g = FOM_g;
            cut_best_g = cut_val;
        }
        graph_g->SetPoint(i-1, cut_val, FOM_g);
    }

    graph_l->SetTitle(Form("FOM vs %s", var_name.Data()));
    graph_l->SetLineColor(kBlue);
    graph_l->SetLineWidth(2);
    graph_l->SetMarkerStyle(20);
    graph_l->SetMarkerSize(1.0);

    graph_g->SetTitle(Form("FOM vs %s", var_name.Data()));
    graph_g->SetLineColor(kBlue);
    graph_g->SetLineWidth(2);
    graph_g->SetMarkerStyle(20);
    graph_g->SetMarkerSize(1.0);
    

    if (FOM_max_g > FOM_max_l){
        std::cout << "The best cut for " << var_name << " is greater than " << cut_best_g << std::endl;
        std::cout << "FOM is " << FOM_max_g << std::endl;
        graph_g->Draw("AL");
        flag = "g";
    }
    else{
        std::cout << "The best cut for " << var_name << " is less than " << cut_best_l << std::endl;
        std::cout << "FOM is " << FOM_max_l << std::endl;
        graph_l->Draw("AL");
        flag = "l";
    }


    c1->SaveAs(Form("./FOM_PSI/%s_%d_v1.pdf", var_name.Data(), nbin));

    std::ofstream fout("FOM_PSI/FOM_summary_v1.txt", std::ios::app); 
    if (fout.is_open()) {
        if (flag == "g"){
            fout
             << "The best cut for " << var_name.Data() << " is greater than " << cut_best_g
             << " FOM_max=" << FOM_max_g << std::endl;    
        }
        else{
            fout
             << "The best cut for " << var_name.Data() << " is less than " << cut_best_l
             << " FOM_max=" << FOM_max_l << std::endl;
        }
        fout << "The range of " << var_name.Data() << "  is [ " << min << " , " << max << " ]" << std::endl;
        fout.close();
    } else {
        std::cerr << "Error opening FOM_summary.txt for writing!" << std::endl;
    }

    file_data->Close();
    file_MC->Close();
}
