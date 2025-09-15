void select() {
    // Open the input ROOT file
    //X(3872)
    //TFile *f_in = TFile::Open("root_files/DATA_XPSI_nonan.root");
    //TFile *f_in = TFile::Open("root_files/MC_X3872.root");
    //TFile *f_in = TFile::Open("root_files/sideband_X3872.root");
    //psi(2S)
    //TFile *f_in = TFile::Open("root_files/DATA_XPSI.root");
    TFile *f_in = TFile::Open("root_files/MC_PSI2S.root");
    //TFile *f_in = TFile::Open("root_files/sideband_PSI.root");

    //TFile *f_in = TFile::Open("root_files/sideband.root");
    TTree *tree = (TTree*)f_in->Get("tree");

    // Get all branches
    TObjArray *branches = tree->GetListOfBranches();

    // Use a map to store variable name -> value
    std::map<std::string, float> vars;

    // Dynamically bind variables for all Float_t scalar branches
    for (int i = 0; i < branches->GetEntries(); ++i) {
        TBranch *br = (TBranch*)branches->At(i);
        TLeaf *leaf = br->GetLeaf(br->GetName());
        if (leaf && strcmp(leaf->GetTypeName(), "Float_t") == 0 && leaf->GetLen() == 1) {
            vars[br->GetName()] = 0.0f;  // initialize
            tree->SetBranchAddress(br->GetName(), &vars[br->GetName()]);
        }
    }

    // Check if B_mass exists
    if (vars.find("B_mass") == vars.end()) {
        Error("select", "Branch 'B_mass' not found or not of Float_t type!");
        return;
    }

    // Create the output file and clone the structure
    //X(3872)
    //TFile *f_out = TFile::Open("DATA_XPSI_v2test1.root", "RECREATE");
    //TFile *f_out = TFile::Open("MC_X3872_trainX1.root", "RECREATE");
    //TFile *f_out = TFile::Open("sideband_X3872_v2test1.root", "RECREATE");
    //psi(2S)
    //TFile *f_out = TFile::Open("DATA_XPSI_cut0.root", "RECREATE");
    TFile *f_out = TFile::Open("MC_PSI2S_trainX1.root", "RECREATE");
    //TFile *f_out = TFile::Open("sideband_PSI_v3test1.root", "RECREATE");

    //TFile *f_out = TFile::Open("sideband.root", "RECREATE");
    //TFile *f_out = TFile::Open("sideband_trainX1.root", "RECREATE");

    TTree *new_tree = tree->CloneTree(0);

    Long64_t nEntries = tree->GetEntries();
    Long64_t kept = 0;

    for (Long64_t i = 0; i < nEntries; ++i) {
        tree->GetEntry(i);
        bool valid = true;

        // B_mass window cut
        //psi(2S)
        //if (valid && (vars["B_mass"] > 3.66 && vars["B_mass"] < 3.72)) valid = false;
        if (valid && (vars["B_mass"] < 3.66 || vars["B_mass"] > 3.72)) valid = false;
        //X(3872)
        //if (valid && (vars["B_mass"] > 3.83 && vars["B_mass"] < 3.91)) valid = false;//sideband
        //if (valid && (vars["B_mass"] < 3.83 || vars["B_mass"] > 3.91)) valid = false;//signal
        //if (valid && (vars["B_mass"] < 3.75 || vars["B_mass"] > 4.0)) valid = false;
        //if (valid && !((vars["B_mass"] > 3.6 && vars["B_mass"] < 3.66)
        //||(vars["B_mass"] > 3.72 && vars["B_mass"] < 3.83)||(vars["B_mass"] > 3.91 && vars["B_mass"] < 4.0))) valid = false;

        // -----------------------
        // Apply pre-cuts (access variables directly by name)
        if (valid && vars["B_chi2cl"] < 0.003) valid = false;
        if (valid && vars["B_Qvalueuj"] > 0.2) valid = false;
        // -----------------------

        if (valid) {
            new_tree->Fill();
            ++kept;
        }

        // Lightweight progress indicator every ~5%
        if ((i % (nEntries / 20 + 1)) == 0) {
            printf("Processed %lld / %lld (%.1f%%)\r",
                   i, nEntries, 100.0 * double(i) / double(nEntries));
            fflush(stdout);
        }
    }

    new_tree->Write();
    f_out->Close();
    f_in->Close();

    printf("Selection completed: original entries = %lld, kept entries = %lld\n", nEntries, kept);
}
