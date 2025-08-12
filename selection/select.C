void select() {
    // Open the input ROOT file
    TFile *f_in = TFile::Open("root_files/MC_PSI2S.root");
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
    TFile *f_out = TFile::Open("MC_PSI2S_test2.root", "RECREATE");
    TTree *new_tree = tree->CloneTree(0);

    Long64_t nEntries = tree->GetEntries();
    Long64_t kept = 0;

    for (Long64_t i = 0; i < nEntries; ++i) {
        tree->GetEntry(i);
        bool valid = true;

        // B_mass window cut
        //if (valid && (vars["B_mass"] > 3.65 && vars["B_mass"] < 3.72)) valid = false;
        if (valid && (vars["B_mass"] < 3.6 || vars["B_mass"] > 3.8)) valid = false;

        // -----------------------
        // Apply pre-cuts (access variables directly by name)
        // Example:
        if (valid && vars["B_chi2cl"] < 0.003) valid = false;
        //if (valid && fabs(vars["B_eta"]) > 2.4) valid = false;
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
