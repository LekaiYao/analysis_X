void new_select() {
    // Open the input ROOT file
    TFile *f_in = TFile::Open("DATA_XPSI.root");
    TTree *tree = (TTree*)f_in->Get("tree");

    // Collect all Float_t branches
    TObjArray *branches = tree->GetListOfBranches();
    std::vector<TBranch*> float_branches;
    std::vector<float*> float_ptrs;

    // Pointer to the B_mass variable
    float *B_mass_ptr = nullptr;

    for (int i = 0; i < branches->GetEntries(); ++i) {
        TBranch *br = (TBranch*)branches->At(i);
        TLeaf *leaf = br->GetLeaf(br->GetName());
        if (leaf && strcmp(leaf->GetTypeName(), "Float_t") == 0 && leaf->GetLen() == 1) {
            float *ptr = new float;
            tree->SetBranchAddress(br->GetName(), ptr);
            float_branches.push_back(br);
            float_ptrs.push_back(ptr);

            // Check if this branch is B_mass
            if (strcmp(br->GetName(), "B_mass") == 0) {
                B_mass_ptr = ptr;
            }
        }
    }

    // Error if B_mass was not found
    if (!B_mass_ptr) {
        Error("select", "Branch 'B_mass' not found or not of Float_t type!");
        return;
    }

    // Create output ROOT file and clone tree structure (no entries yet)
    TFile *f_out = TFile::Open("new_TMVADATA_XPSI.root", "RECREATE");
    TTree *new_tree = tree->CloneTree(0);  // only clone structure

    Long64_t nEntries = tree->GetEntries();
    Long64_t kept = 0;

    for (Long64_t i = 0; i < nEntries; ++i) {
        tree->GetEntry(i);
        bool valid = true;

        // Check for NaN or Inf values in any Float_t branch
        for (float *val : float_ptrs) {
            if (!std::isfinite(*val)) {
                valid = false;
                break;
            }
        }

        // Apply B_mass window cut: must satisfy (B_mass is out of [3.656,3.718], [mean-5sigma,mean+5sigma])
        if (valid && (*B_mass_ptr > 3.656 && *B_mass_ptr <3.718) ) {
            valid = false;
        }
    
        if (valid && (*B_mass_ptr < 3.6 || *B_mass_ptr >3.8) ) {
            valid = false;
        }

        if (valid) {
            new_tree->Fill();  // Add the current entry to the output tree
            ++kept;
        }
    }

    new_tree->Write();
    f_out->Close();
    f_in->Close();

    // Optional: free memory if needed
    // for (float *ptr : float_ptrs) delete ptr;

    printf("Selection completed: original entries = %lld, kept entries = %lld\n", nEntries, kept);
}
