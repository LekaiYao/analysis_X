#include <TFile.h>
#include <TTree.h>
#include <TBranch.h>
#include <TLeaf.h>
#include <TObjArray.h>
#include <vector>
#include <string>
#include <cstdio>
#include <cmath>

void remove_nan() {
    // Open the input ROOT file
    TFile *f_in = TFile::Open("/user/u/u25lekai/work/analysis_B/optimization/data_Bu.root");
    if (!f_in || f_in->IsZombie()) {
        fprintf(stderr, "Error: cannot open input file.\n");
        return;
    }
    TTree *tree = (TTree*)f_in->Get("tree");
    if (!tree) {
        fprintf(stderr, "Error: cannot find TTree 'tree' in input file.\n");
        f_in->Close();
        return;
    }

    // Collect all Float_t scalar leaves ONCE (avoid GetLeaf() inside the event loop)
    TObjArray *branches = tree->GetListOfBranches();
    std::vector<TLeaf*> float_leaves;
    float_leaves.reserve(branches->GetEntries());

    for (int i = 0; i < branches->GetEntries(); ++i) {
        TBranch *br = (TBranch*)branches->At(i);
        if (!br) continue;
        TLeaf *leaf = br->GetLeaf(br->GetName());
        if (leaf && strcmp(leaf->GetTypeName(), "Float_t") == 0 && leaf->GetLen() == 1) {
            float_leaves.push_back(leaf);
        }
    }

    // Create the output file and new tree structure
    TFile *f_out = TFile::Open("data_Bu_nonan.root", "RECREATE");
    if (!f_out || f_out->IsZombie()) {
        fprintf(stderr, "Error: cannot create output file.\n");
        f_in->Close();
        return;
    }
    TTree *new_tree = tree->CloneTree(0);  // clone structure only
    if (!new_tree) {
        fprintf(stderr, "Error: failed to clone tree structure.\n");
        f_out->Close();
        f_in->Close();
        return;
    }
    // Optional: tune autoflush to reduce I/O overhead on large outputs
    // new_tree->SetAutoFlush(100000);

    const Long64_t nEntries = tree->GetEntries();
    Long64_t kept = 0;

    // Loop over events
    for (Long64_t i = 0; i < nEntries; ++i) {
        tree->GetEntry(i);

        bool valid = true;
        // Fast path: iterate cached leaves and check finiteness
        for (TLeaf* leaf : float_leaves) {
            // GetValue() returns a Double_t; convert then check finite
            const double v = leaf->GetValue();
            if (!std::isfinite(v)) {
                valid = false;
                break;
            }
        }

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
    printf("\n");

    // Save
    new_tree->Write();
    f_out->Close();
    f_in->Close();

    printf("Selection completed: original entries = %lld, kept entries = %lld\n", nEntries, kept);
}
