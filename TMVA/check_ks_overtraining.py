#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import sys
import os
import re
import ROOT

"""
Purpose:
  - Read TMVA output ROOT file
  - Locate BDT method directory (e.g. Method_BDT/BDTs)
  - Fetch train/test histograms for signal and background
  - Compute KS p-values (train vs test)
  - Draw overlays with p-values annotated (PNG)
How to run:
  python3 check_ks_overtraining.py /path/to/TMVA_BDT.root
"""

# ---------- helpers ----------

def guess_tmva_path():
    # fallback default if not provided; you can edit this
    return "dataset/dataset_yourlabel/results/rootfiles/TMVA_BDT.root"

def is_hist1d(obj):
    return isinstance(obj, ROOT.TH1) and not isinstance(obj, ROOT.TH2)

def find_method_dir(f, method_regex=r".*Method_BDT.*", classname_regex=r".*BDT.*"):
    """
    Find the first directory that looks like a TMVA BDT method folder,
    typically something like: "dataset/dataset_xxx/Method_BDT/BDTs"
    """
    def walk_dir(d):
        out = []
        for key in d.GetListOfKeys():
            obj = key.ReadObj()
            if isinstance(obj, ROOT.TDirectoryFile):
                out.append(obj)
                out.extend(walk_dir(obj))
        return out

    all_dirs = [obj for obj in walk_dir(f)]
    cand = []
    for d in all_dirs:
        full = d.GetPath()  # e.g. file.root:/dataset/dataset_x/Method_BDT/BDTs
        if re.search(method_regex, full) and re.search(classname_regex, full):
            cand.append(d)
    # pick the deepest path (most specific)
    cand.sort(key=lambda dd: dd.GetPath().count('/'), reverse=True)
    return cand[0] if cand else None

def get_hist_by_patterns(directory, patterns):
    """
    Try several name patterns to robustly fetch a histogram.
    `patterns` is a list of regex strings.
    """
    for key in directory.GetListOfKeys():
        obj = key.ReadObj()
        if is_hist1d(obj):
            name = obj.GetName()
            for pat in patterns:
                if re.search(pat, name):
                    return obj
    return None

def ks_and_flag(htrain, htest):
    """
    Compute Kolmogorov test p-value (ROOT normalizes by default).
    Returns (p, flag_overtrain)
    """
    if not htrain or not htest:
        return None, None
    # Ensure histograms have entries
    if htrain.GetEffectiveEntries() <= 0 or htest.GetEffectiveEntries() <= 0:
        return None, None
    # KS p-value
    p = htrain.KolmogorovTest(htest, "")  # default normalize
    flag = (p < 0.05)
    return p, flag

def ensure_dir(path):
    if not os.path.exists(path):
        os.makedirs(path, exist_ok=True)

def annotate_and_save_overlay(h1, h2, title, pvalue, outpath):
    c = ROOT.TCanvas("c_"+re.sub(r'\W+', '_', title), title, 900, 700)
    c.SetGrid()
    # draw with common style
    h1 = h1.Clone(h1.GetName()+"_cpy")
    h2 = h2.Clone(h2.GetName()+"_cpy")
    h1.SetStats(0)
    h2.SetStats(0)

    # normalize to unit area for visual comparison (KS also does normalization by default)
    if h1.Integral() > 0:
        h1.Scale(1.0 / h1.Integral())
    if h2.Integral() > 0:
        h2.Scale(1.0 / h2.Integral())

    h1.SetLineWidth(2)
    h2.SetLineWidth(2)
    h2.SetLineStyle(7)  # dashed

    # Y range
    ymax = max(h1.GetMaximum(), h2.GetMaximum()) * 1.25
    h1.SetMaximum(ymax)

    h1.SetTitle(title)
    h1.GetXaxis().SetTitle("BDT response")
    h1.GetYaxis().SetTitle("Normalized entries")

    h1.Draw("hist")
    h2.Draw("hist same")

    leg = ROOT.TLegend(0.15, 0.78, 0.55, 0.90)
    leg.SetBorderSize(0)
    leg.SetFillStyle(0)
    leg.AddEntry(h1, "Train", "l")
    leg.AddEntry(h2, "Test", "l")
    leg.Draw()

    pave = ROOT.TPaveText(0.57, 0.78, 0.90, 0.90, "NDC")
    pave.SetFillColor(0)
    pave.SetBorderSize(0)
    pave.AddText("KS p-value = {:.4f}".format(pvalue if pvalue is not None else float('nan')))
    if pvalue is not None and pvalue < 0.05:
        pave.AddText("#color[2]{Possible overtraining}")
    pave.Draw()

    c.SaveAs(outpath)
    c.Close()

# ---------- main ----------

def main():
    # parse input
    if len(sys.argv) >= 2:
        tmva_file = sys.argv[1]
    else:
        tmva_file = guess_tmva_path()

    if not os.path.isfile(tmva_file):
        print(f"[ERROR] TMVA file not found: {tmva_file}")
        sys.exit(1)

    f = ROOT.TFile.Open(tmva_file)
    if not f or f.IsZombie():
        print("[ERROR] Cannot open ROOT file.")
        sys.exit(1)

    method_dir = find_method_dir(f)
    if not method_dir:
        print("[ERROR] Cannot locate BDT method directory (e.g. Method_BDT/BDTs).")
        sys.exit(1)

    print(f"[INFO] Using method directory: {method_dir.GetPath()}")

    # Robust patterns to match TMVA histogram names across versions/configs
    # Typical names: MVA_BDTs_Train_S, MVA_BDTs_Test_S, MVA_BDTs_Train_B, MVA_BDTs_Test_B
    pat_train_S = [r"Train.*_S$", r".*Train.*Signal.*", r".*_Train_S.*"]
    pat_test_S  = [r"Test.*_S$",  r".*Test.*Signal.*",  r".*_Test_S.*"]
    pat_train_B = [r"Train.*_B$", r".*Train.*Back.*",   r".*_Train_B.*"]
    pat_test_B  = [r"Test.*_B$",  r".*Test.*Back.*",    r".*_Test_B.*"]

    h_train_S = get_hist_by_patterns(method_dir, pat_train_S)
    h_test_S  = get_hist_by_patterns(method_dir, pat_test_S)
    h_train_B = get_hist_by_patterns(method_dir, pat_train_B)
    h_test_B  = get_hist_by_patterns(method_dir, pat_test_B)

    # Optional: try fallback by exact common names
    if not h_train_S: h_train_S = method_dir.Get("MVA_BDTs_Train_S")
    if not h_test_S:  h_test_S  = method_dir.Get("MVA_BDTs_Test_S")
    if not h_train_B: h_train_B = method_dir.Get("MVA_BDTs_Train_B")
    if not h_test_B:  h_test_B  = method_dir.Get("MVA_BDTs_Test_B")

    # Compute KS p-values
    pS, flagS = ks_and_flag(h_train_S, h_test_S)
    pB, flagB = ks_and_flag(h_train_B, h_test_B)

    # Print summary
    print("\n========== KS (Train vs Test) ==========")
    print(f"Signal    p-value: {pS if pS is not None else 'N/A'}  "
          f"{'(possible overtraining)' if (pS is not None and flagS) else ''}")
    print(f"Background p-value: {pB if pB is not None else 'N/A'}  "
          f"{'(possible overtraining)' if (pB is not None and flagB) else '''}")
    print("Threshold hint: p < 0.05 may indicate overtraining.\n")

    # Decide output plot directory near the TMVA file (…/plots)
    # For a path like dataset/dataset_x/results/rootfiles/TMVA_BDT.root
    # we go up two levels and into "plots"
    rootfiles_dir = os.path.dirname(os.path.abspath(tmva_file))
    plots_dir = os.path.normpath(os.path.join(rootfiles_dir, "..", "..", "plots"))
    ensure_dir(plots_dir)

    # Draw overlays
    if h_train_S and h_test_S:
        annotate_and_save_overlay(
            h_train_S, h_test_S,
            title="BDT Response (Signal) - Train vs Test",
            pvalue=(pS if pS is not None else float('nan')),
            outpath=os.path.join(plots_dir, "overtraining_BDT_S.png"),
        )

    if h_train_B and h_test_B:
        annotate_and_save_overlay(
            h_train_B, h_test_B,
            title="BDT Response (Background) - Train vs Test",
            pvalue=(pB if pB is not None else float('nan')),
            outpath=os.path.join(plots_dir, "overtraining_BDT_B.png"),
        )

    f.Close()

if __name__ == "__main__":
    ROOT.gROOT.SetBatch(True)
    main()
