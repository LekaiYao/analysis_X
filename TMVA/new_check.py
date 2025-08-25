import ROOT
from array import array
import os

# keep global references to ROOT files so TTree doesn't get invalidated
fdat = None
fmc = None

def ks_on_sidebands(tree, reader, var_arrays, var_names, spec_arrays, thr,
                    mass_var="B_mass", sb_low=(3.6, 3.65), sb_high=(3.72, 3.8), nbins=65, tag=""):
    """
    Run KS test on sidebands. 'tag' is used to give unique histogram names.
    """
    href = ROOT.TH1F(f"h_ref_{tag}", "ref mass", nbins, sb_low[0], sb_high[1])
    hcut = ROOT.TH1F(f"h_cut_{tag}", "cut mass", nbins, sb_low[0], sb_high[1])

    sb_before = 0
    sb_after = 0

    for i in range(tree.GetEntries()):
        tree.GetEntry(i)
        mass_val = getattr(tree, mass_var)
        # check in sideband
        if sb_low[0] <= mass_val <= sb_low[1] or sb_high[0] <= mass_val <= sb_high[1]:
            sb_before += 1
            for v in var_names:
                var_arrays[v][0] = float(getattr(tree, v))
            for s in spec_arrays:
                spec_arrays[s][0] = float(getattr(tree, s))
            score = reader.EvaluateMVA("BDT")
            href.Fill(mass_val)
            if score > thr:
                sb_after += 1
                hcut.Fill(mass_val)

    # normalize histos if nonzero
    if href.Integral() > 0:
        href.Scale(1.0 / href.Integral())
    if hcut.Integral() > 0:
        hcut.Scale(1.0 / hcut.Integral())

    ks_pval = href.KolmogorovTest(hcut) if hcut.Integral() > 0 else -1.0
    # CHANGED: use "M" to get the maximum KS distance D (not pseudo-experiment p-value)
    ks_D    = href.KolmogorovTest(hcut, "M") if hcut.Integral() > 0 else -1.0

    return ks_pval, ks_D, sb_before, sb_after, href, hcut


def build_reader(var_names, xml_path):
    reader = ROOT.TMVA.Reader("!Color:!Silent")
    var_arrays = {}
    for v in var_names:
        var_arrays[v] = array('f', [0.])
        reader.AddVariable(v, var_arrays[v])
    # spectator
    spec_arrays = {}
    spec_name = "B_mass"
    spec_arrays[spec_name] = array('f', [0.])
    reader.AddSpectator(spec_name, spec_arrays[spec_name])
    reader.BookMVA("BDT", xml_path)
    return reader, var_arrays, spec_arrays


def load_tree(file_path, tree_name, run_label):
    """Try to load tree from root file, check both top level and dataset_xxx/"""
    global fdat, fmc
    f = ROOT.TFile.Open(file_path)
    if not f or f.IsZombie():
        raise RuntimeError(f"Cannot open file: {file_path}")

    tree = f.Get(tree_name)
    if not tree:
        tree = f.Get(f"dataset_{run_label}/{tree_name}")

    if not tree:
        raise RuntimeError(f"Cannot find tree '{tree_name}' in {file_path}")

    print(f"[INFO] Loaded {tree_name} from {file_path}, entries = {tree.GetEntries()}")

    # keep file alive
    if "sideband" in file_path:
        fdat = f
    else:
        fmc = f

    return tree


def main():
    run_label = "v3test1"
    weights_xml = f"dataset/dataset_{run_label}/weights/TMVAClassification_BDTs.weights.xml"
    data_file = f"/user/u/u25lekai/work/ppRef/analysis_X/selection/test_root/sideband_PSI_{run_label}.root"
    mc_file = f"/user/u/u25lekai/work/ppRef/analysis_X/selection/test_root/MC_PSI2S_{run_label}.root"
    tree_name = "tree"

    # variables from training
    var_names = ["B_trk1dR", "B_trk2dR", "B_trk1Pt", "B_trk2Pt", "B_trkPtimb"]

    # build reader
    reader, var_arrays, spec_arrays = build_reader(var_names, weights_xml)

    # load trees
    tdat = load_tree(data_file, tree_name, run_label)
    tmc = load_tree(mc_file, tree_name, run_label)

    # get MC BDT scores for thresholds
    scores = []
    for i in range(tmc.GetEntries()):
        tmc.GetEntry(i)
        for v in var_names:
            var_arrays[v][0] = float(getattr(tmc, v))
        for s in spec_arrays:
            spec_arrays[s][0] = float(getattr(tmc, s))
        scores.append(reader.EvaluateMVA("BDT"))
    scores.sort()

    effs = [0.90, 0.70, 0.50, 0.30]
    results = []  # store results for later summary

    for eff in effs:
        idx = int((1.0 - eff) * len(scores))
        thr = scores[idx]
        tag = f"eff{int(eff*100)}"  # unique tag for hist names
        ks_pval, ks_D, sb_before, sb_after, href, hcut = ks_on_sidebands(
            tdat, reader, var_arrays, var_names, spec_arrays, thr, tag=tag
        )

        # save results for later
        results.append((eff, thr, ks_pval, ks_D, sb_before, sb_after))

        # save plot
        c = ROOT.TCanvas()
        href.SetLineColor(ROOT.kBlue)
        hcut.SetLineColor(ROOT.kRed)
        href.Draw("hist")
        hcut.Draw("hist same")
        c.SaveAs(f"dataset/dataset_{run_label}/plots/mass_SB_KS_eff{int(eff*100)}.png")
        c.Close()

    # ---- unified summary output ----
    print("\n======= Summary of KS test results =======")
    print(f"{'Eff':>6} {'Thr':>8} {'KS_pval':>12} {'KS_D':>12} {'SB_before':>12} {'SB_after':>12}")
    for eff, thr, ks_pval, ks_D, sb_before, sb_after in results:
        # CHANGED: scientific notation for p-value; more precision for KS_D
        print(f"{eff:6.2f} {thr:8.3f} {ks_pval:12.3e} {ks_D:12.6f} {sb_before:12d} {sb_after:12d}")
    print("==========================================")



if __name__ == "__main__":
    main()
