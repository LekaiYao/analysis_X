import ROOT
from array import array
import os

# keep global references to ROOT files so TTree doesn't get invalidated
fdat = None
fmc = None

def ks_on_sidebands(tree, reader, var_arrays, var_names, spec_arrays, thr,
                    mass_var="B_mass", sb_low=(3.6, 3.65), sb_high=(3.75, 3.8), nbins=100):
    href = ROOT.TH1F("h_ref", "ref mass", nbins, sb_low[0], sb_high[1])
    hcut = ROOT.TH1F("h_cut", "cut mass", nbins, sb_low[0], sb_high[1])

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

    ks_val = href.KolmogorovTest(hcut) if hcut.Integral() > 0 else -1.0

    return ks_val, sb_before, sb_after, href, hcut


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
    weights_xml = f"dataset_{run_label}/weights/TMVAClassification_BDTs.weights.xml"
    data_file = "/user/u/u25lekai/work/analysis_X/selection/test_root/sideband_PSI_test4.root"
    mc_file = "/user/u/u25lekai/work/analysis_X/selection/test_root/MC_PSI2S_test4.root"
    tree_name = "tree"

    # variables from training
    var_names = ["B_alpha", "B_norm_trk1Dxy", "B_norm_svpvDistance"]

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
    for eff in effs:
        idx = int((1.0 - eff) * len(scores))
        thr = scores[idx]
        ks_val, sb_before, sb_after, href, hcut = ks_on_sidebands(
            tdat, reader, var_arrays, var_names, spec_arrays, thr
        )
        print(f"[RESULT] eff={eff:.2f}, thr={thr:.3f}, KS(sidebands)={ks_val:.4f}, SB before={sb_before}, after={sb_after}")
        # save plot
        c = ROOT.TCanvas()
        href.SetLineColor(ROOT.kBlue)
        hcut.SetLineColor(ROOT.kRed)
        href.Draw("hist")
        hcut.Draw("hist same")
        c.SaveAs(f"dataset_{run_label}/plots/mass_SB_KS_eff{int(eff*100)}.png")


if __name__ == "__main__":
    main()
