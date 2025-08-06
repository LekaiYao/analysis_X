import ROOT
from ROOT import TMVA, TFile, TTree
from array import array

ROOT.gROOT.SetBatch(True)

# -----------------------------
# Configuration parameters
# -----------------------------
model_xml_path = "dataset_test3/weights/TMVAClassification_BDTs.weights.xml"
#input_data_path = "/user/u/u25lekai/work/analysis_B/selection/X_ppRef/root_files/sideband_PSI.root"
#output_root = "sideband_PSI_BDT.root"
#input_data_path = "/user/u/u25lekai/work/analysis_B/selection/X_ppRef/root_files/MC_PSI2S.root"
#output_root = "MC_PSI2S_BDT.root"
input_data_path = "/user/u/u25lekai/work/analysis_B/selection/X_ppRef/root_files/DATA_XPSI_nonan.root"
output_root = "DATA_XPSI_BDT.root"
tree_name = "tree"

# -----------------------------
# Load trained model
# -----------------------------
TMVA.Tools.Instance()
TMVA.PyMethodBase.PyInitialize()

reader = TMVA.Reader("!Color:!Silent")

# -----------------------------
# Define input variables (must match training!)
# -----------------------------
B_alpha = array('f', [0.])
B_trk1dR = array('f', [0.])
B_trk2dR = array('f', [0.])
B_trk1Pt = array('f', [0.])
B_trk2Pt = array('f', [0.])
B_norm_svpvDistance_2D = array('f', [0.])
B_cos_dtheta = array('f', [0.])

reader.AddVariable("B_alpha", B_alpha)
reader.AddVariable("B_trk1dR", B_trk1dR)
reader.AddVariable("B_trk2dR", B_trk2dR)
reader.AddVariable("B_trk1Pt", B_trk1Pt)
reader.AddVariable("B_trk2Pt", B_trk2Pt)
reader.AddVariable("B_norm_svpvDistance_2D", B_norm_svpvDistance_2D)
reader.AddVariable("B_cos_dtheta", B_cos_dtheta)

reader.BookMVA("BDT", model_xml_path)

# -----------------------------
# Load input data
# -----------------------------
input_file = TFile.Open(input_data_path)
tree = input_file.Get(tree_name)

# -----------------------------
# Prepare output file and tree
# -----------------------------
output_file = TFile.Open(output_root, "RECREATE")
output_tree = TTree("tree", "Tree with B_mass and BDT_score")

# Variables to store in output
B_mass = array('f', [0.])
bdt_score = array('f', [0.])

output_tree.Branch("B_mass", B_mass, "B_mass/F")
output_tree.Branch("BDT_score", bdt_score, "BDT_score/F")

# -----------------------------
# Loop over events
# -----------------------------
nentries = tree.GetEntries()
print(f"Processing {nentries} entries...")

for i in range(nentries):
    tree.GetEntry(i)

    # Assign input variable values
    B_alpha[0] = getattr(tree, "B_alpha")
    B_trk1dR[0] = getattr(tree, "B_trk1dR")
    B_trk2dR[0] = getattr(tree, "B_trk2dR")
    B_trk1Pt[0] = getattr(tree, "B_trk1Pt")
    B_trk2Pt[0] = getattr(tree, "B_trk2Pt")
    B_norm_svpvDistance_2D[0] = getattr(tree, "B_norm_svpvDistance_2D")
    B_cos_dtheta[0] = getattr(tree, "B_cos_dtheta")

    # Evaluate BDT
    bdt_score[0] = reader.EvaluateMVA("BDT")

    # Get B_mass
    B_mass[0] = getattr(tree, "B_mass")

    # Fill output tree
    output_tree.Fill()

# -----------------------------
# Save output
# -----------------------------
output_file.cd()
output_tree.Write("", ROOT.TObject.kOverwrite)
output_file.Close()

print(f"Output ROOT file saved: {output_root}")
