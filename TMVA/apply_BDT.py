import ROOT
from ROOT import TMVA, TFile, TTree
from array import array

ROOT.gROOT.SetBatch(True)

# -----------------------------
# Configuration parameters
# -----------------------------
model_xml_path = "dataset_test4/weights/TMVAClassification_BDTs.weights.xml"
#input_data_path = "/user/u/u25lekai/work/analysis_X/selection/test_root/sideband_PSI_test4.root"
#output_root = "sideband_PSI_BDT_test4.root"
#input_data_path = "/user/u/u25lekai/work/analysis_X/selection/test_root/MC_PSI2S_test4.root"
#output_root = "MC_PSI2S_BDT_test4.root"
input_data_path = "/user/u/u25lekai/work/analysis_X/selection/test_root/DATA_XPSI_test4.root"
output_root = "DATA_XPSI_BDT_test4.root"
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
#B_chi2cl = array('f', [0.])
#B_Qvalueuj = array('f', [0.])
B_norm_trk1Dxy = array('f', [0.])
B_norm_svpvDistance = array('f', [0.])

reader.AddVariable("B_alpha", B_alpha)
#reader.AddVariable("B_chi2cl",B_chi2cl)
#reader.AddVariable("B_Qvalueuj",B_Qvalueuj)
reader.AddVariable("B_norm_trk1Dxy",B_norm_trk1Dxy)
reader.AddVariable("B_norm_svpvDistance",B_norm_svpvDistance)

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
    #B_chi2cl[0] = getattr(tree, "B_chi2cl")
    #B_Qvalueuj[0] = getattr(tree, "B_Qvalueuj")
    B_norm_trk1Dxy[0] = getattr(tree, "B_norm_trk1Dxy")
    B_norm_svpvDistance[0] = getattr(tree, "B_norm_svpvDistance")

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
