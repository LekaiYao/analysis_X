import ROOT
from ROOT import TMVA, TFile, TTree, TCut
import os.path
import sys


## run as 'python3 filename.py command_input'
## command_input is a name/label for the train, like train1,test2, etc.

# Get the input argument from the command line
if len(sys.argv) != 2:
    print("Insert a number after the python script")
    sys.exit(1)

input = sys.argv[1]

TMVA.Tools.Instance()               # need to run this two to load up TMVA
TMVA.PyMethodBase.PyInitialize()    # in PyROOT

# input PATH
dir = "/user/l/lekai/work/ppRef/analysis_X/selection/test_root/"
data_file = "sideband_trainX1.root"
mc_file = "MC_X3872_trainX1.root"
 
# Open the ROOT files and access the TTree for data and MC
data = TFile.Open(dir + data_file)
mc = TFile.Open(dir + mc_file)
background = data.Get("tree")
signal = mc.Get("tree")

# directories where the results will be stored
if not os.path.exists("dataset/dataset_" + input + "/results/rootfiles"):
    os.makedirs("dataset/dataset_" + input + "/results/rootfiles")
if not os.path.exists("dataset/dataset_" + input + "/weights"):
    os.makedirs("dataset/dataset_" + input + "/weights")
if not os.path.exists("dataset/dataset_" + input + "/plots"):
    os.makedirs("dataset/dataset_" + input + "/plots")

# Create a ROOT output file where TMVA will store ntuples, histograms, correlationMatrix, etc
outfname='dataset/dataset_' + input + '/results/rootfiles/TMVA_BDT.root' 
output = TFile.Open(outfname, 'RECREATE')

cuts="" # MC
cutb="" # data

mycutS=TCut(cuts)
mycutB=TCut(cutb)

factory = TMVA.Factory('TMVAClassification', output, '!V:!Silent:Color:DrawProgressBar:Transformations=I;D;P;G,D:AnalysisType=Classification')

# loads the outputs to the dir 'dataset'
dataloader = TMVA.DataLoader('dataset/dataset_' + input)

# features to train the BDT (keep as your original)
dataloader.AddVariable("B_trk1dR")
#dataloader.AddVariable("B_chi2cl")
#dataloader.AddVariable("B_cos_dtheta")
#dataloader.AddVariable("B_norm_trk1Dxy")
#dataloader.AddVariable("B_norm_svpvDistance")
dataloader.AddVariable("B_Qvalueuj")


# NEW: add mass as spectator to monitor mass sculpting
dataloader.AddSpectator("B_mass")

signalWeight     = 1.0        # MC/signal
backgroundWeight = 1.0        # dataset sidebands/background

dataloader.AddSignalTree( signal, signalWeight )
dataloader.AddBackgroundTree( background, backgroundWeight )

# Count available events
sigCutEvents = signal.GetEntries(cuts)
bkgCutEvents = background.GetEntries(cutb)

# Train/test split (keep 70% idea but cap training sizes to reduce overfitting/imbalance)
# Caps chosen for your sample sizes: ~10000 signal, ~10000 background for training if available
sigTrainCap = min(10000, int(sigCutEvents * 0.7))
bkgTrainCap = min(10000, int(bkgCutEvents * 0.7))

# Make sure they are not zero and do not exceed totals
sigTrain = max(1, min(sigTrainCap, sigCutEvents))
bkgTrain = max(1, min(bkgTrainCap, bkgCutEvents))

# Prepare training and test trees (explicit training counts; keep NormMode=NumEvents)
dataloader.PrepareTrainingAndTestTree(
    mycutS, mycutB,
    "nTrain_Signal=%i:nTrain_Background=%i:SplitMode=Random:SplitSeed=12345:NormMode=NumEvents:!V" % (sigTrain, bkgTrain)
)


# Book methods
# (keep your BDT method string but with safer hyperparameters and bagging; retain VarTransform=Decorrelate as in your original)
factory.BookMethod(
    dataloader, "BDT", "BDTs",
    "!H:!V:"
    "NTrees=500:"
    "MinNodeSize=5%:"
    "MaxDepth=3:"
    "BoostType=AdaBoost:"
    "VarTransform=Decorrelate:"
    "SeparationType=GiniIndex:"
    "nCuts=30:"
    "UseBaggedBoost=True:"
    "BaggedSampleFraction=0.8"
)

# Run training,BDT test and evaluation
factory.TrainAllMethods()
factory.TestAllMethods()
factory.EvaluateAllMethods() # add full dataset

# Plot ROC Curves AND OTHERS
roc = factory.GetROCCurve(dataloader)
if roc:
    roc.SaveAs('dataset/dataset_' + input + '/plots/ROC_ClassificationBDT.png')

# close the output file
output.Close()

# open the GUI interface
# TMVA.TMVAGui(outfname)
