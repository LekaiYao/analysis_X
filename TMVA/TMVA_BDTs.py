import ROOT
from ROOT import TMVA, TFile, TTree, TCut
import os.path
import sys


##run as 'python3 filename.py command_input'
##command_input is a name/label for the train, like train1,test2, etc.

#Get the input argument from the command line
if len(sys.argv) != 2:
    print("Insert a number after the python script")
    sys.exit(1)

input = sys.argv[1]


TMVA.Tools.Instance()               #need to run this two to load up TMVA
TMVA.PyMethodBase.PyInitialize()    #in PyROOT

#input PATH
dir = "/user/u/u25lekai/work/analysis_X/selection/test_root/"
#data_file = "dataSideband_Bu.root"
data_file = "sideband_PSI_test4.root"
mc_file = "MC_PSI2S_test4.root"
 
# Open the ROOT files and access the TTree for data and MC
data = TFile.Open(dir + data_file)
mc = TFile.Open(dir + mc_file)
background = data.Get("tree")
signal = mc.Get("tree")

#directories where the results will be stored
if not os.path.exists("dataset_" + input + "/results/rootfiles"):
    os.makedirs("dataset_" + input + "/results/rootfiles")
if not os.path.exists("dataset_" + input + "/weights"):
    os.makedirs("dataset_" + input + "/weights")
if not os.path.exists("dataset_" + input + "/plots"):
    os.makedirs("dataset_" + input + "/plots")

# Create a ROOT output file where TMVA will store ntuples, histograms, correlationMatrix, etc
outfname='dataset_' + input + '/results/rootfiles/TMVA_BDT.root' 
output = TFile.Open(outfname, 'RECREATE')


cuts=""#MC
cutb=""#data

mycutS=TCut(cuts)
mycutB=TCut(cutb)


factory = TMVA.Factory('TMVAClassification', output, '!V:!Silent:Color:DrawProgressBar:Transformations=I;D;P;G,D:AnalysisType=Classification')

#loads the outputs to the dir 'dataset'
dataloader = TMVA.DataLoader('dataset_' + input)

# features to train the BDT
# you could add the variables you want here
dataloader.AddVariable("B_alpha")
dataloader.AddVariable("B_chi2cl")
dataloader.AddVariable("B_Qvalueuj")
dataloader.AddVariable("B_norm_trk1Dxy")
dataloader.AddVariable("B_norm_svpvDistance")


signalWeight     = 1.0        #MC/signal
backgroundWeight = 1.0        #dataset sidebands/background

dataloader.AddSignalTree( signal, signalWeight )
dataloader.AddBackgroundTree( background, backgroundWeight )

sigCutEvents = signal.GetEntries(cuts)
bkgCutEvents = background.GetEntries(cutb)

#train 70% of the events, test 30%
sigTrain = int(sigCutEvents * 0.7)   
bkgTrain = int(bkgCutEvents * 0.7) 


dataloader.PrepareTrainingAndTestTree( mycutS, mycutB, "nTrain_Signal=%i:nTrain_Background=%i:SplitMode=Random:NormMode=NumEvents:!V" % (sigTrain, bkgTrain) )


# Book methods
#factory.BookMethod(dataloader, TMVA.Types.kBDT, "BDTs",
#                "!H:!V:NTrees=250:MinNodeSize=2.5%:MaxDepth=5:BoostType=AdaBoost:VarTransform=Decorrelate:SeparationType=GiniIndex:nCuts=30")
factory.BookMethod(dataloader, "BDT", "BDTs",
                "!H:!V:NTrees=250:MinNodeSize=2.5%:MaxDepth=5:BoostType=AdaBoost:VarTransform=Decorrelate:SeparationType=GiniIndex:nCuts=30")


# Run training,BDT test and evaluation
factory.TrainAllMethods()
factory.TestAllMethods()
factory.EvaluateAllMethods() #add full dataset


# Plot ROC Curves AND OTHERS
roc = factory.GetROCCurve(dataloader)
roc.SaveAs('dataset_' + input + '/plots/ROC_ClassificationBDT.png')


#close the output file
output.Close()

#open the GUI interface
#TMVA.TMVAGui(outfname)