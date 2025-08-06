```
This is used to analysis X(3872) in ppRef and PbPb

Function Introduction:
1. selection
   (1) analysis.C: Applies basic selection criteria and flattens TTree
   (2) new_select.C: Removes NaN/INF values and applies further selection criteria to obtain sideband
2. scan
   scan.sh, input.txt, plot_template.C: When range is uncertain, sets range by reading max/min values. Plots distributions of variables of interest for data and MC on the same figure for comparison.
3. fit
   fit_MC_PSI2S.C: Fits MC to confirm mass sideband
   fit_DATA_fsfb.C: Fits data to confirm fs, fb
   fit_data.C: Fits data after applying all cuts
4. TMVA
   Performs BDT on MC and sideband data
   apply_BDT.py: Computes BDT_score for each event in data and MC, retains mass and score
5. optimization
   (1) Optimizes BDT_score for MC and sideband data to find optimal cut
   (2) Run.sh, input.txt, optimization_template.C: Performs single-variable optimization for variables of interest

TMVA Workflow:
    |->scan
selection (MC & data) -> fit MC -> selection (sideband) -> TMVA (train and score) -> fit data fsfb (rough score cut) -> optimization -> fit data
```
