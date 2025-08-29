#!/usr/bin/env python3

import uproot
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt

def apply_cuts_from_accsel(df, is_mc=True, use_sidebands=False):
    """Apply cuts from ACCSEL.h file to dataframe using boolean indexing"""
    df_filtered = df.copy()
    print(f"Starting with {len(df_filtered)} events")

    if is_mc:
        acc_filter = (
            (df_filtered['B_trk1Eta'].abs() < 2.4) & (df_filtered['B_trk2Eta'].abs() < 2.4) &
            (df_filtered['B_trk1Pt'] > 0.5) & (df_filtered['B_trk2Pt'] > 0.5) &
            (df_filtered['B_mu1isAcc'] == 1) & (df_filtered['B_mu2isAcc'] == 1)
        )
        df_filtered = df_filtered[acc_filter]
        print(f"After acceptance cuts: {len(df_filtered)} events")

    fid_filter = ((df_filtered['B_pt'] >= 5) & (df_filtered['B_y'].abs() <= 2.4))
    df_filtered = df_filtered[fid_filter]
    print(f"After fiducial cuts: {len(df_filtered)} events")

    if is_mc:
        sel_filter = (
            ((df_filtered['B_trk1PtErr'] / df_filtered['B_trk1Pt']) < 0.1) &
            ((df_filtered['B_trk2PtErr'] / df_filtered['B_trk2Pt']) < 0.1) &
            ((df_filtered['B_trk1nPixelLayer'] + df_filtered['B_trk1nStripLayer']) > 10) &
            ((df_filtered['B_trk2nPixelLayer'] + df_filtered['B_trk2nStripLayer']) > 10) &
            ((df_filtered['B_trk1Chi2ndf'] / (df_filtered['B_trk1nPixelLayer'] + df_filtered['B_trk1nStripLayer'])) < 0.18) &
            ((df_filtered['B_trk2Chi2ndf'] / (df_filtered['B_trk2nPixelLayer'] + df_filtered['B_trk2nStripLayer'])) < 0.18) &
            (df_filtered['B_trk1highPurity'] == 1) & (df_filtered['B_trk2highPurity'] == 1) &
            (df_filtered['B_mu1SoftMuID'] == 1) & (df_filtered['B_mu2SoftMuID'] == 1) &
            ((df_filtered['B_ujmass'] - 3.096916).abs() < 0.15) &
            (df_filtered['B_ujvProb'] > 0.01)
        )
        df_filtered = df_filtered[sel_filter]
        print(f"After selection cuts: {len(df_filtered)} events")

    trg_filter = ((df_filtered['B_mu1isTriggered'] == 1) & (df_filtered['B_mu2isTriggered'] == 1))
    df_filtered = df_filtered[trg_filter]
    print(f"After trigger matching: {len(df_filtered)} events")

    return df_filtered

# ==========================
# APPLY CUTS
# ==========================
ANALYSIS_TYPE = "MC_signal"  # Options: "MC_signal", "Data_sidebands"

if ANALYSIS_TYPE == "MC_signal":
    root_file = "/user/u/u25lekai/work/ppRef/analysis_X/selection/root_files/MC_X3872.root"
    is_mc = True
    use_sidebands = False
    output_suffix = "MC_signal"
elif ANALYSIS_TYPE == "Data_sidebands":
    root_file = "/user/u/u25lekai/work/ppRef/analysis_X/selection/root_files/sideband.root"
    is_mc = False
    use_sidebands = True
    output_suffix = "Data_sidebands"
else:
    raise ValueError(f"Unknown ANALYSIS_TYPE: {ANALYSIS_TYPE}")

tree_name = "tree"

if is_mc:
    variables = [
        "B_pt", "B_y", "B_mass", "B_alpha", "B_Qvalueuj", "B_Qvalue", "B_cos_dtheta", "B_trkPtimb",
        "B_chi2cl", "B_trk1dR", "B_trk2dR", "B_trk1Pt", "B_trk2Pt",
        "B_norm_svpvDistance_2D", "B_norm_svpvDistance", "B_norm_trk1Dxy", "B_norm_trk2Dxy",
        "B_gen", "B_trk1Eta", "B_trk2Eta", "B_mu1isAcc", "B_mu2isAcc",
        "B_trk1PtErr", "B_trk2PtErr", "B_trk1nPixelLayer", "B_trk1nStripLayer",
        "B_trk2nPixelLayer", "B_trk2nStripLayer", "B_trk1Chi2ndf", "B_trk2Chi2ndf",
        "B_trk1highPurity", "B_trk2highPurity", "B_mu1SoftMuID", "B_mu2SoftMuID",
        "B_ujmass", "B_ujvProb", "B_mu1isTriggered", "B_mu2isTriggered"
    ]
else:
    variables = [
        "B_pt", "B_y", "B_mass", "B_alpha", "B_Qvalueuj", "B_Qvalue", "B_cos_dtheta", "B_trkPtimb",
        "B_chi2cl", "B_trk1dR", "B_trk2dR", "B_trk1Pt", "B_trk2Pt",
        "B_norm_svpvDistance_2D", "B_norm_svpvDistance", "B_norm_trk1Dxy", "B_norm_trk2Dxy",
        "B_mu1isTriggered", "B_mu2isTriggered"
    ]

print(f"Analysis Type: {ANALYSIS_TYPE}")
print(f"Input file: {root_file}")
print(f"Is MC: {is_mc}, Use sidebands: {use_sidebands}")

# ==========================
# LOAD AND FILTER IN CHUNKS (no flatten; data assumed non-nested)
# ==========================
print(f"Opening {root_file} ...")
with uproot.open(root_file) as file:
    if tree_name not in file:
        raise ValueError(f"TTree '{tree_name}' not found in file. Available: {list(file.keys())}")
    tree = file[tree_name]
    n_entries = tree.num_entries
    print(f"Total entries: {n_entries}")
    chunk_size = 500000
    parts = []
    for start in range(0, n_entries, chunk_size):
        end = min(start + chunk_size, n_entries)
        print(f"Chunk {start//chunk_size+1}: entries {start} to {end-1}")

        # Read directly as NumPy arrays (no Awkward, no flatten)
        arrays = tree.arrays(variables, entry_start=start, entry_stop=end, library="np")
        data = {var: arrays[var] for var in variables}
        df_part = pd.DataFrame(data)

        df_part.dropna(inplace=True)

        if is_mc:
            sig_mask = (
                (df_part['B_gen'] == 23333) | (df_part['B_gen'] == 24333) |
                (df_part['B_gen'] == 23433) | (df_part['B_gen'] == 24433)
            )
            df_part = df_part[sig_mask]
        else:
            sb_mask = (
                (df_part['B_mass'] <= 3.65788) |
                ((df_part['B_mass'] >= 3.71698) & (df_part['B_mass'] <= 3.83576)) |
                (df_part['B_mass'] >= 3.91064)
            )
            df_part = df_part[sb_mask]

        df_part = apply_cuts_from_accsel(df_part, is_mc=is_mc, use_sidebands=use_sidebands)
        parts.append(df_part)

    df_filtered = pd.concat(parts, ignore_index=True) if parts else pd.DataFrame(columns=variables)
    print(f"Total retained after all cuts: {len(df_filtered)}")

# ==========================
# PREPARE CORRELATION DATA
# ==========================
correlation_vars = [
    "B_pt", "B_y", "B_mass", "B_alpha", "B_Qvalueuj", "B_Qvalue", "B_cos_dtheta",
    "B_trkPtimb", "B_chi2cl", "B_trk1dR", "B_trk2dR", "B_trk1Pt", "B_trk2Pt",
    "B_norm_svpvDistance_2D", "B_norm_svpvDistance", "B_norm_trk1Dxy", "B_norm_trk2Dxy"
]

if df_filtered.empty:
    print("No events retained after cuts. Exiting without plotting.")
else:
    df_corr = df_filtered[correlation_vars]
    corr_matrix = df_corr.corr(method="pearson")
    print("\nCorrelation matrix:")
    print(corr_matrix)

    plt.figure(figsize=(16, 12))
    sns.heatmap(corr_matrix, annot=True, cmap="coolwarm", vmin=-1, vmax=1,
                fmt='.3f', square=True, cbar_kws={"shrink": .8})

    plot_title = ""
    plt.title(plot_title, fontsize=16)
    plt.xticks(rotation=45, ha='right')
    plt.yticks(rotation=0)
    plt.tight_layout()

    output_file = f"matrix_pdf/correlation_matrix_{output_suffix}_X3872.pdf"
    plt.savefig(output_file, dpi=300, bbox_inches='tight', format='pdf')
    print(f"\nCorrelation matrix plot saved as: {output_file}")
    plt.close()
