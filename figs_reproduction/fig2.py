import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from scipy.stats import ecdf
import seaborn as sns
sns.set_style("white")
sns.despine()


CORRUPT_TAGS = ["2nr1", "1niz", "2l4g", "1v4z", "1yyb"]
df = pd.read_csv("tag_scaffold_results_fix.csv")
# Remove corrupt tags
df = df[~df['tag_id'].isin(CORRUPT_TAGS)]
threshold = 1.0
msa_filter = True

esm_column = "esm_rmsd"

df_no_scaffold = df[df['scaffold_id'] == "no_scaffold"][[f"af2_rmsd", "af3_rmsd", "af2_fix_rmsd", "af3_fix_rmsd", "tag_id", "num_msa_hits", "scaffold_id", "terminus", esm_column]]
# df_no_scaffold = df_no_scaffold.dropna()
good_tags = df_no_scaffold[(df_no_scaffold[f"af3_rmsd"] < threshold) & (df_no_scaffold['num_msa_hits'] > 2)]['tag_id'].unique()
num_good_tags = len(good_tags)
df_no_scaffold = df_no_scaffold[df_no_scaffold['tag_id'].isin(good_tags)]
df_scaffold = df[df['scaffold_id'] != "no_scaffold"][[f"af2_rmsd", "af3_rmsd", "af2_fix_rmsd", "af3_fix_rmsd", "tag_id", "num_msa_hits", "scaffold_id", "terminus", esm_column]]
df_scaffold = df_scaffold[df_scaffold['tag_id'].isin(good_tags)]
print("Num chimeras under consideration: ", len(df_scaffold))

# Define color scheme
af2_colors = {
    'no_scaffold': '#d62728',  # Bright Red
    'scaffold': '#9ecae1',     # Medium Blue
    'fix': '#3182bd',          # Light Blue
}

af3_colors = {
    'no_scaffold': '#ff9898',  # Light Red/Pink
    'scaffold':  '#a1d99b',    # Medium Green
    'fix': '#31a354',          # Light Green
}

# Plot the ECDFs
plt.figure(figsize=(12, 6))

# Remove ESM3 plotting code and just keep AF2/AF3
for algorithm in ["af2", "af3"]:
    # Select colors based on algorithm
    if algorithm == 'af2':
        colors = af2_colors
    else:
        colors = af3_colors
    
    # Calculate ECDFs
    ecdf_no_scaffold = ecdf(df_no_scaffold[f"{algorithm}_rmsd"].dropna()).cdf
    ecdf_scaffold = ecdf(df_scaffold[f"{algorithm}_rmsd"].dropna()).cdf
    ecdf_fix = ecdf(df_scaffold[f"{algorithm}_fix_rmsd"].dropna()).cdf

    # Plot lines with specified colors
    plt.plot(ecdf_no_scaffold.quantiles, ecdf_no_scaffold.probabilities, 
             label=f'AlphaFold-{algorithm[-1]} (No Scaffold)', 
             color=colors['no_scaffold'], linewidth=3)
    plt.plot(ecdf_scaffold.quantiles, ecdf_scaffold.probabilities, 
             label=f'AlphaFold-{algorithm[-1]}', 
             color=colors['scaffold'], linewidth=3)
    plt.plot(ecdf_fix.quantiles, ecdf_fix.probabilities, 
             label=f'AlphaFold-{algorithm[-1]} + Windowed MSA', 
             color=colors['fix'], linewidth=3)

plt.xlim(0, 7)
plt.ylim(0, 1)
plt.xticks(fontsize=18)
plt.yticks(fontsize=18)
plt.xlabel('RMSD (Å)', fontsize=18)
plt.ylabel("Proportion of tags with RMSD < xÅ", fontsize=18)
plt.legend(fontsize=18, ncol=2, loc='lower right')
plt.grid(True, alpha=0.5)
plt.tight_layout()
plt.savefig(f'aggregate_before_vs_after_fix.pdf', dpi=600, bbox_inches='tight')
plt.close()

# Create separate plots for each scaffold + terminus combination
for scaffold_idx, scaffold in enumerate(df_scaffold['scaffold_id'].unique()):
    scaffold_data = df_scaffold[df_scaffold['scaffold_id'] == scaffold]
    for terminus in scaffold_data['terminus'].unique():
        terminus_data = scaffold_data[scaffold_data['terminus'] == terminus]
        
        plt.figure(figsize=(12, 6))
        
        # Plot no_scaffold baseline and other curves (ESM3 removed)
        for algorithm in ["af2", "af3"]:
            ecdf_no = ecdf(df_no_scaffold[f"{algorithm}_rmsd"].dropna()).cdf
            plt.plot(ecdf_no.quantiles, ecdf_no.probabilities,
                     label=f'AlphaFold-{algorithm[-1]} (No Scaffold)',
                     color=af2_colors['no_scaffold'] if algorithm == 'af2' else af3_colors['no_scaffold'],
                     linewidth=3)
        
            # Original scaffold predictions
            ecdf_scaf = ecdf(terminus_data[f"{algorithm}_rmsd"].dropna()).cdf
            plt.plot(ecdf_scaf.quantiles, ecdf_scaf.probabilities,
                     linewidth=3,
                     color=af2_colors['scaffold'] if algorithm == 'af2' else af3_colors['scaffold'],
                     label=f'AlphaFold-{algorithm[-1]}')
            
            # Fixed MSA predictions
            ecdf_fix = ecdf(terminus_data[f"{algorithm}_fix_rmsd"].dropna()).cdf
            plt.plot(ecdf_fix.quantiles, ecdf_fix.probabilities,
                     linewidth=3,
                     color=af2_colors['fix'] if algorithm == 'af2' else af3_colors['fix'],
                     label=f'AlphaFold-{algorithm[-1]} + Windowed MSA')

        plt.xlim(0, 7)
        plt.ylim(0, 1)
        plt.xticks(fontsize=18)
        plt.yticks(fontsize=18)
        plt.xlabel('RMSD (Å)', fontsize=18)
        plt.ylabel("Proportion of tags with RMSD < xÅ", fontsize=18)
        plt.title(f'Scaffold: {scaffold}, Terminus: {terminus}', fontsize=20)
        plt.legend(fontsize=14, ncol=2, loc='lower right')
        plt.grid(True, alpha=0.5)
        plt.tight_layout()
        plt.savefig(f'aggregate_rmsd_{scaffold}_{terminus}.pdf', dpi=600, bbox_inches='tight')
        plt.close()