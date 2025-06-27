import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
sns.set_style("white")
sns.despine()
from scipy.stats import ecdf

df = pd.read_csv("tag_scaffold_results_fix.csv")
# Filter out tags with no_scaffold
df = df[df['scaffold_id'] == "no_scaffold"][["af2_rmsd", "af3_rmsd", "esm_argmax_rmsd"]]
# Drop nan values
df = df.dropna()

# Create a new dataframe in long format for plotting
df_long = pd.melt(df, value_vars=['af3_rmsd', 'af2_rmsd', 'esm_argmax_rmsd'], 
                  var_name='Model', value_name='RMSD')

# Use seaborn to plot violin plots
plt.figure(figsize=(10, 5))
sns.violinplot(data=df_long, x='Model', y='RMSD', bw_adjust=0.25, inner=None,
               palette={'af2_rmsd': 'red', 'af3_rmsd': 'blue', 'esm_argmax_rmsd': 'purple'}, alpha=0.3)
sns.pointplot(x='Model', y='RMSD', data=df_long, color='black', markers='x', scale=1.0, errorbar=None, linestyle="")
# Customize the plot
plt.xticks(['af3_rmsd', 'af2_rmsd', 'esm_argmax_rmsd'], 
           ['AlphaFold-3', 'AlphaFold-2', 'ESM-3'],
           fontsize=18)
plt.yticks(fontsize=18)
# Remove xlabel
plt.xlabel('')
plt.ylabel('RMSD (Å)', fontsize=18)
plt.grid(True)
plt.ylim(0, 10)
plt.tight_layout()
plt.savefig('tag-prediction-violin.pdf', dpi=600, bbox_inches='tight')
plt.close()

# Calculate ECDF for both AF2 and AF3
af2_ecdf = ecdf(df['af2_rmsd']).cdf
af3_ecdf = ecdf(df['af3_rmsd']).cdf
esm_argmax_ecdf = ecdf(df['esm_argmax_rmsd']).cdf

# Create the CDF plot
plt.figure(figsize=(10, 5))
plt.plot(af3_ecdf.quantiles, af3_ecdf.probabilities, label='AlphaFold-3', color='blue', linewidth=3)
plt.plot(af2_ecdf.quantiles, af2_ecdf.probabilities, label='AlphaFold-2', color='red', linewidth=3)
# plt.plot(esm_ecdf.quantiles, esm_ecdf.probabilities, label='ESM', color='green', linewidth=3)
plt.plot(esm_argmax_ecdf.quantiles, esm_argmax_ecdf.probabilities, label='ESM-3', color='purple', linewidth=3)
# Draw a line at 1.5A
plt.axvline(x=1.0, color='gray', linestyle='--', alpha=0.9)
plt.xlim(0, 7)
plt.ylim(0, 1)
plt.xticks(fontsize=18)
plt.yticks(fontsize=18)
plt.xlabel('RMSD (Å)', fontsize=18)
plt.ylabel("Proportion of tags with RMSD < xÅ", fontsize=18)
plt.legend(fontsize=18)
plt.grid(True, alpha=0.3)
plt.tight_layout()
plt.savefig('tag-prediction-cdf.pdf', dpi=600, bbox_inches='tight')
plt.close()

