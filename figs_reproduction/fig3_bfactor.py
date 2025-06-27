import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
import numpy as np
import scipy.stats as st

sns.set_style("white")
sns.despine()

# Load and clean data
CORRUPT_TAGS = ["2nr1", "1niz", "2l4g", "1v4z", "1yyb"]
df = pd.read_csv("tag_scaffold_results_fix_with_bfactors.csv")
df = df[~df['tag_id'].isin(CORRUPT_TAGS)]

# Set analysis parameters
threshold = 1.0
msa_filter = True  # Using MSA filter for good tag selection

# Identify good quality tags
no_scaffold_df = df[df['scaffold_id'] == "no_scaffold"][['af3_rmsd', 'tag_id', 'num_msa_hits']].dropna()
good_tags = no_scaffold_df[(no_scaffold_df['af3_rmsd'] < threshold) & 
                           (no_scaffold_df['num_msa_hits'] > 2)]['tag_id'].unique()

# Prepare scaffold data for visualization
scaffold_data = df[
    df['scaffold_id'].isin(['no_scaffold'] + list(set(df['scaffold_id']) - {'no_scaffold'})) & 
    df['tag_id'].isin(good_tags)
][['tag_id', 'scaffold_id', 'terminus',
   'af2_mean_bfactor', 'af2_fix_mean_bfactor',
   'af3_mean_bfactor', 'af3_fix_mean_bfactor']]

# Add terminus for no_scaffold
scaffold_data.loc[scaffold_data['scaffold_id'] == 'no_scaffold', 'terminus'] = ''
scaffold_data['scaffold_id'] = scaffold_data.apply(
    lambda x: x['scaffold_id'] if x['scaffold_id'] == 'no_scaffold' 
    else x['scaffold_id'] + '-' + x['terminus'], axis=1
)

# Reshape data for plotting
plot_data = pd.melt(
    scaffold_data,
    id_vars=['tag_id', 'scaffold_id'],
    value_vars=['af2_mean_bfactor', 'af2_fix_mean_bfactor', 
                'af3_mean_bfactor', 'af3_fix_mean_bfactor'],
    var_name='metric',
    value_name='bfactor'
)
split_columns = plot_data['metric'].str.split('_', n=1, expand=True)
plot_data['algorithm'] = split_columns[0].str.upper()
plot_data['fix_status'] = split_columns[1].apply(
    lambda x: 'After Fix' if x.startswith('fix') else 'Before Fix'
)
plot_data['category'] = plot_data['algorithm'] + ' ' + plot_data['fix_status']

# Ensure the scaffold_id order is the following:
order_scaffold_id = ["no_scaffold", "GFP-N", "GST-N", "MBP-N", "SUMO-N", "GFP-C", "GST-C", "MBP-C", "SUMO-C"]
plot_data['scaffold_id'] = pd.Categorical(plot_data['scaffold_id'], categories=order_scaffold_id, ordered=True)

# Prepare the figure with increased spacing between scaffold groups.
plt.figure(figsize=(16, 5))

# Define marker style and palette for each category
markers = {
    'AF2 Before Fix': 'o',
    'AF2 After Fix': 'o',
    'AF3 Before Fix': 's',
    'AF3 After Fix': 's'
}
palette = {
    'AF2 Before Fix': '#aec7e8',
    'AF2 After Fix': '#1f77b4',
    'AF3 Before Fix': '#98df8a',
    'AF3 After Fix': '#2ca02c'
}

# Compute mean and 95% CI for each scaffold and category, dropping NaNs.
stats = []
for (scaffold, category), group in plot_data.groupby(['scaffold_id', 'category']):
    data = group['bfactor'].dropna()
    if len(data) == 0:
        continue  # Skip groups with no data
    mean = data.mean()
    # If only one data point, sem becomes NaN; set CI to zero in that case.
    sem = st.sem(data) if len(data) > 1 else 0.0
    ci = 1.96 * sem if sem==sem else 0.0  # Use 0 if sem is nan
    stats.append((scaffold, category, mean, ci))
stats_df = pd.DataFrame(stats, columns=['scaffold_id', 'category', 'mean', 'ci'])

# Multiply scaffold indices for better separation.
scaffold_to_x = {scaffold: idx*0.8 for idx, scaffold in enumerate(order_scaffold_id)}

# Define new offsets with tighter spacing:
category_offsets = {
    'AF2 Before Fix': -0.15,
    'AF2 After Fix': -0.05,
    'AF3 Before Fix':  0.05,
    'AF3 After Fix':   0.15
}
stats_df['x'] = stats_df.apply(
    lambda row: scaffold_to_x[row['scaffold_id']] + category_offsets[row['category']], axis=1
)

# Plot the mean values with error bars using plt.errorbar
for _, row in stats_df.iterrows():
    plt.errorbar(
        row['x'], row['mean'], yerr=row['ci'],
        fmt=markers[row['category']], color=palette[row['category']],
        capsize=4, markersize=10, linestyle='None', zorder=3
    )

# Draw connecting lines for the AF2 and AF3 pairs
for scaffold in order_scaffold_id:
    # Connect AF2 Before Fix to AF2 After Fix
    af2_before = stats_df[(stats_df['scaffold_id'] == scaffold) & (stats_df['category'] == 'AF2 Before Fix')]
    af2_after = stats_df[(stats_df['scaffold_id'] == scaffold) & (stats_df['category'] == 'AF2 After Fix')]
    if not af2_before.empty and not af2_after.empty:
        plt.plot(
            [af2_before['x'].values[0], af2_after['x'].values[0]],
            [af2_before['mean'].values[0], af2_after['mean'].values[0]],
            color=palette['AF2 After Fix'], linestyle='--', alpha=0.5, zorder=2
        )
    # Connect AF3 Before Fix to AF3 After Fix
    af3_before = stats_df[(stats_df['scaffold_id'] == scaffold) & (stats_df['category'] == 'AF3 Before Fix')]
    af3_after = stats_df[(stats_df['scaffold_id'] == scaffold) & (stats_df['category'] == 'AF3 After Fix')]
    if not af3_before.empty and not af3_after.empty:
        plt.plot(
            [af3_before['x'].values[0], af3_after['x'].values[0]],
            [af3_before['mean'].values[0], af3_after['mean'].values[0]],
            color=palette['AF3 After Fix'], linestyle='--', alpha=0.5, zorder=2
        )

# Add shaded regions and vertical separators for terminus groups.
# Adjust the x-range to account for multiplied scaffold indices.
plt.axvspan(scaffold_to_x[order_scaffold_id[0]]-0.6, scaffold_to_x[order_scaffold_id[4]]+0.4, 
            color='#f0f0f0', alpha=0.3, zorder=0)  # N-terminus region
plt.axvspan(scaffold_to_x[order_scaffold_id[4]]+0.4, scaffold_to_x[order_scaffold_id[-1]]+0.6, 
            color='#e0e0e0', alpha=0.3, zorder=0)  # C-terminus region
# Vertical separators for the scaffold groups
for x in scaffold_to_x.values():
    plt.axvline(x=x-0.4, color='gray', linestyle=':', linewidth=1, alpha=0.5)

# Get the actual data range for y-axis
y_max = stats_df['mean'].max() + stats_df['ci'].max()
y_min = max(0, stats_df['mean'].min() - stats_df['ci'].min())  # Don't go below 0
y_padding = (y_max - y_min) * 0.1  # Add 10% padding

# Get the actual x range needed
x_min = min(stats_df['x']) - 0.3  # Add some padding for first group
x_max = max(stats_df['x']) + 0.3  # Add some padding for last group

plt.xlim(x_min, x_max)
plt.ylim(20, 100)

# Set custom xtick labels at the scaffold center positions
xtick_positions = [scaffold_to_x[s] for s in order_scaffold_id]
plt.xticks(xtick_positions, ["No Scaffold"] + order_scaffold_id[1:], fontsize=24)

plt.yticks(fontsize=24)
plt.xlabel('')
plt.ylabel('pLDDT', fontsize=24)
plt.grid(axis='y', alpha=0.3)

# Build a custom legend
from matplotlib.lines import Line2D
legend_elements = [
    Line2D([0], [0], marker=markers['AF2 Before Fix'], color='w', label='AlphaFold-2',
           markerfacecolor=palette['AF2 Before Fix'], markersize=10),
    Line2D([0], [0], marker=markers['AF2 After Fix'], color='w', label='AlphaFold-2 + Windowed MSA',
           markerfacecolor=palette['AF2 After Fix'], markersize=10),
    Line2D([0], [0], marker=markers['AF3 Before Fix'], color='w', label='AlphaFold-3',
           markerfacecolor=palette['AF3 Before Fix'], markersize=10),
    Line2D([0], [0], marker=markers['AF3 After Fix'], color='w', label='AlphaFold-3 + Windowed MSA',
           markerfacecolor=palette['AF3 After Fix'], markersize=10)
]
plt.legend(handles=legend_elements, fontsize=22, loc='upper center', bbox_to_anchor=(0.5, 1.25), ncol=2)

plt.tight_layout()
plt.savefig('scaffold_wise_breakdown_bfactor.pdf', dpi=600, bbox_inches='tight')
plt.show()
