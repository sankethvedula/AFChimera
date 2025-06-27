import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
from matplotlib.lines import Line2D
sns.set_style("white")
sns.despine()

# Load and clean data
CORRUPT_TAGS = ["2nr1", "1niz", "2l4g", "1v4z", "1yyb"]
df = pd.read_csv("tag_scaffold_results_fix.csv")
df = df[~df['tag_id'].isin(CORRUPT_TAGS)]

# Set analysis parameters
threshold = 1.0
msa_filter = True  # Using MSA filter for good tag selection

# Identify good quality tags
no_scaffold_df = df[df['scaffold_id'] == "no_scaffold"][['af3_rmsd', 'tag_id', 'num_msa_hits']].dropna()
good_tags = no_scaffold_df[(no_scaffold_df['af3_rmsd'] < threshold) & 
                          (no_scaffold_df['num_msa_hits'] > 2)]['tag_id'].unique()

# Prepare scaffold data for visualization
scaffold_data = df[df['scaffold_id'].isin(['no_scaffold'] + list(set(df['scaffold_id']) - {'no_scaffold'})) & 
                  df['tag_id'].isin(good_tags)][['tag_id', 'scaffold_id', 'terminus',
                                               'af2_rmsd', 'af2_fix_rmsd',
                                               'af3_rmsd', 'af3_fix_rmsd']]
# Add terminus for no_scaffold
scaffold_data.loc[scaffold_data['scaffold_id'] == 'no_scaffold', 'terminus'] = ''
scaffold_data['scaffold_id'] = scaffold_data.apply(lambda x: x['scaffold_id'] if x['scaffold_id'] == 'no_scaffold' 
                                                 else x['scaffold_id'] + '-' + x['terminus'], axis=1)


scaffold_data['delta_rmsd'] = scaffold_data['af3_fix_rmsd'] - scaffold_data['af3_rmsd']

# Reshape data for plotting
plot_data = pd.melt(
    scaffold_data,
    id_vars=['tag_id', 'scaffold_id'],
    value_vars=['af2_rmsd', 'af2_fix_rmsd', 'af3_rmsd', 'af3_fix_rmsd'],
    var_name='metric',
    value_name='rmsd'
)
split_columns = plot_data['metric'].str.split('_', n=1, expand=True)
plot_data['algorithm'] = split_columns[0].str.upper()
plot_data['fix_status'] = split_columns[1].apply(lambda x: 'After Fix' if x.startswith('fix') else 'Before Fix')
plot_data['category'] = plot_data['algorithm'] + ' ' + plot_data['fix_status']
# Remove no_scaffold entries
plot_data = plot_data[plot_data['scaffold_id'] != 'no_scaffold']

# Load metadata and ensure column names are correct
metadata = pd.read_csv("metadata_tags.csv")
print("Metadata columns:", metadata.columns.tolist())  # Add this line to debug

# Try to merge with available columns
# Comment out the problematic merge line for now
# merged_data = pd.merge(plot_data, metadata[['tag_id', 'num_cases_improved']], on='tag_id')

# Temporarily use plot_data instead of merged_data until we fix the merge
merged_data = plot_data.copy()

tag_deltas = scaffold_data.groupby('tag_id').apply(
    lambda x: (x['af3_fix_rmsd'] - x['af3_rmsd']).mean()
).reset_index(name='delta_rmsd')
selected_tags = tag_deltas.sort_values(by='delta_rmsd')['tag_id'][:50]


# Sort tags without using num_cases_improved for now
merged_data = merged_data[merged_data['tag_id'].isin(selected_tags)]
sorted_tags = merged_data['tag_id'].unique().tolist()

# Get unique tags and sort them
tags = sorted(plot_data['tag_id'].unique())
n_tags = len(tags)
tags_per_row = 6
n_rows = (n_tags + tags_per_row - 1) // tags_per_row  # Number of rows needed

# Create figure with one subplot per row
fig, axes = plt.subplots(n_rows, 1, figsize=(20, 5*n_rows))
if n_rows == 1:
    axes = [axes]  # Ensure axes is always a list

# Define consistent parameters
categories = ['AF2 Before Fix', 'AF2 After Fix', 'AF3 Before Fix', 'AF3 After Fix']
palette = {'AF2 Before Fix': '#a6cee3', 'AF2 After Fix': '#d3e5f3',
           'AF3 Before Fix': '#b2df8a', 'AF3 After Fix': '#d9f0c9'}

# Calculate global y-axis limit
y_max = plot_data['rmsd'].max() * 1.1
# PDF setup parameters
tags_per_page = 10
n_pages = (len(sorted_tags) + tags_per_page - 1) // tags_per_page
pdf_path = "tag_comparison_report.pdf"

# Create PDF with optimized settings
with PdfPages(pdf_path) as pdf:
    for page in range(n_pages):
        # Get tags for this page
        start_idx = page * tags_per_page
        end_idx = start_idx + tags_per_page
        page_tags = sorted_tags[start_idx:end_idx]
        
        # Create figure for this page
        if page == 0:
            figsize = (10.3, 4)
        else:
            figsize = (10, 4)
        fig = plt.figure(figsize=figsize, dpi=300)
        ax = fig.add_subplot(111)
        
        
        # Filter data for current page
        page_data = merged_data[merged_data['tag_id'].isin(page_tags)]
        
        # Add no_scaffold markers with matching colors
        no_scaffold_data = df[df['scaffold_id'] == 'no_scaffold']
        for tag in page_tags:
            tag_data = no_scaffold_data[no_scaffold_data['tag_id'] == tag]
            if not tag_data.empty:
                af2_rmsd = tag_data['af2_rmsd'].values[0]
                af3_rmsd = tag_data['af3_rmsd'].values[0]
                # Get x-coordinate for the tag
                x_pos = page_tags.index(tag)
                # Plot markers using the same colors as 'Before Fix'
                ax.plot(x_pos - 0.2, af2_rmsd, 'x', color='#a6cee3', markersize=8, markeredgewidth=2, label='AlphaFold-2 (No Scaffold)')  # AF2 color
                ax.plot(x_pos + 0.2, af3_rmsd, 'x', color='#b2df8a', markersize=8, markeredgewidth=2, label='AlphaFold-3 (No Scaffold)')  # AF3 color
    
        # Create violin plot
        sns.violinplot(
            data=page_data,
            x='tag_id',
            y='rmsd',
            hue='category',
            order=page_tags,  # Maintain sorted order
            hue_order=categories,
            palette=palette,
            ax=ax,
            inner="points",
            linewidth=0.5,
            dodge=True,
            scale="width",
            cut=0,
            width=0.7,  # Make violins thinner
            inner_kws={"s": 10},  # Increase the size of the dots
            legend=False  # Disable automatic legend
        )
        
        # Add alternating background shading and separators
        for i in range(len(page_tags)):
            # Add light gray background for even-numbered tags
            if i % 2 == 0:
                ax.axvspan(i - 0.5, i + 0.5, color='#f5f5f5', zorder=0)
            # Add vertical separator lines
            if i < len(page_tags) - 1:
                ax.axvline(x=i + 0.5, color='#e0e0e0', linestyle='-', linewidth=0.5, zorder=1)
        
        ax.set_xlim(-0.5, len(page_tags) - 0.5)
        
        ax.set_xlabel('')
        ax.set_ylabel('RMSD (Å)', fontsize=14)
        ax.set_ylim(0, 12)
        ax.grid(axis='y', alpha=0.6)
        ax.yaxis.set_tick_params(labelsize=14)
        page_tags = [tag.upper() for tag in page_tags]
        ax.set_xticklabels(page_tags, ha='center', fontsize=14)
        
        # Create unified legend
        legend_elements = [
            Line2D([0], [0], marker='x', color='#a6cee3', lw=0, markersize=8, markeredgewidth=2, label='AlphaFold-2 (No Scaffold)'),
            Line2D([0], [0], marker='x', color='#b2df8a', lw=0, markersize=8, markeredgewidth=2, label='AlphaFold-3 (No Scaffold)'),
            Line2D([0], [0], color='#a6cee3', lw=2, label='AlphaFold-2'),
            Line2D([0], [0], color='#d3e5f3', lw=2, label='AlphaFold-2 + Windowed MSA'),
            Line2D([0], [0], color='#b2df8a', lw=2, label='AlphaFold-3'),
            Line2D([0], [0], color='#d9f0c9', lw=2, label='AlphaFold-3 + Windowed MSA'),
        ]
        if page == 0:
            ax.legend(handles=legend_elements,
                 loc='upper center', bbox_to_anchor=(0.5, 1.20), ncol=3, fontsize=12)
        
        # Tight layout and save to PDF
        plt.tight_layout()
        pdf.savefig(fig, bbox_inches='tight')
        plt.close()

print(f"Report generated: {pdf_path}")
