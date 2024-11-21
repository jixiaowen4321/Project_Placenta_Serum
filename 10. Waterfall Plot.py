import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from adjustText import adjust_text  # For text label adjustments

# Load the data
df = pd.read_csv("merged_test.csv")
print(df.columns.values)

# Select relevant columns
df = df[['Alignment ID','Adduct type', 'Metabolite name', 'log2fold', 'pvalue']]

# Calculate rank
df['rank'] = -1 * df['log2fold'].rank(method='max')

# Filter top and bottom 5
df_top_5 = df.nlargest(10, 'log2fold')
df_tail_5 = df.nsmallest(10, 'log2fold')

# Set scale factor for marker size
scale_factor = 50

# Create the scatter plot
plt.figure(figsize=(10, 6))
scatter = plt.scatter(
    df['rank'],
    df['log2fold'],
    c=df['pvalue'],  # Color based on pvalue
    s=np.abs(df['log2fold']) * scale_factor,  # Size based on log2FoldChange
    cmap='RdYlBu_r',  # Using a diverging color palette
    alpha=0.7,
    edgecolors='face',  # To avoid marker edge color warning
    vmax=0.05
)

# Add color bar for p-value
cbar = plt.colorbar(scatter)
cbar.set_label('p-value', fontsize=15)
cbar.ax.tick_params(labelsize=14)   

# Add horizontal and vertical lines
plt.axhline(1, color='grey', linestyle='--')
plt.axhline(-1, color='grey', linestyle='--')
#plt.axvline(-500, color='grey', linestyle='-')

# Create custom legend for size
for size in [1, 2, 5]:  # Example log2fold magnitudes
    plt.scatter([], [], s=size * scale_factor, color='gray', alpha=0.6,
                label=f'|log2fold| = {size}')

# Adjust legend to be outside the plot on the right
plt.legend(scatterpoints=1, frameon=True, labelspacing=1, title="log2FoldChange",
           loc="upper left", bbox_to_anchor=(1.2, 1))

# Set axes and titles
plt.xticks(ticks=[-1000, -750, -500, -250, 0], labels=[0, 250, 500, 750, 1000])
plt.xlabel("Rank of differentially expressed chemicals", fontsize=14)
plt.ylabel("Fold Change", fontsize=14)
#plt.grid()
plt.tight_layout()

plt.savefig("Differentially_Expressed_Chemicals-2.png", dpi=600, bbox_inches="tight")
plt.show()

df.to_csv('rank.csv')
