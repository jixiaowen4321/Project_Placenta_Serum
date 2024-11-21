# -*- coding: utf-8 -*-
"""
Created on Thu Oct 24 12:57:51 2024

@author: jix01
"""

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

# Load data
df1 = pd.read_csv('Neg_ttest.csv')
df2 = pd.read_csv('Pos_ttest.csv')

df = pd.concat([df1, df2], axis=0)

# Scatter plot for all points (grey)
plt.scatter(x=df['log2fold'], y=df['P-value'].apply(lambda x: -np.log10(x)),
            s=1, color='grey')

# Highlight down- and up-regulated metabolites
down = df[(df['log2fold'] <= -1) & (df['P-value'] <= 0.05)]
up = df[(df['log2fold'] >= 1) & (df['P-value'] <= 0.05)]

down_mean = np.log(down.loc[:, 'P-Batch-1-Preterm-1':'P-Batch-4-Term-9'].mean(axis=1))
up_mean = np.log(up.loc[:, 'P-Batch-1-Preterm-1':'P-Batch-4-Term-9'].mean(axis=1))

plt.scatter(x=down['log2fold'], y=down['P-value'].apply(lambda x: -np.log10(x)),
            s=1, color="blue")
plt.scatter(x=up['log2fold'], y=up['P-value'].apply(lambda x: -np.log10(x)),
            s=1 , color="red")


# Plot formatting
plt.xlabel("Fold Change")
plt.ylabel("-Log FDR")
plt.axvline(-1, color="black", linewidth=0.5, linestyle="--")
plt.axvline(1, color="black", linewidth=0.5, linestyle="--")
plt.axhline(1.3, color="black", linewidth=0.5, linestyle="--")

# Add legend
#plt.legend()

# Save and show the plot
plt.savefig('VolcanoPlt-3.png', dpi=600, bbox_inches='tight')
plt.show()
