# -*- coding: utf-8 -*-
"""
Created on Tue May 21 19:29:40 2024

@author: Abrahd05
"""

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.cm
import seaborn as sns

# Load the datasets
df1 = pd.read_csv('Pos_Clean.csv')
df2 = pd.read_csv('Neg_Clean.csv')

# Concatenate the datasets
df = pd.concat([df1, df2], axis=0)

# Print the column names for verification
print(df.columns.values)

# Select relevant columns for analysis
dfH = df.loc[:, 'P-Batch-1-lab-control-1':'P-Batch-4-Term-9']

# Apply log10 transformation
dfH = np.log10(dfH)

# Create separate DataFrames for Term and Preterm samples
dfH_T = dfH.loc[:, dfH.columns.str.contains('Term')]
dfH_P = dfH.loc[:, dfH.columns.str.contains('Preterm')]

# Calculate mean logA values for all samples, Term, and Preterm
dfH['meanlogA'] = dfH.mean(axis=1)
dfH['mean_T'] = dfH_T.mean(axis=1)
dfH['mean_P'] = dfH_P.mean(axis=1)

# Define x and y for the joint plot
x = 'mean_T'
y = 'mean_P'

# Set Seaborn aesthetics
sns.set(font_scale=1)
sns.set_style('white')

# Create the joint plot
fig = sns.jointplot(x=x, y=y, data=dfH, color='darkcyan', kind='reg')

# Set labels for the axes
plt.xlabel('Mean logA Term')
plt.ylabel('Mean logA Preterm')

# Show the plot
plt.show()

fig.savefig('regression_before combat.png', dpi=1200, bbox_inches='tight')

from scipy.stats import pearsonr, t as tdist
dfH = dfH.replace([np.nan, -np.inf], 0)
r, p = pearsonr(dfH['meanU'], dfH['meanS'])
