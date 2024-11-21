# -*- coding: utf-8 -*-
"""
Created on Wed May  1 16:12:52 2024

@author: Abrahd05
"""
import numpy as np
import pandas as pd
import seaborn as sns
from matplotlib import pyplot as plt
from matplotlib.patches import Patch

df1 = pd.read_csv('NegPVFStats2c.csv')
df2 = pd.read_csv('PosPVFStats2c.csv')

df = pd.concat([df1, df2], axis=0)

df.columns
df = df[df['BH_sig'] == 1]
df.columns.values
df = df.set_index('chem_id')


dfH = df.loc[:, '0_P-Batch-1-Term-10':'1_P-Batch-4-Preterm-3']
dfH.columns.values

col_names = dfH.columns.values
col_namesDF = pd.DataFrame(col_names, columns=['columns'])

# for term color
col_namesDF['Term'] = np.where(col_namesDF['columns'].str.contains('1_'), 'orangered', 'darkgray')

# for batch color
conditions = [
    col_namesDF['columns'].str.contains('Batch-1'),
    col_namesDF['columns'].str.contains('Batch-2'),
    col_namesDF['columns'].str.contains('Batch-3'),
    col_namesDF['columns'].str.contains('Batch-4')
    ]

choices = ['#4575b4',  # cool blue for Batch-1
           '#91bfdb',  # light blue for Batch-2
           '#fee090',  # light orange for Batch-3
           '#d73027']  # warm red for Batch-4

col_namesDF['Batch'] = np.select(conditions, choices, default='gray')



col_namesDF = col_namesDF.set_index('columns')



print(col_namesDF)

# for row color
df['Annotation'] = np.where(df['Metabolite name'].str.contains('Unknown'), 'Unknown', 'Annotated')
species = df.pop("Annotation")
Row_Name = dict(zip(species.unique(), ['yellow','blue']))
Row_colors = species.map(Row_Name)



# Figure size

sns.set(font_scale=1)
g = sns.clustermap(dfH, cmap='Blues', 
                   col_colors=col_namesDF, 
                   row_colors=Row_colors,
                   colors_ratio=(0.025, 0.025), standard_scale=0, 
                   dendrogram_ratio=0.08, cbar_pos=(-0.07, .45, .03, .2))

g.savefig('Comb_RPL_BH_sig2c_placenta.png', dpi=400)

# Create legend elements
legend_elements = [
    Patch(facecolor='orangered', edgecolor='black', label='Preterm'),
    Patch(facecolor='darkgray', edgecolor='black', label='Term'),
    Patch(facecolor='#4575b4', edgecolor='black', label='Batch-1'),
    Patch(facecolor='#91bfdb', edgecolor='black', label='Batch-2'),
    Patch(facecolor='#fee090', edgecolor='black', label='Batch-3'),
    Patch(facecolor='#d73027', edgecolor='black', label='Batch-4'),
    Patch(facecolor='yellow', edgecolor='black', label='Unknown'),
    Patch(facecolor='blue', edgecolor='black', label='Annotated')
]

# Add legend to plot
plt.legend(handles=legend_elements, loc='upper left', bbox_to_anchor=(1.05, 1), title="Legend")
plt.savefig('legend_only.png', dpi=400, bbox_inches='tight') 


