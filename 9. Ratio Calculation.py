# -*- coding: utf-8 -*-
"""
Created on Fri Nov  8 14:21:47 2024

@author: jix01
"""

import pandas as pd
import numpy as np

df1 = pd.read_csv('NEG.csv')
df1 = df1.drop(columns=df1.columns[df1.columns.str.contains('lab-control')])
df2 = pd.read_excel('Sample Names.xlsx')

df1s = df1.loc[:, df1.columns.str.contains('S-Batch')]
df1p = df1.loc[:, df1.columns.str.contains('P-Batch')]


placenta_columns = df2['Placenta'].values.flatten()
serum_columns = df2['Serum'].values.flatten()

results = pd.DataFrame()

# Divide matching 'S-Batch' and 'P-Batch' columns based on df2 references
for s_column, p_column in zip(serum_columns, placenta_columns):
    if s_column in df1s.columns and p_column in df1p.columns:
        # Add a small constant to avoid division by zero
        results[f'{s_column}_div_{p_column}'] = df1s[s_column] / (df1p[p_column] + 1e-10)

# Replace inf with NaN if any remain
results.replace([np.inf, -np.inf], np.nan, inplace=True)

# Display the resulting DataFrame with divided columns
print(results)

df = pd.concat([df1['Metabolite name'], df1['Level'], results], axis=1)

df.to_csv('NEG_level_1_sdivp_2.csv')

dfL = np.log10(results)
dfLa = pd.concat([df1['Metabolite name'], df1['Level'], dfL], axis=1)
dfLa.to_csv('NEG_level_1_sdivp_Log.csv')
