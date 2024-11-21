# -*- coding: utf-8 -*-
"""
Created on Tue Feb 20 16:57:20 2024

@author: Abrahd05
"""

import pandas as pd
import numpy as np


df = pd.read_csv('POS.csv')
df = df.drop(columns=df.columns[df.columns.str.contains('lab-control')])
df['chem_id'] = (
                np.round(df['Average Mz'], 3).astype(str) + '_@_' 
                + np.round(df['Average Rt(min)'], 3).astype(str)
                )

df.columns.values

loc1a = df['chem_id']
loc1b = df.loc[:, 'Alignment ID':'MS/MS spectrum']
loc1 = pd.concat([loc1a, loc1b], axis=1)
loc1 = loc1.set_index('chem_id')



# 1 for frequency in fullterm in serum
loc2a = df.loc[:, df.columns.str.contains('S-Batch')]
loc3 = pd.concat([loc1a, loc2a], axis=1)
loc3 = loc3.set_index('chem_id')
loc3[loc3 < 10000] = np.nan

loc3_S = loc3.loc[:, loc3.columns.str.contains('Term')]
loc3_S['count'] = loc3_S.count(axis = 1)
loc3_S['freq'] = loc3_S['count']/len(loc3_S.columns)
loc3_S['freq%'] = loc3_S['freq']*100
loc1 = loc1.reset_index()
loc3_S = loc3_S.reset_index()
df1 = pd.merge(loc1, loc3_S, how='inner')
df1 = df1.set_index('chem_id')
df1.to_csv('POS_Serum_Fullterm.csv')


# 2 for frequency in Preterm in serum

loc2a = df.loc[:, df.columns.str.contains('S-Batch')]
loc3 = pd.concat([loc1a, loc2a], axis=1)
loc3 = loc3.set_index('chem_id')
loc3[loc3 < 10000] = np.nan

loc3_S = loc3.loc[:, loc3.columns.str.contains('Preterm')]
loc3_S['count'] = loc3_S.count(axis = 1)
loc3_S['freq'] = loc3_S['count']/len(loc3_S.columns)
loc3_S['freq%'] = loc3_S['freq']*100
loc1 = loc1.reset_index()
loc3_S = loc3_S.reset_index()
df1 = pd.merge(loc1, loc3_S, how='inner')
df1 = df1.set_index('chem_id')
df1.to_csv('POS_Serum_Preterm.csv')



# 3 for frequency in fullterm in placenta
loc2a = df.loc[:, df.columns.str.contains('P-Batch')]
loc3 = pd.concat([loc1a, loc2a], axis=1)
loc3 = loc3.set_index('chem_id')
loc3[loc3 < 10000] = np.nan

loc3_S = loc3.loc[:, loc3.columns.str.contains('Term')]
loc3_S['count'] = loc3_S.count(axis = 1)
loc3_S['freq'] = loc3_S['count']/len(loc3_S.columns)
loc3_S['freq%'] = loc3_S['freq']*100
loc1 = loc1.reset_index()
loc3_S = loc3_S.reset_index()
df1 = pd.merge(loc1, loc3_S, how='inner')
df1 = df1.set_index('chem_id')
df1.to_csv('POS_Placenta_Fullterm.csv')


# 4 for frequency in Preterm in placenta

loc2a = df.loc[:, df.columns.str.contains('P-Batch')]
loc3 = pd.concat([loc1a, loc2a], axis=1)
loc3 = loc3.set_index('chem_id')
loc3[loc3 < 10000] = np.nan

loc3_S = loc3.loc[:, loc3.columns.str.contains('Preterm')]
loc3_S['count'] = loc3_S.count(axis = 1)
loc3_S['freq'] = loc3_S['count']/len(loc3_S.columns)
loc3_S['freq%'] = loc3_S['freq']*100
loc1 = loc1.reset_index()
loc3_S = loc3_S.reset_index()
df1 = pd.merge(loc1, loc3_S, how='inner')
df1 = df1.set_index('chem_id')
df1.to_csv('POS_Placenta_Preterm.csv')



# 5 for frequency in serum
loc2a = df.loc[:, df.columns.str.contains('S-Batch')]
loc3 = pd.concat([loc1a, loc2a], axis=1)
loc3 = loc3.set_index('chem_id')
loc3[loc3 < 10000] = np.nan

loc3_S = loc3.loc[:, loc3.columns.str.contains('S-Batch')]
loc3_S['count'] = loc3_S.count(axis = 1)
loc3_S['freq'] = loc3_S['count']/len(loc3_S.columns)
loc3_S['freq%'] = loc3_S['freq']*100
loc1 = loc1.reset_index()
loc3_S = loc3_S.reset_index()
df1 = pd.merge(loc1, loc3_S, how='inner')
df1 = df1.set_index('chem_id')
df1.to_csv('POS_freq_Serum.csv')

# 6 for frequency in placenta
loc2b = df.loc[:, df.columns.str.contains('P-Batch')]
loc3 = pd.concat([loc1a, loc2b], axis=1)
loc3 = loc3.set_index('chem_id')
loc3[loc3 < 10000] = np.nan
loc3_U = loc3.loc[:, loc3.columns.str.contains('P-Batch')]
loc3_U['count'] = loc3_U.count(axis = 1)
loc3_U['freq'] = loc3_U['count']/len(loc3.columns)
loc3_U['freq%'] = loc3_U['freq']*100
loc3_U = loc3_U.reset_index()
df2 = pd.merge(loc1, loc3_U, on= 'chem_id', how='inner')
df2 = df2.set_index('chem_id')
df2.to_csv('POS_freq_placenta.csv')


# 7 for frequency in urine&serum
loc2a = df.loc[:, df.columns.str.contains('S-Batch')]
loc2b = df.loc[:, df.columns.str.contains('P-Batch')]
loc3 = pd.concat([loc1a, loc2a, loc2b], axis=1)
loc3 = loc3.set_index('chem_id')
loc3[loc3 < 10000] = np.nan
loc3_S_P = loc3.loc[:, loc3.columns.str.contains('S-Batch|P-Batch')]
loc3_S_P['count'] = loc3_S_P.count(axis = 1)
loc3_S_P['freq'] = loc3_S_P['count']/len(loc3.columns)
loc3_S_P['freq%'] = loc3_S_P['freq']*100
loc3_S_P = loc3_S_P.reset_index()
df3 = pd.merge(loc1, loc3_S_P, on='chem_id', how='inner')
#df3 = df.set_index('chem_id')
df3.to_csv('POS_freq_Serum_Placenta.csv')










