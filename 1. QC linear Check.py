# -*- coding: utf-8 -*-
"""
Created on Mon Apr  8 16:14:53 2024

@author: Abrahd05
"""

import numpy as np
import pandas as pd
import scipy
from scipy import stats
from scipy.stats import linregress

df1 = pd.read_csv('POS_EPA.csv')
df2 = pd.read_csv('EPA mixes.csv')

# Subtract the mass of a proton to estimate the M in [M+H]+
df1['Average Mz'] = df1['Average Mz'] - 1.0078250319

df1['Mass_round'] = np.round(df1['Average Mz'], 3)
df2['Mass_round'] = np.round(df2['Monoisotopic mass'], 3)

df = pd.merge(df1, df2, on='Mass_round', how='inner')

df['Mass_diff'] = ((np.absolute(df['Average Mz'] - df['Monoisotopic mass']))/df['Monoisotopic mass'])*10**6
df['Mass_flag'] = np.where(df['Mass_diff'] > 5, 1, 0)

df.columns.values

loc1 = df.loc[:, 'Alignment ID':'MS/MS spectrum']
loc2 = df.loc[:, 'Mass_round':'Mass_flag']
print(df1.columns.values)

dft = df.loc[:, 'EPAmix-1000ppb':'EPAmix-50ppb']
dft = df[['EPAmix-1000ppb', 'EPAmix-200ppb', 'EPAmix-500ppb', 'EPAmix-50ppb']]
dft = dft.replace(0, 0.01)
dft = np.log10(dft)


axisvalues_ = [1000, 500, 200, 50]

dft = dft.astype(float)
def calc_slope(row):
    a = scipy.stats.linregress(axisvalues_ , row)
    return pd.Series(a._asdict())

print (dft.apply(calc_slope,axis=1))

dft = dft.join(dft.apply(calc_slope,axis=1))
dft

dft = pd.concat([loc1, loc2, dft], axis=1)
dft.to_csv('SlopesStats.csv')



dft['Mass_round'] = np.round(dft['Average Mz'], 2)
df2['Mass_round'] = np.round(df2['Monoisotopic mass'], 2)

df3 = pd.merge(dft, df2, on='Mass_round', how='inner')
print(df3.columns.values)

df3['Mass_diff'] = ((np.absolute(df3['Average Mz'] - df3['Monoisotopic mass_x']))/df3['Monoisotopic mass_x'])*10**6
df3['Mass_flag'] = np.where(df3['Mass_diff'] > 5, 1, 0)

print(df3)
df3.to_csv('MassSearch.csv')



