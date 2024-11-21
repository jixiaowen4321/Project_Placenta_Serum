# -*- coding: utf-8 -*-
"""
Created on Tue Feb 20 16:57:20 2024

@author: Abrahd05
"""

import pandas as pd
import numpy as np
from scipy import stats

df = pd.read_csv('NEG_Serum_Placenta_annotated_clean.csv')
df['chem_id'] = (
                np.round(df['Average Mz'], 3).astype(str) + '_@_' 
                + np.round(df['Average Rt(min)'], 3).astype(str)
                )

df.columns.values

loc1a = df['chem_id']
loc1b = df.loc[:, 'Alignment ID':'MS/MS spectrum']
loc1 = pd.concat([loc1a, loc1b], axis=1)
loc1 = loc1.set_index('chem_id')

loc2a = df.loc[:, df.columns.str.contains('S-Batch')]
loc2b = df.loc[:, df.columns.str.contains('P-Batch')]
loc3 = pd.concat([loc1a, loc2a, loc2b], axis=1)
loc3 = loc3.set_index('chem_id')
loc3[loc3 < 10000] = np.nan

loc3['count'] = loc3.count(axis = 1)
loc3['freq'] = loc3['count']/len(loc3.columns)
loc3['freq%'] = loc3['freq']*100
loc3 = loc3[loc3['freq%'] > 60]
loc3 = loc3.loc[:, loc3.columns.str.contains('P-Batch')]

loc3L = np.log10(loc3)

def fillNaN_with_unifrand(df):
    lower, upper = 0, df.min()
    a = df.values
    m = np.isnan(a) # mask of NaNs
    mu, sigma = df.min(), df.std()
    a[m] = stats.truncnorm.rvs(
          (lower-mu)/sigma,(upper-mu)/sigma,loc=mu,scale=sigma,size=m.sum())
    return df


loc3Lmod = loc3L.apply(fillNaN_with_unifrand)
loc3Lmod = 10**loc3Lmod

loc1 = loc1.reset_index()
loc3Lmod = loc3Lmod.reset_index()
df = pd.merge(loc1, loc3Lmod, how='inner')
df = df.set_index('chem_id')
df.to_csv('Neg_Clean.csv')

