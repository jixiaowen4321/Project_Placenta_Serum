#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Apr 20 19:50:34 2024

@author: dabrahamsson
"""

import numpy as np
import pandas as pd
import scipy
from scipy import stats
from scipy.stats import linregress
import seaborn as sns
import matplotlib as mlp
import matplotlib.pyplot as plt

df = pd.read_csv('Neg_CleanFnaBac.csv')
ga = pd.read_csv('GestationalAge.csv')
#ga['ga_weeks'] = ga['ga_delivery']/7
#ga['Term'] = np.where(ga['ga_weeks'] < 37, 'P', 'F')

df.columns.values
ga.columns.values

loc1 = df.loc[:, 'chem_id':'MS/MS spectrum']
loc3 = df.loc[:, df.columns.str.contains('P-Batch')]

def stat(samples):
    df = pd.concat([loc1, samples], axis=1)
    
    print(df.shape)
    
    ids = samples.columns.values.tolist()
    mt = pd.melt(df, id_vars=['Alignment ID', 'Average Mz', 'Average Rt(min)'], value_vars=ids, value_name='Abundance')
    mt.columns = mt.columns.str.replace('variable', 'Sample Name')
    
    mt = pd.merge(mt, ga, on='Sample Name', how='left')
    mt['PVF'] = np.where(mt['Term'] == 'Term', 0, 1)
    mt['logA'] = np.log10(mt['Abundance'])
    
    def get_element(my_list, Position):
        return my_list[Position]
    
    mt['PVF'] = mt['PVF'].astype(str)
    mt['Sample Name'] = mt['Sample Name'].astype(str)
    mt['identifier'] = mt['PVF'] + '_' + mt['Sample Name']
    mt = mt.drop(['PVF','Sample Name'], axis=1)
    mt = mt.pivot_table('logA','Alignment ID','identifier')
    
    axisvalues = mt.columns.values
    axisvalues_mt = pd.DataFrame({'identifier': axisvalues})
    axisvalues_ = axisvalues_mt['identifier'].str.split('_').apply(lambda x: x[0])
    axisvalues_ = axisvalues_.astype(float)
    
    mt = mt.astype(float)
    def calc_slope(row):
        a = scipy.stats.linregress(axisvalues_ , row)
        return pd.Series(a._asdict())
    
    print (mt.apply(calc_slope,axis=1))
    
    mt = mt.join(mt.apply(calc_slope,axis=1))
    return(mt)


mts = stat(loc3)
mts = mts.reset_index()
mts = pd.concat([loc1, mts], axis=1)
mts.to_csv('NegPVFStats1c.csv')

#Fold-change
mts.columns.values
controls = mts.loc[:, mts.columns.str.startswith('0_')]
cases = mts.loc[:, mts.columns.str.startswith('1_')]
controls = 10**controls
cases = 10**cases
conav = controls.mean(axis=1)
casav = cases.mean(axis=1)
fold = casav/conav
mts['log2fold'] = np.log2(fold)
#mts = mts[(mts['log2fold'] > 0.1)|(mts['log2fold'] < -0.1)]
mts.to_csv('NegPVFStats3c.csv')

#Benjamini-Hochberg filtering p-values
mts = mts.sort_values(by='pvalue')

mts['R2'] = mts['rvalue']**2
mts = mts[mts['R2'] > 0.1]
mts['rank'] = mts.reset_index().index + 1
mts['(I/m)Q'] = (mts['rank']/len(mts))*0.05
mts['(I/m)Q - p'] = mts['(I/m)Q'] - mts['pvalue']
mts = mts.sort_values(by='(I/m)Q - p', ascending =False)
mts['BH_sig'] = np.where(mts['(I/m)Q - p'] < 0, 0, 1)


# isolate the features that show a significant relationship with PTB
mts = mts[mts['pvalue'] < 0.05]
mts.to_csv('NegPVFStats2c.csv')

print('Significant after BH:', mts['BH_sig'].sum())

