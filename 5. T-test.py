# -*- coding: utf-8 -*-
"""
Created on Thu Oct 24 09:28:01 2024

@author: jix01
"""

import numpy as np
import pandas as pd
from scipy.stats import ttest_ind

df = pd.read_csv('Neg_CleanFnaBac.csv')
ga = pd.read_csv('GestationalAge.csv')

loc1 = df.loc[:, 'chem_id':'MS/MS spectrum']
loc3 = df.loc[:, df.columns.str.contains('P-Batch')]

def stat(samples):
    df_full = pd.concat([loc1, samples], axis=1)
    
    print(df_full.shape)
    
    ids = samples.columns.values.tolist()
    mt = pd.melt(df_full, id_vars=['Alignment ID', 'Average Mz', 'Average Rt(min)'], value_vars=ids, value_name='Abundance')
    mt.columns = mt.columns.str.replace('variable', 'Sample Name')
    
    mt = pd.merge(mt, ga, on='Sample Name', how='left')
    mt['PVF'] = np.where(mt['Term'] == 'Term', 0, 1)
    mt['logA'] = np.log10(mt['Abundance'])
    
    mt['PVF'] = mt['PVF'].astype(int)  # Ensure PVF is numeric for grouping
    mt['Sample Name'] = mt['Sample Name'].astype(str)
    mt['identifier'] = mt['PVF'].astype(str) + '_' + mt['Sample Name']
    
    # Pivot the data
    mt_pivot = mt.pivot_table('logA', 'Alignment ID', 'identifier')

    # Calculate fold change and log2 fold change
    controls = mt_pivot.loc[:, mt_pivot.columns.str.startswith('0_')]
    cases = mt_pivot.loc[:, mt_pivot.columns.str.startswith('1_')]
    
    # Convert log values back to linear space (since we log-transformed the data earlier)
    controls_linear = 10 ** controls
    cases_linear = 10 ** cases
    
    # Calculate mean abundance for control and case groups
    conav = controls_linear.mean(axis=1)
    casav = cases_linear.mean(axis=1)
    
    # Calculate fold change (case/control) and log2 fold change
    fold = casav / conav
    mt_pivot['log2fold'] = np.log2(fold)
    
    # T-test calculations for each Alignment ID
    ttest_results = []
    for alignment_id, group in mt.groupby('Alignment ID'):
        group_term = group[group['PVF'] == 0]['logA']
        group_preterm = group[group['PVF'] == 1]['logA']
        
        # Perform T-test
        t_stat, p_value = ttest_ind(group_term, group_preterm, nan_policy='omit')
        ttest_results.append([alignment_id, t_stat, p_value])
    
    # Convert results to DataFrame
    ttest_df = pd.DataFrame(ttest_results, columns=['Alignment ID', 'T-stat', 'P-value'])
    
    # Benjamini-Hochberg correction (5% FDR)
    ttest_df['rank'] = ttest_df['P-value'].rank(method='min')  # Rank the p-values
    m = len(ttest_df)  # Total number of tests
    fdr = 0.05  # FDR threshold
    
    # Adjusted p-values using Benjamini-Hochberg formula
    ttest_df['BH-adjusted P-value'] = (ttest_df['P-value'] * m) / ttest_df['rank']
    
    # Clip adjusted p-values to be <= 1
    ttest_df['BH-adjusted P-value'] = ttest_df['BH-adjusted P-value'].clip(upper=1)
    
    # Merge T-test results, Benjamini-Hochberg adjusted p-values, and log2 fold change into the original dataframe
    df_with_ttest = pd.merge(df_full, ttest_df, on='Alignment ID', how='left')
    df_with_ttest = pd.merge(df_with_ttest, mt_pivot[['log2fold']], on='Alignment ID', how='left')
    
    # Return the modified DataFrame with T-test results, BH-adjusted p-values, and log2 fold change
    return df_with_ttest

# Example call to the stat function
samples = loc3  # Assuming 'loc3' contains the samples for analysis
df_with_results = stat(samples)

# Print the dataframe with T-test, fold change, and BH-adjusted p-values
print(df_with_results)



df_with_results.to_csv('ttest.csv')
