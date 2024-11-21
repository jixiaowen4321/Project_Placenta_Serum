# -*- coding: utf-8 -*-
"""
Created on Wed Jul 31 11:48:07 2024

@author: jix01
"""

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from sklearn.cluster import KMeans
from sklearn.decomposition import PCA
from sklearn.preprocessing import StandardScaler

df1 = pd.read_csv('PosPVFserumStats2c.csv')
df2 = pd.read_csv('NegPVFserumStats2c.csv')

df = pd.concat([df1, df2], axis=0)

df['chem_id'] = np.round(df['Average Mz'], 3).astype(str) + '_@_' + np.round(df['Average Rt(min)'], 3).astype(str)

df.columns.values

df = df.set_index('chem_id')
loc1 = df.loc[:, df.columns.str.contains('b1')]
loc2 = df.loc[:, df.columns.str.contains('b2')]
loc3 = df.loc[:, df.columns.str.contains('b3')]
loc4 = df.loc[:, df.columns.str.contains('b4')]
dr = pd.concat([loc1, loc2, loc3, loc4], axis=1)
dr.columns = dr.columns.str.replace('Fullterm', 'F')
dr.columns = dr.columns.str.replace('Preterm', 'P')

dr = dr.T

dr['color_MC'] = np.where(dr.index.str.contains('F'), 'darkgray','dodgerblue')
dr['color_batch'] = np.where(dr.index.str.contains('-b1'), 'darkgray', 
                    np.where(dr.index.str.contains('-b2'),'dodgerblue',
                    np.where(dr.index.str.contains('-b3'), 'red','pink')))

dr12 = dr.loc[:, 'color_MC':'color_batch']
dr11 = dr.drop(['color_MC', 'color_batch'], axis=1)

# Standardize the data to have a mean of ~0 and a variance of 1
X_std = StandardScaler().fit_transform(dr11)

# Create a PCA instance: pca
pca = PCA(n_components=20)
principalComponents = pca.fit_transform(X_std)

pca.explained_variance_ratio_

# Plot the explained variances
features = range(pca.n_components_)

plt.bar(features, pca.explained_variance_ratio_, color='blue')
plt.xlabel('PCA features')
plt.ylabel('Variance %')
plt.xticks(features)
plt.savefig('PCA features and the variance explained.png', dpi=400, bbox_inches = "tight")

plt.show()

# Save components to a DataFrame
PCA_components = pd.DataFrame(principalComponents)

plt.scatter(PCA_components[0], PCA_components[1], alpha=.1, color='blue')
plt.xlabel('PC1')
plt.ylabel('PC2')
plt.savefig('PC1 and PC2 as a scatterplot.png', dpi=400, bbox_inches = "tight")
plt.show()

ks = range(1, 10)
inertias = []
for k in ks:
    # Create a KMeans instance with k clusters: model
    model = KMeans(n_clusters=k, n_init=10)
    
    # Fit model to samples
    model.fit(PCA_components.iloc[:,:3])

    # Append the inertia to the list of inertias
    inertias.append(model.inertia_)
    
plt.plot(ks, inertias, '-o', color='black')
plt.xlabel('number of clusters, k')
plt.ylabel('Inertia')
plt.xticks(ks)
plt.savefig('approximation of the optimal number of clusters in the dataset.png', dpi=400, bbox_inches = "tight")

plt.show()

from sklearn.cluster import KMeans

kmeans4 = KMeans(n_clusters=4)
y_kmeans4 =kmeans4.fit_predict(PCA_components)
print(y_kmeans4)

PCA_components.columns = ['PC'+ str(col) for col in PCA_components.columns]

PCA_components = PCA_components.reset_index(drop=True)
dr12 = dr12.reset_index(drop=True)

PCA_components = pd.concat([PCA_components, dr12], axis=1)
PCA_components.to_csv('PCA_components_POS.csv')

from matplotlib.lines import Line2D
custom_lines = [Line2D([0], [0], color='darkgray', marker='o', lw=0),
                Line2D([0], [0], color='dodgerblue', marker='o', lw=0)]

fig, ax = plt.subplots()
scatter = plt.scatter(PCA_components['PC0'], PCA_components['PC1'], c=PCA_components['color_MC'])
legend1 = ax.legend(custom_lines, ['F', 'P'],
                title="Groups", loc='center left', bbox_to_anchor=(1, 0.5))
plt.xlabel('PC1')
plt.ylabel('PC2')
plt.show()
scatter.figure.savefig('PC1 and PC2 color-coded by sample type.png', dpi=400, bbox_inches = "tight")

from matplotlib.lines import Line2D
custom_lines = [Line2D([0], [0], color='darkgray', marker='o', lw=0),
                Line2D([0], [0], color='dodgerblue', marker='o', lw=0),
                Line2D([0], [0], color='red', marker='o', lw=0),
                Line2D([0], [0], color='pink', marker='o', lw=0)]

fig, ax = plt.subplots()
scatter = plt.scatter(PCA_components['PC0'], PCA_components['PC1'], c=PCA_components['color_batch'])
legend1 = ax.legend(custom_lines, ['1', '2','3','4'],
                title="Batch", loc='center left', bbox_to_anchor=(1, 0.5))
plt.xlabel('PC1')
plt.ylabel('PC2')
plt.show()
scatter.figure.savefig('PC1 and PC2 color-coded by shipment.png', dpi=400, bbox_inches = "tight")

fig, ax = plt.subplots()
scatter = plt.scatter(PCA_components['PC0'], PCA_components['PC1'], c=y_kmeans4, cmap = 'Set2')
legend1 = ax.legend(*scatter.legend_elements(),
                title="Clusters", loc='center left', bbox_to_anchor=(1, 0.5))
plt.xlabel('PC1')
plt.ylabel('PC2')
plt.show()
scatter.figure.savefig('agnostically derived clusters using a k-means algorithm.png', dpi=400, bbox_inches = "tight")

import matplotlib.pyplot as plt
import seaborn as sns
from statannot import add_stat_annotation

df = PCA_components
df['color_batch'] = df['color_batch'].str.replace('darkgray','1')
df['color_batch'] = df['color_batch'].str.replace('dodgerblue','2')
df['color_batch'] = df['color_batch'].str.replace('red','3')
df['color_batch'] = df['color_batch'].str.replace('pink','4')

x = "color_batch"
y = "PC0"
my_pal = {"1": "darkgray", "2": "dodgerblue", "3": "red", "4": "pink"}
order = ['1', '2', '3', '4']
ax = sns.boxplot(data=df, x=x, y=y, order=order, palette=my_pal) 
plt.xlabel('Batch')
plt.ylabel('PC1')
add_stat_annotation(ax, data=df, x=x, y=y, order=order,
                    box_pairs=[("1", "2"),("2", "3"),("2", "4")],
                    test='Mann-Whitney', text_format='star', loc='outside', verbose=2)
plt.savefig('boxplot for PC1 by shipment.png', dpi=400, bbox_inches = "tight")

df = PCA_components
df['color_MC'] = df['color_MC'].str.replace('darkgray','F')
df['color_MC'] = df['color_MC'].str.replace('dodgerblue','P')

x = "color_MC"
y = "PC0"
my_pal = {"F": "darkgray", "P": "orangered"}
order = ['F', 'P']
ax = sns.boxplot(data=df, x=x, y=y, order=order, palette=my_pal) 
plt.xlabel('Sample type')
plt.ylabel('PC1')
add_stat_annotation(ax, data=df, x=x, y=y, order=order,
                    box_pairs=[("F", "P")],
                    test='Mann-Whitney', text_format='star', loc='outside', verbose=2)
plt.savefig('boxplot for PC1 by sample type.png', dpi=400, bbox_inches = "tight")
