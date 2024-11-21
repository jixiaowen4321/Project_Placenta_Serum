import pandas as pd
import numpy as np
from scipy.stats import spearmanr
import networkx as nx
import matplotlib.pyplot as plt
from statsmodels.stats.multitest import multipletests

# Load the data
df1 = pd.read_csv('Endogenous_placenta.csv', index_col=0)
df2 = pd.read_csv('PCP_exogenous_placenta.csv', index_col=0)

# Transpose the data
df1 = df1.T
df2 = df2.T

# Calculate Spearman correlation and p-values
r, p = spearmanr(df1.values, df2.values, axis=0)

# Extract the relevant parts for df1 vs df2 correlation and p-values
r_df1_df2 = r[:df1.shape[1], df1.shape[1]:]
p_values_df1_df2 = p[:df1.shape[1], df1.shape[1]:]

# Create DataFrames for correlation and p-values with appropriate indices and columns
r_df = pd.DataFrame(r_df1_df2, index=df1.columns, columns=df2.columns)
p_df = pd.DataFrame(p_values_df1_df2, index=df1.columns, columns=df2.columns)

# Set correlations to 0 where absolute value is less than 0.9
r_df[abs(r_df) < 0.9] = 0

# Set p-values: p >= 0.05 becomes -1, 0 <= p < 0.05 becomes 1, p == -1 becomes 0
p_df = p_df.applymap(lambda x: -1 if x >= 0.05 else (1 if x < 0.05 and x >= 0 else 0))

# Create matrix z
z_df = r_df * p_df

# Make the correlation matrix symmetric for igraph compatibility
z1 = np.zeros_like(r)
z1[abs(r) != 0] = 0  # Make the non-zero entries 0
for i in range(len(r_df.columns)):
    for j in range(len(r_df.columns)):
        if z_df.iloc[i, j] != 0:
            z1[i, j] = z_df.iloc[i, j]
            z1[j, i] = z_df.iloc[i, j]

# Convert the matrix into a graph using NetworkX
G = nx.from_numpy_array(z1)

# Remove isolated nodes (degree == 0)
isolated_nodes = list(nx.isolates(G))
G.remove_nodes_from(isolated_nodes)

# Add correlation as edge attribute
for u, v, d in G.edges(data=True):
    d['correlation'] = G[u][v]['weight']
    d['weight'] = abs(d['weight'])

# Get the node labels (from the columns of r_df)
node_labels = {i: label for i, label in enumerate(r_df.columns)}

# Plot the graph with node labels
plt.figure(figsize=(10, 10))
pos = nx.spring_layout(G)  # Layout for visualization

# Get the list of nodes in the graph (after removing isolated ones)
nodes = list(G.nodes())

# Create a mapping from the current nodes to the original labels
current_labels = {node: node_labels[node] for node in nodes}

# Draw the graph with node labels
nx.draw(G, pos, with_labels=True, labels=current_labels, 
        node_size=500, font_size=12, node_color='skyblue', font_weight='bold', edge_color='gray')
plt.title('Network')
plt.show()

# Save the graph in GraphML format
nx.write_graphml(G, 'network_90.graphml')

labels_df = pd.DataFrame(list(current_labels.items()), columns=['Node Index', 'Label'])
labels_df.to_csv('node_labels_90.csv', index=False)
