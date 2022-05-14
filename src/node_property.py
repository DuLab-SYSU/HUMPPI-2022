import networkx as nx
import pandas as pd
import os, sys
G = nx.read_weighted_edgelist(sys.argv[1], delimiter='\t')

nodes = sorted(list(G.nodes()))

Degree = nx.degree(G)  # Degree
k = [Degree[node] for node in nodes]

Betweenness_centrality = nx.betweenness_centrality(G)  # Betweenness
BC = [Betweenness_centrality[node] for node in nodes]

Eigenvector_centrality = nx.eigenvector_centrality(G)  # Eigenvector centrality
x = [Eigenvector_centrality[node] for node in nodes]

Clustering_coefficient = nx.clustering(G) # local clustering coefficient
C = [Clustering_coefficient[node] for node in nodes]

Assortativity = nx.average_neighbor_degree(G) # The average neighborhood degree of a node
NC = [Assortativity[node] for node in nodes]

Closeness_centrality = nx.closeness_centrality(G) # reciprocal of the average shortest path distance
SP = [Closeness_centrality[node] for node in nodes]

if sys.argv[2] == 'gene':
    node_info = pd.DataFrame({'Gene': nodes,
                              'Degree': k, 'Betweenness_centrality': BC,
                              'Eigenvector_centrality': x,
                              'Clustering_coefficient': C,
                              'Assortativity': NC,
                              'Closeness_centrality': SP})

    node_info.to_csv(path_or_buf='../intermediate/Gene_centrality.txt', index=False, sep='\t')
else:
    node_info = pd.DataFrame({'Module': nodes,
                              'Degree': k, 'Betweenness_centrality': BC,
                              'Eigenvector_centrality': x,
                              'Clustering_coefficient': C,
                              'Assortativity': NC,
                              'Closeness_centrality': SP})

    node_info.to_csv(path_or_buf='../intermediate/Module_centrality.txt', index = False, sep = '\t')
