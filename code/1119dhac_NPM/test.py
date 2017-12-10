import networkx as nx
G = nx.DiGraph()
G.add_edge(1,2)
print G.number_of_nodes()
for node in G.nodes():
    print node
    print nx.shortest_path_length(G, 1,node)