import networkx as nx
G = nx.DiGraph()
G.add_node(184, annotation=[1, 2], windowsize=[5, 6], weight=0.9)  # generate a term
G.add_node(184, annotation=[1, 2], windowsize=[5, 6,7], weight=0.9)  # generate a term
G.node[184]['windowsize'] = [1,8,6]
G.add_edges_from([(1,184), (2,184)])
print G.node[184]
print G.predecessors(184)
print G.nodes()
print G.edges()
print G.nodes()
print G.edges()
for edge in G.edges():
    print edge[0], edge[1]


# dic = {'a':1,'b':2}
# print dic
# dic2 = dic
# print dic2