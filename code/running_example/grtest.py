#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
@Time    : 2018/3/24 19:49
@Author  : Junya Lu
@Site    :
"""

import networkx as nx
import pandas as pd
import matplotlib.pyplot as plt
G1 = nx.DiGraph()
G1.add_edge(1,2)
G1.add_edge(1,5)
G1.add_edge(5,7)
G1.node[1]['a'] = 9
G1.node[2]['a'] = 4
G1.node[5]['a'] = 19
G1.node[7]['a'] = 5
G1.add_edge(1,7)
G_copy = G1.copy()
print G_copy.node[1]['a']
# nx.draw(G1)
# plt.show()
# nx.draw(G_copy)
# plt.show()
# G_copy.remove_node(5)
# nx.draw(G_copy)
# plt.show()
distance=nx.shortest_path_length(G1,source=1)

print distance

df = pd.DataFrame(distance.items())
print df
print df[1].max()
print df[0].max()
grouped = df.groupby(1, sort=False)

for g in grouped:
    print g[0]
    print g[1], type(g[1])
    for i in g[1][0].tolist():
        print i , type(i)

p=5
print type(p)
# for depth in range (df[1].min(), df[1].max()+1):
#     groupdf[df[1]==depth]
