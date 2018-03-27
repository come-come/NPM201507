#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
@Time    : 2018/3/24 19:49
@Author  : Junya Lu
@Site    :
"""

import networkx as nx
import matplotlib.pyplot as plt
G1 = nx.DiGraph()
G1.add_edge(1,2)
G1.add_edge(1,5)
G1.add_edge(5,7)
distance=nx.shortest_path_length(G1)
print distance

