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
pos = nx.spring_layout(G1)
nx.draw_networkx(G1, pos, with_labels=True, node_size=10)
plt.show()