#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
@Time    : 2018/4/6 15:41
@Author  : Junya Lu
@Site    : 
"""

import networkx as nx
from compiler.ast import flatten
import pandas as pd

import math
if __name__ =='__main__':
    phe1 = pd.read_table('G:\project2\\NPM201507\\data\\IDMapping_consolidated_allPhi2_cleaned_lfc_avg.txt',
                         index_col=0)
    phe2 = pd.read_table('G:\project2\\NPM201507\\data\\IDMapping_consolidated_allQESV_cleaned_LFC_avg.txt',
                         index_col=0)
    phe3 = pd.read_table('G:\project2\\NPM201507\\data\\IDMapping_consolidated_allQI_new_RAW3_adj_LFC_avg.txt',
                         index_col=0)

    phe1.columns = [i for i in range(0, 113)]
    phe2.columns = [i for i in range(0, 113)]
    phe3.columns = [i for i in range(0, 113)]
    p1_list = flatten(phe1.values.tolist())
    p2_list = flatten(phe2.values.tolist())
    p3_list = flatten(phe3.values.tolist())
    #2018.04.06
    p1_list = sorted(list([value for value in p1_list if not math.isnan(value)]), reverse=True)
    p2_list = sorted(list([value for value in p2_list if not math.isnan(value)]), reverse=True)
    p3_list = sorted(list([value for value in p3_list if not math.isnan(value)]), reverse=True)
    print len(p1_list), len(p2_list), len(p3_list)
    print max(p3_list), p3_list[0]

    top_r = int(len(p1_list) * 0.2) - 1
    botoom_r = len(p1_list) - int(len(p1_list) * 0.2)
    u1 = p1_list[top_r]
    b1 = p1_list[botoom_r]

    top_r = int(len(p2_list) * 0.2) - 1
    botoom_r = len(p2_list) - int(len(p2_list) * 0.2)

    u2 = p2_list[top_r]
    b2 = p2_list[botoom_r]

    top_r = int(len(p3_list) * 0.2) - 1
    botoom_r = len(p3_list) - int(len(p3_list) * 0.2)
    u3 = p3_list[top_r]
    b3 = p3_list[botoom_r]

    print 'threshold for phenotype:', u1, b1, u2, b2, u3, b3

G = nx.DiGraph()
G.add_edges_from([(1,2),(1,3),(2,4),(4,6),(1,7)])
print max(nx.shortest_path_length(G, 1).values())


# edgeFile = '0405edges_sign_id1220.txt'
# fr = open(edgeFile, 'r')
# line = fr.readline()
# dic = {}
# for line in fr.readlines():
#     p,c = line.strip().split('\t')
#     if dic.has_key((c,p)):
#         print line
#     else:
#         dic[(p,c)] = 1
# print len(dic)

