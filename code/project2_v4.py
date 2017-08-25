# -*- coding: utf-8 -*-
"""
Created on Mon Jul 31 16:37:37 2017

@author: lu
"""

import pandas as pd
import numpy as np
import networkx as nx
import time
import os
from os.path import join
import matplotlib.pyplot as plt
import gevent
import gevent.pool
import multiprocessing
from multiprocessing import Pool
import os, time, random

filename = 'result_c5_s10_v2_weight.txt'
data = pd.read_csv(filename, index_col=0, sep='\t')
curr_path = os.getcwd() + '/'
fw = open('0824edges.txt', 'w')
fw2 = open('0824terms.txt', 'w')

def tree(weight_value):

    # 2017.8.21

    windowGraph = {}
    cliqueGraph = nx.DiGraph()
    dic_term = {}
    dic_last_time = {}
    dic_temp = {}
    dic_term_num = {}
    term = 183
    w = data.shape[1]
    i = 0
    for window in range(0, w):


        windowGraph[window] = nx.Graph()
        df = data[data[data.columns[window]] >= (weight_value + 0.00001)]
        print weight_value, df.shape
        for edge in range(0, df.shape[0]):
            node_1, node_2 = df.index[edge].split('_')
            windowGraph[window].add_edge(node_1, node_2)
        print windowGraph[window].size()
        if window == 0:
            print 'window:', window, nx.graph_number_of_cliques(windowGraph[window]),len([clique for clique in nx.find_cliques(windowGraph[window]) if len(clique) >5])
            for clique in nx.find_cliques(windowGraph[window]):
                if len(clique) > 5:
                    cliqueGraph.add_node(term, annotation=list(clique), windowsize=[window],
                                         weight=weight_value)  # generate a term
                    dic_term[frozenset(clique)] = [window]  # dic_term 记录 window和clique
                    dic_term_num[frozenset(clique)] = term  # dic_term_num 记录 term 序号和clique
                    dic_last_time[frozenset(clique)] = [window]  # dic_last_time   记录上一时刻生成的交集 用于下一时刻的比较
                    term = term + 1
                else:
                    continue
            print len(dic_last_time), len(dic_term), cliqueGraph.number_of_nodes()
        else:
            print 'window', window, nx.graph_number_of_cliques(windowGraph[window]),len([clique for clique in nx.find_cliques(windowGraph[window]) if len(clique) >5])
            print len(dic_last_time), len(dic_term), cliqueGraph.number_of_nodes()
            for clique in nx.find_cliques(windowGraph[window]):
                if len(clique) > 5:
                    if dic_term.has_key(frozenset(clique)):
                        # 如果下一时刻的clique在上一时刻生成过， 则更新clique的window
                        # dic_term[frozenset(clique)] = dic_term[frozenset(clique)] + [window]
                        # cliqueGraph.node[dic_term_num[frozenset(clique)]]['windowsize'] = dic_term[frozenset(clique)]
                        continue
                    else:
                        # 如果新生成了一些没有出现过的clique，则新生成一个term
                        dic_term[frozenset(clique)] = [window]  # dic_term 记录 clique-window
                        dic_term_num[frozenset(clique)] = term  # dic_term_num 记录 term 序号和clique
                        cliqueGraph.add_node(term, annotation=list(clique), windowsize=[window],
                                             weight=weight_value)  # generate a term
                        term = term + 1
                        # 求交集

                        for key, value in dic_last_time.items():  # key 是clique ,value是 [window]

                            intersect = sorted(set(key).intersection(set(clique)))

                            if len(intersect) >= 5:

                                if dic_term.has_key(frozenset(intersect)):
                                    # 交集已经出现过
                                    parent = cliqueGraph.predecessors(dic_term_num[frozenset(intersect)])
                                    children = cliqueGraph.successors(dic_term_num[frozenset(intersect)])

                                    if len(parent) > 1:
                                        # 是交集生成的term，则重定向
                                        print 'parent:', parent
                                        cliqueGraph.add_node(term, annotation=list(intersect),
                                                             windowsize=value + [window],
                                                             weight=weight_value)
                                        for p in parent:
                                            cliqueGraph.add_edge(p, term)  # 连边


                                        cliqueGraph.remove_node(dic_term_num[frozenset(intersect)])  # 从图中删除冗余结点
                                        i = i + 1
                                        dic_term.pop(frozenset(intersect))  # 字典中删除
                                        dic_term_num.pop(frozenset(intersect))
                                        dic_term[frozenset(intersect)] = value + [window]  # 新节点插入字典
                                        dic_term_num[frozenset(intersect)] = term

                                        dic_temp[frozenset(intersect)] = value + [window]  # 记录到dic_temp里
                                        term = term + 1
                                        continue
                                    else:
                                        # 是window生成的term
                                        continue

                                else:
                                    # 交集没有出现过， 则生成新的term
                                    print 'new term intersect never appear:', term, dic_term_num[frozenset(key)],dic_term_num[frozenset(clique)]
                                    cliqueGraph.add_node(term, annotation=list(intersect), windowsize=value + [window],
                                                         weight=weight_value)  # generate a term
                                    cliqueGraph.add_edges_from([(dic_term_num[frozenset(key)], term),
                                                                (dic_term_num[frozenset(clique)], term)])  # 连边
                                    dic_term[frozenset(intersect)] = value + [window]  # 新节点插入字典
                                    dic_term_num[frozenset(intersect)] = term
                                    dic_temp[frozenset(intersect)] = value + [window]  # 记录到dic_temp里
                                    term = term + 1
                            else:
                                continue
                else:
                    continue
            # print len(dic_last_time), len(dic_temp)
            # dic_last_time.clear()
            # dic_last_time = dic_temp
            # dic_temp.clear()
            # print len(dic_last_time), len(dic_temp)
            dic_last_time.clear()
            for key, value in dic_temp.items():
                dic_last_time[key] = value
            dic_temp.clear()
            print len(dic_last_time), len(dic_temp)
    print cliqueGraph.number_of_nodes(), cliqueGraph.number_of_edges()
    print 'deleted nodes:', i
    fw.write('parent' + '\t' + 'child' + '\n')
    for edge in cliqueGraph.edges():
        fw.write(str(edge[0]) + '\t' + str(edge[1]) + '\n')
    fw.close()
    fw2.write('term_id' + '\t' + 'anno_genes' + '\t' + 'window' + '\n')
    for key, value in dic_term.items():
        fw2.write(str(dic_term_num[key]) + '\t' + str(key) + '\t' + str(value) + '\n')
    fw2.close()

if __name__ == '__main__':
    start = time.clock()
    tree(0.9)
    end = time.clock()
    print 'The function run time is : %.03f seconds' % (end - start)
