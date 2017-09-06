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


def tree0(weight_value, startwindow):

    # 2017.8.21

    windowGraph = {}
    cliqueGraph = nx.DiGraph()
    dic_term = {}
    dic_last_time = {}
    dic_temp = {}
    dic_term_num = {}
    dic_intersect_level = {}
    term = 183
    root = 0
    cliqueGraph.add_node(root, annotation='root')
    w = data.shape[1]
    i = 0
    q = 0
    for window in range(startwindow, w):
        windowGraph[window] = nx.Graph()
        df = data[data[data.columns[window]] >= (weight_value + 0.00001)]
        # print weight_value, df.shape
        for edge in range(0, df.shape[0]):
            node_1, node_2 = df.index[edge].split('_')
            windowGraph[window].add_edge(node_1, node_2)
        # print windowGraph[window].size()
        dic_intersect_level.clear()

        if window == startwindow:
            for clique in nx.find_cliques(windowGraph[window]):
                if len(clique) > 5:
                    cliqueGraph.add_node(term, annotation=list(clique), windowsize=[window],
                                         weight=weight_value)  # generate a term
                    cliqueGraph.add_edge(root, term)
                    dic_term[frozenset(clique)] = [window]  # dic_term 记录 window和clique
                    dic_term_num[frozenset(clique)] = term  # dic_term_num 记录 term 序号和clique
                    dic_last_time[frozenset(clique)] = [window]  # dic_last_time   记录上一时刻生成的交集 用于下一时刻的比较
                    term = term + 1
                else:
                    continue
            # print len(dic_last_time), len(dic_term), cliqueGraph.number_of_nodes()
        else:
            for clique in nx.find_cliques(windowGraph[window]):
                if len(clique) > 5:
                    for key, value in dic_last_time.items():  # key 是clique ,value是 [window]
                        intersect = sorted(set(key).intersection(set(clique)))
                        q = 0
                        if len(intersect) >= 5:
                            # 同一层判断交集之间是否有重复的父子关系。 每生成一个交集， 判断当前层的其他term和交集的关系。
                            for ik, iv in dic_intersect_level.items():
                                if set(intersect) == (set(ik)):  # 生成一模一样的交集
                                    if dic_term_num[frozenset(key)] != dic_term_num[frozenset(ik)]:
                                        cliqueGraph.add_edge(dic_term_num[frozenset(key)], dic_term_num[frozenset(ik)])
                                    q = 1
                                    break
                                elif set(intersect).issuperset(set(ik)):  # 生成了超集
                                    cliqueGraph.remove_node(dic_term_num[frozenset(ik)])
                                    dic_term.pop(frozenset(ik))  # 从四个字典中都删除该节点的信息
                                    dic_term_num.pop(frozenset(ik))
                                    dic_intersect_level.pop(frozenset(ik))
                                    dic_temp.pop(frozenset(ik))
                                elif set(intersect).issubset(set(ik)):  # 生成了子集
                                    q = 1
                                    break
                            if q == 1:
                                continue
                            dic_intersect_level[frozenset(intersect)] = 1

                            if dic_term.has_key(frozenset(intersect)):
                                # 交集已经出现过
                                parent = cliqueGraph.predecessors(dic_term_num[frozenset(intersect)])
                                children = cliqueGraph.successors(dic_term_num[frozenset(intersect)])
                                if len(parent) > 0:
                                    # 是交集生成的term，则重定向
                                    cliqueGraph.add_node(term, annotation=list(intersect),
                                                         windowsize=value + [window],
                                                         weight=weight_value)
                                    for p in parent:
                                        cliqueGraph.add_edge(p, term)  # 连边

                                    for c in children:
                                        cliqueGraph.add_edge(term, c)  # 连边
                                    cliqueGraph.remove_node(dic_term_num[frozenset(intersect)])  # 从图中删除冗余结点

                                    #print 'deleted intersect nodes:',dic_term_num[frozenset(intersect)]
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
                                #print 'new term intersect never appear:', term
                                cliqueGraph.add_node(term, annotation=list(intersect), windowsize=value + [window],
                                                     weight=weight_value)  # generate a term

                                cliqueGraph.add_edge(dic_term_num[frozenset(key)], term)  # 连边，变化：只连接交集作为父亲。
                                dic_term[frozenset(intersect)] = value + [window]  # 新节点插入字典
                                dic_term_num[frozenset(intersect)] = term
                                dic_temp[frozenset(intersect)] = value + [window]  # 记录到dic_temp里
                                term = term + 1
                        else:
                            continue
                else:
                    continue
            dic_last_time.clear()
            for key, value in dic_temp.items():
                dic_last_time[key] = value
            dic_temp.clear()
    print cliqueGraph.number_of_nodes(), cliqueGraph.number_of_edges()
    print 'deleted nodes:', i


# 0905
    for startwindow in range(1, w):
        for window in range(startwindow, w):
            dic_intersect_level.clear()
            if window == startwindow:
                for clique in nx.find_cliques(windowGraph[window]):
                    if len(clique) > 5:
                        if dic_term.has_key(frozenset(clique)) and set([startwindow]).issubset(dic_term[frozenset(clique)]):
                            continue
                        else:
                            p = 0
                            c = 0
                            cliqueGraph.add_node(term, annotation=list(clique), windowsize=[window],
                                                 weight=weight_value)  # generate a term
                            dic_term[frozenset(clique)] = [window]  # dic_term 记录 window和clique
                            dic_term_num[frozenset(clique)] = term  # dic_term_num 记录 term 序号和clique
                            dic_last_time[frozenset(clique)] = [window]  # dic_last_time   记录上一时刻生成的交集 用于下一时刻的比较

                            # find parent and child
                            for geneSet, windowSet in dic_term.items():
                                if set(geneSet).issupset(set(intersect)) and set(windowSet).issubset(set(intersect)):
                                    cliqueGraph.add_edge(dic_term_num[geneSet], term)
                                    p = 1
                                if set(geneSet).issubset(set(intersect)) and set(windowSet).issupset(set(intersect)):
                                    cliqueGraph.add_edge(term, dic_term_num[geneSet])
                            # if no parents, link to the root node.
                            if p == 0:
                                cliqueGraph.add_edge(root, term)
                                print 'link to root :', term
                            term = term + 1
                    else:
                        continue
            else:
                for clique in nx.find_cliques(windowGraph[window]):
                    if len(clique) > 5:
                        for key, value in dic_last_time.items():  # key 是clique ,value是 [window]
                            intersect = sorted(set(key).intersection(set(clique)))
                            q = 0
                            if len(intersect) >= 5:
                                if dic_term.has_key(frozenset(intersect)) and set(value+[window]).issubset(dic_term[frozenset(clique)]):
                                    continue
                                # 同一层判断交集之间是否有重复的父子关系。 每生成一个交集， 判断当前层的其他term和交集的关系。
                                for ik, iv in dic_intersect_level.items():
                                    if set(intersect) == (set(ik)):  # 生成一模一样的交集
                                        if dic_term_num[frozenset(key)] != dic_term_num[frozenset(ik)]:
                                            cliqueGraph.add_edge(dic_term_num[frozenset(key)],
                                                                 dic_term_num[frozenset(ik)])
                                        q = 1
                                        break
                                    elif set(intersect).issuperset(set(ik)):  # 生成了超集
                                        cliqueGraph.remove_node(dic_term_num[frozenset(ik)])
                                        dic_term.pop(frozenset(ik))  # 从四个字典中都删除该节点的信息
                                        dic_term_num.pop(frozenset(ik))
                                        dic_intersect_level.pop(frozenset(ik))
                                        dic_temp.pop(frozenset(ik))
                                    elif set(intersect).issubset(set(ik)):  # 生成了子集
                                        q = 1
                                        break
                                if q == 1:
                                    continue

                                dic_intersect_level[frozenset(intersect)] = 1

                                if dic_term.has_key(frozenset(intersect)):
                                    # 交集已经出现过
                                    parent = cliqueGraph.predecessors(dic_term_num[frozenset(intersect)])
                                    children = cliqueGraph.successors(dic_term_num[frozenset(intersect)])
                                    if len(parent) > 0:
                                        # 是交集生成的term，则重定向
                                        cliqueGraph.add_node(term, annotation=list(intersect),
                                                             windowsize=value + [window],
                                                             weight=weight_value)
                                        for p in parent:
                                            cliqueGraph.add_edge(p, term)  # 连边

                                        for c in children:
                                            cliqueGraph.add_edge(term, c)  # 连边
                                        cliqueGraph.remove_node(dic_term_num[frozenset(intersect)])  # 从图中删除冗余结点

                                        # print 'deleted intersect nodes:',dic_term_num[frozenset(intersect)]
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
                                    # print 'new term intersect never appear:', term
                                    cliqueGraph.add_node(term, annotation=list(intersect),
                                                         windowsize=value + [window],
                                                         weight=weight_value)  # generate a term

                                    cliqueGraph.add_edge(dic_term_num[frozenset(key)], term)  # 连边，变化：只连接交集作为父亲。
                                    dic_term[frozenset(intersect)] = value + [window]  # 新节点插入字典
                                    dic_term_num[frozenset(intersect)] = term
                                    dic_temp[frozenset(intersect)] = value + [window]  # 记录到dic_temp里
                                    term = term + 1
                            else:
                                continue
                    else:
                        continue
                dic_last_time.clear()
                for key, value in dic_temp.items():
                    dic_last_time[key] = value
                dic_temp.clear()
            print cliqueGraph.number_of_nodes(), cliqueGraph.number_of_edges()
            print 'deleted nodes:', i




    # fw = open('0904edges_remove.txt', 'w')
    # fw2 = open('0904terms_remove.txt', 'w')
    # fw.write('parent' + '\t' + 'child' + '\n')
    # for edge in cliqueGraph.edges():
    #     fw.write(str(edge[0]) + '\t' + str(edge[1]) + '\n')
    # fw.close()
    # fw2.write('term_id' + '\t' + 'anno_genes' + '\t' + 'window' + '\t' + 'gene_size' + '\t' + 'window_size' + '\n')
    # for key, value in dic_term.items():
    #     fw2.write(str(dic_term_num[key]) + '\t' + str(key) + '\t' + str(value) + '\t' + str(len(key)) + '\t' + str(len(value)) + '\n')
    # fw2.close()
    # for nodes in cliqueGraph.nodes():
    #     if cliqueGraph.degree(nodes) == 0:
    #         print nodes
    # 20170905
    return cliqueGraph






if __name__ == '__main__':
    start = time.clock()
    s = 0
    cliqueGraph0 = tree0(0.9, 0)
    cliqueGraph1 = tree0(cliqueGraph0, 1)
    for node in cliqueGraph2.nodes():
        if node in cliqueGraph1.nodes():
            s = s + 1
    print s

    end = time.clock()
    print 'The function run time is : %.03f seconds' % (end - start)
