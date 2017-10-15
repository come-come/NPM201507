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
import matplotlib.pyplot as plt
import seaborn as sns
import sys
from scipy import stats
from compiler.ast import flatten





filename = 'result_c5_s10_v2_weight.txt'
data = pd.read_csv(filename, index_col=0, sep='\t')
curr_path = os.getcwd() + '/'

phe1 = pd.read_table('G:\project2\\NPM201507\\data\\IDMapping_consolidated_allPhi2_cleaned_lfc_avg.txt', index_col=0)
phe2 = pd.read_table('G:\project2\\NPM201507\\data\\IDMapping_consolidated_allQESV_cleaned_LFC_avg.txt',index_col=0)
phe3 = pd.read_table('G:\project2\\NPM201507\\data\\IDMapping_consolidated_allQI_new_RAW3_adj_LFC_avg.txt', index_col=0)

phe1.columns = [i for i in range(0, 113)]
phe2.columns = [i for i in range(0, 113)]
phe3.columns = [i for i in range(0, 113)]

def tree0(weight_value, startwindow, term):

    # 2017.8.21

    windowGraph = {}
    cliqueGraph = nx.DiGraph()
    dic_term = {}
    dic_last_time = {}
    dic_temp = {}
    dic_term_num = {}
    dic_intersect_level = {}
    # term = 183
    root = 0
    cliqueGraph.add_node(root, annotation='root', windowsize='root', weight_value='root')
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
    print 'window', startwindow, 'size is', cliqueGraph.number_of_nodes(), cliqueGraph.number_of_edges()
    # print 'deleted nodes:', i
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
    return cliqueGraph,dic_term,dic_term_num, term



def phenotype(gene_set, window_set):
    phe1 = pd.read_table('G:\project2\\NPM201507\\data\\IDMapping_consolidated_allPhi2_cleaned_lfc_avg.txt', index_col=0)
    phe2 = pd.read_table('G:\project2\\NPM201507\\data\\IDMapping_consolidated_allQESV_cleaned_LFC_avg.txt',index_col=0)
    phe3 = pd.read_table('G:\project2\\NPM201507\\data\\IDMapping_consolidated_allQI_new_RAW3_adj_LFC_avg.txt', index_col=0)
    t_min = min(window_set)
    t_max = max(window_set) + 10
    # rename columns index. [0, 1, 2...]
    phe1.columns = [i for i in range(0, 113)]
    phe2.columns = [i for i in range(0, 113)]
    phe3.columns = [i for i in range(0, 113)]
    print phe1.loc[gene_set, window_set]
    print phe2.loc[gene_set, window_set]
    print phe3.loc[gene_set, window_set]
    q = phe1.loc[gene_set, window_set]
    q2 = phe2.loc[gene_set, window_set]
    q3 = phe3.loc[gene_set, window_set]
    phe_avg = (phe1.loc[gene_set, window_set] + phe2.loc[gene_set, window_set] + phe3.loc[gene_set, window_set]).apply(lambda x : x*2)
    print 'avg'
    print phe_avg
    # plt.pcolor(phe_avg)
    # plt.yticks(np.arange(0.5, len(phe_avg.index), 1), phe_avg.index)
    # plt.xticks(np.arange(0.5, len(phe_avg.columns), 1), phe_avg.columns)
    # plt.show()
    sns.set(style="white")
    cmap = sns.diverging_palette(220, 10, as_cmap=True)
    # ax = sns.heatmap(q, vmin=-0.05, vmax=0.05, cmap=cmap, annot=True,center=0, square=True, linewidths=.5, cbar_kws={"shrink": .5})
    # for item in ax.get_yticklabels():
    #     item.set_rotation(0)

    fig = plt.figure()
    plt.subplot(231)
    ax = sns.heatmap(q, vmin=-0.5, vmax=0.5, cmap=cmap,  center=0, square=True, linewidths=.5,
                     cbar_kws={"shrink": .5})
    for item in ax.get_yticklabels():
        item.set_rotation(0)
    plt.title('phenotype 1')
    plt.subplot(232)
    ax = sns.heatmap(q2, vmin=-0.5, vmax=0.5, cmap=cmap, center=0, square=True, linewidths=.5,
                     cbar_kws={"shrink": .5})
    for item in ax.get_yticklabels():
        item.set_rotation(0)
    plt.title('phenotype 2')

    plt.subplot(233)
    ax = sns.heatmap(q3, vmin=0, vmax=1, cmap=cmap, center=0, square=True, linewidths=.5,
                     cbar_kws={"shrink": .5})
    for item in ax.get_yticklabels():
        item.set_rotation(0)
    plt.title('phenotype 3')
    plt.show()

def sign_value(node, gene_set, window_set):
    t_min = min(window_set)
    t_max = max(window_set) + 10
    p1 = phe1.loc[gene_set, window_set]
    p2 = phe2.loc[gene_set, window_set]
    p3 = phe3.loc[gene_set, window_set]
    p1_list = flatten(p1.values.tolist())
    p2_list = flatten(p2.values.tolist())
    p3_list = flatten(p3.values.tolist())
    sign_p1 =  [x for x in p1_list if x<-0.18 or x >0.08]
    sign_p2 =  [x for x in p2_list if x<-0.22 or x >0.27]
    sign_p3 =  [x for x in p3_list if x>0.34]
    sign = len(sign_p1) + len(sign_p2) + len(sign_p3)
    all = len(p1_list) + len(p2_list) + len(p3_list)
    return sign, float(sign)/all

def class_distance(node, gene_set, window_set):
    p1 = phe1.loc[gene_set, window_set].dropna(axis=1,how='any')
    p2 = phe2.loc[gene_set, window_set].dropna(axis=1,how='any')
    p3 = phe3.loc[gene_set, window_set].dropna(axis=1,how='any')
    result = pd.concat([p1, p2, p3], axis=1)
    # print result.head(5)
    distance = np.trace(np.cov(result))
    return distance


if __name__ == '__main__':


    start = time.clock()

    fw1 = open('1015edges_sign.txt', 'w')
    fw2 = open('1015terms_sign.txt', 'w')
    fw3 = open('1015sign_score.txt', 'w')
    s = 0
    # 先产生第一个window的 tree
    term = 183
    cliqueGraph0, dic_term0, dic_term_num0, term = tree0(0.9, 0, 183)
    dic_all = {}
    dic_all = dic_term0.copy()
    copy_clique = cliqueGraph0
    for i in range(1, 50):
        print 'begin term num:', term
        cliqueGraph1, dic_term1, dic_term_num1, term = tree0(0.9, i, term)

        for key in dic_term1.keys():
            if dic_all.has_key(key):
                if set(dic_all[key]).issuperset(set(dic_term1[key])):
                    dic_term1.pop(key)
                    num = dic_term_num1[key]
                    cliqueGraph1.remove_node(num)
                    s = s + 1
                else:
                    dic_all[key] = dic_all[key] + dic_term1[key]
                    dic_term1.pop(key)
                    # num = dic_term_num1[key]
                    # cliqueGraph1.remove_node(num)
                    # else:
                    #     dic_all[key] = dic_term1[key]
        dic_all.update(dic_term1)
        cliqueGraph0 = nx.compose(cliqueGraph0, cliqueGraph1)
        print 'window', i, cliqueGraph0.number_of_nodes(), cliqueGraph0.number_of_edges()
    dic_term_score = {}
    dic_term_distance = {}
    for node in cliqueGraph0.nodes():
        if node == 0:
            continue
        else:
            gene_set = cliqueGraph0.node[node]['annotation']
            window_set = cliqueGraph0.node[node]['windowsize']
            # 判断phenotype是否有意义
            sign, score = sign_value(node, gene_set, window_set)
            dis = class_distance(node, gene_set, window_set)
            if sign == 0:
                parent = cliqueGraph0.predecessors(node)
                child = cliqueGraph0.successors(node)
                for p in parent:
                    for c in child:
                        cliqueGraph0.add_edge(p, c)
                cliqueGraph0.remove_node(node)
                print node
                # 无意义，delete，重定向
            else:
                dic_term_score[node] = score
                dic_term_distance[node] = dis
                continue

    for key, value in sorted(dic_term_distance.items(), key=lambda d: d[1],reverse=True):
        fw3.write(str(key) + '\t' + str(value) + '\n')
    fw3.close()
    fw1.write('parent' + '\t' + 'child' + '\n')
    for edge in cliqueGraph0.edges():
        fw1.write(str(edge[0]) + '\t' + str(edge[1]) + '\n')
    fw2.write(
        'term_id' + '\t' + 'sign_score' + '\t' + 'annotation_gene' + '\t' + 'annotation_window' + '\t' + 'geneSize' + '\t' + 'window_size' + '\n')
    for node, value in sorted(dic_term_score.items(), key=lambda d: d[1],reverse=True):
        fw2.write(str(node) + '\t' + str(round(value,4)) + '\t' + str(cliqueGraph0.node[node]['annotation']) + '\t' + str(
            cliqueGraph0.node[node]['windowsize']) + '\t' +
                  str(len(cliqueGraph0.node[node]['annotation'])) + '\t' + str(
            len(cliqueGraph0.node[node]['windowsize'])) + '\n')
    fw1.close()
    fw2.close()


    '''
    gene = ['AT1G14345', 'AT1G80380', 'AT2G18790', 'AT2G29180', 'AT4G33010', 'AT5G42270']
    window = [19, 20, 21, 22, 23, 24]
    gene1 = ['AT1G03160', 'AT1G14345', 'AT1G80380', 'AT2G18790', 'AT2G29180', 'AT4G33010', 'AT5G42270']
    window1 = [19, 20, 21, 22, 23]
    gene2 = ['AT1G14150', 'AT1G18730', 'AT3G01440', 'AT3G15110', 'AT4G37925']
    window2 = [8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22]
    gene3 = ['AT3G01440', 'AT3G47060', 'AT4G23890', 'AT4G37925', 'AT5G62720']
    window3 = [24, 25, 26, 27, 28, 29, 30, 31, 32, 33, 34, 3                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                   5, 36, 37, 38, 39, 40, 41]


    phenotype(gene, window)
    phenotype(gene1, window1)
    phenotype(gene2, window2)
    phenotype(gene3, window3)
    '''





    end = time.clock()
    print 'The function run time is : %.03f seconds' % (end - start)
