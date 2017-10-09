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





filename = 'result_c5_s10_v2_weight.txt'
data = pd.read_csv(filename, index_col=0, sep='\t')
curr_path = os.getcwd() + '/'


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
    plt.subplot(221)
    ax = sns.heatmap(q, vmin=-0.07, vmax=0.07, cmap=cmap,  center=0, square=True, linewidths=.5,
                     cbar_kws={"shrink": .5})
    for item in ax.get_yticklabels():
        item.set_rotation(0)

    plt.subplot(222)
    ax = sns.heatmap(q2, vmin=-0.07, vmax=0.07, cmap=cmap, center=0, square=True, linewidths=.5,
                     cbar_kws={"shrink": .5})
    for item in ax.get_yticklabels():
        item.set_rotation(0)

    plt.show()


if __name__ == '__main__':


    start = time.clock()

    gene = ['AT1G12250', 'AT1G55370', 'AT1G76760', 'AT1G80030', 'AT2G27290', 'AT2G30950', 'AT4G04020', 'AT4G39960', 'AT5G19370', 'AT5G20140', 'AT5G23890']
    window = [6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19]
    gene1 = ['AT1G31180', 'AT1G54050', 'AT1G66670', 'AT1G74710', 'AT3G26580', 'AT4G01050', 'AT4G01900', 'AT4G31560', 'AT5G42765', 'AT5G62720']
    window1 = [18, 19, 20, 21, 22, 23, 24]
    gene2 = ['AT1G02560', 'AT1G54050', 'AT1G72640', 'AT1G74710', 'AT2G17630', 'AT2G26800', 'AT2G37660', 'AT5G13120', 'AT5G39830', 'AT5G42270']
    window2 = [0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10]
    gene3 = ['AT1G22020', 'AT1G51110', 'AT3G15840', 'AT3G25480', 'AT5G16150']
    window3 = [8, 9, 10, 11, 12, 13, 14]
    gene4 = ['AT1G10170', 'AT1G23740', 'AT1G50250', 'AT3G15360', 'AT4G13010', 'AT4G39460', 'AT5G16660']
    window4 = [8, 9, 10, 11, 12, 13, 14]
    gene5 =['AT1G10170', 'AT1G29920', 'AT1G30510', 'AT1G79230', 'AT3G22690', 'AT4G20820']
    window5 = [14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24, 25, 26, 27, 28, 29, 30, 31, 32, 33, 34, 35, 36, 37]

    phenotype(gene, window)
    phenotype(gene1, window1)
    phenotype(gene2, window2)
    phenotype(gene3, window3)
    phenotype(gene4, window4)
    phenotype(gene5, window5)



    '''
    fw1 = open('0908edges_10.txt', 'w')
    fw2 = open('0908terms_10.txt', 'w')
    s = 0
    # 先产生第一个window的 tree
    term = 183
    cliqueGraph0, dic_term0, dic_term_num0, term = tree0(0.9, 0, 183)
    dic_all = {}
    dic_all = dic_term0.copy()
    copy_clique = cliqueGraph0
    for i in range(1, 10):
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
    fw1.write('parent' + '\t' + 'child' + '\n')
    for edge in cliqueGraph0.edges():
        fw1.write(str(edge[0]) + '\t' + str(edge[1]) + '\n')
    fw2.write('term_id' + '\t' + 'annotation_gene' + '\t' + 'annotation_window' + '\t' + 'geneSize' + '\t' + 'window_size' + '\n')
    for node in cliqueGraph0.nodes():
        fw2.write(str(node) + '\t' + str(cliqueGraph0.node[node]['annotation']) + '\t' + str(cliqueGraph0.node[node]['windowsize']) + '\t' +
                  str(len(cliqueGraph0.node[node]['annotation'])) + '\t' + str(len(cliqueGraph0.node[node]['windowsize'])) +  '\n')
    fw1.close()
    fw2.close()
    '''

    end = time.clock()
    print 'The function run time is : %.03f seconds' % (end - start)
