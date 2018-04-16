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
import operator
from compiler.ast import flatten
import sys, getopt
import math

# python project2.py -i result_c5_s10_20171220weight.txt -d IDMapping_consolidated_allPhi2_cleaned_lfc_avg.txt -f IDMapping_consolidated_allQESV_cleaned_LFC_avg.txt -g IDMapping_consolidated_allQI_new_RAW3_adj_LFC_avg.txt -w 0.9 -s 5 -p 0.5 -m 0.005 -n -0.18 -b 0.08 -v -0.22 -c -100 -z 0.34
phe1 = pd.read_table('G:\project2\\NPM201507\\data\\IDMapping_consolidated_allPhi2_cleaned_lfc_avg.txt', index_col=0)
phe2 = pd.read_table('G:\project2\\NPM201507\\data\\IDMapping_consolidated_allQESV_cleaned_LFC_avg.txt', index_col=0)
phe3 = pd.read_table('G:\project2\\NPM201507\\data\\IDMapping_consolidated_allQI_new_RAW3_adj_LFC_avg.txt', index_col=0)

phe1.columns = [i for i in range(0, 113)]
phe2.columns = [i for i in range(0, 113)]
phe3.columns = [i for i in range(0, 113)]
'''
p1_list = flatten(phe1.values.tolist())
p2_list = flatten(phe2.values.tolist())
p3_list = flatten(phe3.values.tolist())
p1_list = sorted(list([value for value in p1_list if not math.isnan(value)]), reverse=True)
p2_list = sorted(list([value for value in p2_list if not math.isnan(value)]), reverse=True)
p3_list = sorted(list([value for value in p3_list if not math.isnan(value)]), reverse=True)
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
print len(p1_list), len(p2_list), len(p3_list)
print 'threshold for phenotype:', u1, b1, u2, b2, u3, b3
'''

def tree_statis(weight_value, window):
    # output 50 trees
    filename2 = 'edge/edgeInfo' + str(window) + '.txt'
    f_path2 = curr_path + os.path.normpath(filename2)
    f_window = open(f_path2, 'w')
    df = data[data[data.columns[window]] >= (weight_value + 0.00001)]
    print df.shape
    for edge in range(0, df.shape[0]):
        node_1, node_2 = df.index[edge].split('_')
        f_window.write(node_1 + '\t' + node_2 + '\n')


def tree0(weight_value, startwindow, term):
    # 2017.8.21
    print 'start window:', startwindow
    # windowGraph = {}
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
                                    # 判断两个的编号是否一样？
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
    return cliqueGraph, dic_term, dic_term_num, term


def sign_value(node, gene_set, window):
    t_min = min(window) + 49
    t_max = max(window) + 49 + 10
    window_set = [i for i in range(t_min, t_max)]  # 回到最原始的数据上
    p1 = phe1.loc[gene_set, window_set]
    p2 = phe2.loc[gene_set, window_set]
    p3 = phe3.loc[gene_set, window_set]
    p1_list = flatten(p1.values.tolist())
    p2_list = flatten(p2.values.tolist())
    p3_list = flatten(p3.values.tolist())
    #2018.04.06

    # sign_p1 = [x for x in p1_list if x <= b1 or x >= u1]
    # sign_p2 = [x for x in p2_list if x <= b2 or x >= u2]
    # sign_p3 = [x for x in p3_list if x <= b3 or x >= u3]
    sign_p1 = [x for x in p1_list if x < -0.18 or x > 0.08]
    sign_p2 = [x for x in p2_list if x < -0.22 or x > 0.27]
    sign_p3 = [x for x in p3_list if x > 0.34]

    sign = len(sign_p1) + len(sign_p2) + len(sign_p3)
    all = len(p1_list) + len(p2_list) + len(p3_list)
    if p1[t_min].mean() - p1[t_max - 1].mean() < 0:
        trend1 = 1
    else:
        trend1 = 0
    if p2[t_min].mean() - p2[t_max - 1].mean() < 0:
        trend2 = 1
    else:
        trend2 = 0
    if p3[t_min].mean() - p3[t_max - 1].mean() < 0:
        trend3 = 1
    else:
        trend3 = 0
    purity = float(sign) / all
    purity1 = len(sign_p1) / float(len(p1_list))
    purity2 = len(sign_p2) / float(len(p2_list))
    purity3 = len(sign_p3) / float(len(p3_list))
    # if purity1 >= pp or purity2 >= pp or purity3 >= pp:
    #     purity = max(purity1, purity2, purity3)

    return sign, purity, trend1, trend2, trend3

if __name__ == '__main__':

    start = time.clock()

    s = 0

    term = 183
    pp = 0.025 # purity
    filename = 'result_c5_s10_v2_weight.txt'
    # filename = 'result_c5_s10_20171220weight.txt'
    # filename = 'result_c5_s10_20171219weight.txt'
    data = pd.read_csv(filename, index_col=0, sep='\t')
    curr_path = os.getcwd() + '/'
    weight_value = 0.9
    windowGraph = {}

    for window in range(0, 54):
        windowGraph[window] = nx.Graph()
        df = data[data[data.columns[window]] >= (weight_value + 0.00001)]

        for edge in range(0, df.shape[0]):
            node_1, node_2 = df.index[edge].split('_')
            windowGraph[window].add_edge(node_1, node_2)
    print windowGraph[0].number_of_nodes(), windowGraph[0].number_of_edges()
    cliqueGraph0, dic_term0, dic_term_num0, term = tree0(weight_value, 0, 183)


    dic_all = {}
    dic_all_term_num = dic_term_num0.copy()
    dic_all = dic_term0.copy()
    copy_clique = cliqueGraph0



    for i in range(1, 50):
        print 'begin term num:', term

        sum1 = 0
        cliqueGraph1, dic_term1, dic_term_num1, term = tree0(weight_value, i, term)
        # 先合并两个trees
        cliqueGraph0 = nx.compose(cliqueGraph0, cliqueGraph1)
        for key in dic_term1.keys():
            if dic_all.has_key(key):
                if set(dic_all[key]).issuperset(set(dic_term1[key])):
                    dic_term1.pop(key)
                    num = dic_term_num1[key]
                    cliqueGraph0.remove_node(num)
                    dic_term_num1.pop(key)
                    s = s + 1
                else:
                    dic_all[key] = dic_all[key] + dic_term1[key]
                    dic_term1.pop(key)
                    num = dic_term_num1[key]
                    cliqueGraph0.remove_node(num)
                    dic_term_num1.pop(key)
            else:
                # flag 连最上层子结点 子节点的孩子不需要再连接
                dic_this_edge = {}
                for old in dic_all.keys():  # 对之前已经存在的结点
                    old_id = dic_all_term_num[old]
                    this_id = dic_term_num1[key]
                    flag1 = 0
                    flag2 = 0
                    if set(key).issuperset(old) and set(dic_all[old]).issuperset(dic_term1[key]):
                        # 如果当前结点的注释基因是它的父集 注释时间是它的子集
                        # 判断old的parent是否和 key已经连接
                        try:
                            parent = cliqueGraph0.predecessors(old_id)
                        except:
                            continue
                        # parent 按照编号排序
                        dic_parent = {}
                        for p in parent:
                            dic_parent[p] = len(cliqueGraph0.node[p]['annotation'])
                        dic2 = sorted(dic_parent.items(), key=operator.itemgetter(1), reverse=True)
                        parent = []
                        for pair in dic2:
                            parent.append(pair[0])
                        for p_id in parent:
                            try:
                                p_anno = cliqueGraph0.node[p_id]['annotation']
                                p_wind = cliqueGraph0.node[p_id]['windowsize']
                                if set(key).issuperset(set(p_anno)) and set(p_wind).issuperset(dic_term1[key]):
                                    flag1 = 1
                                    break
                                    # print p_id, old, this_id
                                if dic_this_edge[(this_id, p_id)] == 1:
                                    flag1 = 1
                                    break
                                    # print p_id, old, this_id
                                # if flag1 == 1:
                                #     break
                            except:
                                continue
                        if flag1 == 0:
                            # key--> old
                            if this_id!=old_id and not cliqueGraph0.has_edge(this_id, old_id):
                                cliqueGraph0.add_edge(this_id, old_id)
                                dic_this_edge[(this_id, old_id)] = 1
                                sum1 = sum1 + 1
                                break
                            # print 'add edge', this_id, old_id
                        else:
                            # 已经连接到它的父亲上了
                            continue

                    elif set(key).issuperset(old) and set(dic_term1[key]).issuperset(dic_all[old]):
                        # 如果当前结点的注释基因是它的父集 注释时间也是它的父集  delete old

                        child = cliqueGraph0.successors(old_id)
                        for c in child:
                            if not cliqueGraph0.has_edge(this_id, c) and len(cliqueGraph0.node[c]['windowsize'])>len(dic_term1[key]):
                                cliqueGraph0.add_edge(this_id, c)
                        cliqueGraph0.remove_node(old_id)
                        dic_all.pop(old)
                        dic_all_term_num.pop(old)

                        # print 'remove node', old_id
                        ###
                    # 20180404
                    # elif set(old).issuperset(key) and set(dic_term1[key]).issuperset(dic_all[old]):
                    #     # 如果当前结点的注释基因是它的子集 注释时间是它的父集 old--> key
                    #     if this_id != old_id and not cliqueGraph0.has_edge(old_id, this_id):
                    #         cliqueGraph0.add_edge(old_id, this_id)
                    #     print 'add edge', old_id, this_id
                    #20180404
                    elif set(old).issuperset(key) and set(dic_all[old]).issuperset(dic_term1[key]):
                        # 如果当前结点的注释基因是它的子集 注释时间也是它的子集 delete this
                        child = cliqueGraph0.successors(this_id)
                        for c in child:
                            if not cliqueGraph0.has_edge(old_id, c) and len(cliqueGraph0.node[c]['windowsize']) > len(
                                    dic_all[old]):
                                cliqueGraph0.add_edge(old_id, c)

                        cliqueGraph0.remove_node(this_id)
                        dic_term_num1.pop(this_id)
                        dic_term1.pop(key)
                        # print 'remove node', this_id

        # print 'dic_this_edge', dic_this_edge
        # print dic_all_term_num
        # print dic_term_num1
        print 'sum1', sum1
        dic_all.update(dic_term1)
        dic_all_term_num.update(dic_term_num1)
    print 'raw size--------', cliqueGraph0.number_of_nodes(), cliqueGraph0.number_of_edges()


    print 'len(cliqueGraph0.selfloop_edges())-', cliqueGraph0.selfloop_edges(), nx.isolates(cliqueGraph0)

    print 'start calculate purity...'
    dic_term_score = {}
    delete_leave_flag1 = 0
    delete_sum = 0
    for node in cliqueGraph0.nodes():
        if node == 0:
            continue
        else:
            gene_set = cliqueGraph0.node[node]['annotation']
            window_set = cliqueGraph0.node[node]['windowsize']
            # 判断phenotype是否有意义
            sign, score, trend1, trend2, trend3 = sign_value(node, gene_set, window_set)

            if score < pp:
                # 无意义，delete，重定向
                parent = cliqueGraph0.predecessors(node)
                child = cliqueGraph0.successors(node)
                if len(child) == 0:
                    delete_leave_flag1 += 1
                for p in parent:
                    for c in child:
                        if not cliqueGraph0.has_edge(p, c):
                            cliqueGraph0.add_edge(p, c)
                cliqueGraph0.remove_node(node)
                delete_sum += 1
            else:
                dic_term_score[node] = score
                continue
    print 'after purity window', cliqueGraph0.number_of_nodes(), cliqueGraph0.number_of_edges()
    print 'len(cliqueGraph0.selfloop_edges())-----', len(cliqueGraph0.selfloop_edges()), nx.isolates(cliqueGraph0)

    for node in cliqueGraph0.nodes():
        if node == 0:
            continue
        else:
            parent = cliqueGraph0.predecessors(node)
            if len(parent) > 1 and 0 in parent:
                cliqueGraph0.remove_edge(0, node)

    print 'before remove redundant ', cliqueGraph0.number_of_nodes(), cliqueGraph0.number_of_edges()
    print 'len(cliqueGraph0.selfloop_edges())---', len(cliqueGraph0.selfloop_edges()), nx.isolates(cliqueGraph0)


    print 'remove redundant nodes between same level...'
    OntologyGraph = cliqueGraph0
    distance = nx.shortest_path_length(cliqueGraph0, source=0)
    # print distance
    # print len(distance)
    df = pd.DataFrame(distance.items())
    grouped = df.groupby(1, sort=False)
    dic_level = {}
    # 将nodes按照level分层
    for g in grouped:
        l = g[0]
        nodes = g[1][0].tolist()
        dic_level[l] = nodes
    for level in range(1, l):

        for i in range(0, len(dic_level[level])):

            for j in range(i + 1, len(dic_level[level])):
                try:
                    term1 = dic_level[level][i]
                    term2 = dic_level[level][j]
                    gene1 = cliqueGraph0.node[term1]['annotation']
                    gene2 = cliqueGraph0.node[term2]['annotation']
                    time1 = sorted(cliqueGraph0.node[term2]['windowsize'])
                    time2 = sorted(cliqueGraph0.node[term2]['windowsize'])
                except:
                    continue

                if len(set(gene2) | (set(gene1))) - len(set(gene2) & (set(gene1))) <= 2 and time1 == time2:

                    try:
                        # parent2 = cliqueGraph0.predecessors(term2)
                        child2 = cliqueGraph0.successors(term2)

                        cliqueGraph0.node[term1]['geneSet'] = list(set(gene2) | (set(gene1)))
                        # loop
                        # for p in parent2:
                        #     if p!=term1 and not cliqueGraph0.has_edge(p, term1):
                        #         cliqueGraph0.add_edge(p, term1)
                        for c in child2:
                            if c!=term1 and not cliqueGraph0.has_edge(c, term1):
                                cliqueGraph0.add_edge(term1, c)
                        cliqueGraph0.remove_node(term2)
                        # print term2
                        sum = sum + 1
                    except:
                        continue

    print 'after remove redundant the size of OntologyGraph:------', cliqueGraph0.number_of_nodes(), cliqueGraph0.number_of_edges()
    print 'len(cliqueGraph0.selfloop_edges())------', len(cliqueGraph0.selfloop_edges()), nx.isolates(cliqueGraph0)


    print 'start calculate purity...'
    dic_term_score = {}

    for node in cliqueGraph0.nodes():
        if node == 0:
            continue
        else:
            gene_set = cliqueGraph0.node[node]['annotation']
            window_set = cliqueGraph0.node[node]['windowsize']
            # 判断phenotype是否有意义
            sign, score, trend1, trend2, trend3 = sign_value(node, gene_set, window_set)

            if score < pp:
                # 无意义，delete，重定向
                parent = cliqueGraph0.predecessors(node)
                child = cliqueGraph0.successors(node)

                for p in parent:
                    for c in child:
                        if not cliqueGraph0.has_edge(p, c):
                            cliqueGraph0.add_edge(p, c)
                cliqueGraph0.remove_node(node)

            else:
                dic_term_score[node] = score
                continue
    print 'after purity window----------', cliqueGraph0.number_of_nodes(), cliqueGraph0.number_of_edges()


    for node in cliqueGraph0.nodes():
        if node == 0:
            continue
        else:
            parent = cliqueGraph0.predecessors(node)
            if len(parent) > 1 and 0 in parent:
                cliqueGraph0.remove_edge(0, node)

    loop_edges = cliqueGraph0.selfloop_edges()
    cliqueGraph0.remove_edges_from(loop_edges)
    cliqueGraph0 = max(nx.weakly_connected_component_subgraphs(cliqueGraph0), key=len)

    # write obo files


    fw2 = open('C:\Users\lu\Desktop\\0400\\0405annotation1129_025.txt', 'w')

    fw2.write(
        'term_id' + '\t' + 'sign_score' + '\t' + 'level' + '\t' + 'annotation_gene' + '\t' + 'start_time' + '\t'
         + 'end_time' + '\t'+ 'geneSize' + '\t' + 'windowSize' + '\n')

    for node in cliqueGraph0.nodes():
        if node == 0:
            continue
        else:
            try:
                fw2.write(
                    'TPO:' + str(node) + '\t' +
                    str(round(dic_term_score[node], 4)) + '\t' +
                    str(nx.shortest_path_length(cliqueGraph0, 0, node)) + '\t' +
                    ','.join(cliqueGraph0.node[node]['annotation']) + '\t' +
                    str(min(cliqueGraph0.node[node]['windowsize']) + 49) + '\t' +
                    str(max(cliqueGraph0.node[node]['windowsize']) + 58) + '\t' +
                    str(len(cliqueGraph0.node[node]['annotation'])) + '\t' +
                    str(len(cliqueGraph0.node[node]['windowsize']) + 9) + '\n')
            except:
                cliqueGraph0.add_edge(0, node)
                fw2.write(
                    'TPO:' + str(node) + '\t' +
                    str(round(dic_term_score[node], 4)) + '\t' +
                    str(nx.shortest_path_length(cliqueGraph0, 0, node)) + '\t' +
                    ','.join(cliqueGraph0.node[node]['annotation']) + '\t' +
                    str(min(cliqueGraph0.node[node]['windowsize']) + 49) + '\t' +
                    str(max(cliqueGraph0.node[node]['windowsize']) + 58) + '\t' +
                    str(len(cliqueGraph0.node[node]['annotation'])) + '\t' +
                    str(len(cliqueGraph0.node[node]['windowsize']) + 9) + '\n')


    end = time.clock()

    fw1 = open('C:\Users\lu\Desktop\\0400\\0405edges1129_025.obo', 'w')
    for node in cliqueGraph0.nodes():
        if node != 0:
            fw1.write('[Term]' + '\n')
            fw1.write('id:' + 'TPO:' + str(node) + '\n')
            fw1.write('purity:' + str(round(dic_term_score[node], 4)) + '\n')
            fw1.write('number of annotation genes:' + str(len(cliqueGraph0.node[node]['annotation'])) + '\n')
            fw1.write('number of annotation windows:' + str(len(cliqueGraph0.node[node]['windowsize']) + 9) + '\n')
            parent = set(cliqueGraph0.predecessors(node))
            for p in parent:
                fw1.write('is_a:' + 'TPO:' + str(p) + '\n')
            fw1.write('\n')
        else:
            fw1.write('[Term]' + '\n' + 'id:TPO:0' + '\n' + '\n')
    fw3 = open('C:\Users\lu\Desktop\\0400\\0405edges1129_025.txt', 'w')
    for edge in cliqueGraph0.edges():
        fw3.write('TPO:'+str(edge[0]) + '\t' + 'TPO:'+str(edge[1]) + '\n')
    fw1.close()
    fw2.close()
    fw3.close()

    dic = {}
    fr = open('G:\project2\\NPM201507\\code\\term_name_id\\termN_Id.txt', 'r')
    for line in fr:
        term, idd = line.strip().split('\t')
        dic[term] = idd

    fw1 = open('C:\Users\lu\Desktop\\0400\\0405edges_sign_id1129_025.txt', 'w')
    fw2 = open('C:\Users\lu\Desktop\\0400\\0405terms_sign_id1129_025.txt', 'w')
    fw6 = open('C:\Users\lu\Desktop\\0400\\0405terms_sign_list_id1129_025.txt', 'w')

    fw1.write('parent' + '\t' + 'child' + '\n')
    for edge in cliqueGraph0.edges():
        fw1.write(str(edge[0]) + '\t' + str(edge[1]) + '\n')
    fw2.write(
        'term_id' + '\t' + 'sign_score' + '\t' + 'level' + '\t' + 'annotation_gene' + '\t' + 'start_time' + '\t' + 'end_time' + '\t' + 'geneSize' + '\t' + 'time_size' + '\n')

    fw6.write(
        'term_id' + '\t' + 'sign_score' + '\t' + 'level' + '\t' + 'annotation_gene' + '\t' + 'start_time' + '\t' + 'end_time' + '\t' + 'geneSize' + '\t' + 'time_size' + '\n')

    for node, value in sorted(dic_term_score.items(), key=lambda d: d[1],reverse=True):
        try:
            fw2.write(str(node) + '\t' + str(round(value, 4)) + '\t' +
                      str(
                          nx.shortest_path_length(cliqueGraph0, 0, node)) + '\t' +
                      ','.join(cliqueGraph0.node[node]['annotation'])
                      + '\t' + str(
                min(cliqueGraph0.node[node]['windowsize']) + 49) + '\t' +
                str(max(cliqueGraph0.node[node]['windowsize']) + 58) + '\t' +
                      str(len(cliqueGraph0.node[node]['annotation']))+ '\t' +
                 str(len(cliqueGraph0.node[node]['windowsize']) + 9) + '\n')

            fw6.write(
                str(node) + '\t' + str(round(value, 4)) + '\t' + str(
                    nx.shortest_path_length(cliqueGraph0, 0, node)) + '\t' + ','.join(
                    [dic[t] for t in cliqueGraph0.node[node]['annotation']]) + '\t' + str(
                    min(cliqueGraph0.node[node]['windowsize']) + 49) + '\t' + str(
                    max(cliqueGraph0.node[node]['windowsize']) + 58) + '\t' +
                    str(len(cliqueGraph0.node[node]['annotation'])) + '\t' + str(
                    len(cliqueGraph0.node[node]['windowsize']) + 9) + '\n')
        except:
            continue


    fw1.close()
    fw2.close()
    fw6.close()
    print 'FINALLY----------', cliqueGraph0.number_of_nodes(), cliqueGraph0.number_of_edges()

    print 'pp----term_number----level----leaves----alll----percentage'

    print pp,cliqueGraph0.number_of_nodes(), max(nx.shortest_path_length(cliqueGraph0, 0).values()), delete_leave_flag1, delete_sum


    print 'The function run time is : %.03f seconds' % (end - start)