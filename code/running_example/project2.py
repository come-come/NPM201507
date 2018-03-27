#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
@Time    : 2018/3/24 17:07
@Author  : Junya Lu
@Site    : 
"""

import pandas as pd
import numpy as np
import networkx as nx
import os, time, random
import matplotlib.pyplot as plt
import seaborn as sns
import operator
from compiler.ast import flatten
import sys, getopt

# phe1 = pd.read_table('phe1.txt', index_col=0)
# phe2 = pd.read_table('phe2.txt', index_col=0)
# phe3 = pd.read_table('phe3.txt', index_col=0)

curr_path = os.getcwd() + '/'

def tree_statis(weight_value, window):
    # output 50 trees
    # filename2 = 'edge/edgeInfo' + str(window) + '.txt'
    filename2 = str(window) + '.txt'
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
                if len(clique) >size_clique:
                    # print window, 'clique:', clique
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
                if len(clique) > size_clique:
                    # print window, 'clique:', clique

                    for key, value in dic_last_time.items():  # key 是clique ,value是 [window]
                        intersect = sorted(set(key).intersection(set(clique)))
                        q = 0
                        # if len(intersect) >=  size_clique:
                        if len(intersect) > size_clique:
                            # print 'intersect', intersect
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


def phenotype(f_path, termID, gene_set, window):
    phe1 = pd.read_table('G:\project2\\NPM201507\\data\\IDMapping_consolidated_allPhi2_cleaned_lfc_avg.txt',
                         index_col=0)
    phe2 = pd.read_table('G:\project2\\NPM201507\\data\\IDMapping_consolidated_allQESV_cleaned_LFC_avg.txt',
                         index_col=0)
    phe3 = pd.read_table('G:\project2\\NPM201507\\data\\IDMapping_consolidated_allQI_new_RAW3_adj_LFC_avg.txt',
                         index_col=0)
    # if the window is window
    # t_min = min(window) + 49
    # t_max = max(window) + 49 + 10
    # if the window is time point
    t_min = min(window)
    t_max = max(window)
    # rename columns index. [0, 1, 2...]
    window_set = [i for i in range(t_min, t_max + 1)]  # 回到最原始的数据上
    phe1.columns = [i for i in range(0, 113)]
    phe2.columns = [i for i in range(0, 113)]
    phe3.columns = [i for i in range(0, 113)]
    # print phe1.loc[gene_set, window_set]
    # print phe2.loc[gene_set, window_set]
    # print phe3.loc[gene_set, window_set]
    # write the phenotype value into files
    f_fpath1 = f_path + '\\Phi2.txt'
    f_fpath2 = f_path + '\\qE.txt'
    f_fpath3 = f_path + '\\qI.txt'
    phe1.loc[gene_set, window_set].to_csv(f_fpath1, sep='\t')
    phe2.loc[gene_set, window_set].to_csv(f_fpath2, sep='\t')
    phe3.loc[gene_set, window_set].to_csv(f_fpath3, sep='\t')
    q = phe1.loc[gene_set, window_set]
    q2 = phe2.loc[gene_set, window_set]
    q3 = phe3.loc[gene_set, window_set]
    phe_avg = (phe1.loc[gene_set, window_set] + phe2.loc[gene_set, window_set] + phe3.loc[gene_set, window_set]).apply(
        lambda x: x * 2)
    print 'avg'
    # print phe_avg
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
    fig.suptitle(termID, fontsize=16, y=1)
    plt.subplot(231)
    ax = sns.heatmap(q, vmin=-0.5, vmax=0.5, cmap=cmap, center=0, square=True, linewidths=.5,
                     cbar_kws={"shrink": .5})
    for item in ax.get_yticklabels():
        item.set_rotation(0)
    plt.title('Phi2')
    plt.xlabel('time point')

    plt.subplot(232)
    ax = sns.heatmap(q2, vmin=-0.5, vmax=0.5, cmap=cmap, center=0, square=True, linewidths=.5,
                     cbar_kws={"shrink": .5})
    for item in ax.get_yticklabels():
        item.set_rotation(0)
    plt.title('qE')
    plt.xlabel('time point')

    plt.subplot(233)
    ax = sns.heatmap(q3, vmin=0, vmax=1, cmap=cmap, center=0, square=True, linewidths=.5,
                     cbar_kws={"shrink": .5})
    for item in ax.get_yticklabels():
        item.set_rotation(0)
    plt.title('ql')
    plt.xlabel('time point')

    fig.set_size_inches(22, 10.5)
    fig_path = f_path + str(termID) + '.png'
    plt.savefig(fig_path, dpi=100)


def sign_value(node, gene_set, window):

    # window_set = [i for i in range(min(window) + 49, max(window) + 49 + 10)]  # 回到最原始的数据上
    window_set = window
    p1 = phe1.loc[gene_set, window_set]
    p2 = phe2.loc[gene_set, window_set]
    p3 = phe3.loc[gene_set, window_set]

    p1_list = flatten(p1.values.tolist())
    p2_list = flatten(p2.values.tolist())
    p3_list = flatten(p3.values.tolist())
    sign_p1 = [x for x in p1_list if x < phe1_t1 or x > phe1_t2]
    sign_p2 = [x for x in p2_list if x < phe2_t1 or x > phe2_t2]
    sign_p3 = [x for x in p3_list if x < phe3_t1 or x > phe3_t2]
    sign = len(sign_p1) + len(sign_p2) + len(sign_p3)
    all = len(p1_list) + len(p2_list) + len(p3_list)

    return sign, float(sign)/all


def class_distance(node, gene_set, window):
    t_min = min(window) + 49
    t_max = max(window) + 49 + 10
    window_set = [i for i in range(t_min, t_max)]  # 回到最原始的数据上
    p1 = phe1.loc[gene_set, window_set].dropna(axis=1, how='any')
    p2 = phe2.loc[gene_set, window_set].dropna(axis=1, how='any')
    p3 = phe3.loc[gene_set, window_set].dropna(axis=1, how='any')
    result = pd.concat([p1, p2, p3], axis=1)
    # print result.head(5)
    distance = np.trace(np.cov(result))
    return distance


def pearson(node, gene_set, window):
    t_min = min(window) + 49
    t_max = max(window) + 49 + 10
    window_set = [i for i in range(t_min, t_max)]  # 回到最原始的数据上
    p1 = phe1.loc[gene_set, window_set].dropna(axis=1, how='any')
    p2 = phe2.loc[gene_set, window_set].dropna(axis=1, how='any')
    p3 = phe3.loc[gene_set, window_set].dropna(axis=1, how='any')
    result = pd.concat([p1, p2, p3], axis=1)
    p_matrix = np.corrcoef(result)
    num = (p_matrix.shape[0] * p_matrix.shape[0] - p_matrix.shape[0]) / 2
    sum1 = (p_matrix.sum() - p_matrix.shape[0]) / 2
    return sum1 / num


if __name__ == '__main__':

    start = time.clock()

    inputfile = ""

    try:
        opts, args = getopt.getopt(sys.argv[1:], "hi:d:f:g:w:s:p:m:n:b:v:c:z:", ["ifile=", "phe1_file=","phe2_file=","phe3_file=","weight_value=", "size_clique=","purity=",
                                                               "phe1_t1=", "phe1_t2=",
                                                                "phe2_t1=", "phe2_t2=",
                                                                "phe3_t1=", "phe3_t2="])

    except getopt.GetoptError:
        print 'test.py -i <inputfile> -o <outputfile>'
        sys.exit(2)
    print opts
    for opt, arg in opts:
        if opt == '-h':
            print 'project2.py -i <inputfile> -w weight_value -s size_clique -p purity'
            sys.exit()
        elif opt in ("-i", "--ifile"):
            inputfile = arg
        elif opt in ("-d", "--phe1_file"):
            phe1_file = arg
        elif opt in ("-f", "--phe2_file"):
            phe2_file = arg
        elif opt in ("-g", "--phe3_file"):
            phe3_file = arg
        elif opt in ("-w", "--weight_value"):
            weight_value = float(arg)
        elif opt in ("-s", "--size_clique"):
            size_clique = int(arg)-1
        elif opt in ("-p", "--purity"):
            sign_score = float(arg)
        elif opt in ("-m", "--phe1_t1"):
            phe1_t1 = float(arg)
        elif opt in ("-n", "--phe1_t2"):
            phe1_t2 = float(arg)
        elif opt in ("-b", "--phe2_t1"):
            phe2_t1 = float(arg)
        elif opt in ("-v", "--phe2_t2"):
            phe2_t2 = float(arg)
        elif opt in ("-c", "--phe3_t1"):
            phe3_t1 = float(arg)
            print phe3_t1, type(phe3_t1)
        elif opt in ("-z", "--phe3_t2"):
            phe3_t2 = float(arg)

    # python project2.py -i Graph1.txt -w 0.8 -s 1 -p 0.5 -m 0.5 -n 0.7 -b 0.2 -v 0.7 -c 0.1 -z 1

    term = 183
    # filename = 'Graph1.txt'
    data = pd.read_csv(inputfile, index_col=0, sep='\t')
    columns_num = int(data.shape[1])

    phe1 = pd.read_table(phe1_file,index_col=0)
    phe2 = pd.read_table(phe1_file,index_col=0)
    phe3 = pd.read_table(phe1_file,index_col=0)
    phe1.columns = [i for i in range(0, phe1.shape[1])]
    phe2.columns = [i for i in range(0, phe2.shape[1])]
    phe3.columns = [i for i in range(0, phe3.shape[1])]
    '''
    # columns_num 代替所有的50 和 54
    # 定义weight 阈值
    weight_value = 0.8
    # 定义 clique大小 阈值
    size_clique = 1 # >
    # 定义 phenotype value有意义的区间
    phe1_t1 = 0.5
    phe1_t2 = 0.7
    phe2_t1 = 0.2
    phe2_t2 = 0.7
    phe3_t1 = 0.3
    phe3_t2 = 0.9
    # 定义score 阈值
    sign_score = 0.005
    '''
    curr_path = os.getcwd() + '/'
    windowGraph = {}
    # output undirected graph
    # for i in range(0, columns_num):
    #     tree_statis(weight_value, i)

    s = 0

    for window in range(0, columns_num):
        windowGraph[window] = nx.Graph()
        df = data[data[data.columns[window]] >= (weight_value + 0.00001)]
        # print window, weight_value, df.shape
        for edge in range(0, df.shape[0]):
            node_1, node_2 = df.index[edge].split('_')
            windowGraph[window].add_edge(node_1, node_2)
    # 先产生第一个window的 tree
    cliqueGraph0, dic_term0, dic_term_num0, term = tree0(weight_value, 0, 183)
    # print dic_term0
    # print dic_term_num0
    # print term

    # 验证 t=1 时刻下生成的tree 验证无误
    # pos = nx.spring_layout(cliqueGraph0)
    # nx.draw_networkx(cliqueGraph0, pos, with_labels=True, node_size=600)
    # for node in cliqueGraph0.nodes():
    #     x, y = pos[node]
    #     plt.text(x, y + 0.1, s=cliqueGraph0.node[node]['annotation'],
    #              bbox=dict(facecolor='red', alpha=0.5), horizontalalignment='center')
    #     plt.text(x, y + 0.15, s=cliqueGraph0.node[node]['windowsize'],
    #              bbox=dict(facecolor='red', alpha=0.5), horizontalalignment='center')
    #     plt.title('t=1', fontsize=20)
    # plt.show()

    # 生成之后的window tree
    dic_all = {}
    dic_all = dic_term0.copy()
    # dic_all_num = {}
    dic_all_term_num = dic_term_num0.copy()
    copy_clique = cliqueGraph0
    for i in range(1, columns_num):
        print 'begin term num:', term
        cliqueGraph1, dic_term1, dic_term_num1, term = tree0(weight_value, i, term)

        # 验证 t>1 时刻下生成的tree 验证无误
        # print 'show 当前时刻下的tree t>1 :'
        # print dic_term1
        # print dic_term_num1
        # print term
        # pos = nx.spring_layout(cliqueGraph1)
        # nx.draw_networkx(cliqueGraph1, pos, with_labels=True, node_size=600)
        # for node in cliqueGraph1.nodes():
        #     x, y = pos[node]
        #     plt.text(x, y + 0.1, s=cliqueGraph1.node[node]['annotation'],
        #              bbox=dict(facecolor='red', alpha=0.5), horizontalalignment='center')
        #     plt.text(x, y + 0.15, s=cliqueGraph1.node[node]['windowsize'],
        #              bbox=dict(facecolor='red', alpha=0.5), horizontalalignment='center')
        # title= 't='+str(i+1)
        # plt.title(title, fontsize=20)
        # plt.show()
        # 先合并两个trees
        cliqueGraph0 = nx.compose(cliqueGraph0, cliqueGraph1)
        # 判断冗余信息
        for key in dic_term1.keys(): # dic_term1 当前时间生成的结点信息 geneList:windowList
            if dic_all.has_key(key): # dic_all 之前时间点生成的结点信息 累计
                # print 'exit key',key
                # 如果之前 已经生成 过这样的geneList
                if set(dic_all[key]).issuperset(set(dic_term1[key])):
                    # 如果之前的时间片段较长 就舍弃当前这个term
                    dic_term1.pop(key)
                    num = dic_term_num1[key]
                    cliqueGraph0.remove_node(num)
                    dic_term_num1.pop(key)
                    s = s + 1
                    # print 'remove node', num
                else:
                    #如果现在的时间片段长或者是两者没有子集关系 则 合并时间片段 并舍弃当前这个term
                    dic_all[key] = dic_all[key] + dic_term1[key]
                    dic_term1.pop(key)
                    num = dic_term_num1[key]
                    cliqueGraph0.remove_node(num)
                    dic_term_num1.pop(key)
                    # print 'update node',dic_all_term_num[key], 'and remove node', num
            # 2018.3.23
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
                        parent = cliqueGraph0.predecessors(old_id)
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
                                    # print p_id, old, this_id
                                if dic_this_edge[(this_id, p_id)] == 1:
                                    flag1 = 1
                                    # print p_id, old, this_id

                                if flag1 == 1:
                                    break
                            except:
                                continue
                        if flag1 == 0:
                            # key--> old
                            cliqueGraph0.add_edge(this_id, old_id)
                            dic_this_edge[(this_id, old_id)] = 1
                            # print 'add edge', this_id, old_id
                        else:
                            # 已经连接到它的父亲上了
                            continue

                    elif set(key).issuperset(old) and  set(dic_term1[key]).issuperset(dic_all[old]):
                        # 如果当前结点的注释基因是它的父集 注释时间也是它的父集  delete old
                        cliqueGraph0.remove_node(old_id)
                        dic_all.pop(old)
                        dic_all_term_num.pop(old)

                        # print 'remove node', old_id
                        ###
                    elif set(old).issuperset(key) and  set(dic_term1[key]).issuperset(dic_all[old]):
                        # 如果当前结点的注释基因是它的子集 注释时间是它的父集 old--> key
                        cliqueGraph0.add_edge(old_id, this_id)
                        print 'add edge', old_id, this_id
                    elif set(old).issuperset(key) and set(dic_all[old]).issuperset(dic_term1[key]):
                        # 如果当前结点的注释基因是它的子集 注释时间也是它的子集 delete this
                        cliqueGraph0.remove_node(this_id)
                        dic_term_num1.pop(this_id)
                        dic_term1.pop(key)
                        # print 'remove node', this_id
                        #########
            # 2018.3.23
        # print 'dic_this_edge', dic_this_edge
        # print dic_all_term_num
        # print dic_term_num1
        dic_all.update(dic_term1)
        dic_all_term_num.update(dic_term_num1)



        # print '融合之后的tree'
        # print dic_all
        # print cliqueGraph0.edges()
        # print 'edge:', cliqueGraph0.number_of_edges(), 'node:', cliqueGraph0.number_of_nodes()
        # pos = nx.spring_layout(cliqueGraph0)
        # nx.draw_networkx(cliqueGraph0, pos, with_labels=True, node_size=600)
        # for node in cliqueGraph0.nodes():
        #     x, y = pos[node]
        #     plt.text(x, y + 0.1, s=cliqueGraph0.node[node]['annotation'],
        #              bbox=dict(facecolor='red', alpha=0.5), horizontalalignment='center')
        #     plt.text(x, y + 0.15, s=cliqueGraph0.node[node]['windowsize'],
        #              bbox=dict(facecolor='red', alpha=0.5), horizontalalignment='center')
        # title= 'fusing' + str(i)
        # plt.title(title, fontsize=20)
        # plt.show()
    # remove link 0


    print 'start calculate purity...'
    dic_term_score = {}
    for node in cliqueGraph0.nodes():
        if node == 0:
            continue
        else:
            gene_set = cliqueGraph0.node[node]['annotation']
            window_set = cliqueGraph0.node[node]['windowsize']
            # 判断phenotype是否有意义

            sign, score = sign_value(node, gene_set, window_set)

            if score < sign_score:
                # 无意义，delete，重定向
                parent = cliqueGraph0.predecessors(node)
                child = cliqueGraph0.successors(node)
                for p in parent:
                    for c in child:
                        cliqueGraph0.add_edge(p, c)
                cliqueGraph0.remove_node(node)
                # print 'purity is low...then remove this node:', node
            else:
                dic_term_score[node] = score
                continue
    print 'after purity window', i, cliqueGraph0.number_of_nodes(), cliqueGraph0.number_of_edges()

    for node in cliqueGraph0.nodes():
        if node == 0:
            continue
        else:
            parent = cliqueGraph0.predecessors(node)
            if len(parent) > 1 and 0 in parent:
                cliqueGraph0.remove_edge(0, node)

    # pos = nx.spring_layout(cliqueGraph0)
    # nx.draw_networkx(cliqueGraph0, pos, with_labels=True, node_size=600)
    # for node in cliqueGraph0.nodes():
    #     x, y = pos[node]
    #     plt.text(x, y + 0.1, s=cliqueGraph0.node[node]['annotation'],
    #              bbox=dict(facecolor='red', alpha=0.5), horizontalalignment='center')
    #     plt.text(x, y + 0.15, s=cliqueGraph0.node[node]['windowsize'],
    #              bbox=dict(facecolor='red', alpha=0.5), horizontalalignment='center')
    #
    # plt.title('After purity option', fontsize=20)
    # plt.show()


    # redudant
    print 'remove redundant nodes between same level...'
    OntologyGraph = cliqueGraph0
    distance = nx.shortest_path_length(cliqueGraph0, source=0)
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
                term1 = dic_level[level][i]
                term2 = dic_level[level][j]
                gene1 = cliqueGraph0.node[term1]['annotation']
                gene2 = cliqueGraph0.node[term2]['annotation']
                time1 = sorted(cliqueGraph0.node[term2]['windowsize'])
                time2 = sorted(cliqueGraph0.node[term2]['windowsize'])
                if len(set(gene2) | (set(gene1))) - len(set(gene2) & (set(gene1))) == 2 and time1 == time2:

                    try:
                        parent2 = OntologyGraph.predecessors(term2)
                        child2 = OntologyGraph.successors(term2)

                        OntologyGraph.node[term1]['geneSet'] = list(set(gene2) | (set(gene1)))
                        for p in parent2:
                            OntologyGraph.add_edge(p, term1)
                        for c in child2:
                            OntologyGraph.add_edge(term1, c)
                        OntologyGraph.remove_node(term2)
                        # print term2
                        sum = sum + 1
                    except:
                        continue

    print 'after remove redundant the size of OntologyGraph:', OntologyGraph.number_of_nodes(), OntologyGraph.number_of_edges()

    # pos = nx.spring_layout(cliqueGraph0)
    # nx.draw_networkx(cliqueGraph0, pos, with_labels=True, node_size=600)
    # for node in cliqueGraph0.nodes():
    #     x, y = pos[node]
    #     plt.text(x, y + 0.1, s=cliqueGraph0.node[node]['annotation'],
    #              bbox=dict(facecolor='red', alpha=0.5), horizontalalignment='center')
    #     plt.text(x, y + 0.15, s=cliqueGraph0.node[node]['windowsize'],
    #              bbox=dict(facecolor='red', alpha=0.5), horizontalalignment='center')
    #
    # plt.title('After redundant option', fontsize=20)
    # plt.show()


    print 'start calculate purity...'
    dic_term_score = {}
    for node in cliqueGraph0.nodes():
        if node == 0:
            continue
        else:
            gene_set = cliqueGraph0.node[node]['annotation']
            window_set = cliqueGraph0.node[node]['windowsize']
            # 判断phenotype是否有意义

            sign, score = sign_value(node, gene_set, window_set)

            if score < sign_score:
                # 无意义，delete，重定向
                parent = cliqueGraph0.predecessors(node)
                child = cliqueGraph0.successors(node)
                for p in parent:
                    for c in child:
                        cliqueGraph0.add_edge(p, c)
                cliqueGraph0.remove_node(node)
                # print 'purity is low...then remove this node:', node
            else:
                dic_term_score[node] = score
                continue
    print 'after purity window', i, cliqueGraph0.number_of_nodes(), cliqueGraph0.number_of_edges()


    # draw picture
    pos = nx.spring_layout(cliqueGraph0)
    nx.draw_networkx(cliqueGraph0, pos, with_labels=True, node_size=600)
    for node in cliqueGraph0.nodes():
        x, y = pos[node]
        plt.text(x, y + 0.1, s=cliqueGraph0.node[node]['annotation'],
                 bbox=dict(facecolor='red', alpha=0.5), horizontalalignment='center')
        plt.text(x, y + 0.15, s=cliqueGraph0.node[node]['windowsize'],
                 bbox=dict(facecolor='red', alpha=0.5), horizontalalignment='center')

    plt.title('After purity option', fontsize=20)
    plt.show()

    # write files
    fw1 = open('edges.obo', 'w')

    for node in cliqueGraph0.nodes():
        if node != 0:
            fw1.write('[Term]' + '\n')
            fw1.write('id:' + str(node) + '\n')
            fw1.write('purity:' + str(round(dic_term_score[node], 4)) + '\n')
            fw1.write('number of annotation genes:' + str(len(cliqueGraph0.node[node]['annotation'])) + '\n')
            fw1.write('number of annotation windows:' + str(len(cliqueGraph0.node[node]['windowsize'])) + '\n')
            parent = cliqueGraph0.predecessors(node)
            for p in parent:
                fw1.write('is_a:' + str(p) + '\n')
            fw1.write('\n')
        else:
            fw1.write('[Term]' + '\n' + 'id:0' + '\n' + '\n')

    fw2 = open('annotation.txt', 'w')

    fw2.write(
        'term_id' + '\t' + 'sign_score' + '\t' + 'level' + '\t' + 'annotation_gene' + '\t' + 'start_time' + '\t' + 'end_time' + '\t' + 'geneSize' + '\t' + 'time_size' + '\n')

    for node in cliqueGraph0.nodes():
        if node == 0:
            continue
        else:
            fw2.write(
                str(node) + '\t' +
                str(round(dic_term_score[node], 4)) + '\t' +
                str(nx.shortest_path_length(cliqueGraph0, 0, node)) +'\t'+
                ','.join(cliqueGraph0.node[node]['annotation']) + '\t' +
                str(cliqueGraph0.node[node]['windowsize']) + '\t' +
                str(len(cliqueGraph0.node[node]['annotation'])) + '\t' +
                str(len(cliqueGraph0.node[node]['windowsize'])) + '\n')
    end = time.clock()
    fw3 = open('edges.txt', 'w')
    for edge in cliqueGraph0.edges():
        fw3.write(str(edge[0]) + '\t' + str(edge[1]) + '\n')
    fw1.close()
    fw2.close()
    fw3.close()

    print 'The function run time is : %.03f seconds' % (end - start)
