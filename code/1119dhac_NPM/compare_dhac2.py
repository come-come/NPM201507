# -*- coding: utf-8 -*-
import pandas as pd
import numpy as np
import os
from os.path import join
from compiler.ast import flatten
import networkx as nx
dic_result = {}  # dhac result
phe1 = pd.read_table('G:\project2\\NPM201507\\data\\IDMapping_consolidated_allPhi2_cleaned_lfc_avg.txt', index_col=0)
phe2 = pd.read_table('G:\project2\\NPM201507\\data\\IDMapping_consolidated_allQESV_cleaned_LFC_avg.txt', index_col=0)
phe3 = pd.read_table('G:\project2\\NPM201507\\data\\IDMapping_consolidated_allQI_new_RAW3_adj_LFC_avg.txt', index_col=0)

phe1.columns = [i for i in range(0, 113)]
phe2.columns = [i for i in range(0, 113)]
phe3.columns = [i for i in range(0, 113)]

dic = {}
fr = open('G:\project2\\NPM201507\\code\\term_name_id\\termN_Id.txt', 'r')
for line in fr:
    term, idd = line.strip().split('\t')
    dic[idd] = term
print dic

def read_dhac_result(dest):
    for root, dirs, files in os.walk(dest):
        for OneFileName in files:
            if OneFileName.find('.group') == -1:
                continue
            OneFullFileName = join(root, OneFileName)
            # print OneFullFileName
            # print OneFileName.strip().split('edgeInfo')[1].split('.group')[0]
            dic_result[int(OneFileName.strip().split('edgeInfo')[1].split('.group')[0])] = []
            for line in open(OneFullFileName, 'r'):
                b = line.strip().split('\t')
                # print b  # ['821137', '829490']
                geneName = [dic[i] for i in b]
                if len(geneName) > 5:
                    dic_result[int(OneFileName.strip().split('edgeInfo')[1].split('.group')[0])].append(geneName)
def read_NPM_result(dest):
    print 'df'
    for root, dirs, files in os.walk(dest):

        for OneFileName in files:

            if OneFileName.find('.txt') == -1:
                continue
            OneFullFileName = join(root, OneFileName)

            dic_result[int(OneFileName.strip().split('.txt')[0])] = []
            print int(OneFileName.strip().split('.txt')[0])
            for line in open(OneFullFileName, 'r'):
                b = line.strip().split('\t')
                # print b  # ['AT656', 'AT43']
                if len(b) > 5:
                    dic_result[int(OneFileName.strip().split('.txt')[0])].append(b)

def tree0(weight_value, startwindow, term):

    # 2017.12.4

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

    i = 0
    q = 0
    for window in range(startwindow, 50):
        windowGraph[window] = nx.Graph()

        dic_intersect_level.clear()

        if window == startwindow:
            for clique in dic_result[window]:
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
            for clique in dic_result[window]:
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

def tree_old(weight_value, startwindow, term):


    # 2017.11.28
    dic_term = {}
    dic_last_time = {}
    dic_term_num = {}
    dic_intersect_level = {}

    for window in range(startwindow, 50):
        if window == startwindow:
            for clique in dic_result[window]:
                    dic_term[frozenset(clique)] = [window]  # dic_term 记录 window和clique
                    dic_term_num[frozenset(clique)] = term  # dic_term_num
                    dic_last_time[frozenset(clique)] = [window]  # dic_last_time
                    term = term + 1
        else:
            for clique in dic_result[window]: # 当前window下直接生成的
                for key, value in dic_last_time.items():  # key 是clique ,value是 [window]
                    intersect = sorted(set(key).intersection(set(clique)))
                    q = 0
                    if len(intersect) >= 5:
                        if dic_term.has_key(frozenset(intersect)):
                            dic_term[frozenset(intersect)] = value + [window]
                            dic_intersect_level[frozenset(intersect)] = dic_term[frozenset(intersect)]
                            break
                        else:
                            dic_term[frozenset(intersect)] = value + [window]
                            dic_intersect_level[frozenset(intersect)] = dic_term[frozenset(intersect)]
                            q = 1

                            try:
                                dic_term.pop(key)
                            except:
                                continue
                            break
                    else:
                        continue
            dic_last_time.clear()
            for key, value in dic_intersect_level.items():
                dic_last_time[key] = value
            dic_intersect_level.clear()
    return dic_term

def statis():
    # 画直方图
    dic_statis_my = {}
    dic_statis = {}
    dic_statis_dhac = {}
    fr1 = open('G:\\project2\\NPM201507\\code\\1117terms_sign_list_id.txt','r')
    fr2 = open('dhac_term.txt', 'r')
    fr = open('super_than_dhac1128.txt', 'r')

    dic_statis[1] = 0
    dic_statis_my[1] = 0
    dic_statis_dhac[1] = 0
    title = fr.readline()
    title2 = fr1.readline()
    title3 = fr2.readline()
    dica_statis_my = {}
    dica_statis = {}
    dica_statis_dhac = {}
    dica_statis[1] = 0
    dica_statis_my[1] = 0
    dica_statis_dhac[1] = 0

    for line in fr1:
        line_str = line.strip().split('\t')
        geneSize = int(line_str[5])
        sign_score = float(line_str[1])
        if  geneSize > 40:
            dic_statis_my[1] = dic_statis_my[1] + 1
        else:
            try:
                dic_statis_my[geneSize] = dic_statis_my[geneSize] + 1
                dica_statis_my[geneSize] = dica_statis_my[geneSize] + sign_score
            except:
                dic_statis_my[geneSize] = 1
                dica_statis_my[geneSize] = sign_score
    for line in fr2:
        line_str = line.strip().split('\t')
        geneSize = int(line_str[4])
        sign_score = float(line_str[0])
        if  geneSize > 40:
            dic_statis_dhac[1] = dic_statis[1] + 1
        else:
            try:
                dic_statis_dhac[geneSize] = dic_statis_dhac[geneSize] + 1
                dica_statis_dhac[geneSize] = dica_statis_dhac[geneSize] + sign_score
            except:
                dic_statis_dhac[geneSize] = 1
                dica_statis_dhac[geneSize] = sign_score

    for line in fr:
        line_str = line.strip().split('\t')
        geneSize = int(line_str[2])
        if  geneSize > 40:
            dic_statis[1] = dic_statis[1] + 1
        else:
            try:
                dic_statis[geneSize] = dic_statis[geneSize] + 1
            except:
                dic_statis[geneSize] = 1

    fw = open('statis40_1201_a.txt', 'a')
    for key, value in dic_statis_my.items():
        fw.write(str(key) + '\t')
    fw.write('\n')
    n = 0
    for key, value in dic_statis_my.items():
        fw.write(str(value) + '\t')
        n = n + value
    fw.write('\n')

    for key, value in dic_statis_my.items():
        try:
            fw.write(str(dic_statis_dhac[key]) + '\t')
        except:
            fw.write('0' + '\t')
    fw.write('\n')

    for key, value in dic_statis_my.items():
        try:
            fw.write(str(dic_statis[key]) + '\t')
        except:
            fw.write('0' + '\t')
    fw.write('\n')

    for key, value in dic_statis_my.items():
        try:
            fw.write(str(dica_statis_dhac[key]/dic_statis_dhac[key]) + '\t')
        except:
            fw.write('0' + '\t')
    fw.write('\n')
    for key, value in dic_statis_my.items():
        try:
            fw.write(str(dica_statis_my[key]/dic_statis_my[key]) + '\t')
        except:
            fw.write('0' + '\t')
    fw.write('\n')

    fw.close()

    print dic_statis
    print dic_statis_my
    print dica_statis_my
    print dica_statis_dhac
    print n
def pearson(node, gene_set, window):
    t_min = min(window) + 49
    t_max = max(window) + 49 + 10
    window_set = [i for i in range(t_min, t_max)]
    p1 = phe1.loc[gene_set, window_set].dropna(axis=1,how='any')
    p2 = phe2.loc[gene_set, window_set].dropna(axis=1,how='any')
    p3 = phe3.loc[gene_set, window_set].dropna(axis=1,how='any')
    result = pd.concat([p1, p2, p3], axis=1)
    p_matrix = np.corrcoef(result)
    num = (p_matrix.shape[0] * p_matrix.shape[0] - p_matrix.shape[0]) / 2
    sum1 = (p_matrix.sum() - p_matrix.shape[0]) / 2
    return sum1/num
def class_distance(node, gene_set, window):
    t_min = min(window) + 49
    t_max = max(window) + 49 + 10
    window_set = [i for i in range(t_min, t_max)]  # 回到最原始的数据上
    p1 = phe1.loc[gene_set, window_set].dropna(axis=1,how='any')
    p2 = phe2.loc[gene_set, window_set].dropna(axis=1,how='any')
    p3 = phe3.loc[gene_set, window_set].dropna(axis=1,how='any')
    result = pd.concat([p1, p2, p3], axis=1)
    # print result.head(5)
    distance = np.trace(np.cov(result))
    return distance
def sign_value(node, gene_set, window):
    t_min = min(window) + 49
    t_max = max(window) + 49 + 10
    window_set = [i for i in range(t_min,t_max)] # 回到最原始的数据上
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
    if p1[t_min].mean()-p1[t_max - 1].mean()<0:
        trend1 = 1
    else :
        trend1 = 0
    if p2[t_min].mean()-p2[t_max - 1].mean()<0:
        trend2 = 1
    else :
        trend2 = 0
    if p3[t_min].mean()-p3[t_max - 1].mean()<0:
        trend3 = 1
    else :
        trend3 = 0
    return sign, float(sign)/all, trend1, trend2, trend3
    # return float(sign)/all
def average_compare():
    dhac = pd.read_table('dhac_term.txt')
    our_pearson = pd.read_table('G:\\project2\\NPM201507\\code\\1117term_pearson_id.txt')
    our_distance = pd.read_table('G:\\project2\\NPM201507\\code\\1117sign_distance_id.txt')
    our_sign = pd.read_table('G:\\project2\\NPM201507\\code\\1117terms_sign_list_id.txt')
    print 'dhac_pearson:', dhac['pearson'].mean()
    print 'dhac_purity:', dhac['purity'].mean()
    print 'dhac_class_distance:', dhac['class_distance'].mean()
    print 'our_pearson:', our_pearson['pearson'].mean()
    print 'our_purity:', our_sign['sign_score'].mean()
    print 'our_class_distance:', our_distance['distance'].mean()


if __name__ == "__main__":
    # statis()
    # average_compare()
    # 2017.12.4
    # DHAC
    # read_dhac_result('G:\\project2\\NPM201507\\code\\edge\\20171113dhac')
    # NPM  5 clusters
    # read_NPM_result('G:\\project2\\NPM201507\\code\\edge\\NPM')
    # NPM  10 clusters
    read_NPM_result('G:\\project2\\NPM201507\\code\\edge\\NPM_cluster10_(1)')
    # read_NPM_result('G:\\project2\\NPM201507\\code\\edge\\NPM_cluster10_(2)')
    # read_NPM_result('G:\\project2\\NPM201507\\code\\edge\\NPM_cluster10_1222')
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
    dic_vector = {}
    dic_pearson = {}
    dic = {}
    fr = open('G:\project2\\NPM201507\\code\\term_name_id\\termN_Id.txt', 'r')
    for line in fr:
        term, idd = line.strip().split('\t')
        dic[term] = idd
    for node in cliqueGraph0.nodes():
        if node == 0:
            continue
        else:
            gene_set = cliqueGraph0.node[node]['annotation']
            window_set = cliqueGraph0.node[node]['windowsize']
            # 判断phenotype是否有意义

            sign, score, trend1, trend2, trend3 = sign_value(node, gene_set, window_set)
            dis = class_distance(node, gene_set, window_set)
            pearson_coff = pearson(node, gene_set, window_set)
            if score < 0.005:
                # 无意义，delete，重定向
                parent = cliqueGraph0.predecessors(node)
                child = cliqueGraph0.successors(node)
                for p in parent:
                    for c in child:
                        cliqueGraph0.add_edge(p, c)
                cliqueGraph0.remove_node(node)
            else:
                dic_term_score[node] = score
                dic_term_distance[node] = dis
                dic_pearson[node] = pearson_coff
                dic_vector[node] = []
                dic_vector[node].append(len(cliqueGraph0.node[node]['annotation']))
                dic_vector[node].append(len(cliqueGraph0.node[node]['windowsize']) + 9)
                if len(cliqueGraph0.successors(node)) == 0:
                    dic_vector[node].append(1)
                else:
                    dic_vector[node].append(0)
                dic_vector[node].append(trend1)
                dic_vector[node].append(trend2)
                dic_vector[node].append(trend3)
                continue
    print 'finally:', cliqueGraph0.number_of_nodes(), cliqueGraph0.number_of_edges()


    # write into files
    fw1 = open('1211_NPM10_edges_sign_id.txt', 'w')
    fw2 = open('1211_NPM10_terms_sign_id.txt', 'w')
    fw3 = open('1211_NPM10_sign_distance_id.txt', 'w')
    fw4 = open('1211_NPM10_term_vector_id.txt', 'w')
    fw5 = open('1211_NPM10_term_pearson_id.txt', 'w')
    fw6 = open('1211_NPM10_terms_sign_list_id.txt', 'w')
    # fw1 = open('1208_DHAC_edges_sign_id.txt', 'w')
    # fw2 = open('1208_DHAC_terms_sign_id.txt', 'w')
    # fw3 = open('1208_DHAC_sign_distance_id.txt', 'w')
    # fw4 = open('1208_DHAC_term_vector_id.txt', 'w')
    # fw5 = open('1208_DHAC_term_pearson_id.txt', 'w')
    # fw6 = open('1208_DHAC_terms_sign_list_id.txt', 'w')

    fw1.write('parent' + '\t' + 'child' + '\n')
    for edge in cliqueGraph0.edges():
        fw1.write(str(edge[0]) + '\t' + str(edge[1]) + '\n')
    fw2.write(
        'term_id' + '\t' + 'sign_score' + '\t' + 'level' + '\t' + 'annotation_gene' + '\t' + 'start_time' + '\t' + 'end_time' + '\t' + 'geneSize' + '\t' + 'time_size' + '\n')
    fw3.write('Term_Id' + '\t' + 'distance' + '\n')
    fw4.write('Term_id' + '\t' + 'Sign_score' +  '\t' + 'Gene_number' + '\t' + 'Time_point_length' + '\t' + 'leaf' + '\t' + 'Trend1' + '\t' + 'Trend2' + '\t' + 'Trend3' + '\n')
    fw5.write('Term_id' + '\t' + 'pearson' + '\n')
    fw6.write(
        'term_id' + '\t' + 'sign_score' + '\t' + 'level' + '\t' + 'annotation_gene' + '\t' + 'start_time' + '\t' + 'end_time' + '\t' + 'geneSize' + '\t' + 'time_size' + '\n')

    for node, value in sorted(dic_term_score.items(), key=lambda d: d[1],reverse=True):
        fw2.write(
            str(node) + '\t' + str(round(value, 4)) + '\t' + str(
                nx.shortest_path_length(cliqueGraph0, 0, node)) + '\t' + ','.join(
                [t for t in cliqueGraph0.node[node]['annotation']]) + '\t' + str(
                min(cliqueGraph0.node[node]['windowsize']) + 49) + '\t' + str(
                max(cliqueGraph0.node[node]['windowsize']) + 58) + '\t' +
            str(len(cliqueGraph0.node[node]['annotation'])) + '\t' + str(
                len(cliqueGraph0.node[node]['windowsize']) + 9) + '\n')
        fw3.write(str(node) + '\t' + str(dic_term_distance[node]) + '\n')
        fw4.write(str(node) + '\t' + str(round(value,4)) + '\t' + '\t'.join(str(i) for i in dic_vector[node]) + '\n')
        fw5.write(str(node) + '\t' + str(dic_pearson[node]) + '\n')
        fw6.write(
            str(node) + '\t' + str(round(value, 4)) + '\t' + str(nx.shortest_path_length(cliqueGraph0, 0, node)) + '\t' + ','.join([dic[t] for t in cliqueGraph0.node[node]['annotation']]) + '\t' + str(
                min(cliqueGraph0.node[node]['windowsize']) + 49) + '\t' + str(
                max(cliqueGraph0.node[node]['windowsize']) + 58) + '\t' +
            str(len(cliqueGraph0.node[node]['annotation'])) + '\t' + str(
                len(cliqueGraph0.node[node]['windowsize']) + 9) + '\n')

    fw1.close()
    fw2.close()
    fw3.close()
    fw4.close()
    fw5.close()
    fw6.close()

    '''
    fr_leaf = open('my_method_leaf_nodes.txt', 'r')
    fr_super = open('super_than_dhac.txt', 'r')
    title = fr_super.readline()
    dic = {}
    leaf = 0
    s = 0
    for line in fr_leaf:
        dic[line.strip().split('\n')[0]] = 1
        leaf = leaf + 1
    for line in fr_super:
        if dic.has_key(line.strip().split('\t')[0]):
            s = s + 1
    print 'total leaf:',leaf
    print 'leaf term that super than dhac', s
    # total leaf: 4186
    # leaf term that super than dhac 339
    '''
    # 得到442个dhac term
    '''
    geneList = []
    windowList = []
    dict = {}
    # geneName-geneId
    dic = {}
    fr = open('G:\project2\\NPM201507\\code\\term_name_id\\termN_Id.txt', 'r')
    for line in fr:
        term, idd = line.strip().split('\t')
        dic[idd] = term

    read_dhac_result('G:\\project2\\NPM201507\\code\\edge\\20171113dhac')
    term0 = tree0(0.9, 0, 1)

    fw = open('dhac_term_addGeneSize.txt', 'a')
    fw.write('purity' + '\t' + 'class_distance' + '\t' + 'pearson' + '\t' + 'geneSet' + '\t' + 'geneSize' + '\t' + 'min_TimePoint' + '\t' + 'max_TimePoint' + '\n')
    for i in range(1,50):
        term = tree0(0.9, i, 1)
        for key, value in term.items():
            if term0.has_key(key):
                if set(term0[key]).issuperset(set(value)):
                    continue
                else:
                    term0[key]  = term0[key] + value
            else:
                term0[key] = value
    print len(term0)
    for key, value in term0.items():
        genes_name = [dic[i] for i in key]
        fw.write(str(sign_value(genes_name, value)) + '\t' + str(class_distance(genes_name, value)) + '\t' + str(pearson(genes_name, value)) + '\t' + ','.join(t for t in key) + '\t' + str(len(key)) + '\t' + str(min(value) + 49) + '\t' + str(max(value) + 58) + '\n')
    '''

    # 比较dhac 和 my method   11/13
    '''
    fr = open("dhac_term.txt", 'r')s
    title = fr.readline()
    dic_dhac_term = {}
    for line in fr:
        line_str = line.strip().split('\t')
        geneSet = line_str[0].strip().split(',')
        windowSet = range(int(line_str[1]), int(line_str[2]))
        dic_dhac_term[frozenset(geneSet)] = windowSet
    print len(dic_dhac_term)
    fr2 = open('G:\project2\\NPM201507\\code\\1117terms_sign_list_id.txt', 'r')
    title2 = fr2.readline()
    dic_my_term = {}
    for line in fr2:
        line_str = line.strip().split('\t')
        geneSet = line_str[2].strip().split(',')
        windowSet = range(int(line_str[3]), int(line_str[4]))
        dic_my_term[frozenset(geneSet)] = windowSet
    print len(dic_my_term)
    n = 0
    s = 0
    for key,value in dic_dhac_term.items():
        if dic_my_term.has_key(key):
            n = n + 1
            print value,dic_my_term[key]
            if set(dic_my_term[key]).issuperset(set(value)):
                print 'ok'
                s = s +1
    print n
    print 'my is super than dhac:', s
    '''



