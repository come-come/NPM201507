# -*- coding: utf-8 -*-
import pandas as pd
import numpy as np
import os
from os.path import join
from compiler.ast import flatten
dic_result = {}  # dhac result
phe1 = pd.read_table('G:\project2\\NPM201507\\data\\IDMapping_consolidated_allPhi2_cleaned_lfc_avg.txt', index_col=0)
phe2 = pd.read_table('G:\project2\\NPM201507\\data\\IDMapping_consolidated_allQESV_cleaned_LFC_avg.txt', index_col=0)
phe3 = pd.read_table('G:\project2\\NPM201507\\data\\IDMapping_consolidated_allQI_new_RAW3_adj_LFC_avg.txt', index_col=0)

phe1.columns = [i for i in range(0, 113)]
phe2.columns = [i for i in range(0, 113)]
phe3.columns = [i for i in range(0, 113)]
def read_dhac_result(dest):
    for root, dirs, files in os.walk(dest):
        for OneFileName in files:
            if OneFileName.find('.group') == -1:
                continue
            OneFullFileName = join(root, OneFileName)
            # print OneFullFileName
            # print OneFileName.strip().split('edgeInfo')[1].split('.group')[0]
            dic_result [int(OneFileName.strip().split('edgeInfo')[1].split('.group')[0])] = []
            for line in open(OneFullFileName, 'r'):
                b = line.strip().split('\t')
                # print b  # ['821137', '829490']
                if len(b) > 4:
                    dic_result[int(OneFileName.strip().split('edgeInfo')[1].split('.group')[0])].append(b)
term = 1

def tree0(weight_value, startwindow, term):
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
def pearson(gene_set, window):
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
def class_distance(gene_set, window):
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
def sign_value(gene_set, window):
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
    # return sign, float(sign)/all, trend1, trend2, trend3
    return float(sign)/all
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
    statis()
    # average_compare()

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



