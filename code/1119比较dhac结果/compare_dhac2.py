import pandas as pd
import numpy as np
import os
from os.path import join

dic_result = {}  # dhac result
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
cliqueGraph = nx.DiGraph()
def tree0(weight_value, startwindow, term):
    dic_term = {}
    dic_last_time = {}
    dic_temp = {}
    dic_term_num = {}
    dic_intersect_level = {}

    for window in range(startwindow, 50):
        dic_intersect_level.clear()
        if window == startwindow:
            for clique in dic_result[window]:
                    cliqueGraph.add_node(term, annotation=list(clique), windowsize=[window])  # generate a term
                    dic_term[frozenset(clique)] = [window]  # dic_term 记录 window和clique
                    dic_term_num[frozenset(clique)] = term  # dic_term_num
                    dic_last_time[frozenset(clique)] = [window]  # dic_last_time
                    term = term + 1
        else:
            for clique in dic_result[window]:
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


if __name__ == "__main__":
    geneList = []
    windowList = []
    dict = {}
    read_dhac_result('G:\\project2\\NPM201507\\code\\edge\\20171113dhac')
    # print dic_result[1]
    print len(dic_result)
    for i in range(0,49):
        print len(dic_result[i])


            # clique_result = open('G:\\project2\\NPM201507\\code\\1117terms_sign_list_id.txt','r')
    # title = clique_result.readline()
    #
    # dic_dhac_term = {}
    # dic_my_term= {}
    #
    # for line in clique_result:
    #     s = line.strip().split('\t')
    #     genes = s[2].strip().split(',')
    #     windows = [i for i in range(int(s[3]) - 49, int(s[4]) - 58 + 1)]
    #
    #     term = int(s[0])
    #     dic_my_term[term] = windows
    #     print type(term)
    #     dic_dhac_term[term] = []
    #     for window in range(0, 50):
    #         for value in dic_result[window]:
    #             if set(genes).issubset(set(value)):
    #                 dic_dhac_term[term].append(window)
    #                 break
    #     print dic_my_term[term], dic_dhac_term[term]
