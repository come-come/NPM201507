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

def clique3(startwindow):
    print ('start window %d:' % startwindow)
    # 自底向上       2017.7.19
    filename = 'result_c5_s10_v2_weight.txt'
    data = pd.read_csv(filename, index_col=0, sep='\t')
    windowGraph = {}
    cliqueGraph = nx.DiGraph()
    dic_term = {}
    dic_wind = {}  
    term = 183
    for window in range(0, data.shape[1]):
        windowGraph[window] = nx.Graph()
    #weight_value = np.arange(1.0, -0.01, -0.01)
    weight_value = np.arange(1.0, -0.05, -0.05)
    # weight_value.sort(reverse=True) #[1, 0.99, 0.98,...0]
    dic_now = {}
    dic_now_temp = {}
    dic_compare = {}
    dic_dup = {}
    w = data.shape[1]
    filename1 = 'termInfo' + str(startwindow) + '.txt'
    filename2 = 'edgeInfo' + str(startwindow) + '.txt'
    fw = open(filename1, 'w')
    fw.write('threshold' + '\t' + 'term' + '\t' + 'GeneSize' + '\t' + 'WindowSize' + '\t' + 'key' + '\t' + 'value' + '\n')
    
    for t in weight_value:
        # print 't',t,'-------------------------------','w',w
        flag = startwindow - 1
        for window in range(startwindow, w):
            if flag == window-1:
                df = data[(data[data.columns[window]] >= (t - 0.05 - 0.00001)) & (data[data.columns[window]] <= (t + 0.00001))]
                for edge in range(0, df.shape[0]):  # get each row(gene)
                    node_1, node_2 = df.index[edge].split('_')
                    windowGraph[window].add_edge(node_1, node_2)  # generate WindowGraph
                if window == startwindow:
                    for clique in nx.find_cliques(windowGraph[window]):
                        if len(clique) > 4:
                            dic_now[frozenset(clique)] = [window]  # 用于比较的当前t的clique key--基因集和 value--window区间
                            flag = window
                else:
                    for key, value in dic_now.items():  # key 是基因集合 value是时间区间
                        # maxx = 0
                        # tag = []
                        for clique in nx.find_cliques(windowGraph[window]):
                            intersect = sorted(set(key).intersection(set(clique)))  # 求交集
                            # if dic_now_temp.has_key(frozenset(intersect):
                            if not dic_now_temp.has_key(frozenset(intersect)) and len(intersect) > 4 and window == max(value)+1:  # 交集基因大于4个
                                dic_compare[frozenset(intersect)] = value + [window]                                                                         
                            else:
                                continue
                        if len(dic_compare) > 0:  # 有交集,找出最大的交集
                            max_key = sorted(dic_compare.items(), key=lambda d: len(d[0]), reverse=True)[0][0]
                            max_value = sorted(dic_compare.items(), key=lambda d: len(d[0]), reverse=True)[0][1]
                            # for key1, value1 in dic_compare.items():
                            #     if len(key1) > maxx:
                            #         maxx = len(key1)
                            #         tag = key1
                            dic_now_temp[max_key] = max_value
                            flag = window                            
                        else:  # 无交集
                            dic_now_temp[key] = value  
                        dic_compare.clear()                           
                    # print 'len dic_now \dic_now_temp',len(dic_now), len(dic_now_temp)
                    if len(dic_now_temp) > 0:
                        dic_now.clear()
                        
                        for key1, value1 in dic_now_temp.items():
                            for key, value in dic_now_temp.items():
                                if key1 != key and set(key1).issubset(set(key)) and set(value1).issubset(value):
                                    dic_dup[frozenset(key1)] = value1

                        #print 'before pop', len(dic_now_temp)
                        for key, value in dic_dup.items():
                            dic_now_temp.pop(key)
                        #print 'after pop', len(dic_now_temp)
                        dic_dup.clear()
                        
                        for key, value in dic_now_temp.items():
                            dic_now[key] = value   
                            
                        dic_now_temp.clear()                                   
            else:
                w = flag
                break

        for key, value in dic_now.items():
            if startwindow in value and len(value) > 1 and len(key) > 4 : #包含t=startwindow的 时间区间大于1并且geneSize大于4
                if sorted(key) in dic_term.values() and sorted(value) == dic_wind[list(dic_term.keys())[list(dic_term.values()).index(sorted(key))]]:
                    continue
                else:
                    dic_term[term] = sorted(key)
                    dic_wind[term] = sorted(value)                    
                    cliqueGraph.add_node(term, annotation=list(key), windowsize=value, weight=t)# generate a term
                    genes = ','.join(key)
                    values = ','.join(str(i) for i in value)
                    fw.write(str(t) + '\t' + str(term) +'\t' + str(len(key)) + '\t' + str(len(value)) + '\t' + genes + '\t' + values + '\n')
                    term = term + 1 
        dic_now.clear()

    fw1 = open(filename2, 'w')
    fw1.write('Parent' + '\t' + 'Child' + '\n')
    print ('window %d term number %d.' % (startwindow, (term-183)))

    for key, value in sorted(dic_term.items(), key=lambda d:d[0]):  # sorted by term id 183,184...
        for parent in range(key+1, term): # 从下一个term开始寻找父亲结点
            # print key,i,value,dic_term[i],dic_wind[i],dic_wind[i],dic_wind[key]
            if set(value).issubset(set(dic_term[parent])) and set(dic_wind[parent]).issubset(set(dic_wind[key])):
                cliqueGraph.add_edge(parent, key)
                fw1.write(str(parent)+'\t'+str(key)+'\n')
                break
    fw1.close()          
    fw.close()
    print ('ending window %d ' % startwindow)

 
def asy():
    pool = gevent.pool.Pool(20)
    for window in xrange(4):
        pool.add(gevent.spawn(clique3,window))
    pool.join()
    print ('end')
    

    
if __name__ == '__main__':   
    start = time.clock()
    ###########################running time: 734.632
    jobs = []
    for i in xrange(53):
        p = multiprocessing.Process(target=clique3, args=(i,))
        jobs.append(p)
        p.start()
    for p in jobs:
        p.join()
    ###########################running time: 734.632

    # p = Pool(4)
    # for i in xrange(50):
    #     p.apply_async(clique3, args=(i,))
    # p.close()
    # p.join()
    # print('All subprocesses done.')

    end = time.clock()
    print 'The function run time is : %.03f seconds' % (end-start)
