# -*- coding: utf-8 -*-
import pandas as pd
import numpy as np
import os
from os.path import join
from compiler.ast import flatten
import networkx as nx

fr = open('G:\project2\\NPM201507\\code\\term_name_id\\termN_Id.txt', 'r')
dic = {}
for line in fr:
    term, idd = line.strip().split('\t')
    dic[idd] = term
def read_tree(filename,termFile):
    G = nx.DiGraph()
    fr = open(filename, 'r')
    title = fr.readline()
    fr2 = open(termFile, 'r')
    title2 = fr2.readline()
    for line in fr:
        G.add_edge(line.strip().split('\t')[0], line.strip().split('\t')[1])
    for line in fr2:
        line_str = line.strip().split('\t')
        termId = line_str[0]
        annotation_gene = line_str[3].strip().split(',')
        annotation_time = [int(line_str[4]), int(line_str[5])]
        G.node[termId]['geneSet'] = annotation_gene
        G.node[termId]['windowSet'] = annotation_time
    return G

def read_level_dic(termFile, level, OntologyGraph):

    fr = open(termFile, 'r')
    title = fr.readline()
    dic_gene = {}
    dic_time = {}
    dic_level = {}
    for line in fr:
        line_str = line.strip().split('\t')
        term = str(line_str[0])
        dic_gene[term] = line_str[3].strip().split(',')
        dic_time[term] = [int(line_str[4]), int(line_str[5])]
        try:
            dic_level[int(line_str[2])].append(term)
        except:
            dic_level[int(line_str[2])] = []
            dic_level[int(line_str[2])].append(term)
    sum = 0
    print 'size of OntologyGraph', OntologyGraph.number_of_nodes(), OntologyGraph.number_of_edges()
    for i in range(0, len(dic_level[level])):
        for j in range(i+1, len(dic_level[level])):
            term1 = dic_level[level][i]
            term2 = dic_level[level][j]
            if len(set(dic_gene[term2])|(set(dic_gene[term1])))-len(set(dic_gene[term2])&(set(dic_gene[term1]))) == 2 and dic_time[term1] == dic_time[term2]:

                try:
                    parent2 = OntologyGraph.predecessors(term2)
                    child2 = OntologyGraph.successors(term2)

                    OntologyGraph.node[term1]['geneSet'] = list(set(dic_gene[term2]) | (set(dic_gene[term1])))
                    for p in parent2:
                        OntologyGraph.add_edge(p, term1)
                    for c in child2:
                        OntologyGraph.add_edge(term1, c)
                    OntologyGraph.remove_node(term2)
                    # print term2
                    sum = sum + 1
                except:
                    continue

    print 'after,,,size of OntologyGraph', OntologyGraph.number_of_nodes(), OntologyGraph.number_of_edges()
    print 'redundant for level (difference==1): ',level, sum,'/',len(dic_level[level])

def read_level(termFile, level,OntologyGraph):
    data = pd.read_csv(termFile, sep='\t')
    level_data = data[data['level'].isin([level])].reset_index(drop=True)
    remove = 0
    # print level_data['annotation_gene'].loc[0]
    # print level_data.loc[0]['annotation_gene'].strip().split(',')
    for term1 in range(0, level_data.shape[0]):
        for term2 in range(term1+1, level_data.shape[0]):
            geneSet1 = level_data.loc[term1]['annotation_gene'].strip().split(',')
            geneStartTime1 = level_data.loc[term1]['start_time']
            geneEndTime1 = level_data.loc[term1]['end_time']
            geneSet2 = level_data.loc[term2]['annotation_gene'].strip().split(',')
            geneStartTime2 = level_data.loc[term2]['start_time']
            geneEndTime2 = level_data.loc[term2]['end_time']
            if len(set(geneSet2).difference(
                    set(geneSet1))) < 2 and geneEndTime1 == geneEndTime2 and geneStartTime1 == geneStartTime2:
                # print term1, geneSet1
                # print term2, geneSet2, 'q'
                remove = remove + 1
    print remove

def sign_value(node, geneIdset, window):
    t_min = min(window)
    t_max = max(window) + 1
    window_set = [i for i in range(t_min,t_max)] # 回到最原始的数据上

    gene_set = [dic[i] for i in geneIdset]
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
    score = float(sign)/all
    return score


if __name__ =="__main__":
    # edge_filename = 'G:\\project2\\NPM201507\\code\\1129edges_sign_id.txt'
    # node_filename = 'G:\\project2\\NPM201507\\code\\1129terms_sign_list_id.txt'
    # Run our method one more time 1218
    edge_filename = 'G:\\project2\\NPM201507\\code\\0519edges_sign_id.txt'
    node_filename = 'G:\\project2\\NPM201507\\code\\0519terms_sign_list_id.txt'

    phe1 = pd.read_table('G:\project2\\NPM201507\\data\\IDMapping_consolidated_allPhi2_cleaned_lfc_avg.txt',
                         index_col=0)
    phe2 = pd.read_table('G:\project2\\NPM201507\\data\\IDMapping_consolidated_allQESV_cleaned_LFC_avg.txt',
                         index_col=0)
    phe3 = pd.read_table('G:\project2\\NPM201507\\data\\IDMapping_consolidated_allQI_new_RAW3_adj_LFC_avg.txt',
                         index_col=0)
    phe1.columns = [i for i in range(0, 113)]
    phe2.columns = [i for i in range(0, 113)]
    phe3.columns = [i for i in range(0, 113)]
    OntologyGraph = read_tree(edge_filename,node_filename)
    # for i in range(1, 17):
    for i in range(1, 17):
        read_level_dic(node_filename, i, OntologyGraph)
    print OntologyGraph.number_of_nodes()

    fw = open('0519terms.txt', 'w')
    fw2 = open('0519edges.txt', 'w')

    fw.write('term_id' + '\t' + 'sign_score' + '\t' + 'level' + '\t' + 'annotation_gene' + '\t' + 'start_time' + '\t' + 'end_time' + '\t' + 'geneSize' + '\t' + 'time_size' + '\n')
    for node in OntologyGraph.nodes():
        if node == '0':
            continue
        else:
            fw.write(node + '\t' +
                     str(round(sign_value(node, OntologyGraph.node[node]['geneSet'],
                                    OntologyGraph.node[node]['windowSet']),5))+ '\t' +
                     str(nx.shortest_path_length(OntologyGraph, '0', node)) + '\t' +
                     ','.join([dic[t] for t in OntologyGraph.node[node]['geneSet']]) + '\t' +
                     str(min(OntologyGraph.node[node]['windowSet'])) + '\t' +
                     str(max(OntologyGraph.node[node]['windowSet'])) + '\t' +
                     str(len(OntologyGraph.node[node]['geneSet'])) + '\t' +
                     str(max(OntologyGraph.node[node]['windowSet']) - min(
                         OntologyGraph.node[node]['windowSet']) + 1) + '\n')

    fw.close()
    fw2.write('parent' + '\t' + 'child' + '\n')
    for edge in OntologyGraph.edges():
        fw2.write(edge[0] + '\t' + edge[1] + '\n')
    fw2.close()




