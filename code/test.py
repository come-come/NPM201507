# -*- coding: utf-8 -*-
import networkx as nx
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
filename = 'G:\\project2\\NPM201507\\code\\0908terms(2).txt'
data = pd.read_table(filename,  index_col=0)

data_gene = data['geneSize']
data_window = data['window_size']
# 根据value_counts 的index进行排序
data_gene = data_gene.value_counts().sort_index()
data_window = data_window.value_counts().sort_index()
print 'term注释基因个数的不同有几种:', data_gene.value_counts().shape, 'window的不同长度有几种:', data_window.value_counts().shape

# data_gene.to_csv('geneSize.csv')
# data_window.to_csv('windowSize.csv')
plt.hist(np.array(data_gene.index), data_gene.values, label ='gene size')
plt.xticks(data_gene.index, data_gene.index, rotation=0)
plt.xlabel('gene size')
plt.ylabel('term number')
plt.legend()
plt.show()
plt.plot(np.array(data_window.index),data_window.values, label='window size')
plt.xticks(data_window.index, data_window.index, rotation=0)
plt.xlabel('window size')
plt.ylabel('term number')
plt.legend()
plt.show()

# plt.show()
# try:
#     f1 = open(filename, 'r')
# finally:
#     f1.close()

# G = nx.DiGraph()
# G.add_node(184, annotation=[1, 2], windowsize=[5, 6], weight=0.9)  # generate a term
# G.add_node(184, annotation=[1, 2], windowsize=[5, 6,7], weight=0.9)  # generate a term
# G.node[184]['windowsize'] = [1,8,6]
# G.add_edges_from([(1,184), (2,184)])
# print G.node[184]
# print G.predecessors(184)
# print G.nodes()
# print G.edges()
# print G.nodes()
# print G.edges()
# for edge in G.edges():
#     print edge[0], edge[1]
#
# for node in G.nodes():
#     print G.degree(node), G.in_degree(node), G.out_degree(node)



# dic = {'a':1,'b':2}
# print dic
# dic2 = dic
# print dic2

# 2017.8.25
# windowGraph = {}
# window = 0
# windowGraph[window] = nx.Graph()
#
# fw = open('countWin1.txt', 'w')
# filename = 'result_c5_s10_v2_weight.txt'
# data = pd.read_csv(filename, index_col=0, sep='\t')
# weight_value = np.arange(1.0, -0.1, -0.1)
# for  weight in weight_value:
#     df = data[data[data.columns[0]] >= (weight + 0.00001)]
#     for edge in range(0, df.shape[0]):
#         node_1, node_2 = df.index[edge].split('_')
#         # print node_1, node_2
#         windowGraph[window].add_edge(node_1, node_2)
#     print weight, df.shape[0]
#     if df.shape[0] >0 :
#
#         fw.write(str(weight) + '\t' + str(df.shape[0]) + '\t' + str(windowGraph[window].number_of_nodes() )+ '\n')
# df = data[data[data.columns[0]] == 0]
# print  '0',df.shape[0]
# for edge in range(0, df.shape[0]):
#     node_1, node_2 = df.index[edge].split('_')
#     windowGraph[0].add_edge(node_1, node_2)
# fw.write(str(0) + '\t' + str(df.shape[0]) + '\t' + str(windowGraph[window].number_of_nodes()) + '\n')
# fw.close()