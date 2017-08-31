import networkx as nx
import pandas as pd
import numpy as np


a = [1,2,3,4]
b = [1,2,3,4]
print  set(b).issubset(set(a))
print  set(b).issuperset(set(a))
t = ['AT2G17630', 'AT4G38380', 'AT3G10470', 'AT3G46780', 'AT5G53170', 'AT4G31560', 'AT1G07510', 'AT1G75690', 'AT3G26580', 'AT1G42960', 'AT1G14150', 'AT5G52970', 'AT1G12250', 'AT5G19500', 'AT5G42070', 'AT1G02560', 'AT5G13120', 'AT5G16660', 'AT4G04020', 'AT2G43010', 'AT5G39830', 'AT1G74710', 'AT3G21055', 'AT1G03160', 'AT1G66670', 'AT3G60370', 'AT3G28850', 'AT1G54500', 'AT2G40100', 'AT1G54050', 'AT3G22690', 'AT1G51350', 'AT4G39460', 'AT5G59250', 'AT5G12470', 'AT1G14345', 'AT1G80030', 'AT1G56500', 'AT1G50460', 'AT4G19830', 'AT1G17850', 'AT1G29920', 'AT2G35040', 'AT5G20140', 'AT5G15250', 'AT4G13010', 'AT2G24020', 'AT5G42270', 'AT1G15820', 'AT2G37660', 'AT4G18370', 'AT1G67700', 'AT2G29500', 'AT1G67280', 'AT2G27290', 'AT2G40760', 'AT1G72640', 'AT3G58010', 'AT4G13670', 'AT5G01500', 'AT1G06690', 'AT5G19370', 'AT2G30950', 'AT1G15140', 'AT1G55370', 'AT2G20260', 'AT3G55330', 'AT1G52220', 'AT3G47060', 'AT1G79230', 'AT4G16155', 'AT2G20270', 'AT1G51110', 'AT5G16150', 'AT5G08540', 'AT5G03455', 'AT3G57680', 'AT2G44990', 'AT4G20820', 'AT4G35770', 'AT2G41330', 'AT2G26800']
m = ['AT2G17630', 'AT4G38380', 'AT3G10470', 'AT3G46780', 'AT5G53170', 'AT4G31560', 'AT1G07510', 'AT1G75690', 'AT3G26580', 'AT1G42960', 'AT1G14150', 'AT5G52970', 'AT1G12250', 'AT5G19500', 'AT5G42070', 'AT1G02560', 'AT5G13120', 'AT5G16660', 'AT4G04020', 'AT2G43010', 'AT5G39830', 'AT1G74710', 'AT4G16155', 'AT1G03160', 'AT1G66670', 'AT3G60370', 'AT3G28850', 'AT3G58010', 'AT2G40100', 'AT1G54050', 'AT3G22690', 'AT1G51350', 'AT4G39460', 'AT5G59250', 'AT5G12470', 'AT1G14345', 'AT1G80030', 'AT1G56500', 'AT1G50460', 'AT4G19830', 'AT1G17850', 'AT1G29920', 'AT2G35040', 'AT5G15250', 'AT4G13010', 'AT2G24020', 'AT5G42270', 'AT1G15820', 'AT2G37660', 'AT4G18370', 'AT1G67700', 'AT4G20820', 'AT2G29500', 'AT4G33010', 'AT1G67280', 'AT2G27290', 'AT2G40760', 'AT1G72640', 'AT5G20140', 'AT4G13670', 'AT5G01500', 'AT1G06690', 'AT5G19370', 'AT2G30950', 'AT1G15140', 'AT1G55370', 'AT2G20260', 'AT3G55330', 'AT1G52220', 'AT3G08920', 'AT1G79230', 'AT2G20270', 'AT1G51110', 'AT5G16150', 'AT5G08540', 'AT5G03455', 'AT3G57680', 'AT2G44990', 'AT3G47060', 'AT4G35770', 'AT2G41330', 'AT2G26800']
x = ['AT2G17630', 'AT4G38380', 'AT3G10470', 'AT3G46780', 'AT5G53170', 'AT4G31560', 'AT1G07510', 'AT1G75690', 'AT3G26580', 'AT1G42960', 'AT1G14150', 'AT5G52970', 'AT1G12250', 'AT5G19500', 'AT5G42070', 'AT1G02560', 'AT5G13120', 'AT5G16660', 'AT4G04020', 'AT2G43010', 'AT5G39830', 'AT1G74710', 'AT4G16155', 'AT1G03160', 'AT1G66670', 'AT3G60370', 'AT3G28850', 'AT3G58010', 'AT2G40100', 'AT1G54050', 'AT3G22690', 'AT1G51350', 'AT4G39460', 'AT5G59250', 'AT5G12470', 'AT1G14345', 'AT1G80030', 'AT1G56500', 'AT1G50460', 'AT4G19830', 'AT1G17850', 'AT1G29920', 'AT2G35040', 'AT5G15250', 'AT4G13010', 'AT2G24020', 'AT5G42270', 'AT1G15820', 'AT2G37660', 'AT4G18370', 'AT1G67700', 'AT4G20820', 'AT2G29500', 'AT4G33010', 'AT1G67280', 'AT2G27290', 'AT2G40760', 'AT1G72640', 'AT5G20140', 'AT4G13670', 'AT5G01500', 'AT1G06690', 'AT5G19370', 'AT2G30950', 'AT1G15140', 'AT1G55370', 'AT2G20260', 'AT3G55330', 'AT1G52220', 'AT3G08920', 'AT1G79230', 'AT2G20270', 'AT1G51110', 'AT5G16150', 'AT5G08540', 'AT5G03455', 'AT3G57680', 'AT2G44990', 'AT3G47060', 'AT4G35770', 'AT2G41330', 'AT2G26800']
print set(t)-set(m)
print set(m) - set(t)
print set(x) - set(t)
print set(x) - set(m)
print set(x)==set(m), set(x) == set(t)
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