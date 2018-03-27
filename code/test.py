# -*- coding: utf-8 -*-

import networkx as nx
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
import statsmodels.api as sm
import pylab
import sys, getopt
from scipy import stats
from scipy.stats import kstest
from scipy.stats import anderson
from compiler.ast import flatten

import matplotlib.mlab as mlab

import numpy as np
import pandas as pd
from mpl_toolkits.mplot3d import Axes3D
import matplotlib.pyplot as plt
import random
from mpl_toolkits.mplot3d import Axes3D
from matplotlib import cm
from matplotlib.ticker import LinearLocator, FormatStrFormatter
import matplotlib.pyplot as plt
import numpy as np
from scipy.interpolate import griddata
from mpl_toolkits.mplot3d import Axes3D
from matplotlib import cm
from matplotlib.colors import LightSource
from sys import argv
from scipy.stats import gaussian_kde
from pykrige.ok import OrdinaryKriging
import math

filename = 'result_c5_s10_v2_weight.txt'
phe1 = pd.read_table('G:\project2\\NPM201507\\data\\IDMapping_consolidated_allPhi2_cleaned_lfc_avg.txt', index_col=0)
phe2 = pd.read_table('G:\project2\\NPM201507\\data\\IDMapping_consolidated_allQESV_cleaned_LFC_avg.txt', index_col=0)
phe3 = pd.read_table('G:\project2\\NPM201507\\data\\IDMapping_consolidated_allQI_new_RAW3_adj_LFC_avg.txt', index_col=0)
phe1.columns = [i for i in range(0, 113)]
phe2.columns = [i for i in range(0, 113)]
phe3.columns = [i for i in range(0, 113)]
xmin = 0
xmax=1
ymin = 0
ymax=1
kriging_resolution = 200
X = [-0.105625455, -0.174359608, -0.105873346, -0.20576805399999998, -0.120525806, -0.28215594, -0.093060075, -0.14982578900000001, -0.119486435, -0.151960108, -0.10104339999999999, -0.23913942, -0.106686264, -0.16127649300000002]
Y  = [0.177165099, -0.029949987, 0.145072897, -0.038238821, 0.157384682, -0.035901472000000004, 0.138161626, -0.037059772000000005, 0.151773712, -0.04529447, 0.14621367300000002, -0.036671638, 0.16619656800000002, -0.032974679]
Z = [0.166806127, 0.20898661899999998, 0.187877354, 0.226367115, 0.18765406699999998, 0.22494885899999997, 0.19316028699999999, 0.212110438, 0.204398195, 0.22953668100000002, 0.205943133, 0.231321515, 0.194322528, 0.246532118]
x = [-0.105625455, -0.174359608, -0.105873346, -0.20576805399999998, -0.120525806, -0.28215594, -0.093060075, -0.14982578900000001, -0.119486435, -0.151960108, -0.10104339999999999, -0.23913942, -0.106686264, -0.16127649300000002]
y  = [0.177165099, -0.029949987, 0.145072897, -0.038238821, 0.157384682, -0.035901472000000004, 0.138161626, -0.037059772000000005, 0.151773712, -0.04529447, 0.14621367300000002, -0.036671638, 0.16619656800000002, -0.032974679]
z = [0.166806127, 0.20898661899999998, 0.187877354, 0.226367115, 0.18765406699999998, 0.22494885899999997, 0.19316028699999999, 0.212110438, 0.204398195, 0.22953668100000002, 0.205943133, 0.231321515, 0.194322528, 0.246532118]
xy = np.vstack([x,y])
xi = np.linspace(xmin,xmax,kriging_resolution)
yi = np.linspace(ymin,ymax,kriging_resolution)
zi = griddata((x, y), z, (xi[None,:], yi[:,None]), method='cubic')
OK = OrdinaryKriging(x, y, z, variogram_model='linear')
mask = np.array([[math.isnan(x) for x in row] for row in zi])
zs, ss = OK.execute('masked', xi, yi, mask=mask, backend="loop" )
zz=np.ma.filled(zs, fill_value=0)

X, Y = np.meshgrid(xi, yi)
positions = np.vstack([X.ravel(), Y.ravel()])
values = xy
kernel = gaussian_kde(values)
zs_density = np.reshape(kernel(positions).T, X.shape)

# create 3D figure
Z = zs_density  # zs #griddata(x, y, z, xi, yi)

# generate colors
zs_cor = zs - np.amin(zs)
zs_cor = zs_cor / np.amax(zs_cor)
colors = cm.jet(zs_cor)

fig = plt.figure()
ax = Axes3D(fig)
# ax.plot_surface(X, Y, Z, facecolors=colors, rstride=1, cstride=1, linewidth=1, antialiased=True)
# ax.plot_surface(X, Y, Z, facecolors=colors, rstride=1, cstride=1, linewidth=1, antialiased=True)
surf = ax.plot_surface(X, Y, Z, facecolors=colors, rstride=1, cstride=1, linewidth=0, antialiased=False)
# surf=ax.plot_wireframe(X, Y, Z, rstride=5, cstride=5)



# Customize the z axis.
# ax.set_zlim(0, 18.5)
ax.set_xlim(xmin, xmax)
ax.set_ylim(ymin, ymax)
ax.view_init(40, 30)
# ax.view_init(90, -90)
ax.set_xlabel('rtwert')
ax.set_ylabel('rgg')
ax.set_zlabel("density")

# ax.zaxis.set_major_locator(LinearLocator(10))
# ax.zaxis.set_major_formatter(FormatStrFormatter('%.02f'))
m = cm.ScalarMappable(cmap=cm.jet)

# plt.colorbar(m, shrink=0.5, aspect=5)

# Add a color bar which maps values to colors.
# fig.colorbar(surf, shrink=0.5, aspect=5)
plt.savefig('last_krig_fig')
plt.show()
#plt.colorbar(m, shrink=0.5, aspect=5)
# Add a color bar which maps values to colors.
# fig.colorbar(surf, shrink=0.5, aspect=5)
plt.savefig('last_krig_fig')
plt.show()
print zi


def sign_value(phe1, gene_set, window):
    t_min = min(window) + 49
    t_max = max(window) + 58
    window_set = [i for i in range(t_min,t_max)] # 回到最原始的数据上

    p1 = phe1.loc[gene_set, window_set].dropna(axis=1, how='any')
    p2 = phe2.loc[gene_set, window_set].dropna(axis=1, how='any')
    p3 = phe3.loc[gene_set, window_set].dropna(axis=1, how='any')
    result = pd.concat([p1, p2, p3], axis=1)
    p_matrix = np.corrcoef(result) # 行与行之间的相关系数
    print p_matrix
    print 'tttt'
    print  p_matrix.sum(),p_matrix.shape,p_matrix.shape[0]
    print p_matrix.sum(axis=1),sum(p_matrix.sum(axis=1))
    num = (p_matrix.shape[0] * p_matrix.shape[0]-p_matrix.shape[0])/2
    print num
    sum1= ( p_matrix.sum() - p_matrix.shape[0] )/2
    print sum1
    print sum1/num
    # for i in range(0, p1.shape[0]):
    #     print p1[t_min][i], p1[t_max][i]




# plot the distrubution of phenotypes  [subplot (1,3)]
columns = phe1.shape[1]
phe1.columns = [i for i in range(0, columns)]
print phe1[0][0], phe1[columns-1][0]
data = pd.read_csv(filename, index_col=0, sep='\t')
data.columns = [i for i in range(0, 54)]

gene = ['AT1G14345', 'AT1G80380', 'AT2G18790', 'AT2G29180', 'AT4G33010']
window = [19, 20, 21, 22, 23, 24]
sign_value(phe1,gene,window)
# print kstest(phe1[0], 'norm')
# print stats.normaltest(phe1[0], axis=0)
# print anderson(phe1[0])


phe_data_raw = flatten(phe1.values.tolist())
phe_data_raw2 = flatten(phe2.values.tolist())
phe_data_raw3 = flatten(phe3.values.tolist())
phe_data = [round(x, 3) for x in phe_data_raw if str(x) != 'nan']
phe_data2 = [round(x, 3) for x in phe_data_raw2 if str(x) != 'nan']
phe_data3 = [round(x, 3) for x in phe_data_raw3 if str(x) != 'nan']

phe_data_p = [round(x, 3) for x in phe_data if not float(x) < 0]
phe_data_n = [round(x, 3) for x in phe_data if float(x) < 0]
print len(phe_data), len(phe_data_p), len(phe_data_n), phe_data[0], phe_data[20204],type(phe_data[0])

print kstest(phe_data_p, 'norm',  alternative='greater')
print stats.normaltest(phe_data, axis=0)
print anderson(phe_data)

sum1 = 0
sum2 = 0
sum3 = 0
for i in phe_data:
    if i<0.08 and i > -0.18:
        sum1 = sum1 + 1

for i in phe_data2:
    if i<0.27 and i > -0.22:
        sum2 = sum2 + 1

for i in phe_data3:
    if i<0.34:
        sum3 = sum3 + 1
print sum1, sum2, sum3


plt.subplot(231)
n, bins, patches = plt.hist(phe_data, 20)
print n
print bins
plt.title('phenotype 1')
plt.subplot(232)
n, bins, patches = plt.hist(phe_data2, 20)
print n
print bins
plt.title('phenotype 2')
plt.subplot(233)

n, bins, patches  = plt.hist(phe_data3,20)
print n
print bins
plt.title('phenotype 3')
# plt.subplot(224)
# n, bins, patches = plt.hist(phe_data_p, 150, normed=1)
# mu = np.mean(phe_data_p)
# sigma = np.std(phe_data_p)
# plt.plot(bins, mlab.normpdf(bins, mu, sigma))
plt.show()




'''
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
'''

'''
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
'''

'''
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
'''