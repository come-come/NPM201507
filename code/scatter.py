# -*- coding: utf-8 -*-
"""
@author: JunYa Lu
"""
'''import scipy.io as sio
from mpl_toolkits.mplot3d import Axes3D
import matplotlib.pyplot as plt
import pandas as pd
import numpy as np'''

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

phe1 = pd.read_table('G:\project2\\NPM201507\\data\\IDMapping_consolidated_allPhi2_cleaned_lfc_avg.txt', index_col=0)
phe2 = pd.read_table('G:\project2\\NPM201507\\data\\IDMapping_consolidated_allQESV_cleaned_LFC_avg.txt', index_col=0)
phe3 = pd.read_table('G:\project2\\NPM201507\\data\\IDMapping_consolidated_allQI_new_RAW3_adj_LFC_avg.txt', index_col=0)
phe1.columns = [i for i in range(0, 113)]
phe2.columns = [i for i in range(0, 113)]
phe3.columns = [i for i in range(0, 113)]

def draw(term, gene_set, start_time, end_time):
    window_set = [i for i in range(start_time, end_time + 1)]  # 回到最原始的数据上
    q = phe1.loc[gene_set, window_set]
    q2 = phe2.loc[gene_set, window_set]
    q3 = phe3.loc[gene_set, window_set]
    sns.set(style="white")
    cmap = sns.diverging_palette(220, 10, as_cmap=True)
    # ax = sns.heatmap(q, vmin=-0.05, vmax=0.05, cmap=cmap, annot=True,center=0, square=True, linewidths=.5, cbar_kws={"shrink": .5})
    # for item in ax.get_yticklabels():
    #     item.set_rotation(0)
    fig = plt.figure()
    fig.suptitle(term)
    plt.subplot(231)
    ax = sns.heatmap(q, vmin=-0.5, vmax=0.5, cmap=cmap,  center=0, square=True, linewidths=.5,
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
    plt.title('QE')
    plt.xlabel('time point')

    plt.subplot(233)
    ax = sns.heatmap(q3, vmin=0, vmax=1, cmap=cmap, center=0, square=True, linewidths=.5,
                     cbar_kws={"shrink": .5})
    for item in ax.get_yticklabels():
        item.set_rotation(0)
    plt.title('QI')
    plt.xlabel('time point')
    plt.show()

def phenotype_draw(term, gene_set, start_time, end_time):
    window_set = [i for i in range(start_time, end_time + 1)] # 回到最原始的数据上
    q = phe1.loc[gene_set, window_set]
    q2 = phe2.loc[gene_set, window_set]
    q3 = phe3.loc[gene_set, window_set]
    phe_avg = (phe1.loc[gene_set, window_set] + phe2.loc[gene_set, window_set] + phe3.loc[gene_set, window_set]).apply(lambda x : x*2)
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
    fig.suptitle(term)
    plt.subplot(131)
    ax = sns.heatmap(q, vmin=-0.5, vmax=0.5, cmap=cmap,  center=0, square=True, linewidths=.5,
                     cbar_kws={"shrink": .5})
    for item in ax.get_yticklabels():
        item.set_rotation(0)
    plt.title('Phi2')
    plt.xlabel('time point',fontsize=8)
    plt.xticks(fontsize=8)  # 对坐标的值数值，大小限制
    plt.yticks(fontsize=8)

    plt.subplot(132)
    ax = sns.heatmap(q2, vmin=-0.5, vmax=0.5, cmap=cmap, center=0, square=True, linewidths=.5,
                     cbar_kws={"shrink": .5})
    for item in ax.get_yticklabels():
        item.set_rotation(0)
    plt.title('QE')
    #plt.xlabel('time point', fontsize=8)
    plt.xticks(fontsize=8)  # 对坐标的值数值，大小限制
    plt.yticks(fontsize=8)

    plt.subplot(133)
    ax = sns.heatmap(q3, vmin=0, vmax=1, cmap=cmap, center=0, square=True, linewidths=.5,
                     cbar_kws={"shrink": .5})
    for item in ax.get_yticklabels():
        item.set_rotation(0)
    plt.title('QI')
    plt.xlabel('time point', fontsize=8)
    plt.xticks(fontsize=8)  # 对坐标的值数值，大小限制
    plt.yticks(fontsize=8)
    filename = str(term) + '.png'
    plt.savefig(filename)
    plt.show()

def phenotype(term, gene_set, start_time, end_time):
    window_set = [i for i in range(start_time,end_time+1)]
    q = phe1.loc[gene_set, window_set]
    q2 = phe2.loc[gene_set, window_set]
    q3 = phe3.loc[gene_set, window_set]

    #fw.write('term_id:' + term + '\t' + 'annotation_gene:' + str(gene_set) + '\t' + 'time_point:' + str(start_time) + '--' + str(end_time) + '\n')
    # filename = str(term) + '.csv'
    # q.to_csv(filename, sep=',', header=True, index=True, mode='a')
    #
    # q2.to_csv(filename, sep=',', header=True, index=True, mode='a')
    #
    # q3.to_csv(filename, sep=',', header=True, index=True, mode='a')

    # colors = ['b', 'c', 'y', 'm', 'r', 'g', 'k', 'gray', 'cyan']
    # fig = plt.figure()
    # ax = fig.gca(projection='3d')
    # #ax = Axes3D(fig)
    # c = 0
    print term
    for i in gene_set:
        X = np.array(q.loc[i]).tolist()
        Y = np.array(q2.loc[i]).tolist()
        Z = np.array(q3.loc[i]).tolist()
        # l = ax.scatter(X, Y, Z, c=colors[c], label=i)  # 绘制数据点
        print i
        print X
        print Y
        print Z
        print '----------------------------------------'
        # c = c + 1
    # ax.set_zlabel('Z')  # 坐标轴
    # ax.set_ylabel('Y')
    # ax.set_xlabel('X')
    # ax.set_title('term id: ' + str(term))
    # # print min(q),min(q2),max(q3)
    # plt.legend(loc='upper left', numpoints=1, ncol=3, fontsize=8, bbox_to_anchor=(0,0))
    # ax.text(min(q)+ 0.2,min(q2)+0.2,max(q3)+0.2,'%s' % (window_set), size=8, zorder=1)
    #filename = 'G:\project2\\NPM201507\\code\\figure\\' + str(term) + '.png'
    #plt.savefig(filename)
    #plt.show()

#phenotype(88001,['AT1G14345', 'AT1G80380', 'AT2G18790', 'AT2G29180', 'AT4G33010', 'AT5G42270'], 68, 82)
# fw = open('top11.txt', 'a')

file = open('G:\project2\\NPM201507\\code\\1102terms_sign_list.txt', 'r')
title = file.readline()
for line in file.readlines()[0:12]:
    line_sp = line.strip().split('\t')
    term = line_sp[0]
    gene_set = line_sp[2]
    start_time = line_sp[3]
    end_time = line_sp[4]
    result = gene_set.strip().split(',')
    phenotype(term, result, int(start_time), int(end_time))

