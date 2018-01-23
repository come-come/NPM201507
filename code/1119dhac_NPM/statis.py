# -*- coding: utf-8 -*-
import networkx as nx
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib
import numpy as np


#  get the color bar
color= []
for name, hex in matplotlib.colors.cnames.iteritems():
    color.append(str(name))

def draw (termData, a, colorId, legend):
    # 散点图 --count2
    b = [int(i) for i in a]  # index--int
    geneSize = []
    timeSize = []
    for group in termData.loc[b].groupby('geneSize'):
        geneSize.append(group[0])
        timeSize.append(group[1]['time_size'].mean())
        p = plt.plot(group[0], group[1]['time_size'].mean(), marker='o', c=color[colorId])
        # print group[0], group[1].shape[0], group[1]['time_size'].mean()
    p = plt.plot(group[0], group[1]['time_size'].mean(), marker='o', c=color[colorId], label=legend)

def count2(edge_filename, term_filename,method):
    # 读取的文件中有一列是结点所在的level
    print method
    fr = open(edge_filename, 'r')
    termData = pd.read_csv(term_filename, sep='\t', index_col=0)
    title = fr.readline()
    colorId = 0
    for group in termData.groupby('level'):
        # 每一层的平均gene size 和 平均time
        print group[0], '\t', round(group[1]['geneSize'].mean(), 3), '\t', round(group[1]['time_size'].mean(),3)
        a = list(group[1].index)
        # print  group[0], len(a) # 统计每一层有多少个结点
        legend = 'level' + str(colorId + 1)
        draw(termData, a, colorId, legend)
        colorId = colorId + 1
    # plt.xlabel("geneSize")
    plt.axis([0, 100, 10, 45])
    # plt.legend(bbox_to_anchor=(1.02, 3), loc=2, borderaxespad=0.)
    # plt.legend(loc=0, numpoints=1)
    plt.title(method, fontsize=16)
    # plt.show()



    #     plt.plot(group[1]['geneSize'].mean(), group[1]['time_size'].mean(), marker='o', c=color[0])
    # plt.plot(group[1]['geneSize'].mean(), group[1]['time_size'].mean(), marker='o', c=color[0], label = 'our_method')
    # plt.xlabel("geneSize")
    # plt.ylabel("average_timeSize")
    # plt.title(method)
    # plt.legend(loc=0, numpoints=1)
    # plt.show()
def everylevel(purity_my_terms):
    plt.figure(1)
    termData = pd.read_csv(purity_my_terms, sep='\t', index_col=0)
    colorId = 0
    row = 1
    for group in termData.groupby('level'):
        if group[0] < 9:
            ax = plt.subplot(8, 1, row)
        else:
            continue

        plt.sca(ax)
        a = list(group[1].index)
        legend = 'level' + str(colorId + 1)
        b = [int(i) for i in a]  # index--int
        geneSize = []
        timeSize = []
        for group1 in termData.loc[b].groupby('geneSize'):
            geneSize.append(group1[0])
            timeSize.append(group1[1]['time_size'].mean())
            p = plt.plot(group1[0], group1[1]['time_size'].mean(), marker='o', c=color[colorId])
        p = plt.plot(group1[0], group1[1]['time_size'].mean(), marker='o', c=color[colorId], label=legend)
        row = row + 1
        colorId = colorId + 1
        plt.legend(bbox_to_anchor=(1.05, 0.9), loc=2, borderaxespad=0., numpoints=1, handlelength=0)
        plt.axis([5, 100, 10, 45])
    plt.ylabel("average length of time point", position=(1.8, 4.5), fontsize=20)
    plt.xlabel("geneSize", fontsize=20)
    plt.show()

def plot(edge_filename, term_filename,color_bar, fig, label):
    fr = open(edge_filename, 'r')
    termData = pd.read_csv(term_filename, sep='\t', index_col=0)
    title = fr.readline()
    ave_gene = []; ave_time = []; level = []; numbers = []
    for group in termData.groupby('level'):
        # 每一层的平均gene size 和 平均time length
        print group[0], '\t', round(group[1]['geneSize'].mean(), 3), '\t', round(group[1]['time_size'].mean(),3),group[1].shape[0]
        if group[0]< 11:
            ave_gene.append(round(group[1]['geneSize'].mean(), 3))
            ave_time.append(round(group[1]['time_size'].mean(), 3))
            level.append(group[0])
            numbers.append(group[1].shape[0]/2.0)


    # 设置标题
    ax1.set_title('Annotation Time',size=16)
    # ax1.set_title('Annotation Gene',size=16)
    # 设置X轴标签
    plt.xlabel('level',size=13)
    plt.xticks(np.arange(20))
    # 设置Y轴标签
    plt.ylabel('average length of time point',size=13)
    # plt.ylabel('average size of annotation genes',size=13)
    # 设置点的大小
    sValue = numbers
    # 画散点图
    ax1.scatter(level, ave_time, s=sValue, c=color_bar, marker='o', label=label,alpha=0.5)
    # ax1.scatter(level, ave_gene, s=sValue, c=color_bar, marker='o', label=label, alpha=0.5)
    # 设置图标

    plt.legend(scatterpoints=1,fontsize = 'large',loc=4)

    # 显示所画的图
    return fig

if __name__ == "__main__":
    edge_dhac_filename = 'G:\project2\\NPM201507\\code\\1119dhac_NPM\\1208_DHAC_edges_sign_id.txt'
    term_dhac_filename = 'G:\project2\\NPM201507\\code\\1119dhac_NPM\\1208_DHAC_terms_sign_list_id.txt'
    edge_NPM_filename3 = 'G:\project2\\NPM201507\\code\\1119dhac_NPM\\1211_NPM10_edges_sign_id.txt'
    term_NPM_filename3 = 'G:\project2\\NPM201507\\code\\1119dhac_NPM\\1211_NPM10_terms_sign_list_id.txt'
    edge_NPM_filename2 = 'G:\project2\\NPM201507\\code\\1119dhac_NPM\\1211_NPM10_edges_sign_id(2).txt'
    term_NPM_filename2 = 'G:\project2\\NPM201507\\code\\1119dhac_NPM\\1211_NPM10_terms_sign_list_id(2).txt'
    edge_NPM_filename1 = 'G:\project2\\NPM201507\\code\\1119dhac_NPM\\1222_NPM_edges_sign_id(3).txt'
    term_NPM_filename1 = 'G:\project2\\NPM201507\\code\\1119dhac_NPM\\1222_NPM_terms_sign_list_id(3).txt'
    purity_my_terms = '1129terms.txt'
    purity_my_edges = '1129edges.txt'
    purity_my_terms1 = '1218terms.txt'
    purity_my_edges1 = '1218edges.txt'
    purity_my_terms2 = '1219terms.txt'
    purity_my_edges2 = '1219edges.txt'

    fig = plt.figure()
    ax1 = fig.add_subplot(111)
    color_dic = {1:'b', 2: 'r', 3:'y'}
    l_plot = []
    fig = plot(purity_my_edges, purity_my_terms, color_dic[1],fig, 'Our Method(1)')
    fig = plot(purity_my_edges1, purity_my_terms1, color_dic[2], fig, 'Our Method(2)')
    fig = plot(purity_my_edges2, purity_my_terms2, color_dic[3], fig, 'Our Method(3)')
    # fig = plot(edge_NPM_filename1, term_NPM_filename1, color_dic[1], fig, 'NPM(1)')
    # fig = plot(edge_NPM_filename2, term_NPM_filename2, color_dic[2], fig, 'NPM(2)')
    # fig = plot(edge_NPM_filename3, term_NPM_filename3, color_dic[3], fig, 'NPM(3)')

    plt.show()

    '''
    # everylevel(purity_my_terms)
    plt.figure(1)
    ax1 = plt.subplot(411)  # 在图表2中创建子图1
    ax2 = plt.subplot(412)
    ax3 = plt.subplot(413)
    ax4 = plt.subplot(414)
    plt.sca(ax1)
    count2(purity_my_edges2, purity_my_terms2, 'Our_Method')
    plt.sca(ax2)
    count2(edge_NPM_filename1, term_NPM_filename1, 'NPM(1)')
    # count2(purity_my_edges1, purity_my_terms1, 'Our_Method(1)')
    plt.sca(ax3)
    count2(edge_NPM_filename2, term_NPM_filename2, 'NPM(2)')
    # count2(purity_my_edges2, purity_my_terms2, 'Our_Method(2)')
    plt.sca(ax4)
    count2(edge_NPM_filename3, term_NPM_filename3, 'NPM(3)')
    # count2(purity_my_edges, purity_my_terms, 'Our_Method(3)')
    plt.ylabel("average length of time point", position=(1.8,2.4), fontsize=20)
    plt.xlabel("geneSize", fontsize=20)
    plt.legend(bbox_to_anchor=(1.04, 3.5), loc=2, borderaxespad=0., numpoints=1,  handlelength=0)
    plt.show()
    '''






