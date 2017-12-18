# -*- coding: utf-8 -*-
import networkx as nx
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib
#  get the color bar
color= []
for name, hex in matplotlib.colors.cnames.iteritems():
    color.append(str(name))

def draw (termData, a, colorId, legend):
    b = [int(i) for i in a]  # index--int
    geneSize = []
    timeSize = []
    for group in termData.loc[b].groupby('geneSize'):
        geneSize.append(group[0])
        timeSize.append(group[1]['time_size'].mean())
        p = plt.plot(group[0], group[1]['time_size'].mean(), marker='o', c=color[colorId])
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
        print group[0], group[1]['geneSize'].mean(), group[1]['time_size'].mean()
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

if __name__ == "__main__":
    edge_dhac_filename = 'G:\project2\\NPM201507\\code\\1119dhac_NPM\\1208_DHAC_edges_sign_id.txt'
    term_dhac_filename = 'G:\project2\\NPM201507\\code\\1119dhac_NPM\\1208_DHAC_terms_sign_list_id.txt'

    edge_NPM_filename = 'G:\project2\\NPM201507\\code\\1119dhac_NPM\\1211_NPM10_edges_sign_id.txt'
    term_NPM_filename = 'G:\project2\\NPM201507\\code\\1119dhac_NPM\\1211_NPM10_terms_sign_list_id.txt'
    edge_NPM_filename2 = 'G:\project2\\NPM201507\\code\\1119dhac_NPM\\1211_NPM10_edges_sign_id(2).txt'
    term_NPM_filename2 = 'G:\project2\\NPM201507\\code\\1119dhac_NPM\\1211_NPM10_terms_sign_list_id(2).txt'

    purity_my_terms = '1129terms.txt'
    purity_my_edges = '1129edges.txt'

    plt.figure(1)
    ax1 = plt.subplot(411)  # 在图表2中创建子图1
    ax2 = plt.subplot(412)
    ax3 = plt.subplot(413)
    ax4 = plt.subplot(414)
    plt.sca(ax1)
    count2(edge_dhac_filename, term_dhac_filename, 'DHAC')

    plt.sca(ax2)
    count2(edge_NPM_filename, term_NPM_filename, 'NPM(1)')

    plt.sca(ax3)
    count2(edge_NPM_filename2, term_NPM_filename2, 'NPM(2)')

    plt.sca(ax4)
    count2(purity_my_edges, purity_my_terms, 'Our_Method')
    plt.ylabel("average length of time point", position=(1.8,2.4), fontsize=20)
    plt.xlabel("geneSize", fontsize=20)
    plt.legend(bbox_to_anchor=(1.04, 3.7), loc=2, borderaxespad=0., numpoints=1,  handlelength=0)
    plt.show()



    # count2(edge_NPM_filename, term_NPM_filename, 'NPM1')
    # count2(edge_NPM_filename2, term_NPM_filename2, 'NPM2')
    # count2(purity_my_edges, purity_my_terms,'MY')

