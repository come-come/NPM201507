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
        p = plt.plot(group[0], group[1]['time_size'].mean(), marker='o', c = color[colorId])
    p = plt.plot(group[0], group[1]['time_size'].mean(), marker='o', c=color[colorId], label=legend)

def count(edge_filename, term_filename,method):
    fr = open(edge_filename, 'r')
    termData = pd.read_csv(term_filename, sep='\t', index_col=0)
    title = fr.readline()
    G = nx.DiGraph()
    for line in fr:
        parent, child = line.strip().split('\t')
        G.add_edge(parent, child)
    print 'G.number_of_nodes()', G.number_of_nodes()
    a = G.successors('0')
    colorId = 0
    count_list = []
    while len(a) > 0:
        print len(a)
        count_list.append(len(a))
        next_level = []
        legend = 'level' + str(colorId + 1)
        draw(termData, a, colorId, legend)
        for node in a:
            next_level.extend(G.successors(node))
        a = next_level
        colorId = colorId + 1
    print method, count_list
    plt.xlabel("geneSize")
    plt.ylabel("average_timeSize")
    # plt.legend(loc=0, numpoints=1)
    plt.title(method)
    plt.legend(loc=0, numpoints=1)
    plt.show()

def count2(edge_filename, term_filename,method):
    # 读取的文件中有一列是结点所在的level
    fr = open(edge_filename, 'r')
    termData = pd.read_csv(term_filename, sep='\t', index_col=0)
    title = fr.readline()
    colorId = 0
    for group in termData.groupby('level'):
        print group[0], group[1]['geneSize'].mean(), group[1]['time_size'].mean(), group[1].shape[0]
        a = list(group[1].index)
        print len(a)
        legend = 'level' + str(colorId + 1)
        draw(termData, a, colorId, legend)
        colorId = colorId + 1

    print method
    plt.xlabel("geneSize")
    plt.ylabel("average_timeSize")
    # plt.legend(bbox_to_anchor=(1.02, 3), loc=2, borderaxespad=0.)
    plt.legend(loc=0, numpoints=1)
    plt.title(method)
    plt.show()



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
    edge_NPM_filename = 'G:\project2\\NPM201507\\code\\1119dhac_NPM\\1208_NPM_edges_sign_id.txt'
    term_NPM_filename = 'G:\project2\\NPM201507\\code\\1119dhac_NPM\\1208_NPM_terms_sign_list_id.txt'
    edge_our_filename = 'G:\\project2\\NPM201507\\code\\1208edges_sign_id.txt'
    term_our_filename = 'G:\\project2\\NPM201507\\code\\1208terms_sign_list_id.txt'
    # plt.figure(1)
    # ax1 = plt.subplot(311)  # 在图表2中创建子图1
    # ax2 = plt.subplot(312)
    # ax3 = plt.subplot(313)
    # plt.sca(ax1)
    # count(edge_dhac_filename, term_dhac_filename, 'DHAC')
    # plt.sca(ax2)
    # count(edge_NPM_filename, term_NPM_filename, 'NPM')
    # plt.sca(ax3)
    # plt.legend(loc=0, numpoints=1)
    # plt.show()


    print 'DHAC'
    count2(edge_dhac_filename, term_dhac_filename, 'DHAC')
    print 'NPM：'
    count2(edge_NPM_filename, term_NPM_filename, 'NPM')
    print 'Our:'
    count2(edge_our_filename, term_our_filename, 'Our_Method')

