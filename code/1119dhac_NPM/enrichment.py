import os
from  os.path import join
import pandas as pd
import numpy as np
from matplotlib import pyplot as plt
from matplotlib import cm
from matplotlib import axes

def draw_heatmap(figure, data,xlabels,ylabels,i,title):
    cmap = cm.Blues

    ax=figure.add_subplot(3,1,i,position=[0.15,0.15,0.8,0.8])
    ax.set_yticks(range(len(ylabels)))
    ax.set_yticklabels(ylabels)
    ax.set_xticks(range(len(xlabels)))
    ax.set_xticklabels(xlabels)
    plt.title(title)
    vmax=data[0][0]
    vmin=data[0][0]
    for i in data:
        for j in i:
            if j>vmax:
                vmax=j
            if j<vmin:
                vmin=j
    map=ax.imshow(data,interpolation='nearest',cmap=cmap,aspect='auto',vmin=0,vmax=32)
    cb=plt.colorbar(mappable=map,cax=None,ax=None,shrink=0.5)

def sub_heatmap_numbers():
    termFile = {1: 'G:\\project2\\NPM201507\\code\\1119dhac_NPM\\1129terms.txt',
                2: 'G:\\project2\\NPM201507\\code\\1119dhac_NPM\\1208_DHAC_terms_sign_list_id.txt',
                3: 'G:\\project2\\NPM201507\\code\\1119dhac_NPM\\1211_NPM10_terms_sign_list_id.txt'}
    title = {1: 'Our Method', 2: 'DHAC', 3: 'NPM'}

    figure = plt.figure(facecolor='w')
    for tf in range(1, 4):
        triDict = {}
        rowSet = set([])
        colSet = set([])
        termData = pd.read_csv(termFile[tf], sep='\t', index_col=0)
        for group in termData.groupby(['geneSize', 'time_size']):
            triDict[(group[0][0], group[0][1])] = group[1].shape[0]
            rowSet.add(group[0][0])
            colSet.add(group[0][1])
            print group[0], group[1].shape[0]
        m = len(rowSet);
        n = len(colSet)
        print colSet
        print rowSet
        print min(colSet), max(colSet), min(rowSet), max(rowSet)
        print n, m
        rowSet = sorted(list(rowSet));
        colSet = sorted(list(colSet))
        res = np.zeros((m, n))
        dictKeys = triDict.keys()
        for i, rowNum in enumerate(rowSet):
            for j, colNum in enumerate(colSet):
                if (rowNum, colNum) in dictKeys:
                    res[i][j] = triDict[(rowNum, colNum)]
        # np.savetxt('res.txt', res, fmt='%d')
        xlabels = list(colSet)
        ylabels = list(rowSet)
        draw_heatmap(figure, res, xlabels, ylabels, tf, title[tf])
    plt.show()

def sub_heatmap_percentage():
    print 'percentage'


if __name__ == "__main__":
    # write the number of enrichment GO terms into dictionary
    dic_go = {}
    # Our Method
    path = 'G:\\project2\\NPM201507\\code\\luguilin\\2\\data_need\\'
    # DHAC
    # path = 'G:\\project2\\NPM201507\\code\\luguilin\\3\\data_need\\'
    # NPM
    # path = 'G:\\project2\\NPM201507\\code\\luguilin\\4\\data_need\\'
    for root, dirs, files in os.walk(path):
        for f in files:
            filename = join(root, f)
            dic_go[int(f.strip().split('.')[0])] = len(open(filename, 'r').readlines())-1
    data = pd.DataFrame.from_dict(dic_go, orient='index')
    data.columns = ['enrichment']
    # read level information
    termFile = 'G:\\project2\\NPM201507\\code\\1119dhac_NPM\\1129terms.txt'
    # DHAC
    # termFile = 'G:\\project2\\NPM201507\\code\\1119dhac_NPM\\1208_DHAC_terms_sign_list_id.txt'
    # NPM
    # termFile = 'G:\\project2\\NPM201507\\code\\1119dhac_NPM\\1211_NPM10_terms_sign_list_id.txt'
    termData = pd.read_csv(termFile, sep='\t',index_col=0)
    # concat two data
    result = pd.concat([data, termData], axis=1)
    for group in result.groupby('level'):
        df1 = group[1]
        zeroNum = df1[df1['enrichment'] == 0].shape[0]
        total = group[1]['enrichment'].shape[0]
        nonzero = total - zeroNum
        percentage = nonzero/float(total)
        # print group[0], '\t', nonzero, '\t', total, '\t', percentage
    sub_heatmap_numbers()
    sub_heatmap_pecentage()






