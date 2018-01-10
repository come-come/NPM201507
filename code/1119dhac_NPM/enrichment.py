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
    # ax.set_yticks(range(len(ylabels)))
    # ax.set_yticklabels(ylabels)
    # ax.set_xticks(range(len(xlabels)))
    # ax.set_xticklabels(xlabels)
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
    # map = ax.imshow(data, interpolation='nearest', cmap=cmap, aspect='auto', vmin=0, vmax=1)
    cb=plt.colorbar(mappable=map,cax=None,ax=None,shrink=0.5)


def sub_heatmap_numbers():
    termFile = {1: 'G:\\project2\\NPM201507\\code\\1119dhac_NPM\\1129terms.txt',
                2: 'G:\\project2\\NPM201507\\code\\1119dhac_NPM\\1219terms.txt',
                3: 'G:\\project2\\NPM201507\\code\\1119dhac_NPM\\1218terms.txt',
                4: 'G:\project2\\NPM201507\\code\\1119dhac_NPM\\1222_NPM_terms_sign_list_id(3).txt',
                5: 'G:\\project2\\NPM201507\\code\\1119dhac_NPM\\1211_NPM10_terms_sign_list_id(2).txt',
                6: 'G:\\project2\\NPM201507\\code\\1119dhac_NPM\\1211_NPM10_terms_sign_list_id.txt',
                7: 'G:\\project2\\NPM201507\\code\\1119dhac_NPM\\1208_DHAC_terms_sign_list_id.txt',}
    title = {1: 'Our Method(1)',2: 'Our Method(2)',3: 'Our Method(3)', 4: 'NPM(1)', 5: 'NPM(2)', 6: 'NPM(3)', 7: 'DHAC'}

    figure = plt.figure(facecolor='w')
    for tf in range(4, 7):
        triDict = {}
        rowSet = set([])
        colSet = set([])
        termData = pd.read_csv(termFile[tf], sep='\t', index_col=0)
        for group in termData.groupby(['geneSize', 'time_size']):
            print group[0]
            triDict[(group[0][0], group[0][1])] = group[1].shape[0]
            rowSet.add(group[0][0])
            colSet.add(group[0][1])
            # print group[0], group[1].shape[0]
        m = len(rowSet);
        n = len(colSet)
        print 'rowSet(gene) min:', min(rowSet), 'max', max(rowSet), 'rowSet', 'length:', m, rowSet
        print 'colSet(time) min:', min(colSet), 'max', max(colSet), 'colSet', 'length:', n, colSet

        rowSet = sorted(list(rowSet));
        colSet = sorted(list(colSet))
        res = np.zeros((m, n))
        dictKeys = triDict.keys()
        for i, rowNum in enumerate(rowSet):
            for j, colNum in enumerate(colSet):
                if (rowNum, colNum) in dictKeys:
                    res[i][j] = triDict[(rowNum, colNum)]
        # np.savetxt('res.txt', res, fmt='%d')
        df = pd.DataFrame(res)

        df.to_csv('res2.txt',sep = '\t')
        xlabels = list(colSet)
        ylabels = list(rowSet)
        i = tf % 3
        if i ==0: i = i + 3
        draw_heatmap(figure, res, xlabels, ylabels, i, title[tf])
    plt.xlabel("annotation time point")
    plt.ylabel("annotation gene Size (5-134) ",position=(1.1,1.6))

    plt.show()

def sub_heatmap_percentage():
    termFiles = {1: 'G:\\project2\\NPM201507\\code\\1119dhac_NPM\\1129terms.txt',
                2: 'G:\\project2\\NPM201507\\code\\1119dhac_NPM\\1208_DHAC_terms_sign_list_id.txt',
                3: 'G:\\project2\\NPM201507\\code\\1119dhac_NPM\\1211_NPM10_terms_sign_list_id.txt'}
    title = {1: 'Our Method', 2: 'DHAC', 3: 'NPM'}
    paths = {1: 'G:\\project2\\NPM201507\\code\\luguilin\\2\\data_need\\',
             2:'G:\\project2\\NPM201507\\code\\luguilin\\3\\data_need\\',
             3: 'G:\\project2\\NPM201507\\code\\luguilin\\4\\data_need\\'}
    figure = plt.figure(1)
    for tf in range(1,4):
        dic_go = {}
        paths = path[tf]
        for root, dirs, files in os.walk(path):
            for f in files:
                filename = join(root, f)
                dic_go[int(f.strip().split('.')[0])] = len(open(filename, 'r').readlines()) - 1
        data = pd.DataFrame.from_dict(dic_go, orient='index')
        data.columns = ['enrichment']
        # read level information
        termFile = termFiles[tf]
        termData = pd.read_csv(termFile, sep='\t', index_col=0)
        # concat two data
        termData = pd.concat([data, termData], axis=1)
        triDict = {}
        rowSet = set([])
        colSet = set([])
        for group in termData.groupby(['geneSize', 'time_size']):
            df1 = group[1]
            zeroNum = df1[df1['enrichment'] == 0].shape[0]
            total = group[1]['enrichment'].shape[0]
            nonzero = total - zeroNum
            percentage = round(nonzero / float(total), 3)
            triDict[(group[0][0], group[0][1])] = percentage
            rowSet.add(group[0][0])
            colSet.add(group[0][1])
            # print group[0], group[1].shape[0]
        m = len(rowSet);
        n = len(colSet)
        print 'rowSet(gene) min:', min(rowSet), 'max', max(rowSet), 'rowSet', 'length:', m, rowSet
        print 'colSet(time) min:', min(colSet), 'max', max(colSet), 'colSet', 'length:', n, colSet

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

    plt.xlabel("annotation time point")
    plt.ylabel("enrichment percentage ", position=(1.1, 1.6))
    plt.show()



if __name__ == "__main__":
    # write the number of enrichment GO terms into dictionary
    dic_go = {}
    # Our Method
    # path = 'G:\\project2\\NPM201507\\code\\luguilin\\2\\data_need\\'
    # path = 'G:\\project2\\enrichment\\enrich20180103_1129my\\'
    # DHAC
    # path = 'G:\\project2\\NPM201507\\code\\luguilin\\3\\data_need\\'
    # path = 'G:\\project2\\enrichment\\enrich20180104_1208DHAC\\'
    # NPM
    # path = 'G:\\project2\\NPM201507\\code\\luguilin\\4\\data_need\\'
    path = 'G:\\project2\\enrichment\\enrich20180103_1219my\\'
    for root, dirs, files in os.walk(path):
        for f in files:
            filename = join(root, f)
            dic_go[int(f.strip().split('.')[0])] = len(open(filename, 'r').readlines())-1
    data = pd.DataFrame.from_dict(dic_go, orient='index')
    data.columns = ['enrichment']
    # read level information
    # termFile = 'G:\\project2\\NPM201507\\code\\1119dhac_NPM\\1218terms.txt'
    termFile = 'G:\\project2\\NPM201507\\code\\1119dhac_NPM\\1219terms.txt'
    # DHAC
    # termFile = 'G:\\project2\\NPM201507\\code\\1119dhac_NPM\\1208_DHAC_terms_sign_list_id.txt'
    # NPM
    # termFile = 'G:\\project2\\NPM201507\\code\\1119dhac_NPM\\1211_NPM10_terms_sign_id(2).txt'
    termData = pd.read_csv(termFile, sep='\t',index_col=0)
    # concat two data
    result = pd.concat([data, termData], axis=1)
    print 'group', '\t', 'nonzero', '\t', 'total', '\t', 'percentage', '\t', 'avg'
    for group in result.groupby('level'):
        df1 = group[1]
        zeroNum = df1[df1['enrichment'] == 0].shape[0]
        nz = df1['enrichment' ].mean()
        total = group[1]['enrichment'].shape[0]
        nonzero = total - zeroNum
        percentage = round(nonzero/float(total),3)
        print group[0], '\t', nonzero, '\t', total, '\t', percentage, '\t', round(nz,3)
        if group[0] == 14:
            print group[1]

    # sub_heatmap_numbers()
    # sub_heatmap_percentage()




