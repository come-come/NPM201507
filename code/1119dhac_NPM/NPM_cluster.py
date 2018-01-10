# -*- coding: utf-8 -*-
import pandas as pd
import os
from os.path import join

# c5 s10
def comparison_npm_file(desk):
    print desk

    for root1, dirs1, files1 in os.walk(desk):
        for dirs2 in dirs1:
            root2 = os.path.join('%s%s' % (root1, dirs2))
            child = os.path.join('%s%s' % (root2, '\\0.txt'))
            # print child
            print dirs2
            try:
                data = pd.read_csv(child, sep = '\t')
            except:
                print child
            group = data.groupby('1')
            filename = str(int(dirs2) - 49) + '.txt'
            print filename
            try:
                fw = open(filename, 'a')

            except:
                print filename
            for i in group:
                clusters = list(i[1]['Gene'])
                # print clusters
                # for cluster in clusters:
                fw.write('\t'.join(gene for gene in clusters))
                fw.write('\n')
            fw.close()

if __name__ == "__main__":

    # desk = 'G:\\project2\\NPM201507\\clusterResultAverageDHACc10_step10_ljy(2)\\'
    desk = 'G:\\project2\\NPM201507\\clusterResultAverageDHACc10_step10_ljy_1222\\'
    # remove writen file into G:\project2\NPM201507\code\edge\NPM_cluster10_(1)

    comparison_npm_file(desk)
