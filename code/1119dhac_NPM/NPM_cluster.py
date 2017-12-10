# -*- coding: utf-8 -*-
import pandas as pd
import os
from os.path import join

# c5 s10
def comparison_npm_file(desk):
    n = 0
    for root1, dirs1, files1 in os.walk(desk):
        for dirs2 in dirs1:
            root2 = os.path.join('%s%s' % (root1, dirs2))
            child = os.path.join('%s%s' % (root2, '\\1.txt'))
            print child
            data = pd.read_csv(child, sep = '\t')
            group = data.groupby('1')
            filename = str(n) + '.txt'
            fw = open(filename, 'a')
            for i in group:
                clusters = list(i[1]['name[flat]'])
                print clusters
                # for cluster in clusters:
                fw.write('\t'.join(gene for gene in clusters))
                fw.write('\n')
            fw.close()
            n = n + 1
if __name__ == "__main__":
    desk = 'G:\\project2\\NPM201507\\clusterResult0516\\'
    comparison_npm_file(desk)
