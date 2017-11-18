import pandas as pd
import numpy as np
import os
from os.path import join

dic_result = {}
def read_dhac_result(dest):
    for root, dirs, files in os.walk(dest):
        for OneFileName in files:
            if OneFileName.find('.group') == -1:
                continue
            OneFullFileName = join(root, OneFileName)
            # print OneFullFileName
            # print OneFileName.strip().split('edgeInfo')[1].split('.group')[0]
            dic_result [int(OneFileName.strip().split('edgeInfo')[1].split('.group')[0])] = []
            for line in open(OneFullFileName, 'r'):
                b = line.strip().split('\t')
                # print b  # ['821137', '829490']
                dic_result[int(OneFileName.strip().split('edgeInfo')[1].split('.group')[0])].append(b)
if __name__ == "__main__" :
    geneList = []
    windowList = []
    dict = {}
    read_dhac_result('G:\\project2\\NPM201507\\code\\edge\\20171113dhac')
    # print dic_result[1]
    clique_result = open('G:\\project2\\NPM201507\\code\\1117terms_sign_list_id.txt','r')
    title = clique_result.readline()
    '''
    first_row = clique_result.readline() 
    s = first_row.strip().split('\t')
    genes = s[2].strip().split(',')
    windows = [i for i in range(int(s[3])-49, int(s[4])-58 + 1)]
    print 'MY METHOD:', genes, windows
    for window in range(0, 50):
        for value in dic_result[window]:
            if set(genes).issubset(set(value)):
                print window
                break
    '''
    dic_dhac_term = {}
    dic_my_term= {}
    fw = open('super_than_dhac.txt', 'w')
    fw2 = open('sub_than_dhac.txt', 'w')
    fw3 = open('sssub_than_dhac.txt', 'w')
    fw.write('term_id' + '\t' + 'annotation_genes' + '\t' + 'my_method' + '\t' + 'dhac' + '\n')
    fw2.write('term_id' + '\t' + 'annotation_genes' + '\t' + 'my_method' + '\t' + 'dhac' + '\n')

    for line in clique_result:

        s = line.strip().split('\t')
        genes = s[2].strip().split(',')
        windows = [i for i in range(int(s[3]) - 49, int(s[4]) - 58 + 1)]

        term = int(s[0])
        dic_my_term[term] = windows
        print type(term)
        dic_dhac_term[term] = []
        for window in range(0, 50):
            for value in dic_result[window]:
                if set(genes).issubset(set(value)):
                    dic_dhac_term[term].append(window)
                    break
        print dic_my_term[term], dic_dhac_term[term]


        if len(dic_dhac_term[term])>0 and set(dic_my_term[term]).issuperset(set(dic_dhac_term[term])):
            fw.write(str(term) + '\t' + str(genes) + '\t' + str(dic_my_term[term]) + '\t' + str(dic_dhac_term[term]) + '\n')
        elif len(dic_dhac_term[term])>0 and set(dic_my_term[term]).issubset(set(dic_dhac_term[term])):
            fw2.write(str(term) + '\t' + str(genes) + '\t' + str(dic_my_term[term]) + '\t' + str(dic_dhac_term[term]) + '\n')

        if len(dic_dhac_term[term])>0 and set(dic_my_term[term]).issubset(set(dic_dhac_term[term])):
            fw3.write(str(term) + '\t' + str(genes) + '\t' + str(dic_my_term[term]) + '\t' + str(dic_dhac_term[term]) + '\n')
    fw.close()
    fw2.close()







