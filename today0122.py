# -*- coding: utf-8 -*-
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import association

if __name__ == "__main__":
    # read_files()
    # chose_name()
    # genes_not_in_pathway()
    # background_genes()
    #split_background_gene()
    # split_by()
    # count_overlap()
    # output_file()
    dic_files = {}
    path = 'G:\\project2\\enrichment\\enrich20180103_1219my_BP\\'
    termFile = 'G:\\project2\\NPM201507\\code\\1119dhac_NPM\\1219terms.txt'
    dic_files[path] = termFile
    path2 = 'G:\\project2\\enrichment\\enrich20180103_1218my_BP\\'
    termFile2 = 'G:\\project2\\NPM201507\\code\\1119dhac_NPM\\1218terms.txt'
    dic_files[path2] = termFile2
    path3 = 'G:\\project2\\enrichment\\enrich20180103_1129my_BP\\'
    termFile3 = 'G:\\project2\\NPM201507\\code\\1119dhac_NPM\\1129terms.txt'
    dic_files[path3] = termFile3

    dic_files2 = {}
    path = 'G:\\project2\\enrichment\\enrich20180104_1211NPM_BP\\'
    termFile = 'G:\\project2\\NPM201507\\code\\1119dhac_NPM\\1211_NPM10_terms_sign_id.txt'
    dic_files2[path] = termFile
    path = 'G:\\project2\\enrichment\\enrich20180104_1211NPM(2)_BP\\'
    termFile = 'G:\\project2\\NPM201507\\code\\1119dhac_NPM\\1211_NPM10_terms_sign_id(2).txt'
    dic_files2[path] = termFile
    path = 'G:\\project2\\enrichment\\enrich20180104_1211NPM(3)_BP\\'
    termFile = 'G:\\project2\\NPM201507\\code\\1119dhac_NPM\\1222_NPM_terms_sign_id(3).txt'
    dic_files2[path] = termFile

    dic_files3 = {}
    path = 'G:\\project2\\enrichment\\enrich20180104_1208DHAC_BP\\'
    termFile = 'G:\\project2\\NPM201507\\code\\1119dhac_NPM\\1208_DHAC_terms_sign_list_id.txt'
    dic_files3[path] = termFile
    dic_depth, dic_depth_term = association.BP_tree()

    dic_files4 = {}
    path = 'G:\\project2\\enrichment\\enrich20180103_1219my_MF\\'
    termFile = 'G:\\project2\\NPM201507\\code\\1119dhac_NPM\\1219terms.txt'
    dic_files4[path] = termFile
    path2 = 'G:\\project2\\enrichment\\enrich20180103_1218my_MF\\'
    termFile2 = 'G:\\project2\\NPM201507\\code\\1119dhac_NPM\\1218terms.txt'
    dic_files4[path2] = termFile2
    path3 = 'G:\\project2\\enrichment\\enrich20180103_1129my_MF\\'
    termFile3 = 'G:\\project2\\NPM201507\\code\\1119dhac_NPM\\1129terms.txt'
    dic_files4[path3] = termFile3

    dic_files5 = {}
    path = 'G:\\project2\\enrichment\\enrich20180104_1211NPM_MF\\'
    termFile = 'G:\\project2\\NPM201507\\code\\1119dhac_NPM\\1211_NPM10_terms_sign_id.txt'
    dic_files5[path] = termFile
    path = 'G:\\project2\\enrichment\\enrich20180104_1211NPM(2)_MF\\'
    termFile = 'G:\\project2\\NPM201507\\code\\1119dhac_NPM\\1211_NPM10_terms_sign_id(2).txt'
    dic_files5[path] = termFile
    path = 'G:\\project2\\enrichment\\enrich20180104_1211NPM(3)_MF\\'
    termFile = 'G:\\project2\\NPM201507\\code\\1119dhac_NPM\\1222_NPM_terms_sign_id(3).txt'
    dic_files5[path] = termFile

    dic_files6 = {}
    path = 'G:\\project2\\enrichment\\enrich20180104_1208DHAC_MF\\'
    termFile = 'G:\\project2\\NPM201507\\code\\1119dhac_NPM\\1208_DHAC_terms_sign_list_id.txt'
    dic_files6[path] = termFile

    '''
    # BP分支上 富集到最深的term代表node富集到的深度，分析前三个level的nodes的深度
    '''
    print '.......step 1 : generate depth dictionary......'
    dic_depth, dic_depth_MFterm = association.BP_tree()
    fig = plt.figure(figsize=(12,5))
    i = 0
    print '......step 2 : start drawing......'
    print '2.1 Our method'
    Our_DataFrame1 = pd.DataFrame()
    Our_DataFrame2 = pd.DataFrame()
    Our_DataFrame3 = pd.DataFrame()
    for key, value in dic_files.items():
        # dic_output,x, y = avg_depth_of_enrich_restult(key, value, dic_depth_MFterm, i,'Our Method')
        x, level1, level2, level3 = association.avg_depth_of_enrich_restult3(key, value, dic_depth_MFterm, i,'Our Method')
        Our_DataFrame1[i] = level1
        Our_DataFrame2[i] = level2
        Our_DataFrame3[i] = level3
        print x,level1,level2,level3
        i = i + 1
    print x
    # Our_DataFrame1['depth'] = x
    # Our_DataFrame2['depth'] = x
    # Our_DataFrame3['depth'] = x

    level1_avg = Our_DataFrame1.mean(axis=1)
    std1 = Our_DataFrame1.std(axis=1)

    level2_avg = Our_DataFrame2.mean(axis=1)    # level2 行均值
    std2 = Our_DataFrame2.std(axis=1)

    level3_avg = Our_DataFrame3.mean(axis=1)
    std3 = Our_DataFrame3.std(axis=1)

    # tmpArray_Avg = level1_avg.values.reshape((-1,1))
    # tmpArray_Avg = np.concatenate((tmpArray_Avg, level2_avg.values.reshape((-1,1)), level3_avg.values.reshape((-1,1))), axis=1)
    # Our_Avg_DataFrame = pd.DataFrame(tmpArray_Avg,index=x,columns=('level1','level2','level3'))
    # print Our_Avg_DataFrame

    Our_Avg_DataFrame = pd.DataFrame()
    Our_Avg_DataFrame['level0'] = level1_avg
    Our_Avg_DataFrame['level1'] = level2_avg
    Our_Avg_DataFrame['level2'] = level3_avg
    Our_Avg_DataFrame.index = x
    print Our_Avg_DataFrame

    Our_Std_DataFrame = pd.DataFrame()
    Our_Std_DataFrame['level0'] = std1
    Our_Std_DataFrame['level1'] = std2
    Our_Std_DataFrame['level2'] = std3
    Our_Std_DataFrame.index = x
    print Our_Std_DataFrame
    # tmpArray_Std = std1.values.reshape((-1, 1))
    # tmpArray_Std = np.concatenate(
    #     (tmpArray_Std, std2.values.reshape((-1, 1)), std3.values.reshape((-1, 1))), axis=1)
    # Our_Std_DataFrame = pd.DataFrame(tmpArray_Std, index=x, columns=('level1', 'level2', 'level3'))

    fig, ax = plt.subplots()
    Our_Avg_DataFrame.plot.bar(yerr=Our_Std_DataFrame, ax=ax,title='Our Method')
    # Our_Avg_DataFrame.plot(kind='bar')
    #plt.xticks(x, x, rotation=0)
    plt.xlabel('depth')
    plt.ylabel('percentage')
    plt.yticks(np.arange(0, 0.6,0.1),fontsize=10)
    plt.show()
    print '2.2 NPM'
    i = 0
    for key, value in dic_files2.items():
        # dic_output, x, y = avg_depth_of_enrich_restult(key, value, dic_depth_MFterm, i, 'NPM')
        x, level1, level2, level3 = association.avg_depth_of_enrich_restult3(key, value, dic_depth_MFterm, i, 'NPM')
        i = i + 1


    print '2.3 DHAC'

    i = 0
    for key, value in dic_files3.items():
        # dic_output, x, y= avg_depth_of_enrich_restult(key, value, dic_depth_MFterm, i, 'DHAC')
        x, level1, level2, level3 = association.avg_depth_of_enrich_restult3(key, value, dic_depth_MFterm, i, 'DHAC')
        i = i + 1


    print '-----------MF---------'
    dic_depth, dic_depth_MFterm = association.MF_tree()

    i = 0
    for key, value in dic_files4.items():
        x, level1, level2, level3 = association.avg_depth_of_enrich_restult3(key, value, dic_depth_MFterm, i, 'Our Method')
        i = i + 1

    i = 0
    for key, value in dic_files5.items():
        x, level1, level2, level3 = association.avg_depth_of_enrich_restult3(key, value, dic_depth_MFterm, i, 'NPM')
        i = i + 1

    i = 0
    for key, value in dic_files6.items():
        x, level1, level2, level3 = association.avg_depth_of_enrich_restult3(key, value, dic_depth_MFterm, i, 'DHAC')
        i = i + 1











