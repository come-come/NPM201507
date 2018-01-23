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


    print '......step 2 : start drawing......'
    print '2.1 Our method'
    i = 0
    Our_DataFrame1 = pd.DataFrame()
    Our_DataFrame2 = pd.DataFrame()
    Our_DataFrame3 = pd.DataFrame()
    for key, value in dic_files.items():
        x, level1, level2, level3 = association.avg_depth_of_enrich_restult3(key, value, dic_depth_MFterm, i,'Our Method')
        Our_DataFrame1[i] = level1
        Our_DataFrame2[i] = level2
        Our_DataFrame3[i] = level3
        i = i + 1
    print x

    level1_avg = Our_DataFrame1.mean(axis=1)
    std1 = Our_DataFrame1.std(axis=1)

    level2_avg = Our_DataFrame2.mean(axis=1)
    std2 = Our_DataFrame2.std(axis=1)

    level3_avg = Our_DataFrame3.mean(axis=1)
    std3 = Our_DataFrame3.std(axis=1)

    Our_Avg_DataFrame = pd.DataFrame()
    Our_Avg_DataFrame['level1'] = level1_avg
    Our_Avg_DataFrame['level2'] = level2_avg
    Our_Avg_DataFrame['level3'] = level3_avg
    Our_Avg_DataFrame.index = x
    print Our_Avg_DataFrame

    Our_Std_DataFrame = pd.DataFrame()
    Our_Std_DataFrame['level1'] = std1
    Our_Std_DataFrame['level2'] = std2
    Our_Std_DataFrame['level3'] = std3
    Our_Std_DataFrame.index = x
    print Our_Std_DataFrame

    fig, (ax1, ax2, ax3) = plt.subplots(1, 3, sharex=True)
    Our_Avg_DataFrame.plot.bar(yerr=Our_Std_DataFrame, ax=ax1,title='Our Method',yticks=np.arange(-0.1, 0.6, 0.1))
    plt.xlabel('depth')
    plt.ylabel('percentage')
    plt.yticks(np.arange(0, 0.6,0.1),fontsize=10)


    print '2.2 NPM'
    print '-----------'
    i = 0
    Our_DataFrame1 = pd.DataFrame()
    Our_DataFrame2 = pd.DataFrame()
    Our_DataFrame3 = pd.DataFrame()
    x = [0, 2, 3, 4, 5, 6, 7]
    Our_DataFrame1[0] = [0.212, 0.011, 0.076, 0.061, 0.252, 0.129, 0.176]
    Our_DataFrame2[0] = [0.223, 0.01, 0.076, 0.056, 0.213, 0.127, 0.178]
    Our_DataFrame3[0] = [0.318, 0.009, 0.036, 0.073, 0.136, 0.164, 0.164]
    Our_DataFrame1[1] = [0.255, 0.004, 0.066, 0.099, 0.252, 0.084, 0.186]
    Our_DataFrame2[1] = [0.234, 0, 0.069, 0.096, 0.213, 0.085, 0.213]
    Our_DataFrame3[1] = [0.214, 0, 0.068, 0.137, 0.179, 0.068, 0.231]
    Our_DataFrame1[2] = [0.231, 0.008, 0.049, 0.057, 0.231, 0.148, 0.212]
    Our_DataFrame2[2] = [0.269, 0.01, 0.021, 0.047, 0.233, 0.104, 0.254]
    Our_DataFrame3[3] = [0.26, 0, 0.065, 0.098, 0.211, 0.114, 0.203]
    print Our_DataFrame1
    print Our_DataFrame2
    print Our_DataFrame3
    level1_avg = Our_DataFrame1.mean(axis=1)
    std1 = Our_DataFrame1.std(axis=1)

    level2_avg = Our_DataFrame2.mean(axis=1)
    std2 = Our_DataFrame2.std(axis=1)

    level3_avg = Our_DataFrame3.mean(axis=1)
    std3 = Our_DataFrame3.std(axis=1)

    Our_Avg_DataFrame = pd.DataFrame()
    Our_Avg_DataFrame['level1'] = level1_avg
    Our_Avg_DataFrame['level2'] = level2_avg
    Our_Avg_DataFrame['level3'] = level3_avg
    Our_Avg_DataFrame.index = x
    print Our_Avg_DataFrame

    Our_Std_DataFrame = pd.DataFrame()
    Our_Std_DataFrame['level1'] = std1
    Our_Std_DataFrame['level2'] = std2
    Our_Std_DataFrame['level3'] = std3
    Our_Std_DataFrame.index = x
    print Our_Std_DataFrame
    Our_Avg_DataFrame.plot.bar(yerr=Our_Std_DataFrame, ax=ax2, title='NPM', yticks=np.arange(-0.1, 0.6, 0.1))
    plt.xlabel('depth')
    plt.ylabel('percentage')

    print '2.3 DHAC'
    i = 0
    Our_DataFrame1 = pd.DataFrame()
    Our_DataFrame2 = pd.DataFrame()
    Our_DataFrame3 = pd.DataFrame()
    x = [0, 2, 3, 4, 5, 6, 7]
    Our_DataFrame1[0] = [0.256, 0.01, 0.044, 0.128, 0.207, 0.138, 0.172]
    Our_DataFrame2[0] = [0.25, 0, 0.083, 0.135, 0.208, 0.135, 0.125]
    Our_DataFrame3[0] = [0.293, 0, 0.049, 0.122, 0.146, 0.244, 0.122]
    level1_avg = Our_DataFrame1.mean(axis=1)
    std1 = Our_DataFrame1.std(axis=1)

    level2_avg = Our_DataFrame2.mean(axis=1)
    std2 = Our_DataFrame2.std(axis=1)

    level3_avg = Our_DataFrame3.mean(axis=1)
    std3 = Our_DataFrame3.std(axis=1)

    Our_Avg_DataFrame = pd.DataFrame()
    Our_Avg_DataFrame['level1'] = level1_avg
    Our_Avg_DataFrame['level2'] = level2_avg
    Our_Avg_DataFrame['level3'] = level3_avg
    Our_Avg_DataFrame.index = x
    print Our_Avg_DataFrame

    Our_Std_DataFrame = pd.DataFrame()
    Our_Std_DataFrame['level1'] = std1
    Our_Std_DataFrame['level2'] = std2
    Our_Std_DataFrame['level3'] = std3
    Our_Std_DataFrame.index = x
    print Our_Std_DataFrame
    Our_Avg_DataFrame.plot.bar(yerr=Our_Std_DataFrame,ax=ax3, title='DHAC',yticks=np.arange(-0.1, 0.6, 0.1))
    plt.xlabel('depth')
    plt.ylabel('percentage')
    plt.yticks(np.arange(0, 0.6, 0.1), fontsize=10)
    plt.show()


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











