# -*- coding: utf-8 -*-
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import association
import os
from os.path import join

def enrich_percentage(path, termFile, dic_depth_term):
    # 每一层 有多少node 可以富集到 GO term
    dic_go = {}
    print 'len(dic_depth_term):', len(dic_depth_term)
    for root, dirs, files in os.walk(path):
        for f in files:
            filename = join(root, f)
            sum = 0
            fr = open(filename, 'r')
            title = fr.readline()
            for line in fr.readlines():
                line_str = line.strip().split(',')
                dep = dic_depth_term[line_str[0][1:-1].strip()]
                # if dep < sum:
                if dep < 3:
                    continue
                else:
                    # sum = dep
                    sum = sum + 1
            dic_go[int(f.strip().split('.')[0])] = sum

    print len(dic_go)

    data = pd.DataFrame.from_dict(dic_go, orient='index')
    data.columns = ['enrichment']
    termData = pd.read_csv(termFile, sep='\t', index_col=0)
    # concat two data
    result = pd.concat([data, termData], axis=1)
    dic_output = {}
    print 'group', '\t', 'nonzero', '\t', 'total', '\t', 'percentage', '\t', 'avg'
    color = {1: 'r', 2: 'g', 3: 'b'}
    # fig = plt.figure()
    # fig, axes = plt.subplots(1, 3, figsize=(6, 6))
    # fig.suptitle('Our Method')
    for group in result.groupby('level'):
        if group[0] < 8:
            df1 = group[1]
            zeroNum = df1[df1['enrichment'] == 0].shape[0]
            nz = df1['enrichment'].mean()
            total = group[1]['enrichment'].shape[0]
            nonzero = total - zeroNum
            percentage = round(nonzero / float(total), 3)
            avg = round(nz, 3)
            print group[0], '\t', nonzero, '\t', total, '\t', percentage, '\t', avg
            dic_output[group[0]] = [int(nonzero), int(total), percentage, avg]
    # title_name = 'level' + str(group[0])
    #         group_Data = group[1].reset_index()
    #         ax=group_Data.plot(x='index', y='enrichment', kind='scatter', color=color[group[0]], ax=axes[group[0]-1], title=title_name,
    #                         yticks=range(0,10))
    #         ax.set_ylabel('depth')
    #
    # print 'plt.show'
    #
    # plt.show()

    return dic_output

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
    # The Phenotypic Ontology Construction Based on clique
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
    # 图片 1行3列
    fig, (ax1, ax2, ax3) = plt.subplots(1, 3, sharex=True)
    # 第一个子图
    Our_Avg_DataFrame.plot.bar(yerr=Our_Std_DataFrame, ax=ax1,title='Our Method',yticks=np.arange(-0.1, 0.6, 0.1))
    ax1.set_xlabel('depth')
    ax1.set_ylabel('percentage')



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
    Our_DataFrame3[2] = [0.26, 0, 0.065, 0.098, 0.211, 0.114, 0.203]
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
    ax2.set_xlabel('depth')
    ax2.set_ylabel('percentage')


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
    ax3.set_xlabel('depth')
    ax3.set_ylabel('percentage')

    plt.show()


    print '-----------MF---------'
    dic_depth, dic_depth_MFterm = association.MF_tree()

    print '......step 2 : start drawing......MF'
    print '2.1 Our method'
    # The Phenotypic Ontology Construction Based on clique
    i = 0
    Our_DataFrame1 = pd.DataFrame()
    Our_DataFrame2 = pd.DataFrame()
    Our_DataFrame3 = pd.DataFrame()
    x = [0, 1, 2, 3, 4, 5, 6, 7]
    Our_DataFrame1[0] = [0.222, 0.003, 0.06, 0.02, 0.052, 0.159, 0.247, 0.238]
    Our_DataFrame2[0] = [0.269, 0.001, 0.059, 0.019, 0.059, 0.16, 0.228, 0.204]
    Our_DataFrame3[0] = [0.255, 0, 0.083, 0.032, 0.078, 0.156, 0.209, 0.188]
    Our_DataFrame1[1] = [0.193, 0, 0.074, 0.015, 0.041, 0.2, 0.193, 0.284]
    Our_DataFrame2[1] = [0.259, 0, 0.062, 0.017, 0.044, 0.155, 0.214, 0.248]
    Our_DataFrame3[1] = [0.243, 0.001, 0.075, 0.02, 0.057, 0.149, 0.218, 0.236]
    Our_DataFrame1[2] = [0.226, 0.001, 0.036, 0.021, 0.068, 0.179, 0.224, 0.244]
    Our_DataFrame2[2] = [0.247, 0, 0.054, 0.022, 0.087, 0.195, 0.17, 0.225]
    Our_DataFrame3[2] = [0.243, 0, 0.072, 0.027, 0.065, 0.189, 0.199, 0.206]


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
    # 图片 1行3列
    fig, (ax1, ax2, ax3) = plt.subplots(1, 3, sharex=True)
    # 第一个子图
    Our_Avg_DataFrame.plot.bar(yerr=Our_Std_DataFrame, ax=ax1, title='Our Method', yticks=np.arange(-0.1, 0.6, 0.1))
    ax1.set_xlabel('depth')
    ax1.set_ylabel('percentage')

    print '2.2 NPM'
    print '-----------'

    Our_DataFrame1 = pd.DataFrame()
    Our_DataFrame2 = pd.DataFrame()
    Our_DataFrame3 = pd.DataFrame()
    x = [0, 1, 2, 3, 4, 5, 6, 7]
    Our_DataFrame1[0] = [0.396, 0.014, 0.029, 0.04, 0.043, 0.14, 0.101, 0.187]
    Our_DataFrame2[0] = [0.299, 0.005, 0.051, 0.091, 0.036, 0.223, 0.071, 0.188]
    Our_DataFrame3[0] = [0.391, 0, 0.073, 0.045, 0.018, 0.164, 0.064, 0.245]
    Our_DataFrame1[1] = [0.386, 0, 0.03, 0.049, 0.038, 0.148, 0.125, 0.17]
    Our_DataFrame2[1] = [0.404, 0, 0.036, 0.067, 0.005, 0.155, 0.171, 0.15]
    Our_DataFrame3[1] = [0.423, 0, 0.041, 0.024, 0.008, 0.195, 0.163, 0.146]
    Our_DataFrame1[2] = [0.42, 0.007, 0.026, 0.04, 0.022, 0.131, 0.113, 0.179]
    Our_DataFrame2[2] = [0.394, 0, 0.048, 0.064, 0.011, 0.181, 0.117, 0.176]
    Our_DataFrame3[2] = [0.393, 0, 0.051, 0.051, 0.009, 0.197, 0.111, 0.188]

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
    ax2.set_xlabel('depth')
    ax2.set_ylabel('percentage')

    print '2.3 DHAC'
    i = 0
    Our_DataFrame1 = pd.DataFrame()
    Our_DataFrame2 = pd.DataFrame()
    Our_DataFrame3 = pd.DataFrame()
    x = [0, 1, 2, 3, 4, 5, 6, 7]
    Our_DataFrame1[0] = [0.379, 0.005, 0.054, 0.059, 0.064, 0.202, 0.079, 0.158]
    Our_DataFrame2[0] = [0.396, 0, 0.042, 0.01, 0.052, 0.219, 0.104, 0.177]
    Our_DataFrame3[0] = [0.39, 0, 0.024, 0, 0.049, 0.22, 0.122, 0.195]

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
    Our_Avg_DataFrame.plot.bar(yerr=Our_Std_DataFrame, ax=ax3, title='DHAC', yticks=np.arange(-0.1, 0.6, 0.1))
    ax3.set_xlabel('depth')
    ax3.set_ylabel('percentage')

    plt.show()

    print '--------percentage----------MF'

    dic_depth, dic_depth_MFterm = association.MF_tree()

    o_data = pd.DataFrame()
    n_data = pd.DataFrame()
    d_data = pd.DataFrame()
    i=0
    for key, value in dic_files4.items():
        dic_output = enrich_percentage(key, value, dic_depth_MFterm)
        o_data[i] = pd.DataFrame(dic_output).loc[2]
        i = i + 1
    i = 0
    for key, value in dic_files5.items():
        dic_output = enrich_percentage(key, value, dic_depth_MFterm)
        n_data[i] = pd.DataFrame(dic_output).loc[2]
        i = i + 1
    i = 0
    for key, value in dic_files6.items():
        dic_output = enrich_percentage(key, value, dic_depth_MFterm)
        d_data[i] = pd.DataFrame(dic_output).loc[2]
        i = i + 1
    print 'our',o_data.mean(axis=1)
    print 'npm',n_data.mean(axis=1)
    print 'dhac',d_data.mean(axis=1)
    fig = plt.figure()
    plt.plot(o_data.mean(axis=1), marker='*',  ms=10, label='Our Method')
    plt.plot(n_data.mean(axis=1), marker='*', ms=10, label='NPM')
    plt.plot(d_data.mean(axis=1), marker='*', ms=10, label='DHAC')
    plt.legend(loc=0, numpoints=1)
    plt.title('MF Enrichment Results',fontsize=20)
    plt.xticks(range(0,9))
    plt.yticks(np.arange(0, 1.1, 0.1))
    # plt.yticks(range(0,4))
    plt.xlabel('Level', fontsize=15)
    plt.ylabel('Percentage', fontsize=15)
    plt.show()

    print '--------percentage----------BP'

    dic_depth, dic_depth_MFterm = association.BP_tree()

    o_data = pd.DataFrame()
    n_data = pd.DataFrame()
    d_data = pd.DataFrame()
    i = 0
    for key, value in dic_files.items():
        dic_output = enrich_percentage(key, value, dic_depth_MFterm)
        o_data[i] = pd.DataFrame(dic_output).loc[2]
        i = i + 1
    i = 0
    for key, value in dic_files2.items():
        dic_output = enrich_percentage(key, value, dic_depth_MFterm)
        n_data[i] = pd.DataFrame(dic_output).loc[2]
        i = i + 1
    i = 0
    for key, value in dic_files3.items():
        dic_output = enrich_percentage(key, value, dic_depth_MFterm)
        d_data[i] = pd.DataFrame(dic_output).loc[2]
        i = i + 1
    print 'our', o_data.mean(axis=1)
    print 'npm', n_data.mean(axis=1)
    print 'dhac', d_data.mean(axis=1)
    fig = plt.figure()
    plt.plot(o_data.mean(axis=1), marker='*', ms=10, label='Our Method')
    plt.plot(n_data.mean(axis=1), marker='*', ms=10, label='NPM')
    plt.plot(d_data.mean(axis=1), marker='*', ms=10, label='DHAC')
    plt.legend(loc=0, numpoints=1)
    plt.title('BP Enrichment Results', fontsize=20)
    plt.xticks(range(0, 9))
    plt.yticks(np.arange(0, 1.1, 0.1))
    # plt.yticks(range(0,4))
    plt.xlabel('Level', fontsize=15)
    plt.ylabel('Percentage', fontsize=15)
    plt.show()