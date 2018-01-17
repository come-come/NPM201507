# -*- coding: utf-8 -*-
import pandas as pd
import numpy as np
import pylab as plt
from matplotlib_venn import venn3, venn3_circles,venn2
import networkx as nx
import os
from os.path import join

def read_files():

    data = pd.read_table('aracyc_pathways.20140902.txt', sep = '\t')
    a = data[data['Gene-id'] == 'unknown']
    print a.shape
    association = data[['Pathway-id', 'Gene-id']].drop_duplicates()


    print 'Before remove the redundant rows:', data.shape
    print 'After remove the redundant rows:', association.shape

    # print data.head(5)
    association.to_csv('association1228.txt', sep='\t', index='False')
    df = association[association['Gene-id'] == 'unknown']
    print df.shape
    # Gene id -- Gene name
    # unknown

def chose_name():
    # Pathway2Gene
    data = pd.DataFrame(columns=['Pathway-id', 'Gene-name'])
    fr = open('aracyc_pathways.20140902.txt', 'r')
    title = fr.readline()
    i = 0
    sums = 0
    for line in fr.readlines():

        line_str = line.strip().split('\t')

        # Gene-Id is 'AT'
        if len(line_str[6].strip().split('AT')) > 1:
            data.loc[i] = [line_str[0], line_str[6]]

        # Gene-Name is 'AT.3'
        elif len(line_str[7].strip().split('AT')) > 1:
            data.loc[i] = [line_str[0], line_str[7].split('.')[0]]

        # Gene-Name==Gene-Id=='unknown'
        elif line_str[7] == 'unknown' and line_str[6] == 'unknown':
            continue


        else:
        # G-15187	cMDH
            data.loc[i] = [line_str[0], line_str[7]]

        i = i + 1
    print i, sums
    print data.shape
    data.drop_duplicates().to_csv('TERM2PATHWAY.txt', sep='\t', header=None,index=None)
    print data.drop_duplicates().shape

def genes_not_in_pathway():
    dic_genes = pd.read_table('termN_Id.txt', header=None, index_col=0).to_dict('index')
    pathways = pd.read_table('TERM2PATHWAY.txt',header=None)
    not_in_pathway = pd.DataFrame(columns=['Pathway-id', 'Gene-name'])
    print pathways.head(5)
    print pathways.shape
    j = 0
    gene = set([])
    for i in range(0, pathways.shape[0]):
        if pathways.loc[i][1] not in dic_genes.keys():
            not_in_pathway.loc[j] = [pathways.loc[i][0], pathways.loc[i][1]]
            j += 1
            pathways = pathways.drop([i])
        else :
            gene.add(pathways.loc[i][1])
    # not_in_pathway.loc[0] = [pathways.loc[0][0], pathways.loc[0][1]]
    # dd = pathways.drop([0])
    print pathways.head(5)
    print 'in the 182 gene list:', pathways.shape
    print not_in_pathway.head(5)
    print 'not in the 182 gene list:', not_in_pathway.shape
    print len(gene)
def background_genes():
    data = pd.read_table("association_new.txt",header=None).drop_duplicates()
    print data.shape
    data.to_csv('association_new2.txt',index=None,header=None)
def split_background_gene():
    data = pd.read_table("gene_association.tair", header=None, sep='\t')
    print data.head(5)
    # for i in range(0,data.shape[0]):
    #
    #     data.loc[i,10] = data.loc[i][10].split('|')[0]

    # print data.loc[1][10],type(data[1][10]),data.loc[1][10].split('|')[0]
    # print data.loc[1][10]
    # data.loc[1,3] = 25
    #
    # print data.loc[1][3]
    print data.head(5)
    for group in data.groupby(8):
        print group[0],group[1].shape
        filename = group[0] + '.txt'
        group[1][[4,10]].to_csv(filename, sep = '\t', header=None, index=None)
        print group[1][[4,10]].drop_duplicates().shape
def split_by():
    # 细分GO的三个分支
    fr = open('F.txt', 'r')
    fw = open('F_new.txt', 'w')
    for line in fr.readlines():
        line_str = line.strip().split('\t')
        line_ss= line_str[1].strip().split('|')[0].split('"')
        if len(line_ss) > 1:
            fw.write(line_str[0] + '\t' + line_ss[1] + '\n')
        else:
            fw.write(line_str[0] + '\t' + line_ss[0] + '\n')
    fw.close()
    fr.close()
    data = pd.read_table('F_new.txt', sep='\t', header=None)
    print data.shape
    data.drop_duplicates().to_csv('F_new2.txt', sep='\t', header=None, index=None)
    print data.drop_duplicates().shape
def count_overlap():
    background = pd.read_table('background_arab.txt', header=None)
    # background = pd.read_table('termN_Id.txt', header=None)
    print 'The number of background genes:', background.shape
    data = set(np.array(background[0]).tolist())

    anno_Gene = pd.read_table('termN_Id.txt', header=None)
    print 'The number of annotation genes:', anno_Gene.shape
    anno_Genes = set(np.array(anno_Gene[0]).tolist())

    all_genes = pd.read_table('F_new2.txt', header=None)
    # all_genes = pd.read_table('TERM2PATHWAY.txt', header=None)

    print all_genes[1].drop_duplicates().shape
    MF_genes = set(np.array(all_genes[1].tolist()))
    jiao = MF_genes & data
    print jiao
    fw = open('background_jiao_MF.txt', 'w')
    for gene in jiao:
        fw.write(gene+ '\n')
    fw.close()
    print 'The number of BP genes:', len(MF_genes)
    print 'The number of intersect genes:', len(jiao)
    venn2([data, MF_genes], ('data', 'MF_Genes'))
    # plt.figure()
    # ax1 = plt.subplot(221)
    # plt.sca(ax1)
    # venn2([data, anno_Genes], ('BG_Genes','Anno_Genes'))
    # ax2 = plt.subplot(223)
    # plt.sca(ax2)
    # venn2([data, BP_genes], ('BG_Genes', 'BP_Genes'))
    # ax3 = plt.subplot(222)
    # plt.sca(ax3)
    # venn2([anno_Genes, BP_genes], ('Anno_Genes', 'BP_Genes'))
    # ax3 = plt.subplot(224)
    # plt.sca(ax3)
    # venn3([data, anno_Genes, BP_genes], ('BG_Genes', 'Anno_Genes', 'Pathway_Genes'))
    plt.show()

def BP_tree():
    # root GO:0008150
    # 统计前三层的term 返回字典 term-depth
    G = nx.DiGraph()
    filename = 'G:\project1\GO_tree\BP_edges.txt'
    depth = 'G:\\project1\\GO_tree\\BP_Depth_NumInheritGene_Id.txt'
    fr2 = open(depth, 'r')
    dic_depth = {}
    fr = open(filename, 'r')
    title = fr2.readline()
    for line in fr.readlines():
        parent, child = line.strip().split('\t')
        G.add_edge(parent, child)
    dic_depth_term = {}
    for node in G.nodes():
        dic_depth_term[node] = nx.shortest_path_length(G, 'GO:0008150',node)
    for line in fr2.readlines():
        line_str= line.strip().split('\t')
        dep = int(line_str[0])
        term = line_str[1]
        dic_depth[term] = dep
    print 'BP_term (has annotation gene)', len(dic_depth)
    print G.number_of_nodes(),G.number_of_edges()
    all_genes = pd.read_table('P_new2.txt', header=None)
    print all_genes[0].drop_duplicates().shape
    BP_genes = set(np.array(all_genes[0].tolist()))
    jiao = BP_genes & set(dic_depth_term.keys())
    print 'The number of BP terms in tair files:', len(BP_genes)
    print 'The number of BP terms in ontology:', len(dic_depth_term)
    print 'intersection:' , len(jiao)
    return dic_depth,dic_depth_term

def MF_tree():
    # root GO:0008150
    # 统计前三层的term 返回字典 term-depth
    G = nx.DiGraph()
    filename = 'G:\project1\GO_tree\MF_edges.txt'
    depth = 'G:\\project1\\GO_tree\\MF_Depth_NumInheritGene_Id.txt'
    fr2 = open(depth, 'r')
    dic_depth = {}
    fr = open(filename, 'r')
    title = fr2.readline()
    for line in fr.readlines():
        parent, child = line.strip().split('\t')
        G.add_edge(parent, child)
    dic_depth_term = {}
    for node in G.nodes():
        dic_depth_term[node] = nx.shortest_path_length(G, 'GO:0003674',node)
    for line in fr2.readlines():
        line_str= line.strip().split('\t')
        dep = int(line_str[0])
        term = line_str[1]
        dic_depth[term] = dep
    print 'MF_term (has annotation gene)', len(dic_depth)
    print G.number_of_nodes(),G.number_of_edges()
    all_genes = pd.read_table('F_new2.txt', header=None)
    print all_genes[0].drop_duplicates().shape
    BP_genes = set(np.array(all_genes[0].tolist()))
    jiao = BP_genes & set(dic_depth_term.keys())
    print 'The number of MF terms in tair files:', len(BP_genes)
    print 'The number of MF terms in ontology:', len(dic_depth_term)
    print 'intersection:' , len(jiao)
    return dic_depth,dic_depth_term

#2018.1.4
def output_file():
    # 统计每个term注释基因可以被富集的比率是多少
    # 注释的182个基因有11个不在BP的term里
    annoGeneNoBP = set(['AT1G20810', 'AT4G16155', 'AT4G19830', 'AT4G20820', 'AT1G14345', 'AT4G38380', 'AT1G51110', 'AT4G38100', 'AT5G35170', 'AT1G72640', 'AT5G23890'])
    # 注释的182个基因有8个不在背景基因（5109）里
    annoGeneNoBack = set(['AT3G47450', 'AT1G03090', 'AT5G03455', 'AT1G54790', 'AT2G04030', 'AT4G37925', 'AT2G43010', 'AT2G26800'])
    print 'annoGeneNoBP & annoGeneNoBack:', len(annoGeneNoBP & annoGeneNoBack)
    # 1129terms.txt  0
    data = pd.read_table("G:\\project2\\NPM201507\\code\\1119dhac_NPM\\1211_NPM10_terms_sign_list_id.txt")
    # background = pd.read_table('background_arab.txt', header=None)
    # # background = pd.read_table('termN_Id.txt', header=None)
    # print 'The number of backgrount genes:', background.shape
    # background_gene = set(np.array(background[0]).tolist())
    t = 0
    for i in range(0,data.shape[0]):
        annoGene = set(data.loc[i]['annotation_gene'].split(','))
        # print annoGene
        if len(annoGene & annoGeneNoBack)/float(len(annoGene)) > 0.5:
            t = t + 1
            print annoGene
            # print data.loc[i]
    print t


    print data.loc[0]['annotation_gene']

def analysis_enrichment(path, termFile, dic_depth_term):
    # 平均能富集到多少个GO BP/MF term
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
    termData = pd.read_csv(termFile, sep='\t',index_col=0)
    # concat two data
    result = pd.concat([data, termData], axis=1)
    dic_output = {}
    print 'group', '\t', 'nonzero', '\t', 'total', '\t', 'percentage', '\t', 'avg'
    color = {1: 'r', 2:'g', 3:'b'}
    # fig = plt.figure()
    # fig, axes = plt.subplots(1, 3, figsize=(6, 6))
    # fig.suptitle('Our Method')
    for group in result.groupby('level'):
        if group[0]<11:
            df1 = group[1]
            zeroNum = df1[df1['enrichment'] == 0].shape[0]
            nz = df1['enrichment'].mean()
            total = group[1]['enrichment'].shape[0]
            nonzero = total - zeroNum
            percentage = round(nonzero / float(total), 3)
            avg = round(nz, 3)
            print group[0], '\t', nonzero, '\t', total, '\t', percentage, '\t', avg
            dic_output[group[0]] = [int(nonzero), int(total), percentage, avg]
    #         title_name = 'level' + str(group[0])
    #         group_Data = group[1].reset_index()
    #         ax=group_Data.plot(x='index', y='enrichment', kind='scatter', color=color[group[0]], ax=axes[group[0]-1], title=title_name,
    #                         yticks=range(0,10))
    #         ax.set_ylabel('depth')
    #
    # print 'plt.show'
    #
    # plt.show()
    return dic_output

def avg_depth_of_enrich_restult(path, termFile, dic_depth_term, figure_i, figure_title):
    # 富集的深度
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
                # 富集到的最深的term是第几层
                if dep < sum:
                    continue
                else:
                    sum = dep
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

    fig.suptitle(figure_title)
    plt.subplot(130 + figure_i +1)
    for group in result.groupby('level'):
        dic_this_level = {}
        if group[0] < 4:
            df1 = group[1]
            num_level_node = df1.shape[0]
            group2 = group[1].groupby('enrichment')
            print 'group2'
            for group3 in group2:
                print group3[0], group3[1].shape[0]
                dic_this_level[group3[0]] = round(group3[1].shape[0] /float(num_level_node),3)

            print dic_this_level
            x = sorted(dic_this_level.keys())
            y = [dic_this_level[key] for key in x]
            plt.plot(x, y, marker='*', label='level'+ str(group[0]))
            # zeroNum = df1[df1['enrichment'] == 0].shape[0]
            # nz = df1['enrichment'].mean()
            # total = group[1]['enrichment'].shape[0]
            # nonzero = total - zeroNum
            # percentage = round(nonzero / float(total), 3)
            # avg = round(nz, 3)
            # print group[0], '\t', nonzero, '\t', total, '\t', percentage, '\t', avg
            # dic_output[group[0]] = [int(nonzero), int(total), percentage, avg]
            # title_name = 'level' + str(group[0])
            # group_Data = group[1].reset_index()
            # ax=group_Data.plot(x='index', y='enrichment', kind='scatter', color=color[group[0]], ax=axes[group[0]-1], title=title_name,
            #                 yticks=range(0,10))
            # ax.set_ylabel('depth')
    plt.legend(numpoints=1)
    plt.xticks(range(0,10))
    plt.yticks(np.arange(0, 1.1,0.1))
    plt.xlabel('depth')
    plt.ylabel('percentage')

    print 'plt.show'

    return dic_output


# 2018.1.11
def avg_num_of_BP_enrich_restult():
    # 平均每层每个node可以富集到多少个GO term

    avg_dataFrame_NPM = pd.DataFrame()
    avg_dataFrame_OurMethod = pd.DataFrame()
    print 'NPM'
    i = 0
    for key, value in dic_files2.items():
        dic_output = analysis_enrichment(key, value, dic_depth_term)
        label = 'NPM' + str(i + 1)
        d = pd.DataFrame.from_dict(dic_output, orient='index')
        #        plt.plot(list(d.index), d[3].tolist(), color='green',label=label)
        avg_dataFrame_NPM[i] = d[3]
        i = i + 1
    print avg_dataFrame_NPM.mean(1)
    plt.title('NPM')
    # plt.show()
    j = 0
    print 'Our method'
    for key, value in dic_files.items():
        dic_output = analysis_enrichment(key, value, dic_depth_term)
        d = pd.DataFrame.from_dict(dic_output, orient='index')
        label = 'OurMethod' + str(j + 1)
        #         plt.plot(list(d.index), d[3].tolist(),  color='red',label=label)
        avg_dataFrame_OurMethod[j] = d[3]
        j = j + 1
    print avg_dataFrame_OurMethod.mean(1)
    plt.title('Our method')

    print 'DHAC'
    plt.plot(list(avg_dataFrame_OurMethod.index), avg_dataFrame_OurMethod.mean(1).tolist(), color='red',
             label='OurMethod')
    plt.plot(list(avg_dataFrame_NPM.index), avg_dataFrame_NPM.mean(1).tolist(), color='green', label='NPM')

    for key, value in dic_files3.items():
        dic_output = analysis_enrichment(key, value, dic_depth_term)
        d = pd.DataFrame.from_dict(dic_output, orient='index')
        label = 'DHAC'
        plt.plot(list(d.index), d[3].tolist(), color='black', label=label)

    plt.title('BP Enrichment Results')
    plt.xticks(range(0, 12, 1))
    plt.yticks(range(0, 15, 1))
    plt.xlabel('level')
    plt.ylabel('average number of enricment results')
    plt.legend()
    plt.show()

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
    path = 'G:\\project2\\enrichment\\enrich20180103_1219my\\'
    termFile = 'G:\\project2\\NPM201507\\code\\1119dhac_NPM\\1219terms.txt'
    dic_files[path] = termFile 
    path2 = 'G:\\project2\\enrichment\\enrich20180103_1218my\\'
    termFile2 = 'G:\\project2\\NPM201507\\code\\1119dhac_NPM\\1218terms.txt'
    dic_files[path2] = termFile2
    path3 = 'G:\\project2\\enrichment\\enrich20180103_1129my\\'
    termFile3 = 'G:\\project2\\NPM201507\\code\\1119dhac_NPM\\1129terms.txt'
    dic_files[path3] = termFile3

    dic_files2 = {}
    path = 'G:\\project2\\enrichment\\enrich20180104_1211NPM\\'
    termFile = 'G:\\project2\\NPM201507\\code\\1119dhac_NPM\\1211_NPM10_terms_sign_id.txt'
    dic_files2[path] = termFile
    path = 'G:\\project2\\enrichment\\enrich20180104_1211NPM(2)\\'
    termFile = 'G:\\project2\\NPM201507\\code\\1119dhac_NPM\\1211_NPM10_terms_sign_id(2).txt'
    dic_files2[path] = termFile
    path = 'G:\\project2\\enrichment\\enrich20180104_1211NPM(3)\\'
    termFile = 'G:\\project2\\NPM201507\\code\\1119dhac_NPM\\1222_NPM_terms_sign_id(3).txt'
    dic_files2[path] = termFile

    dic_files3 = {}
    path = 'G:\\project2\\enrichment\\enrich20180104_1208DHAC\\'
    termFile = 'G:\\project2\\NPM201507\\code\\1119dhac_NPM\\1208_DHAC_terms_sign_list_id.txt'
    dic_files3[path] = termFile
    dic_depth, dic_depth_term = BP_tree()

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

    dic_depth, dic_depth_MFterm = BP_tree()

    fig = plt.figure()
    i = 0
    for key, value in dic_files.items():
        dic_output = avg_depth_of_enrich_restult(key, value, dic_depth_MFterm, i,'Our Method')
        i = i + 1
    plt.show()

    fig = plt.figure()
    i = 0
    for key, value in dic_files2.items():
        dic_output = avg_depth_of_enrich_restult(key, value, dic_depth_MFterm, i, 'NPM')
        i = i + 1
    plt.show()

    fig = plt.figure()
    i = 0
    for key, value in dic_files3.items():
        dic_output = avg_depth_of_enrich_restult(key, value, dic_depth_MFterm, i, 'DHAC')
        i = i + 1
    plt.show()


    # avg_num_of_BP_enrich_restult()

    #dic_depth, dic_depth_MFterm = MF_tree()

    '''
    print 'Our method'
    o_data = pd.DataFrame()
    n_data = pd.DataFrame()
    d_data = pd.DataFrame()
    i=0
    for key, value in dic_files.items():
        dic_output = analysis_enrichment(key, value, dic_depth_MFterm)
        o_data[i] = pd.DataFrame(dic_output).loc[3]
        i = i + 1
    i = 0
    for key, value in dic_files2.items():
        dic_output = analysis_enrichment(key, value, dic_depth_MFterm)
        n_data[i] = pd.DataFrame(dic_output).loc[3]
        i = i + 1
    i = 0
    for key, value in dic_files3.items():
        dic_output = analysis_enrichment(key, value, dic_depth_MFterm)
        d_data[i] = pd.DataFrame(dic_output).loc[3]
        i = i + 1
    print 'our',o_data.mean(axis=1)
    print 'npm',n_data.mean(axis=1)
    print 'dhac',d_data.mean(axis=1)
    fig = plt.figure()
    plt.plot(o_data.mean(axis=1), marker='*', ms=10, label='Our Method')
    plt.plot(n_data.mean(axis=1), marker='*', ms=10, label='NPM')
    plt.plot(d_data.mean(axis=1), marker='*', ms=10, label='DHAC')
    plt.legend(loc=0, numpoints=1)
    plt.title('BP Enrichment Results',fontsize=20)
    plt.xticks(range(0,12))
    plt.yticks(range(0,14))
    plt.xlabel('Level', fontsize=15)
    plt.ylabel('Average number of enrichment results', fontsize=15)
    plt.show()
    '''










