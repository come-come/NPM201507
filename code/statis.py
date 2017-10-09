import pandas as pd

import networkx as nx

dic1 = {}
dic2 = {}
dic3 = {}
dic = {}
dic11 = {}
dic22 = {}
dic33 = {}

def readf():
    data = pd.read_table('gene_association.tair', header=None)
    fr = open('geneName.txt', 'r')
    for line in fr:
        dic[line.strip()] = []
        dic11[line.strip()] = []
        dic22[line.strip()] = []
        dic33[line.strip()] = []
        dic1[line.strip()] = [] # Molecular Function
        dic2[line.strip()] = [] # Biological Process
        dic3[line.strip()] = [] # Cellular Component

    fr2 = open('gene_association.tair', 'r')
    for line in fr2:
        line_str = line.strip().split('\t')
        if dic.has_key(line_str[9].strip()):
            print line_str[8]
            if line_str[8].strip() == 'C':
                dic3[line_str[9].strip()].append(line_str[4].strip())
            if line_str[8].strip() == 'P':
                dic2[line_str[9].strip()].append(line_str[4].strip())
            if line_str[8].strip() == 'F':
                dic1[line_str[9].strip()].append(line_str[4].strip())
        #print line_str[4].strip(), line_str[8].strip(), line_str[9].strip()
    # fw1 = open('F.txt', 'w')
    # fw2 = open('P.txt', 'w')
    # fw3 = open('C.txt', 'w')
    # for key in dic1:
    #     fw1.write(key + '\t' + str(dic1[key]) + '\n')
    #     fw2.write(key + '\t' + str(dic2[key]) + '\n')
    #     fw3.write(key + '\t' + str(dic3[key]) + '\n')
    # fw1.close()
    # fw2.close()
    # fw3.close()
    print 'ok'
    return dic1,dic2,dic3,dic11,dic22,dic33

def go_tree(dic1,dic2,dic3,dic11,dic22,dic33):
    fr = open('go_go_one.txt', 'r')
    G = nx.DiGraph()
    for line in fr :
        q = line.strip().split('\t')
        if len(q) > 1:
            G.add_edge(q[0],q[1])
        else:
            G.add_node(q[0])
    print G.number_of_nodes(),G.number_of_edges()
    length = nx.all_pairs_shortest_path_length(G)
    for key,value in dic1.items():  # Molecular Function
        if len(value) > 0:
            for i in value:
                try:
                    dic11[key].append(length['GO:0003674'][i])
                except:
                    continue
    for key,value in dic2.items():  # Biological Process
        if len(value) > 0:
            for i in value:
                try:
                    dic22[key].append(length['GO:0008150'][i])
                except:
                    continue
    for key,value in dic3.items():  # Cellular Component
        if len(value) > 0:
            for i in value:
                try:
                    dic33[key].append(length['GO:0005575'][i])
                except:
                    continue
    print  length['GO:0003674']['GO:0003676']
    print  length['GO:0003674']['GO:0016787']
    fw1 = open('MF.txt', 'w')
    fw2 = open('BP.txt', 'w')
    fw3 = open('CC.txt', 'w')
    for key in dic11:
        fw1.write(key + '\t' + str(dic11[key]) + '\n')
        fw2.write(key + '\t' + str(dic22[key]) + '\n')
        fw3.write(key + '\t' + str(dic33[key]) + '\n')
    fw1.close()
    fw2.close()
    fw3.close()

if __name__ == '__main__':
    dic1, dic2, dic3, dic11, dic22, dic33 = readf()
    go_tree(dic1, dic2, dic3, dic11, dic22, dic33)
