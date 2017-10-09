# -*- coding: utf-8 -*-
"""
Created on Wed Mar 08 14:47:06 2017

@author: hww
"""
import networkx as nx

from sklearn.metrics import roc_curve
import numpy as np



#this function is used to find a depth for a unique term in the HPO tree.return a dictionary.
def get_term_depth(G2) :
    dic_term = {}
    set1 = set()
    '''set2 = set()
    set3 = set()'''
    length = nx.all_pairs_shortest_path_length(G2)
    err1=0
    for i in G2.nodes():
        try :
            dic_term[i] =length['HP:0000001'][i]
            set1.add(i)
        except :
            err1= err1+1
    
    return dic_term



#this function is used to 
def get_go_graph() :
    G2=nx.DiGraph()
    #a = []
    fr = open(r'J:\data\huiweiwei\gene similarity_from\hpo_hpo.txt',"r")
    #fw = open ('edges_error.txt','w')
    for line in fr :
        line_str = line.strip().split("\t")
        
        if len(line_str)==1 :
            G2.add_node(line_str[0])
        else :
            G2.add_node(line_str[0])
            G2.add_node(line_str[1])
            G2.add_edge(line_str[0],line_str[1])

    return G2
    
    
    
    
'''def get_gene_graph(filename) :
    G1=nx.Graph()
    fr = open(filename,"r")
    #fw = open ('edges_error.txt','w')
    for line in fr :
        line_str = line.strip().split("\t")
        G1.add_edge(line_str[0],line_str[1])
    sub_gragh = {}
    for i in nx.connected_component_subgraphs(G1):
        sub_gragh[n] = i
        n=n+1
    #print 'gene_graph子网数',len(sub_gragh)
    #print [len(c) for c in sorted(nx.connected_component_subgraphs(G1), key=len, reverse=True)]

    return G1'''
if not __name__ != '__main__':
    g=get_go_graph()
    d=get_term_depth(g)# d stores the information of terms' depth
    #print d['HP:0010628']
    #print max(d.values())
                    
    with open(r"J:\data\huiweiwei\gene similarity_from\selected1_term_depth.txt",'a') as ff:
        with open(r"J:\data\huiweiwei\gene similarity_from\selected1.txt",'r')as f:
            for line in f:
                line_data=line.strip().split('\t')
                #a.append(line_data[0])
                first_depth=d["HP:"+line_data[0].zfill(7)]
                second_depth=d["HP:"+line_data[1].zfill(7)]
                result="%s\t%s\t%s\t%s\t%s\n" % (line_data[0].zfill(7),line_data[1].zfill(7),line_data[2],first_depth,second_depth)
                ff.write(result)
                
                
    
    
    
    
    
    
    
    
    
    
    
    
    
    

                
                
            
