import pandas as pd
import numpy as np
import networkx as nx

if __name__ == '__main__':


    s = 0
    term = 183
    filename = 'result_c5_s10_20171217weight.txt'
    # filename = 'result_c5_s10_v2_weight.txt'
    data = pd.read_csv(filename, index_col=0, sep='\t')
    print data.head(5), type(data.loc['AT4G33520_AT1G67840']['49'])
    weight_value = 0.9
    windowGraph = {}
    for window in range(0, 54):
        windowGraph[window] = nx.Graph()
        df = data[data[data.columns[window]] >= (weight_value + 0.00001)]
        # print window, weight_value, df.shape
        for edge in range(0, df.shape[0]):
            node_1, node_2 = df.index[edge].split('_')
            windowGraph[window].add_edge(node_1, node_2)
        cliques = [i for i in nx.find_cliques(windowGraph[window]) if len(i)>5]
        print len(cliques)
        print window, windowGraph[window].number_of_nodes(), windowGraph[window].number_of_edges()