# -*- coding: utf-8 -*-

import pandas as pd
import numpy as np
import glob


def term():
    filename = 'result_c5_s10_v2_weight.txt'
    data = pd.read_csv(filename, index_col=0, sep='\t')
    weight_value = np.arange(1, -0.05, -0.05)
    # 注意此处，生成的 weight_value 比正常值少0.0000001，所以上下线都要+0.000001
    f_cout = open('count.txt', 'w')
    for weight in weight_value:
        f_cout.write('(' + str(weight - 0.05) + ',' + str(weight) + ']' + '\t')
        for window in range(49, 103):
            df = data[str(window)][
                (data[str(window)] > (weight - 0.05 + 0.00001)) & (data[str(window)] <= weight + 0.0001)]
            f_cout.write(str(df.shape[0]) + '\t')
        f_cout.write('\n')
    f_cout.close()


def count_edge():
    f_count = open('count_edge2.txt', 'w')
    txt_filename = glob.glob('G:\project2\\NPM201507\\code\\0802files\\edge_0802aft\\*.txt')
    print 'file number', len(txt_filename)
    for filename in txt_filename:
        window = filename.split('edgeInfo')[1].split('.')[0]
        count = len(open(filename, 'r').readlines())
        f_count.write(str(window) + '\t' + str(count - 1) + '\n')
    f_count.close()


def count_term():
    f_count = open('count_term2.txt', 'w')
    txt_filename = glob.glob('G:\project2\\NPM201507\\code\\0802files\\term_0802aft\\*.txt')
    print 'file number', len(txt_filename)
    for filename in txt_filename:
        window = filename.split('termInfo')[1].split('.')[0]
        count = len(open(filename, 'r').readlines())
        f_count.write(str(window) + '\t' + str(count - 1) + '\n')
    f_count.close()

if __name__ == '__main__':
    count_term()
    # count_edge()