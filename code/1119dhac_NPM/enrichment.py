import os
from  os.path import join
import pandas as pd



if __name__ == "__main__":
    # write the number of enrichment GO terms into dictionary
    dic_go = {}
    # Our Method
    path = 'G:\\project2\\NPM201507\\code\\luguilin\\2\\data_need\\'
    # DHAC
    # path = 'G:\\project2\\NPM201507\\code\\luguilin\\3\\data_need\\'
    # NPM
    # path = 'G:\\project2\\NPM201507\\code\\luguilin\\4\\data_need\\'
    for root, dirs, files in os.walk(path):
        for f in files:
            filename = join(root, f)
            dic_go[int(f.strip().split('.')[0])] = len(open(filename, 'r').readlines())-1
    data = pd.DataFrame.from_dict(dic_go, orient='index')
    data.columns = ['enrichment']
    # read level information
    termFile = 'G:\\project2\\NPM201507\\code\\1119dhac_NPM\\1129terms.txt'
    # DHAC
    # termFile = 'G:\\project2\\NPM201507\\code\\1119dhac_NPM\\1208_DHAC_terms_sign_list_id.txt'
    # NPM
    # termFile = 'G:\\project2\\NPM201507\\code\\1119dhac_NPM\\1211_NPM10_terms_sign_list_id.txt'
    termData = pd.read_csv(termFile, sep='\t',index_col=0)
    # concat two data
    result = pd.concat([data, termData], axis=1)
    for group in result.groupby('level'):
        df1 = group[1]
        zeroNum = df1[df1['enrichment'] == 0].shape[0]
        total = group[1]['enrichment'].shape[0]
        nonzero = total - zeroNum
        percentage = nonzero/float(total)
        # print group[0], '\t', nonzero, '\t', total, '\t', percentage
    for group in termData.groupby(['geneSize','time_size']):
        print group[0],group[1].shape[0]
