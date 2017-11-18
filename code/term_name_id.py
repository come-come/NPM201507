import os
def mapping():
    filename = 'G:\project2\\NPM201507\\code\\term_name_id\\termN_Id.txt'
    term = 'G:\project2\\NPM201507\\code\\term_name_id\\terms.txt'

    dic = {}
    fr = open(term, 'r')
    fr2 = open(filename, 'r')
    fw = open('G:\project2\\NPM201507\\code\\term_name_id\\none.txt', 'w')
    for line in fr2.readlines():
        term, idd = line.strip().split('\t')
        dic[term] = idd
    for line in fr.readlines():
        if dic.has_key(line.strip()):
            continue
        else:
            fw.write(line)

def translate(filename, dic):

    fw_name = 'edge\edge_id\edgeInfo' + str(window) + '.txt'
    fw_path = curr_path + os.path.normpath(fw_name)
    fw = open(fw_path, 'w')
    for line in open(filename):
        gn1, gn2 = line.strip().split('\t')
        fw.write(dic[gn1] + '\t' + dic[gn2] + '\n')
    fw.close()


if __name__ == '__main__':
    curr_path = os.getcwd() + '/'
    dic = {}
    fr = open('G:\project2\\NPM201507\\code\\term_name_id\\termN_Id.txt', 'r')
    for line in fr:
        term, idd = line.strip().split('\t')
        dic[term] = idd


    fr2 = open('G:\project2\\NPM201507\\code\\result_c5_s10_v2_weight.txt', 'r')
    title = fr2.readline()
    fw = open('qqqq.txt','w')
    for line in fr2:
        names = line.strip().split('\t')[0]
        name1, name2 = names.strip().split('_')
        id1 = dic[name1]
        id2 = dic[name2]
        ids = id1 + '_' + id2
        fw.write(ids+'\n')

    '''
    # translate(filename, dic)
    for window in range(0, 50):
        filename = 'edge/edgeInfo' + str(window) + '.txt'
        f_path= curr_path + os.path.normpath(filename)
        translate(f_path, dic)
        print f_path
    '''