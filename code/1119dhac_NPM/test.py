import os
from os.path import join

def read_NPM_result(dest):
    print 'df'
    for root, dirs, files in os.walk(dest):

        for OneFileName in files:

            if OneFileName.find('.txt') == -1:
                continue
            OneFullFileName = join(root, OneFileName)

            dic_result[int(OneFileName.strip().split('.txt')[0])] = []
            print int(OneFileName.strip().split('.txt')[0])
            for line in open(OneFullFileName, 'r'):
                b = line.strip().split('\t')
                # print b  # ['AT656', 'AT43']
                if len(b) > 5:
                    dic_result[int(OneFileName.strip().split('.txt')[0])].append(b)

def read_dhac_result(dest):
    n = 0
    sum = 0
    for root, dirs, files in os.walk(dest):
        for OneFileName in files:
            count = 0
            if OneFileName.find('.group') == -1:
                continue
            OneFullFileName = join(root, OneFileName)
            # print OneFullFileName
            # print OneFileName.strip().split('edgeInfo')[1].split('.group')[0]
            for line in open(OneFullFileName, 'r'):
                b = line.strip().split('\t')
                if len(b) > 5:
                    count = count + 1


            # count = len(open(OneFullFileName, 'r').readlines())
            n = n + 1
            sum = sum + count
            print OneFileName.strip().split('edgeInfo')[1].split('.group')[0], count

    print n, sum
    print sum/float(n)
if __name__ == "__main__":
    # get average cluster number of dhac 38
    # get average cluster number of dhac>5 genes  10
    read_dhac_result('G:\\project2\\NPM201507\\code\\edge\\20171113dhac')