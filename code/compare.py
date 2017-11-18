#!/usr/bin/env python
# -*- coding: utf-8 -*-
# __author__ == Sherry

import os


def terms_sign_pro():
    cwd1 = os.getcwd()
    if os.path.exists("1016terms_sign.txt"):
        pass
    f1 = open("1016terms_sign.txt", "r")
    datas = f1.readlines()

    del datas[0:2]
    annotation = []
    annotation_gene = []
    for line in datas:
        annotation.append(line.split("\t"))
    for i in range(100):
        annotation_gene.append(annotation[i][1][1:-1].split(","))

    os.chdir("./data")
    if os.path.exists("annotation_list"):
        pass
    else:
        os.mkdir("annotation_list")
    os.chdir("./annotation_list")
    for i in range(100):
        f = open(str(i + 1)+".txt", "w")
        f.writelines(annotation_gene[i][0][1:-1] + "\n")
        for j in range(len(annotation_gene[i]) - 1):
            f.writelines(annotation_gene[i][j+1][2:-1]+ "\n")
    os.chdir(cwd1)

def GOanalysis():
    annotation_dict = datapro.upgrate_annotation_dict()
    




if __name__ == "__main__":
    terms_sign_pro()
