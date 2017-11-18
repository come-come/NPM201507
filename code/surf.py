import numpy as np
import pandas as pd
from mpl_toolkits.mplot3d import Axes3D
import matplotlib.pyplot as plt
import random
from mpl_toolkits.mplot3d import Axes3D
from matplotlib import cm
from matplotlib.ticker import LinearLocator, FormatStrFormatter
import matplotlib.pyplot as plt
import numpy as np
from scipy.interpolate import griddata
from mpl_toolkits.mplot3d import Axes3D
from matplotlib import cm
from matplotlib.colors import LightSource
from sys import argv

phe1 = pd.read_table('G:\project2\\NPM201507\\data\\IDMapping_consolidated_allPhi2_cleaned_lfc_avg.txt', index_col=0)
phe2 = pd.read_table('G:\project2\\NPM201507\\data\\IDMapping_consolidated_allQESV_cleaned_LFC_avg.txt', index_col=0)
phe3 = pd.read_table('G:\project2\\NPM201507\\data\\IDMapping_consolidated_allQI_new_RAW3_adj_LFC_avg.txt', index_col=0)
phe1.columns = [i for i in range(0, 113)]
phe2.columns = [i for i in range(0, 113)]
phe3.columns = [i for i in range(0, 113)]

def draw_3D (term, gene_set, start_time, end_time):
    window_set = [i for i in range(start_time,end_time+1)]
    q = phe1.loc[gene_set, window_set]
    q2 = phe2.loc[gene_set, window_set]
    q3 = phe3.loc[gene_set, window_set]
    print term
    for i in gene_set:
        X = np.array(q.loc[i]).tolist()
        Y = np.array(q2.loc[i]).tolist()
        Z = np.array(q3.loc[i]).tolist()
        xyz = {'x': X, 'y': Y, 'z': Z}
        print X
        print Y
        print Z
        df = pd.DataFrame(xyz, index=range(len(xyz['x'])))
        x1 = np.linspace(df['x'].min(), df['x'].max(), len(df['x'].unique()))
        y1 = np.linspace(df['y'].min(), df['y'].max(), len(df['y'].unique()))
        x2, y2 = np.meshgrid(x1, y1)
        z2 = griddata((df['x'], df['y']), df['z'], (x2, y2), method='cubic')
        fig = plt.figure()
        # ax = fig.gca(projection='3d')
        ax = Axes3D(fig)
        surf = ax.plot_surface(x2, y2, z2, rstride=1, cstride=1, cmap=cm.coolwarm,
                               linewidth=0, antialiased=False)

        ax.zaxis.set_major_locator(LinearLocator(10))
        ax.zaxis.set_major_formatter(FormatStrFormatter('%.02f'))
        fig.colorbar(surf, shrink=0.5, aspect=5)
        plt.show()
        fig = plt.figure()
        ax = Axes3D(fig)


file = open('G:\project2\\NPM201507\\code\\1102terms_sign_list.txt', 'r')
title = file.readline()
for line in file.readlines()[0:12]:
    line_sp = line.strip().split('\t')
    term = line_sp[0]
    gene_set = line_sp[2]
    start_time = line_sp[3]
    end_time = line_sp[4]
    result = gene_set.strip().split(',')
    draw_3D(term, result, int(start_time), int(end_time))