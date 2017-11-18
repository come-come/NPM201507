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

X = [-0.105625455, -0.174359608, -0.105873346, -0.20576805399999998, -0.120525806, -0.28215594, -0.093060075, -0.14982578900000001, -0.119486435, -0.151960108, -0.10104339999999999, -0.23913942, -0.106686264, -0.16127649300000002]
Y  = [0.177165099, -0.029949987, 0.145072897, -0.038238821, 0.157384682, -0.035901472000000004, 0.138161626, -0.037059772000000005, 0.151773712, -0.04529447, 0.14621367300000002, -0.036671638, 0.16619656800000002, -0.032974679]
Z = [0.166806127, 0.20898661899999998, 0.187877354, 0.226367115, 0.18765406699999998, 0.22494885899999997, 0.19316028699999999, 0.212110438, 0.204398195, 0.22953668100000002, 0.205943133, 0.231321515, 0.194322528, 0.246532118]
x = [-0.105625455, -0.174359608, -0.105873346, -0.20576805399999998, -0.120525806, -0.28215594, -0.093060075, -0.14982578900000001, -0.119486435, -0.151960108, -0.10104339999999999, -0.23913942, -0.106686264, -0.16127649300000002]
y  = [0.177165099, -0.029949987, 0.145072897, -0.038238821, 0.157384682, -0.035901472000000004, 0.138161626, -0.037059772000000005, 0.151773712, -0.04529447, 0.14621367300000002, -0.036671638, 0.16619656800000002, -0.032974679]
z = [0.166806127, 0.20898661899999998, 0.187877354, 0.226367115, 0.18765406699999998, 0.22494885899999997, 0.19316028699999999, 0.212110438, 0.204398195, 0.22953668100000002, 0.205943133, 0.231321515, 0.194322528, 0.246532118]

'''
[X1,Y1,Z1]=griddata(x,y,z, method = 'cubic')
print X1, Y1, Z1

print Z1
# Z=griddata((x,y),z,np.linspace(min(x),max(x)),np.linspace(min(y),max(y)))
fig = plt.figure()
ax = Axes3D(fig)
surf = ax.plot_surface(
    X1, Y1, Z1 )
plt.show()
'''

'''
xyz = {'x': X, 'y': Y, 'z':Z}
df = pd.DataFrame(xyz, index=range(len(xyz['x'])))
fig = plt.figure()
ax = Axes3D(fig)
surf = ax.plot_trisurf(df.x, df.y, df.z, cmap=cm.jet, linewidth=0.1)
fig.colorbar(surf, shrink=0.5, aspect=5)
plt.savefig('teste.pdf')
plt.show()
'''



# xyz = {'x': X, 'y': Y, 'z':Z}
# df = pd.DataFrame(xyz, index=range(len(xyz['x'])))
# x1 = np.linspace(df['x'].min(), df['x'].max())
# y1 = np.linspace(df['y'].min(), df['y'].max())
x1 = np.linspace(min(X), max(X))
y1 = np.linspace(min(Y), max(Y))
x2, y2 = np.meshgrid(x1, y1)
z2 = griddata((X, Y), Z, (x2, y2), method='cubic')

# z2 = griddata((df['x'], df['y']), df['z'], (x2, y2), method='cubic')
fig = plt.figure()
# ax = fig.gca(projection='3d')
ax = Axes3D(fig)
surf = ax.plot_surface(x2, y2, z2, rstride=1, cstride=1, cmap=cm.coolwarm, linewidth=0, antialiased=False)
ax.zaxis.set_major_locator(LinearLocator(10))
ax.zaxis.set_major_formatter(FormatStrFormatter('%.02f'))
fig.colorbar(surf, shrink=0.5, aspect=5)
plt.show()


def phenotype(gene_set, start_time, end_time):
    phe1 = pd.read_table('G:\project2\\NPM201507\\data\\IDMapping_consolidated_allPhi2_cleaned_lfc_avg.txt', index_col=0)
    phe2 = pd.read_table('G:\project2\\NPM201507\\data\\IDMapping_consolidated_allQESV_cleaned_LFC_avg.txt',index_col=0)
    phe3 = pd.read_table('G:\project2\\NPM201507\\data\\IDMapping_consolidated_allQI_new_RAW3_adj_LFC_avg.txt', index_col=0)

    window_set = [i for i in range(start_time,end_time)]
    phe1.columns = [i for i in range(0, 113)]
    phe2.columns = [i for i in range(0, 113)]
    phe3.columns = [i for i in range(0, 113)]

    q = phe1.loc[gene_set, window_set]
    q2 = phe2.loc[gene_set, window_set]
    q3 = phe3.loc[gene_set, window_set]

    # c = ['red', 'blue', 'yellow', 'black', 'purple','orange']
    # j = 0
    # for i in q.index:
    #     color = plt.cm.Set2(random.choice(xrange(plt.cm.Set2.N)))
    #     color_list = plt.cm.Set3(np.linspace(0, 1, 12))
    #     ax.scatter(q.loc[i], q2.loc[i],q3.loc[i], c=c[j],s = 40)
    #     j = j + 1
    #
    # ax.set_xlabel('Phi')
    # ax.set_ylabel('QESV')
    # ax.set_zlabel('QI')
    #
    # plt.show()

'''
    # 88001
    gene = ['AT1G14345', 'AT1G80380', 'AT2G18790', 'AT2G29180', 'AT4G33010', 'AT5G42270']
    window = [19, 20, 21, 22, 23, 24]
    # 87603
    gene1 = ['AT1G03160', 'AT1G14345', 'AT1G80380', 'AT2G18790', 'AT2G29180', 'AT4G33010', 'AT5G42270']
    window1 = [19, 20, 21, 22, 23]
    # 125845
    gene2 = ['AT1G22450', 'AT1G51350', 'AT3G15840', 'AT3G28850', 'AT4G16155', 'AT4G27700', 'AT5G53170']
    window2 = [25, 26, 27, 28, 29, 30, 31, 32, 33]
    # 141747
    gene3 = ['AT1G50450', 'AT1G67840', 'AT1G74880', 'AT2G21530', 'AT4G19830', 'AT4G38100']
    window3 = [29, 30, 31, 32, 33, 34, 35, 36]
    # 121926
    gene4 = ['AT1G22450', 'AT1G51350', 'AT3G15840', 'AT3G28850', 'AT4G16155', 'AT5G53170']
    window4 = [24, 25, 26, 27, 28, 29, 30, 31, 32, 33, 34, 35, 36, 37, 38, 39]

    gene5 = ['AT1G03160', 'AT2G18790', 'AT3G08920', 'AT4G33010', 'AT5G53170']
    window5 = [0, 1, 2]
    # 5594
    gene6 = ['AT1G03160', 'AT1G14345', 'AT1G67700', 'AT1G72640', 'AT2G20260', 'AT2G29180', 'AT3G01440', 'AT3G08920', 'AT3G14420', 'AT4G27700', 'AT4G33010']
    window6 = [3, 4]
    # 77463
    gene7 = ['AT1G22450', 'AT1G51350', 'AT1G75690', 'AT3G15840', 'AT3G28850', 'AT4G16155', 'AT4G27700']
    window7 = [17, 18, 19, 20, 21, 22, 23, 24, 25, 26, 27, 28, 29, 30, 31, 32, 33]
    # 23520
    gene8 = ['AT1G03090', 'AT1G50450', 'AT1G54500', 'AT1G67840', 'AT1G74880', 'AT2G21530', 'AT4G38100']
    window8 = [8, 9, 10, 11]
    phenotype(gene, window)
    phenotype(gene1, window1)
    phenotype(gene2, window2)
    phenotype(gene3, window3)
    phenotype(gene4, window4)
    phenotype(gene5, window5)
    phenotype(gene6, window6)
    phenotype(gene7, window7)
    phenotype(gene8, window8)
    '''



phenotype(['AT1G14345', 'AT1G80380', 'AT2G18790', 'AT2G29180', 'AT4G33010', 'AT5G42270'], 68, 82)
