#!/usr/bin/env python

import numpy as np
import matplotlib.pyplot as plt
# import matplotlib.ticker as ticker
import seaborn as sns

'''

 Objects and functions for generation of heat maps: 
    for MO occupations vs perturbation field parameters 

'''
def get_xyz_from_dat_file(dat_file_path):
    '''
    get x, y, z value from dat file.
    dat file format:
        x0 y0 z0
        x1 y1 z1
        ...
    '''
    x, y, map_value = [], [], {}
    for line in open(dat_file_path):
        lst = line.split()
        temp_x = float(lst[0])*100
        temp_y = float(lst[1])
        temp_z = float(lst[2])
        x.append(temp_x)
        y.append(temp_y)
        map_value[(temp_x, temp_y)] = temp_z
    x = np.asarray(list(set(x)))
    y = np.asarray(list(set(y)))[::-1]
    return x, y, map_value

def draw_heatmap(x, y, map_value, title, orb, mp=False, colorbar=False, sq_shape=False,savefig=False):
    #
    plt_z = np.zeros(shape = (len(x), len(y)))
    for i in range(len(x)):
        for j in range(len(y)):
            if map_value.__contains__((x.item(i), y.item(j))):
                plt_z[i,j] = map_value[(x.item(i), y.item(j))]
    plt_z = np.transpose(plt_z)
    z_min = plt_z.min()
    z_max = plt_z.max()
    #
    img_name = 'avg_occup_orb_{}.png'.format(orb)
    #
    # def myfmt(x, pos):
    #     return '{0:.2f}'.format(x)
    #
    plt.clf()
    # 
    ax = plt.gca()
    ax.set_xticks(list(x))
    ax.set_yticks(list(y))
    xlabels = list(np.around(x,decimals=2))
    ylabels = list(np.around(y,decimals=2))
    #
    if mp:
        #
        color_map = plt.cm.rainbow # cm.rainbow, cm.hot, cm.gist_heat
        y_diff = abs(y[0] - y[1])/2
        x_diff = abs(x[0] - x[1])/2
        plt.pcolor(x, y, plt_z, cmap=color_map, vmin=z_min, vmax=z_max)
        plt.axis([x.min()-x_diff, x.max()+x_diff, y.min()-y_diff, y.max()+y_diff])
        # colorbar(format=ticker.FuncFormatter(myfmt)).set_label(y=0.45)
        if colorbar:
            plt.colorbar(format='%.2f').set_label('MO occupation number', rotation=270, labelpad=20)
        #
        ax.set_xticklabels(xlabels)
        ax.set_yticklabels(ylabels) 
        ax.set_aspect('auto') # set_aspect('equal')
        plt.grid()
    else:
        #
        color_bar = {'format':'%.2f',}#'label':'MO occupation number',
        sns.heatmap(plt_z,vmin=0.0, vmax=2.0,cmap='coolwarm',square=sq_shape,cbar=colorbar,cbar_kws=color_bar,xticklabels=xlabels,yticklabels=ylabels,)
        #
        cbar = ax.collections[0].colorbar
        cbar.set_label(label='MO occupation number', rotation=270, labelpad=20)
    #
    plt.title(title)
    plt.xlabel(r"field-strength ($\times 10^2$) (au)")
    plt.ylabel("field-on time (fs)")    
    #
    figure = plt.gcf()
    # plt.show()
    #
    if savefig:
        plt.savefig(img_name,bbox_inches='tight',)
    return figure