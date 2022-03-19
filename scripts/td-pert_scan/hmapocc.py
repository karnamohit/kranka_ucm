#!/usr/bin/env python

import numpy as np
import matplotlib.pyplot as plt
import matplotlib.ticker as ticker

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
    x = []
    y = []
    z = []
    map_value = {}
    for line in open(dat_file_path):
        list = line.split()
        temp_x = float(list[0])*100
        temp_y = float(list[1])
        temp_z = 2*float(list[2])
        x.append(temp_x)
        y.append(temp_y)
        z.append(temp_z)
        map_value[(temp_x, temp_y)] = temp_z
    return x, y, map_value

def draw_heatmap(x, y, map_value, title, orb):
    #
    plt_x = np.asarray(list(set(x)))
    plt_y = np.asarray(list(set(y)))
    plt_z = np.zeros(shape = (len(plt_x), len(plt_y)))
    for i in range(len(plt_x)):
        for j in range(len(plt_y)):
            if map_value.__contains__((plt_x.item(i), plt_y.item(j))):
                plt_z[i][j] = map_value[(plt_x.item(i), plt_y.item(j))]
    z_min = plt_z.min()
    z_max = plt_z.max()
    plt_z = np.transpose(plt_z)
    #
    plot_name = 'avg_occup_orb_{}'.format(orb)
    img_name = plot_name+'.png'
    #
    def myfmt(x, pos):
        return '{0:.2f}'.format(x)
    #
    color_map = plt.cm.rainbow # cm.rainbow, cm.hot, cm.gist_heat
    plt.clf()
    y_diff = abs(plt_y[0] - plt_y[1])/2
    x_diff = abs(plt_x[0] - plt_x[1])/2
    plt.pcolor(plt_x, plt_y, plt_z, cmap=color_map, vmin=z_min, vmax=z_max)
    plt.axis([plt_x.min()-x_diff, plt_x.max()+x_diff, plt_y.min()-y_diff, plt_y.max()+y_diff])
    plt.title(title)
    plt.xlabel(r"field-strength ($\times 10^2$) (au)")
    plt.ylabel("field-on time (fs)")    
    plt.colorbar(format=ticker.FuncFormatter(myfmt)).set_label('MO occupation number', rotation=270, labelpad=20) # set_label(y=0.45)
    #
    ax = plt.gca()
    ax.set_xticks(list(plt_x))
    labels = []
    for i in plt_x:
        labels.append('{:.1f}'.format(i))
    ax.set_xticklabels(labels)
    ax.set_yticks(list(plt_y))
    labels = []
    for i in plt_y:
        labels.append('{:.2f}'.format(i))
    ax.set_yticklabels(labels) 
    ax.set_aspect('auto') # set_aspect('equal')
    figure = plt.gcf()
    plt.grid()
    # plt.show()
    #
    plt.savefig(img_name)  
    return figure