#!/usr/bin/env python

import sys
import os
import numpy as np
import subprocess
import shutil
import matplotlib.pyplot as plt

'''

 Objects and methods for generation of heat maps

'''
def get_xyz_from_dat_file(dat_file_path):
    '''
    get x, y, z value from dat file
    dat file format: x0 y0 z0
    '''
    x = []
    y = []
    z = []
    map_value = {}

    for line in open(dat_file_path):
        list = line.split()
        temp_x = float(list[0]) # -ve for GDV i14+
        temp_y = float(list[1])
        temp_z = float(list[2])
        x.append(temp_x) 
        y.append(temp_y)
        z.append(temp_z)
        map_value[(temp_x, temp_y)] = temp_z

    return x, y, map_value

def draw_heatmap(x, y, map_value, mo, title):

    plt_x = np.asarray(list(set(x)))
    plt_y = np.asarray(list(set(y)))
    plt_z = np.zeros(shape = (len(plt_x), len(plt_y)))

    for i in range(len(plt_x)):
        for j in range(len(plt_y)):
            if map_value.has_key((plt_x.item(i), plt_y.item(j))):
                plt_z[i][j] = map_value[(plt_x.item(i), plt_y.item(j))]

    z_min = 0.0 #plt_z.min()
    z_max = 2.0 #plt_z.max()
    plt_z = np.transpose(plt_z)

    if (title == 'avg_occup_orb_20'):
        plot_name = 'HOMO'
    elif (title == 'avg_occup_orb_21'):
        plot_name = 'LUMO'
    else:
        plot_name = mo
    img_name = title+'.png'


    color_map = plt.cm.gist_heat #plt.cm.rainbow #plt.cm.hot #plt.cm.gist_heat
    plt.clf()
    plt.pcolor(plt_x, plt_y, plt_z, cmap=color_map, vmin=z_min, vmax=z_max)
    plt.axis([plt_x.min(), plt_x.max(), plt_y.min(), plt_y.max()])
    plt.title(plot_name)
    plt.colorbar().set_label('', rotation=270)
    ax = plt.gca()
    ax.set_aspect('auto')
    figure = plt.gcf()
    plt.savefig(img_name)
#    plt.show()
    return figure

'''

 End of the heat map part of the code

'''

