#!/usr/bin/env python

import sys
import os
import numpy as np
import subprocess
import shutil
import matplotlib.pyplot as plt

# Set 'flag' to 0, for preparing and submitting input files for RT-TD Gaussian calculations.
# Set --"--  to 1, for creating the data files for analyses (currently just occupation).
# Set --"--  to 2, for getting the avg. nos. for calculations with different parameters.

flag = int(float(sys.argv[1]))


#
# Objects and methods for generation of heat maps
#
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
        temp_x = float(list[0])*100
        temp_y = float(list[1])
        temp_z = float(list[2])
        x.append(temp_x)
        y.append(temp_y)
        z.append(temp_z)
        map_value[(temp_x, temp_y)] = temp_z

    return x, y, map_value

def draw_heatmap(x, y, map_value, title):

    plt_x = np.asarray(list(set(x)))
    plt_y = np.asarray(list(set(y)))
    plt_z = np.zeros(shape = (len(plt_x), len(plt_y)))

    for i in range(len(plt_x)):
        for j in range(len(plt_y)):
            if map_value.has_key((plt_x.item(i), plt_y.item(j))):
                plt_z[i][j] = map_value[(plt_x.item(i), plt_y.item(j))]

    z_min = plt_z.min()
    z_max = plt_z.max()
    plt_z = np.transpose(plt_z)

    plot_name = title
    img_name = title+'.png'


    color_map = plt.cm.gist_heat #plt.cm.rainbow #plt.cm.hot #plt.cm.gist_heat
    plt.clf()
    plt.pcolor(plt_x, plt_y, plt_z, cmap=color_map, vmin=z_min, vmax=z_max)
    plt.axis([plt_x.min(), plt_x.max(), plt_y.min(), plt_y.max()])
    plt.title(plot_name)
    plt.colorbar().set_label(plot_name, rotation=270)
    ax = plt.gca()
    ax.set_aspect('equal')
    figure = plt.gcf()
    plt.savefig(img_name)
#    plt.show()
    return figure

#
# End of the heat map part of the code
#


h = 6.626e-34

NMOs = 34			# no. of MOs

Nsimsteps = 10000		# no. of time-steps to simulate post-field
Nsteps = 10			# no. of points for field time-steps
Nfield = 15			# no. of points for field strength
buffsteps = 3000		# no. of time-steps (after switching off the field) to avg. quantities AFTER
avgwindow = 6000		# no. of time-steps to avg. quantities over

energyeV = 8.95
stepsizefs = 0.002

dens_file = 'hf_h_s1_density.txt'

inp_key = 'main.txt'
key_file = open(inp_key, 'r')

key_lines = key_file.readlines()

key_file.close()

energyau = energyeV*0.0367493	# S0-S1 resonant frequency in Ha/a.u.

periodfs = 1.0e15*h/(energyeV*1.60218e-19)

tonsteps = int(periodfs/stepsizefs)	# total time-steps per cycle
#timesteps = Nsimsteps + tonsteps

steps = []
time = []

for i in range(1,Nsteps+1):
    steps.append(i*tonsteps+Nsimsteps)
    time.append('%3.3f' % float(i*periodfs))

#print 'size of steps, time = ',len(steps),',',len(time)

#print 'ton_steps: \n', steps
#print 'ton: \n', time

strength = 0.005		# basic unit value of field strength

field = []

for i in range(1,Nfield+1):
    field.append('%3.3f' % float(i*strength))

#print 'field_strength: \n', field

for (n, line) in enumerate(key_lines):
    if ('%chk=' in line):
        elements = line.split('=')
        elements = elements[1].split('.')
        chkfile = elements[0]
#        print 'chkpt file: ', chkfile
    if (' 5/33' in line):
        elements = []
#        print 'line: ',line
        elements = line.split("=")
#        print 'elements: \n', elements
        elements = elements[4].split(',')
#        print 'elements[4].split(','): ',elements
        timesteps_tot = elements[0]
    if ('! 1   2' in line):
        elements = []
        elements = key_lines[n+1].split()
#        print 'elements: \n',elements
#        print 'elements[1/7]: ',elements[1],' ',elements[7]
        fieldstrengthx = elements[1]
        fieldtime = elements[7]

tmp = 'tmp.tmp'
tmp2 = 'tmp2.tmp'
new_inp_key = 'main.tmp'

name1 = 'sys1_h_field_'
name2 = '_cyc_'
fileformat = '.com'

if (flag == 0):
    for i in range(Nsteps):
        for j in range(Nfield):
            filename = name1+str(field[j])+name2+str(i+1)
            comfile = filename+fileformat
            dumstr1 = 's/'+str(timesteps_tot)+'/'+str(steps[i])+'/g;'
            dumstr2 = 's/'+str(fieldstrengthx)+'/'+str(field[j])+'/g;'
            dumstr3 = 's/'+str(fieldtime)+'/'+str(time[i])+'/g;'
            dumstr4 = 's/'+str(chkfile)+'/'+filename+'/g'
            dumstr = dumstr1+dumstr2+dumstr3+dumstr4
#        print type(dumstr1)
#
            new_key_file = open(new_inp_key, 'w')
            tmp_file = open(tmp, 'w')
#
            sub = subprocess.call(['sed', dumstr, inp_key], stdout=new_key_file)
            shutil.copy(new_inp_key, comfile)
            dumstr = 'cat '+dens_file+' >> '+comfile
            sub = subprocess.call([dumstr], shell=True, stdout=tmp_file)
            dumstr = 'gsub '+comfile
            sub = subprocess.call([dumstr], shell=True, stdout=tmp_file)
#
            new_key_file.close()
            tmp_file.close()
            
elif (flag == 1):
    for i in range(Nsteps):
        for j in range(Nfield):
            filename = name1+str(field[j])+name2+str(i+1)
            comfile = filename+fileformat
#
            tmp_file = open(tmp, 'w')
#
            dumstr = './rttddft_extractor.py '+filename+' '+str(NMOs)+' '+str(steps[i])
            sub = subprocess.call([dumstr], shell=True, stdout=tmp_file)
#
            tmp_file.close()
            
elif (flag == 2):
    orb = int(float(sys.argv[2]))   # MO averaged over
    dumstr = 'avg_occup_orb_'+str(orb)+'.dat'
    occup_file = open(dumstr, 'w')
    for i in range(Nsteps):
        for j in range(Nfield):
            tocc = np.zeros([avgwindow+1,NMOs], np.float64)
            occup = 'time_occup.'+name1+str(field[j])+name2+str(i+1)+'.txt'
#            field = 'time_field.'+name1+str(field[j])+name2+str(i+1)+'.txt'
#            dipole = 'time_dipole.'+name1+str(field[j])+name2+str(i+1)+'.txt'
            file_occup = open(occup, 'r')
            occup_lines = file_occup.readlines()
            avg = 0.0
            k=0
            #n,line = enumerate(occup_lines)
#            elements = np.zeros([NMOs+1], np.float64)
            for o in range(1,avgwindow+1):
                elements = occup_lines[buffsteps+o].split()
#                elements = float(elements)
                tocc[o][:] = elements[1:]
#                if (o == 1 or o == avgwindow):
#                    print 'tocc ',o,' \n', tocc[o][:]
                avg += tocc[o][orb-1]
                k += 1
            avg /= k
            occup_file.write('%2.6f %3.6f %3.6f \n' % (float(field[j]), float(time[i]), float(avg)))
    occup_file.close()
elif (flag == 3):
    gp = 'occnum_hf.gp'
    gp_file = open(gp, 'r')
    for i in range(Nsteps):
        for j in range(Nfield):
            dumstr1 = name1+str(field[j])+name2+str(i+1)
            dumstr = 's/'+'sys1_h_field_0.015_cyc_4'+'/'+dumstr1+'/g'
            tmp = 'occ_plot.gp'
#
            tmp_file = open(tmp, 'w')
            tmp_file2 = open(tmp2, 'w')
#
            sub = subprocess.call(['sed', dumstr, gp], stdout=tmp_file)
            sub = subprocess.call(['gnuplot', tmp], stdout=tmp_file2)
#
            tmp_file.close()
            tmp_file2.close()
    gp_file.close()
elif (flag == 4):
    orb = int(float(sys.argv[2]))
    orb_name = 'avg_occup_orb_'+str(orb)
    dat_file_name = orb_name+'.dat'
    x, y, map_value = get_xyz_from_dat_file(dat_file_name)
    draw_heatmap(x, y, map_value, orb_name)
