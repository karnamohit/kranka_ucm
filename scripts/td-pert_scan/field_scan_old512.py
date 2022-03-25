#!/usr/bin/env python

import sys
import os
import numpy as np

'''
    running field_scan_old512.py:

        $ python field_scan_old512.py <flag>
            OR
        $ ./field_scan_old512.py <flag>


    Set "flag" to:
        0, for preparing and submitting input files for RT-TDHF Gaussian calculations.
        1, for creating the TD data files for analyses (MO occupation, field, dipole).
        2, for getting the avg. MO occupations for different field perturbation 
               parameters.
        3, for plotting time-dependent MO occupations.
        4, for plotting avg MO occupation numbers as a function of field-on time (fs)
               and field-strength (au)
    Files needed:
        <flag> = 0
            in:     "main_old512.txt", density-matrix file (cf. "hf_h_gs_density.txt"),
                    job submission file (generated by `gsub` here)
            out:    "*.log"
        <flag> = 1
            in:     "*.log"
            out:    "time_*.*.txt"
        <flag> = 2
            in:     "time_occup.*.txt"
            out:    "avg_occup_orb_*.dat"
        <flag> = 3
            in:     "time_occup.*.txt", "occnum.gp"
            out:    "time_occup.*.eps"
        <flag> = 4
            in:     "avg_occup_orb_*.dat"
            out:    "avg_occup_orb_*.png"
'''

flag = int(float(sys.argv[1]))      # get the flag value
read_densmat_from_txt = False 
dens_file = 'hf_h_s1_density.txt'

# Setting up parameters and defining constants
h = 6.626e-34                   # Planck's constant in SI units
NMOs = 58                       # total no. of MOs
HOMO = 34	        		    # HOMO MO no.
LUMO = 35			            # LUMO MO no.
Nsimsteps = 10000               # no. of time-steps to simulate post-field
Nsteps = 10                     # no. of points for field time-steps
Nfield = 15                     # no. of points for field strength
buffsteps = 3000                # no. of time-steps (after switching off the field) to avg. quantities AFTER
avgwindow = 6000                # no. of time-steps to avg. quantities over
energyeV = 5.989                # S0-S1 resonant frequency in eV
stepsizefs = 0.002              # step-size of electron-dynamics in fs
strength = 0.005		        # field-strength minimum, and increment, in au
energyau = energyeV*0.0367493
periodfs = 0.690 #1.0e15*h/(energyeV*1.60218e-19)  # time-period in fs, with freq. "energyeV" (in eV)
tonsteps = 345 #int(periodfs/stepsizefs)	        # total time-steps per cycle

steps = []
time = []
for i in range(1,Nsteps+1):
    steps.append(i*tonsteps+Nsimsteps)
    time.append('%3.3f' % float(i*periodfs))
field = []
for i in range(1,Nfield+1):
    field.append('%3.3f' % float(i*strength))

if (flag != 4):
    inp_key = 'main_old512.txt'
    key_file = open(inp_key, 'r')
    key_lines = key_file.readlines()
    key_file.close()
    for (n, line) in enumerate(key_lines):
        if ('%chk=' in line):
            elements = line.split('=')
            elements = elements[1].split('.chk')
            chkfile = elements[0]
        if (' 5/33' in line):
            elements = []
            elements = line.split("=")
            elements = elements[4].split(',')
            timesteps_tot = elements[0]
        if ('! 1   2' in line):
            elements = []
            elements = key_lines[n+1].split()
            fieldstrengthx = elements[1]
            fieldtime = elements[7]
    tmp = 'tmp.tmp'
    tmp2 = 'tmp2.tmp'
    new_inp_key = 'main.tmp'
    name1 = 'sys3_h_field_'
    name2 = '_cyc_'
    fileformat = '.com'
    #
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
                # 
                os.system("sed -e '{}' {} > {}".format(dumstr,inp_key,new_inp_key))
                os.system('cp {} {}'.format(new_inp_key,comfile))
                if read_densmat_from_txt:
                    dumstr = 'cat '+dens_file+' >> '+comfile
                    os.system(dumstr)
                # queue submission command
                dumstr = 'gsub -gdvi14+ '+comfile
                os.system(dumstr)
        # 
    elif (flag == 1):
        for i in range(Nsteps):
            for j in range(Nfield):
                filename = name1+str(field[j])+name2+str(i+1)
                comfile = filename+fileformat
                # 
                dumstr = 'python rttddft_extractor_old512.py {} {} {}'.format(filename,NMOs,steps[i])
                os.system(dumstr)
        # 
    elif (flag == 2):
        for orb in range(HOMO-4,LUMO+5):
            dumstr = 'avg_occup_orb_'+str(orb)+'.dat'
            occup_file = open(dumstr, 'w')
            for i in range(Nsteps):
                for j in range(Nfield):
                    tocc = np.zeros([avgwindow+1,NMOs], np.float64)
                    occup = 'time_occup.'+name1+str(field[j])+name2+str(i+1)+'.txt'
                    file_occup = open(occup, 'r')
                    occup_lines = file_occup.readlines()
                    file_occup.close()
                    avg = 0.0
                    k=0
                    for o in range(1,avgwindow+1):
                        elements = occup_lines[buffsteps+o].split()
                        tocc[o][:] = float(elements[1:])
                        avg += tocc[o][orb-1]
                        k += 1
                    avg /= k
                    occup_file.write('{:2.6f} {:3.6f} {:3.6f} \n'.format(float(field[j]), float(time[i]), float(avg)))
            occup_file.close()
        #
    elif (flag == 3):
        gp = 'occnum_hf.gp'
        for i in range(Nsteps):
            for j in range(Nfield):
                dumstr1 = name1+str(field[j])+name2+str(i+1)
                dumstr = 's/'+'dens_nofield_h_gs_20fs_y_rhf'+'/'+dumstr1+'/g'
                tmp = 'occ_plot.gp'
                # 
                os.system("sed -e '{}' {} > {}".format(dumstr,gp,tmp))
                os.system('gnuplot {}'.format(tmp))
    #
elif (flag == 4):
    # sys.path.append('/path/to/kranka_ucm/sripts/td-pert_scan/')
    from hmapocc import *
    for orb in range(HOMO-4,LUMO+5):
        orb_name = 'avg_occup_orb_'+str(orb)
        dat_file_name = orb_name+'.dat'
        orb_title = 'avg. occupation of MO #{}'.format(orb)
        x, y, map_value = get_xyz_from_dat_file(dat_file_name)
        draw_heatmap(x, y, map_value, orb_title, orb)