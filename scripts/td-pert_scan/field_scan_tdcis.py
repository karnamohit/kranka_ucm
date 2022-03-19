#!/usr/bin/env python

import sys
import os
import numpy as np

'''
    running field_scan_tdcis.py:

        $ python field_scan_tdcis.py <flag>
            OR
        $ ./field_scan_tdcis.py <flag>


    Set "flag" to:
        0, for preparing and submitting input files for RT-TDCIS calculations.
        2, for getting the avg. MO occupations for different field perturbation 
               parameters.
        3, for plotting time-dependent MO occupations.
        4, for plotting avg MO occupation numbers as a function of field-on time (fs)
               and field-strength (au)
    Files needed:
        <flag> = 0
            in:     "TDCI_1pulse_scan.py", "tdm_tdci.txt", job submission file ("pyjob.sub")
            out:    "time_occup.*.txt"
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
NMOs = 34                       # total no. of MOs
HOMO = 20        		        # HOMO MO no.
LUMO = 21                       # LUMO MO no.
Nsimsteps = 10000               # no. of time-steps to simulate post-field
Nsteps = 10                     # no. of points for field time-steps
Nfield = 15                     # no. of points for field strength
buffsteps = 3000                # no. of time-steps (after switching off the field) to avg. quantities AFTER
avgwindow = 6000                # no. of time-steps to avg. quantities over
energyeV = 9.4668               # S0-S1 resonant frequency in eV
stepsizefs = 0.002              # step-size of electron-dynamics in fs
strength = 0.005		        # field-strength minimum, and increment, in au
energyau = energyeV*0.0367493
periodfs = 1.0e15*h/(energyeV*1.60218e-19)  # time-period in fs, with freq. "energyeV" (in eV)
tonsteps = int(periodfs/stepsizefs)	        # total time-steps per cycle

steps = []
time = []
for i in range(1,Nsteps+1):
    # steps.append(i*tonsteps+Nsimsteps)
    time.append('%3.3f' % float(i*periodfs))
    steps.append(i)
field = []
for i in range(1,Nfield+1):
    field.append('%3.3f' % float(i*strength))

if (flag != 4):
    inp_key = 'TDCI_1pulse_scan.py'
    key_file = open(inp_key, 'r')
    key_lines = key_file.readlines()
    key_file.close()
    for (n, line) in enumerate(key_lines):
        if ('REPLACE_CYC' in line):
            elements = line.split('=')[1].split()
            ncyc = elements[0]
        if ('REPLACE_AMP' in line):
            elements = line.split("=")[1].split()
            emax = elements[0]
        if ('REPLACE_TTL' in line):
            elements = line.split("=")[1].split(' ')
            title = elements[2]
    tmp = 'tmp.tmp'
    tmp2 = 'tmp2.tmp'
    tmp3 = 'tmp3.sub'
    name1 = 'sys1_h_field_'
    name2 = '_cyc_'
    name3 = 'temp_'
    name4 = 'pyjob.sub' # queue-job submission file
    name5 = 'calc_tdm_cis'
    #
    if (flag == 0):
        for i in range(Nsteps):
            for j in range(Nfield):
                filename = name1+str(field[j])+name2+str(i+1)
                dumstr1 = 's/'+ncyc+'/'+str(steps[i])+'/g;'
                dumstr2 = 's/'+emax+'/'+str(field[j])+'/g;'
                dumstr3 = 's/'+title+'/'+filename+'/g'
                dumstr = dumstr1+dumstr2+dumstr3
                #
                new_inp_key = name3+str(j)+'_'+str(i)
                os.system("sed -e '{}' {} > {}.py".format(dumstr,inp_key,new_inp_key))
                dumstr = 's/'+str(name5)+'/'+str(new_inp_key)+'/g'
                os.system("sed -e '{}' {} > {}".format(dumstr,name4,tmp3))
                # queue-job submission command
                dumstr = 'qsub '+tmp3
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
                        tocc[o][:] = float(elements[1:-1])
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