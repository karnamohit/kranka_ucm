#!/usr/bin/env python
"""
program to run several RT-TD calculations over varying field strengths and cycles.
useful for MO occupation vs field time/strength analysis.

Running the script:
$ ./field_scan_gdvi14p.py <flag>

 - Karnamohit, Jan 2019.
"""

import sys
import os
import numpy as np
import subprocess
import shutil
#import hmapocc 

'''

 Set 'flag' to 0, for preparing and submitting input files for RT-TD Gaussian calculations.
 Set --"--  to 1, for creating the data files for analyses (currently just occupation).
 Set --"--  to 2, for getting the avg. nos. for calculations with different parameters.
 Set --"--  to 3, for getting MO occupation numbers vs. time plots.
 Set --"--  to 4, for getting heat maps of occupation number averages for different MO's,
                      plotted on a field time vs field intensity graph.

'''

flag = int(float(sys.argv[1]))




'''

parameters and constants

'''

h = 6.626e-34

NMOs = 34			# no. of MOs

Nsimsteps = 15000		# no. of time-steps to simulate post-field
Nsteps = 10			# no. of field cycles
Nfield = 15			# no. of points for field strength
buffsteps = 3000		# no. of time-steps (after switching off the field) to avg. quantities AFTER
avgwindow = 11000		# no. of time-steps to avg. quantities over

energyeV = 8.95			# excitation energy in eV
stepsizefs = 0.0012		# step-size of propagation in femtoseconds
strength = 0.005		# value of field strength in arbitrary units


inp_key = 'main.txt'
key_file = open(inp_key, 'r')

key_lines = key_file.readlines()

key_file.close()

energyau = energyeV*0.0367493		# S0-S1 resonant frequency in Ha/a.u.

periodfs = 1.0e15*h/(energyeV*1.60218e-19)
periodau = periodfs*41.34

tonsteps = int(periodfs/stepsizefs)	# total time-steps per field-cycle

steps = []
time = []

for i in range(1,Nsteps+1):
    steps.append(i*tonsteps+Nsimsteps)
    time.append('%3.3f' % float(i*periodau))

#print 'size of steps, time = ',len(steps),',',len(time)

#print 'ton_steps: \n', steps
#print 'ton: \n', time

field = []

for i in range(1,Nfield+1):
    field.append('%3.3f' % float(-i*strength))

#print 'field_strength: \n', field

for (n, line) in enumerate(key_lines):
    if ('%chk=' in line and '_lr' in line):
        elements = line.split('=')
        lrchkfile = elements[1].rstrip()
#        print chkfile
    elif ('%chk=' in line):
        elements = line.split('=')
        chkfile = elements[1].rstrip()
#        print lrchkfile
    if ('maxpoints' in line):
        elements = []
#        print 'line: ',line
        elements = line.split("=")
#        print 'elements: \n', elements
        elements = elements[3].split(')')
#        print 'elements[4].split(','): ',elements
        timesteps_tot = elements[0]
    if ('ElectricLength' in line):
        elements = []
        elements = key_lines[n].split('=')
#        print 'elements: \n',elements
#        print 'elements[1/7]: ',elements[1],' ',elements[7]
        dumstr = elements[5].split('\n')
        fieldstrengthx = dumstr[0]
        elements = elements[2].split()
        fieldtime = elements[0]

tmp = 'tmp.tmp'
tmp2 = 'tmp2.tmp'
new_inp_key = 'main.tmp'

name1 = 'sys1_h_gs_field_'
name2 = '_cyc_'
name3 = '_lr'
comfileformat = '.com'

if (flag == 0):
    for i in range(Nsteps):
        for j in range(Nfield):
            filename = name1+str(field[j])+name2+str(i+1)
            comfile = filename+comfileformat
            dumstr1 = 's/'+str(timesteps_tot)+'/'+str(steps[i])+'/g; '
            dumstr2 = 's/'+str(fieldstrengthx)+'/'+str(field[j])+'/g; '
            dumstr3 = 's/'+str(fieldtime)+'/'+str(time[i])+'/g; '
            dumstr4 = 's/'+str(chkfile)+'/'+filename+'/g; '
            dumstr5 = 's/'+str(lrchkfile)+'/'+filename+name3+'/g; '
            dumstr = dumstr1+dumstr2+dumstr3+dumstr4+dumstr5
#            print dumstr
#
            new_key_file = open(new_inp_key, 'w')
            tmp_file = open(tmp, 'w')
#
#            sedcmd = 'sed '+dumstr
#            print sedcmd
            sub = subprocess.call(['sed', dumstr, inp_key], stdout=new_key_file)
            shutil.copy(new_inp_key, comfile)
#
            dumstr = 'gsub -gdvi14+ '+comfile
            sub = subprocess.call([dumstr], shell=True, stdout=tmp_file)
#
            new_key_file.close()
            tmp_file.close()
            
elif (flag == 1):
    for i in range(Nsteps):
        for j in range(Nfield):
            filename = name1+str(field[j])+name2+str(i+1)
            comfile = filename+comfileformat
#
            tmp_file = open(tmp, 'w')
#
            dumstr = './rttddft_extractor_gdvi14p.py '+filename+' '+str(NMOs)+' '+str(steps[i])
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
                tocc[o][:] = elements[1:]
                avg += tocc[o][orb-1]
                k += 1
            avg /= k
            occup_file.write('%2.6f %3.6f %3.6f \n' % (float(field[j]), float(time[i]), float(avg)))
    occup_file.close()
elif (flag == 3):
    gp = 'occnum_hf_h.gp'
    gp_file = open(gp, 'r')
    for i in range(Nsteps):
        for j in range(Nfield):
            dumstr1 = name1+str(field[j])+name2+str(i+1)
            dumstr = 's/'+'field_h_gs_negdot025_2cyc_8p95eV'+'/'+dumstr1+'/g'
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
#elif (flag == 4):
#    orb = int(float(sys.argv[2]))
#    orb_name = 'avg_occup_orb_'+str(orb)
#    dat_file_name = orb_name+'.dat'
#    x, y, map_value = hmapocc.get_xyz_from_dat_file(dat_file_name)
#    draw_heatmap(x, y, map_value, orb_name)
