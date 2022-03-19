#!/usr/bin/env python
"""
author(s): Brad Habenicht (bhf), Karnamohit Ranka (KR)
program to extract time, MO populations, instantaneous dipole moment and electric field info from 
an rttddft simulation using the older modified link 512 and the Gaussian development version (gdv)

Running the script:
$ python rttddft_extractor_old512.py <.LOG file-name, w/o the extension> <# of MOs> <# of time-steps>
    OR
$ ./rttddft_extractor_old512.py <.LOG> <# MOs> <# time-steps>
"""
#
# imports
 
import sys
import numpy as np
import subprocess

def find_linenum(s,txtfile): # locate matching string's line number(s) using `grep`
    strg = subprocess.Popen(['grep','-n', s, txtfile], stdout=subprocess.PIPE).stdout.read().decode("utf-8")
    lst = strg.split("\n")[:-1]
    lst = [int(lst[i].split(":")[0]) - 1 for i in range(len(lst))]
    return lst

#
# parameters

logfile_name = sys.argv[1]
logfile = logfile_name + '.log'			    # 'koop_nofield_rhf.log' the gdv-generated .log file
NMOs = int(float(sys.argv[2])) 			    # of MOs
sim_time = int(float(sys.argv[3]))
sim_time = sim_time + (sim_time % 2) + 2	# of time-steps
# tstep = 0.02

#
# open logfile, read all lines into memory and close logfile
infile = open(logfile, 'r')
lines = infile.readlines()
infile.close

# initialize time array and extract system time
time = np.zeros([sim_time], np.float64)
ct = 0
lst = find_linenum('Ehrenfest time :',logfile)
for n in lst:
    elements = lines[n].split()
    time[ct] = float(elements[5])
    ct = ct + 1 
ntime = ct 
print('ntime = ', ntime)

# '''
# initialize population array and extract MO populations
# RTTDDFT program only outputs MOs every other time step, 
# so we'll just read in t(i) and t(i+1) as the same population
mo_pop = np.zeros([sim_time,NMOs], np.float64)
ct = 0
loops = int(NMOs / 6) + 1
last = NMOs % 6
lst = find_linenum('Alpha orbital occupation numbers',logfile)
occ = 2
if len(lst) == 0:
    lst = find_linenum('Orbital occupation numbers',logfile)
    occ = 1
for n in lst:
    for k in range(loops):
        if k == loops-1:
            end = last
        else:
            end = 6
        for j in range(end):
            elements = lines[n+k+1].split()
            s = k*6 + j
            mo_pop[ct][s] = occ*float(elements[j])
            mo_pop[ct+1][s] = occ*float(elements[j])
    ct = ct + 2
mo_time = ct - 2
print('mo_time = ', mo_time)
# '''

# '''
# output MO occupation data to file
popfile_name = 'time_occup.' + logfile_name + '.txt'
tpop = open(popfile_name, 'w')

for i in range(mo_time):
    tpop.write('{:14.8f}'.format(time[i]))
    for j in range(NMOs):
        tpop.write("{:14.8f}".format(mo_pop[i][j]))
    tpop.write(" \n")

tpop.close()  
# '''

# '''
# initialize dipole array and extract instanteous dipole moment
idipx = np.zeros([sim_time+1], np.float64)
idipy = np.zeros([sim_time+1], np.float64)
idipz = np.zeros([sim_time+1], np.float64)
ct = 0
lst = find_linenum('Real-Time Dipole Moment',logfile)
for n in lst:
    elements = lines[n].split()
    idipx[ct] = float(elements[7])
    idipy[ct] = float(elements[11])
    idipz[ct] = float(elements[15])
    ct = ct + 1
idip_time = ct - 1
print('idip_time = ', idip_time)

# initialize electric field array and extract x, y, z components
efield = np.zeros([sim_time, 3], np.float64)
ener = np.zeros([sim_time], np.float64)
ct = 0
totener = 0
lst = find_linenum('Current electric field:',logfile)
for n in lst:
    elements = lines[n+1].split()
    efield[ct][0] = float(elements[2])
    efield[ct][1] = float(elements[6])
    efield[ct][2] = float(elements[10])
    etmp = np.linalg.norm(efield[ct][:])
    totener = totener + etmp
    ener[ct] = totener
    ct = ct + 1
print()
field_time = ct - 1
print('field_time = ', field_time)
print('total_energy = ', totener)

#
# output data to files
fieldfile_name = 'time_field.' + logfile_name + '.txt'
dipolefile_name = 'time_dipole.' + logfile_name + '.txt'
field = open(fieldfile_name, 'w')
dipole = open(dipolefile_name, 'w')

try:
    for i in range(field_time):
        dipole.write("{:14.8f} {:14.8f} {:14.8f} {:14.8f}\n".format(time[i], idipx[i], idipy[i], idipz[i]))    
        field.write("{:14.8f} {:14.8f} {:14.8f} {:14.8f} {:14.8f}\n".format(time[i], efield[i][0], efield[i][1], efield[i][2], ener[i]))
except Exception as e:
    print(e)

field.close()
dipole.close()
# '''