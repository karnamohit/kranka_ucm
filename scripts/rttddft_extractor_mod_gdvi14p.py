#!/usr/bin/env python
"""
program to extract time, MO populations, instantaneous dipole moment and electric field info from 
an rttddft simulation using the older modified link 512 and the Gaussian development version (gdv)

Running the script:
$ ./rttddft_extractor_gdvi14p.py <.LOG file-name, w/o the extension> <# of MOs> <# of time-steps>
"""
#
# imports
 
import sys
import numpy
import math
import shutil
import subprocess
import os

#
# parameters

logfile_name = sys.argv[1]
logfile = logfile_name + '.log' # 'koop_nofield_rhf.log' the gdv-generated .log file
outp = 'logfile.tmp'

shutil.copy(logfile, outp)         # Copying the orig output file to 'logfile.tmp', without having to do this manually as previously required - KR

out_file = open(outp, 'w')

# Gaussian prints exponentials with the form 1.0D+00, which python doesn't recognize, so we
# go through the log file and change it to 1.0E+00
sub = subprocess.call(['sed', 's/D/E/g', logfile], stdout=out_file )

#infile.close
out_file.close

NMOs = int(float(sys.argv[2])) # # of MOs
sim_time = int(float(sys.argv[3])) + 1 # # of time-steps
test_steps = 1500
#tstep = 0.02

#
# open logfile, read all lines into memory and close logfile

infile = open(outp, 'r')
lines = infile.readlines()
infile.close

# initialize time array and extract system time
time = numpy.zeros([sim_time], numpy.float)
ct = 0
for (n,line) in enumerate(lines):
    if (' Time = ' in line and ct < sim_time):
#    if (' Time = ' in line and ct < test_steps):
        elements = lines[n].split()
#        print elements
        time[ct] = elements[4]
        ct += 1 
ntime = ct 
print 'ntime = ', ntime

# initialize population array and extract MO populations
# RTTDDFT program only outputs MOs every other time step, 
# so we'll just read in t(i) and t(i+1) as the same population
#
mo_pop = numpy.zeros([sim_time,NMOs], numpy.float64)
ct = 0
loops = NMOs / 6 + 1
last = NMOs % 6
print 'loops, last = ', loops, last
for (n,line) in enumerate(lines):
    if (' occupation numbers:' in line and ct < (sim_time-1)):
#    if (' occupation numbers:' in line and ct < test_steps):
        for k in range(loops):
           if k == loops-1:
               end = last
           else:
               end = 6
           for j in range(end):
              elements = lines[n+k+1].split()
              s = k*6 + j
              mo_pop[ct][s] = elements[j]
              mo_pop[ct+1][s] = elements[j]
              mo_pop[ct][s] *= 1.
              mo_pop[ct+1][s] *= 1. 
        ct += 1 
mo_time = ct        
print 'mo_time = ', mo_time

#initialize dipole array and extract instanteous dipole moment
idipx = numpy.zeros([sim_time+1], numpy.float64)
idipy = numpy.zeros([sim_time+1], numpy.float64)
idipz = numpy.zeros([sim_time+1], numpy.float64)
totdip = numpy.zeros([sim_time+1], numpy.float64)
ct = 0
for (n,line) in enumerate(lines):
    if ('Eipole Moment (Eebye):' in line and ct < sim_time):
#    if ('Eipole Moment (Eebye):' in line and ct < test_steps):
        elements = lines[n+1].split()
        idipx[ct] = elements[1]
        idipy[ct] = elements[3]
        idipz[ct] = elements[5]
        totdip[ct] = elements[7]
        ct += 1
idip_time = ct
print 'idip_time = ', idip_time 

#initialize electric field array and extract x, y, z components
efield = numpy.zeros([sim_time, 3], numpy.float64)
ener = numpy.zeros([sim_time, 1], numpy.float64)
ct = 0
totener = 0
etmp = 0
for (n,line) in enumerate(lines):
    if ('A dipole field of (real)' in line and ct < sim_time):
#    if ('A dipole field of (real)' in line and ct < test_steps):
        elements = lines[n].split()
        efield[ct][0] = elements[5]
        efield[ct][1] = elements[6]
        efield[ct][2] = elements[7]
        if (ct > 0):
            etmp = idipz[ct-1]*efield[ct][2]
        ener[ct] = etmp
        totener += etmp
#        ener[ct] = totener
        ct += 1
#ct = 0
#totener = 0
#for (n,line) in enumerate(lines):
#    if (' Energy =  ' in line):
#    if (' Energy =  ' in line and ct < test_steps):
#        elements = lines[n].split()
#        print elements[2]
#        ener[ct] = elements[2]
#        totener += ener[ct]
#        ct += 1
field_time = ct
print 'field_time = ', field_time
print 'total_energy = ', totener

# have all the data we need. now output to files for plotting
#

fieldfile_name = 'time_field.' + logfile_name + '.txt'
dipolefile_name = 'time_dipole.' + logfile_name + '.txt'
popfile_name = 'time_occup.' + logfile_name + '.txt'

field = open(fieldfile_name, 'w')
dipole = open(dipolefile_name, 'w')
tpop = open(popfile_name, 'w') 

for i in range(field_time):
#    dipole.write('%14.8f %14.8f %14.8f %14.8f' % (time[i], idipx[i], idipy[i], idipz[i]) + "\n")    
    field.write('%14.8f %14.8f %14.8f %14.8f %14.8f' % (time[i], efield[i][0], efield[i][1], efield[i][2], ener[i]) + "\n")

for i in range(mo_time):
    dipole.write('%14.8f %14.8f %14.8f %14.8f %14.8f' % (time[i], idipx[i], idipy[i], idipz[i], totdip[i]) + "\n")
    tpop.write('%14.8f' % time[i])
    for j in range(NMOs):
        tpop.write('%14.8f' % mo_pop[i][j])
    tpop.write(" \n")

field.close()
dipole.close()
tpop.close()  















