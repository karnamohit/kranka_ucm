#!/usr/bin/env python
"""
April 17, 2014
@author B.F. Habenicht
"""
"""
script to extract MO coefficients, dipole matrix, and CI coefficients from Gaussian simulation
and use these to calculate the transition dipole moments between all pairs of states. these will 
later be used by the python script TD_CIS.py to do excited state electron dynamics in an electric 
field. we first find the data in the gaussian log file. we then rotate the dipole matrix into the 
MO basis (from AO) using the MO coefficient matrix. once we have this, we then can multiply the 
appropriate CI coefficients together and multiply those terms times the dipole matrix in the MO basis.
these terms are collected and summed to calculate the transition dipole moments between states.

NOTE!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
In a CASSCF calculation in Gaussian, there are usually several iterations to optimize the CI coefficients
and MO coefficients. Unfortunately, these have the same label 'FINAL COEFFICIENT MATRIX' in the log file.
You will need to go into the file by hand and change all of them to something else (say, 'INAL COEFFICIENT')
except for the last one, so that the program only finds that one. otherwise python responds with- 

Traceback (most recent call last):
  File "calc_tran_dip.py", line 133, in <module>
    MO[i][j] = float(elements[k])
IndexError: list index out of range

NOTE!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

"""

#===============
# IMPORTS
#===============

import sys
import numpy
import os
import subprocess
import shutil

#===============
# PARAMETERS
#===============

Logfile = 'sa-casscf.log'   # Gaussian output file,
#p CASSCF(2,2,nroot=3)/sto-3g scf(tight) nosymm IOp(3/33=3,3/36=1,4/33=3,5/33=3,6/33=3,9/33=3)
NMOs = 2           # total number of MOs 
NStates = 3        # total number of states (CSFs)
NCAS=2              # size of CAS active space  
NOCC=1             # size of CAS space plus orbitals with frozen occupations
NFRZ=0              # number of orbitals with frozen occupations
NELECT=2            # number of electrons in CAS space
NAtoms=2            # number of atoms in molecule 

#===============
debug = 0      # my debug flag set = 1 or 2 to print stuff
#===============

# open log file for reading
#infile = open(Logfile, 'r')


outp = 'logfile.tmp'

shutil.copy(Logfile, outp)         # Copying the orig output file to 'logfile.tmp', without having to do this manually as previously required - KR

out_file = open(outp, "w")

# Gaussian prints exponentials with the form 1.0D+00, which python doesn't recognize, so we 
# go through the log file and change it to 1.0E+00
sub = subprocess.call(['sed', 's/D/E/g', Logfile], stdout=out_file )

#infile.close
out_file.close

infile = open(outp, 'r')

# read entire file 
lines = infile.readlines()

# close file
infile.close()

# define square root of 2
rt2 = numpy.sqrt(2)

# initialize array for molecular dipole matrices
dipX = numpy.zeros([NMOs,NMOs], numpy.float64)
dipY = numpy.zeros([NMOs,NMOs], numpy.float64)
dipZ = numpy.zeros([NMOs,NMOs], numpy.float64)

# find dipoles and fill arrays 
# Gaussian prints the dipole matrices as only lower triangle and the matrix is
# 5 columns across and NMOs rows down. If have less than 6 atoms, don't need to 
# worry about sections of matrix being printed later. we'll fix that part later.

# dipole X

if (NMOs < 6):
    for (n,line) in enumerate(lines):
        if ('Multipole matrices IBuc=  518 IX=    1' in line):
            for i in range (NMOs):
                for j in range (i):
                    elements = lines[n+i+1].split()
                    dipX[i][j] = float(elements[j+1])
                    if ( i != j):
                        dipX[j][i] = dipX[i][j]
else:
    loops = NMOs / 5 + 1
#    print 'loops = ', loops
    last = NMOs % 5
#    print 'last = ', last
    shift = 0
    for (n,line) in enumerate(lines):
        if ('Multipole matrices IBuc=  518 IX=    1' in line):
            for k in range(loops):
                irange = NMOs - k*5
                for i in range(irange):
                    if k == loops:
                        if i < last:
                            end = i + 1
                        else:
                            end = last
                    else:
                        if i <= 4:
                            end = i + 1
                        else:
                            end = 5
#                    print 'k, i, end', k, i, end
                    for j in range(end):
                        elements=lines[n+shift+k*1+i+2].split()
                        s = k*5 + j
                        m = i + k*5
#                        print k, s, m, elements
                        dipX[m][s] = float(elements[j+1])
                        if i != j:
                            dipX[s][m] = dipX[m][s]
                shift = shift + irange

# dipole Y
if (NMOs) < 6:
    for (n,line) in enumerate(lines):
        if ('Multipole matrices IBuc=  518 IX=    2' in line):
            for i in range (NMOs):
                for j in range (i):
                    elements = lines[n+i+1].split()
                    dipY[i][j] = float(elements[j+1])
                    if ( i != j):
                        dipY[j][i] = dipY[i][j]
else:
    loops = NMOs / 5 + 1
#    print 'loops = ', loops
    last = NMOs % 5
#    print 'last = ', last
    shift = 0
    for (n,line) in enumerate(lines):
        if ('Multipole matrices IBuc=  518 IX=    2' in line):
            for k in range(loops):
                irange = NMOs - k*5
                for i in range(irange):
                    if k == loops:
                        if i < last:
                            end = i + 1
                        else:
                            end = last
                    else:
                        if i <= 4:
                            end = i + 1
                        else:
                            end = 5
#                    print 'k, i, end', k, i, end
                    for j in range(end):
                        elements=lines[n+shift+k*1+i+2].split()
                        s = k*5 + j
                        m = i + k*5
#                        print k, s, m, elements
                        dipY[m][s] = float(elements[j+1])
                        if i != j:
                            dipY[s][m] = dipY[m][s]
                shift = shift + irange

if (NMOs < 6):
#dipole Z
    for (n,line) in enumerate(lines):
        if ('Multipole matrices IBuc=  518 IX=    3' in line):
            for i in range (NMOs):
                j = 0
                while (j <= i): 
                    elements = lines[n+i+2].split()
                    k = j + 1
                    dipZ[i][j] = float(elements[k])
                    if ( i != j):
                        dipZ[j][i] = dipZ[i][j]
                    j=j+1
else:
    loops = NMOs / 5 + 1
#    print 'loops = ', loops
    last = NMOs % 5
#    print 'last = ', last
    shift = 0
    for (n,line) in enumerate(lines):
        if ('Multipole matrices IBuc=  518 IX=    3' in line):
            for k in range(loops):
                irange = NMOs - k*5
                for i in range(irange):
                    if k == loops:
                        if i < last:
                            end = i + 1
                        else:
                            end = last
                    else:
                        if i <= 4:
                            end = i + 1
                        else:
                            end = 5
#                    print 'k, i, end', k, i, end
                    for j in range(end):
                        elements=lines[n+shift+k*1+i+2].split()
                        s = k*5 + j
                        m = i + k*5
#                        print k, s, m, elements
                        dipZ[m][s] = float(elements[j+1])
                        if i != j:
                            dipZ[s][m] = dipZ[m][s]
                shift = shift + irange

#print 'dipole moment along z-direction: \n', dipZ

# debug print statement
if debug == 1:
  print 'Dipole X'
  for i in range(NMOs):
      for j in range(NMOs):
          print "%d  %d  %f" % (i,j, dipX[i][j])

  print 'Dipole Y'
  for i in range(NMOs):
      for j in range(NMOs):
         print "%d  %d  %f" % (i,j, dipY[i][j])

  print 'Dipole Z'
  for i in range(NMOs):
      for j in range(NMOs):
         print "%d  %d  %f" % (i,j, dipZ[i][j])

#initialize MO matrix natoms x natoms
MO = numpy.zeros([NMOs,NMOs], numpy.float64)

# pulling out the line number of the last instance of "FINAL COEFFICIENT MATRIX" in the log file to extract the MO coefficients - KR

#print "entering dummy phase..."

dummy1 = []

for (n, line) in enumerate(lines):
    if ('FINAL COEFFICIENT MATRIX' in line):
        dummy1.append(n)
        #print "n ", n, "; line ", line

nn1 = dummy1[-1]

#print "nn1:",nn1

#"it is done."

n = nn1		# dummy variable for the line number

# again less than 10 MOs is easy
if (NMOs < 11):

# need to get MO coefs- in Gaussian stored as rows
# this is the part that requires changing the log file,     << can be
# only one instance of 'FINAL COEFFICIENT MATRIX' can exist << ignored - KR
#    for (n, line) in enumerate(lines):
#        if ('FINAL COEFFICIENT MATRIX' in line):
#    print "n ", n
    for i in range(NMOs):
        for j in range(NMOs):
            elements = lines[n+i*2+2].split()                     
            MO[j][i] = float(elements[j])
else:
    loops = NMOs / 10 + 1
    last = NMOs % 10
#    ii=0
#    for (n,line) in enumerate(lines):
#        if ('FINAL COEFFICIENT MATRIX' in line):
#    print "n ",n
    for i in range(NMOs):
        for j in range(loops):
            #ii += 1
            if j == loops-1:
                end = last
            else:
                end = 10
                for k in range(end):
                    s = k+10*j
                    elements=lines[n+i+2+j+loops*i].split()
                    #print "s",s,"i",i,"j",j,"k",k,"end",end,"line",(n+i+j+ii+2)
                    #print elements[k] 
                    MO[s][i] = elements[k]

# need to truncate MOs so only those in occupied space (NOCC) have values

#if (NOCC < NMOs):
#    for i in range(NOCC,NMOs):
#        for j in range(NOCC,NMOs):
#            MO[i][j] = 0.0
#            MO[j][i] = 0.0 
#            dipX[i][j] = 0.0
#            dipX[j][i] = 0.0
#            dipY[i][j] = 0.0
#            dipY[j][i] = 0.0
#            dipZ[i][j] = 0.0
#            dipZ[j][i] = 0.0


# debug print statement
if debug == 1:
  print 'MO Coefs'
  for i in range(NMOs):
      for j in range(NMOs):
          print i, j, MO[i][j]

# initialize new matrix for dipole in MO basis
dipXMO = numpy.zeros([NMOs,NMOs], numpy.float64)
dipYMO = numpy.zeros([NMOs,NMOs], numpy.float64)
dipZMO = numpy.zeros([NMOs,NMOs], numpy.float64)

# rotate dipole matrix into MO basis with C'muC
dipXMO = MO.T.dot(dipX).dot(MO)
dipYMO = MO.T.dot(dipY).dot(MO)
dipZMO = MO.T.dot(dipZ).dot(MO)

if (debug == 2): 
    print 'MO Dipole x,y,z'
    print dipXMO
    print dipYMO
    print dipZMO

#print 'dipole moment along z-direction in MO basis: \n', dipZMO, '\n', dipXMO.shape

#initialize matrix for CI vectors
CI_vect = numpy.zeros([NStates,NStates], numpy.float64)
CI_E = numpy.zeros([NStates], numpy.float64)

# less than 6 MOs as a first go
if (NStates < 6):
    for (n, line) in enumerate(lines):
# note that AND is now spelled ANE due to sed call at beginning
        if ('FINAL EIGENVALUES ANE EIGENVECTORS' in line):
            for i in range(NStates):
                for j in range(NStates):
                    elements = lines[n+i*2+4].split()
                    k = j+2
                    #print elements
                    CI_E[i] = elements[1]
                    CI_vect[i,j] = elements[k]
else:
    loops = NStates / 5 + 1
    last = NStates % 5
    for (n,line) in enumerate(lines):
      if ('FINAL EIGENVALUES ANE EIGENVECTORS' in line):
        for i in range(NStates):
            #print "i = ",i
            for j in range(loops):
                #print "j = ",j
                if j == 0:
                    elements = lines[n+4+i*(loops+1)].split()
                    #print "n",n,"i",i,"j",j,"loops",loops,"skipped lines",(4+i*(loops))
                    #print elements
                    CI_E[i] = elements[1]
                    for k in range(5):
                         #print "k = ",k  
#                         print 'i, j, k, elements', i,j,k,elements
                         s = k + 2
                         #print "elements","\n",elements,"\n","s","\n",s 
                         CI_vect[i,k] = elements[s]
                else: 
                    if j == loops - 1:
                        end = last
                    else:
                        end = 5
                    for k in range(end):
                        #print "k = ",k
                        elements = lines[n+4+i*(loops+1)+j].split()
#                        print 'i,j,k,s,elements', i,j,k,s,elements
                        #print "n",n,"i",i,"j",j,"loops",loops,"skipped lines",(4+i*(loops)+j),"end",end
                        #print elements
                        s = k+5*j
                        CI_vect[i,s] = elements[k]


# debug print statement
if debug == 2:
  print 'CI Eigenvalues'
  for i in range(NStates):
     print i, CI_E[i]
  print 'CI Eigenvectors'
  for i in range(NStates):
     for j in range(NStates):
         print i,j,CI_vect[i,j]

# initialize arrays to read and analyze MO configuration for each configuration
tmp = numpy.zeros((NStates), dtype='object')
config = numpy.zeros((NStates,NCAS), dtype='object')
multi = numpy.zeros((NStates,NStates), dtype='int')
alpha = numpy.zeros((NStates), dtype='int')

for (n, line) in enumerate(lines):
    if ('BOTTOM WEIGHT' in line):
        for i in range(NStates):
            elements = lines[n + i + 1].split()
            tmp = list(elements[4])
            for j in range(NCAS):
               if tmp[j] == '1':
                   tmp[j] = '2'
               if tmp[j] == 'a' or tmp[j] == 'b':
                   tmp[j] = '1'
               config[i][j] = tmp[j]
               if (tmp[j] == 'a'):
                  alpha[i] = alpha[i] + 1

#print 'config = ', config 

coupling = numpy.zeros([NStates,NStates], dtype='int')
if (debug == 1):
  print 'config = ', 
  for i in range(NStates):
    for j in range(NCAS):
       print i,j,config[i][j]

ct_trans = 0
tmp_conf = numpy.zeros([NCAS],dtype='object')
tmp_conf2 = numpy.zeros([NCAS],dtype='object')
tmp_conf3 = numpy.zeros([NCAS],dtype='object')
keep_conf = numpy.zeros([NCAS],dtype='object')
orig_conf = numpy.zeros([NCAS],dtype='object')
follow = numpy.zeros([NStates,NCAS], dtype='int')
leng_list = NStates*3  
orb_info = numpy.zeros([NCAS,NCAS,leng_list], dtype='int')
ct_orb = numpy.zeros([NCAS,NCAS], dtype='int')

# need to cycle through the configurations and electron occupations of the MOs
# to figure out which CI coefficients should be multiplied by which MOs and
# what the prefactor needs to be
#
# do diagonal terms first, MO#1 = MO#2 

for i in range(NStates):
    for j in range(NCAS):
        if config[i][j] == '1':
            sp = ct_orb[j][j]
#            print 'i, j, sp = ', i, j, sp
            orb_info[j][j][sp] = i
            orb_info[j][j][sp+1] = i 
            orb_info[j][j][sp+2] = 1
	    ct_orb[j][j] = sp + 3
	elif config[i][j] == '2':     
            sp = ct_orb[j][j]
#            print 'i, j, sp = ', i, j, sp
            orb_info[j][j][sp] = i
            orb_info[j][j][sp+1] = i 
            orb_info[j][j][sp+2] = 2
            ct_orb[j][j] = sp + 3


# get off-diagonal MOs and CI coefficients
for i in range(NStates):
    for k in range(NCAS):
        orig_conf[k] = config[i][k]
    for k in range(NCAS):
        for p in range(NCAS):
            keep_conf[p] = orig_conf[p] 
        if config[i][k] != '0':
            if config[i][k] == '2':
		keep_conf[k] = '1'
#                print 'keep_conf = ', keep_conf
                for j in range(k+1,NCAS):
                    for s in range(NCAS):
                        tmp_conf[s] = keep_conf[s]
#                    print 'i, k, j = ', i, k, j, tmp_conf, keep_conf, tmp_conf[j] 
                    if tmp_conf[j] == '1':
                        tmp_conf[j] = '2'
                        ct = -1
                        for n in range(NStates):
                            check = 1
                            for m in range(NCAS):
                               tmp_conf2[m] = config[i][m]
                               if tmp_conf[m] != config[n][m]:
                                  check = 0
                            if check == 1: 
                               ct = j
#                               print 'coupling CIs ', i, n, 'with MOs ', k, ct
#                               print tmp_conf2
#                               print tmp_conf
                               sp = ct_orb[k][ct]
#                               print 'sp = ', sp
                               orb_info[k][ct][sp] = i
                               orb_info[ct][k][sp] = i
                               orb_info[k][ct][sp+1] = n
                               orb_info[ct][k][sp+1] = n
                               if config[i][k] == '2' or config[i][ct] == '2' or config[n][k] == '2' or config[n][ct] == '2':
                                   ifactor = 5  # note that 5 = square root of two later in the code, integer placeholder here 
                               else:
                                   ifactor = 1 
                               orb_info[ct][k][sp+2] = ifactor
                               orb_info[k][ct][sp+2] = ifactor
                               ct_orb[k][ct] = sp + 3
                               ct_orb[ct][k] = sp + 3
                               ct_trans = ct_trans + 1
                    elif tmp_conf[j] == '0':
                        tmp_conf[j] = '1'
                        ct = -1
#                        print 'in loop i, k, j = ', i, k, j, tmp_conf  
                        for n in range(NStates):
                            check = 1
                            for m in range(NCAS):
                               tmp_conf2[m] = config[i][m]
                               tmp_conf3[m] = config[n][m]
                               if tmp_conf[m] != config[n][m]:
                                  check = 0
#                            print 'weird ', tmp_conf, tmp_conf3, check
                            if check == 1:                              
                               ct = j
#                               print 'coupling CIs ', i, n, 'with MOs ', k, ct
#                               print tmp_conf2
#                               print tmp_conf
                               sp = ct_orb[k][ct]
#                               print 'sp = ', sp 
                               orb_info[k][ct][sp] = i
                               orb_info[ct][k][sp] = i
                               orb_info[k][ct][sp+1] = n
                               orb_info[ct][k][sp+1] = n
                               if config[i][k] == '2' or config[i][ct] == '2' or config[n][k] == '2' or config[n][ct] == '2':
                                   ifactor = 5  # note that 5 = square root of two later in the code, integer placeholder here 
                               else:
                                   ifactor = 1
                               orb_info[ct][k][sp+2] = ifactor
                               orb_info[k][ct][sp+2] = ifactor
                               ct_orb[k][ct] = sp + 3
                               ct_orb[ct][k] = sp + 3
                               ct_trans = ct_trans + 1
            elif config[i][k] == '1':
                keep_conf[k] = '0'
#                print 'keep_conf = ', keep_conf
                for j in range(k+1,NCAS):
                    for s in range(NCAS):
                        tmp_conf[s] = keep_conf[s]
#                    print 'i, j, k = ', i, j, k, tmp_conf, keep_conf
                    if tmp_conf[j] == '1':
                        tmp_conf[j] = '2'
                        ct = -1
                        for n in range(NStates):
                            check = 1
                            for m in range(NCAS):
                               tmp_conf2[m] = config[i][m]
                               if tmp_conf[m] != config[n][m]:
                                  check = 0
                            if check == 1:
                               ct = j
#                               print 'coupling CIs ', i, n, 'with MOs ', k, ct
#                               print tmp_conf2
#                               print tmp_conf
                               sp = ct_orb[k][ct]
#                               print 'sp = ', sp
                               orb_info[k][ct][sp] = i
                               orb_info[ct][k][sp] = i
                               orb_info[k][ct][sp+1] = n
                               orb_info[ct][k][sp+1] = n
                               if config[i][k] == '2' or config[i][ct] == '2' or config[n][k] == '2' or config[n][ct] == '2':
                                   ifactor = 5  # note that 5 = square root of two later in the code, integer placeholder here 
                               else:
                                   ifactor = 1
                               orb_info[ct][k][sp+2] = ifactor
                               orb_info[k][ct][sp+2] = ifactor
                               ct_orb[k][ct] = sp + 3
                               ct_orb[ct][k] = sp + 3
                               ct_trans = ct_trans + 1
                    elif tmp_conf[j] == '0':
                        tmp_conf[j] = '1'
                        ct = -1
                        for n in range(NStates):
                            check = 1
                            for m in range(NCAS):
                               tmp_conf2[m] = config[i][m]
                               if tmp_conf[m] != config[n][m]:
                                  check = 0
                            if check == 1:
                               ct = j
#                               print 'coupling CIs ', i, n, 'with MOs ', k, ct
#                               print tmp_conf2
#                               print tmp_conf
#                               print 'sp = ', sp 
                               sp = ct_orb[k][ct]
                               orb_info[k][ct][sp] = i
                               orb_info[ct][k][sp] = i
                               orb_info[k][ct][sp+1] = n
                               orb_info[ct][k][sp+1] = n
                               if config[i][k] == '2' or config[i][ct] == '2' or config[n][k] == '2' or config[n][ct] == '2':
                                   ifactor = 5  # note that 5 = square root of two later in the code, integer placeholder here 
                               else:
                                   ifactor = 1
                               orb_info[ct][k][sp+2] = ifactor
                               orb_info[k][ct][sp+2] = ifactor
                               ct_orb[k][ct] = sp + 3
                               ct_orb[ct][k] = sp + 3
                               ct_trans = ct_trans + 1


#print 'orb_info: \n', orb_info
#print 'ct_orb: \n', ct_orb
#print 'orb_info'
#for i in range(NCAS):
#    for j in range(NCAS):
#        print i, j, orb_info[i][j][:]

if debug == 1:
    print 'ct_orb' 
    print ct_orb

    print 'orb_info'
    for i in range(NCAS):
        for j in range(NCAS):
            print i, j, orb_info[i][j][:]
 
num_trans = 0
for i in range(NStates):
  num_trans = num_trans + i

#print 'num_trans = ', num_trans 
mem_trans = num_trans*2*NCAS*NCAS
trans_den = numpy.zeros([mem_trans,NMOs+1,NMOs+1], numpy.float64) 
index = numpy.zeros([4,mem_trans])


# now we calculate the transition dipole moment bfh
# open output files 
tdm = open('transition_dipole.txt', 'w')
db = open('debug.out', 'w')
tdip = numpy.zeros([NStates,NStates,3], numpy.float64)

for n in range(NStates):
    for k in range(n,NStates):
#for n in range(1):
#    for k in range(1,2):
#     if n == 0 and k == 3:
        tsumij = 0.0
        tsumji = 0.0
        left = n 
        right = k
        for i in range(NCAS):
            for j in range(i,NCAS):
              for s in range(NStates):
                ct = s*3 
                ci1 = orb_info[i][j][ct]
                ci2 = orb_info[i][j][ct+1]
                fact = orb_info[i][j][ct+2]
#                print 'n, k, i, j, ci1, ci2, fact = ', n, k, i+NFRZ, j+NFRZ, ci1, ci2, fact
#                print 'ci coefs = ', CI_vect[n,ci1], CI_vect[k,ci2], CI_vect[k,ci1], CI_vect[n,ci2]
                if fact == 1:
                    factor = 1.0
                elif fact == 2:
                    factor = 2.0
                elif fact == 5:
                    factor = rt2
                elif fact == 0:
                    factor = 0
                else:
                    print 'problem with factors'                
#               adjust MOs to CASSCF space
                r = i + NFRZ
                t = j + NFRZ
                if i == j:
                    tsumij = CI_vect[left,ci1]*CI_vect[right,ci2]
                    tsumji = 0.0
                    tdip[n][k][0] = tdip[n][k][0] + (tsumij*dipXMO[r][t] + tsumji*dipXMO[t][r])*factor
                    tdip[n][k][1] = tdip[n][k][1] + (tsumij*dipYMO[r][t] + tsumji*dipYMO[t][r])*factor
                    tdip[n][k][2] = tdip[n][k][2] + (tsumij*dipZMO[r][t] + tsumji*dipZMO[t][r])*factor
#                    print 'tsumij, MO dipz, tdipz = ', tsumij, dipZMO[r][t], tdip[n][k][2]
                    if (left == 2 and left == right and ci1 == 1 and ci2 == ci1):
                        print '\n \n i, j = ', r,', ',t
                        print 'c_01(S2)**2 = ', tsumij
                        print 'normalization factor = ', factor
                        print 'MO dipole moment matrix element (ij, ji)= ',dipZMO[r,t],', ',dipZMO[t,r] 
                        print '\n \n '
                else:
                    tsumij = CI_vect[left,ci1]*CI_vect[right,ci2] 
                    tsumji = CI_vect[right,ci1]*CI_vect[left,ci2]
                    tdip[n][k][0] = tdip[n][k][0] + (tsumij*dipXMO[t][r] + tsumji*dipXMO[t][r])*factor
                    tdip[n][k][1] = tdip[n][k][1] + (tsumij*dipYMO[t][r] + tsumji*dipYMO[t][r])*factor
                    tdip[n][k][2] = tdip[n][k][2] + (tsumij*dipZMO[t][r] + tsumji*dipZMO[t][r])*factor
#                    print 'tsumij, tsumji, MO dipz, tdipz = ', tsumij, tsumji, dipZMO[r][t], tdip[n][k][2]
                    if (left == 2 and right == left and ci1 == 1 and ci2 == 2):
                        print '\n \n i, j = ', r,', ', t
                        print 'c_01(S2)*c_0011(S2) = ', tsumij
                        print 'c_0011(S2)*c_01(S2) = ', tsumji
                        print 'normalization factor = ', factor
                        print 'MO dipole moment matrix element (ij, ji)= ',dipZMO[r,t],', ',dipZMO[t,r]
                        print '\n \n '

print 'Transition Dipole Moment X \n', tdip[:,:,0]
print 'Transition Dipole Moment Y \n', tdip[:,:,1]
print 'Transition Dipole Moment Z \n', tdip[:,:,2]
 
#for n in range(NStates):
#    for k in range(n+1,NStates):
#        print n, k, tdip[n][k][2]*(-1.0) 

# have now calculated the Transition dipole moments. time to get atom coordinates and dump data
# to file so TD_CIS.py can read it in.

tdm = open('tdm_tdcis.txt', 'w')
tdm.write('This file contains the information one needs to run a time-dependent CASSCF simulation \n')
tdm.write('An initial Gaussian calculation was used to calculate the CI eigenvalues and eigenvectors \n')
tdm.write('as well as the CI state configurations. The program calc_tdm.py was then used to \n')
tdm.write('calculate the transition dipole moments between the CI states as well as the project of \n ')
tdm.write('the CI states onto the molecular orbitals. This file was generated by calc_tdm.py as \n')
tdm.write('input to TD_CIS.py. -BFH 9 May 2014 \n')
tdm.write(' \n ')
tdm.write(' \n')
tdm.write('First, we set up the atomic coordinates \n')
tdm.write(' \n')
tdm.write('Standard orientation: \n ')
tdm.write('------------------------------------------------------------------- \n')
tdm.write('Center   Atomic     Atomic            Coordinates (Angstroms)  \n')
tdm.write('Number   Number      Type            X            Y          Z   \n')
tdm.write('------------------------------------------------------------------- \n')


coords = numpy.zeros([NAtoms,3], numpy.float64)
atom_info = numpy.zeros([NAtoms,3], numpy.float64)

for (n, line) in enumerate(lines):
    if ('Standard basis:' in line):
        for i in range(NAtoms):
            elements = lines[n + i + 6].split()
            atom_info[i][0] = float(elements[0])
            atom_info[i][1] = float(elements[1])
            atom_info[i][2] = float(elements[2])
            coords[i][0] = float(elements[3])
            coords[i][1] = float(elements[4])
            coords[i][2] = float(elements[5])
            tdm.write('%-10s %-10s %-10s %-10s %-10s %-10s' % (atom_info[i][0], atom_info[i][1], atom_info[i][2], coords[i][0], coords[i][1], coords[i][2]) + "\n") 
#            tdm.write('  ' + elements[0] + str(atom_info[i][0]) + '   ' + str(atom_info[i][1]) + '   ' + str(atom_info[i][2]) + '    '  + str(coords[i][0]) + '   ' + str(coords[i][1]) + '    ' + str(coords[i][2]) + "\n")  

# input SCF energy- assume energy of first NState
tdm.write(' \n ')
tdm.write('SCF Energy \n')
tdm.write('SCF Done: E(RHF) = ' + str(CI_E[0]) + "\n")
tdm.write(' \n ')

# Excited State energies
tdm.write('Excited state energies \n')
tdm.write(' \n')
for i in range(1,NStates):
    prt_E = (CI_E[i] - CI_E[0])*27.211396 
    tdm.write('Excited State   ' + str(i) + '     Singlet-BU     ' + str(prt_E) + ' eV  \n ')

# Ground state to excited state transition dipole moments
tdm.write(' \n')
tdm.write('\n')
tdm.write('Ground state to excited state transition dipole moments \n')
tdm.write(' \n ')
tdm.write('Ground to excited state transition electric dipole moments \n ')
tdm.write('     state      X       Y     Z     Dip. S.      Osc.     \n ' )
for i in range(1,NStates):
    tdm.write('%-10d %-10.6f %-10.6f %-10.6f %-10.6f %-10.6f' % (i, tdip[0][i][0]*(-1.0), tdip[0][i][1]*(-1.0), tdip[0][i][2]*(-1.0), 0, 0 ) + "\n")

# Excited to excited state transition dipole moments
tdm.write(' \n')
tdm.write('\n')
tdm.write('Excited state to excited state transition dipole moments \n')
tdm.write(' \n ')
tdm.write('Excited to excited state transition electric dipole moments \n ')
tdm.write('  state I    state J      X       Y     Z     Dip. S.      Osc.     \n ' )
for i in range(2,NStates):
    for j in range(1,i):
        tdm.write('%10d %10d %-10.6f %-10.6f %-10.6f %-10.6f %-10.6f' % (i, j, tdip[j][i][0]*(-1.0), tdip[j][i][1]*(-1.0), tdip[j][i][2]*(-1.0), 0, 0 ) + "\n")
 
tdm.write(' \n ')
tdm.write(' \n ')
 
# Finally, need to write out information on CI configs and MOs occupations
tdm.write('Output CI Vectors \n')
tdm.write('CI State Number, CI coefs 1...N \n')
for i in range(NStates):
    tdm.write('%12d' % i )
    for j in range(NStates):
        tdm.write('%14.8f' % CI_vect[i][j])
    tdm.write('\n')
tdm.write('\n')
tdm.write('\n')
    

#note, want state # then number for each orbital. something like
#   state       MO 1     MO 2    MO 3     MO 4
#     1         1.95     1.95    0.05     0.05

mo_weight = numpy.zeros([NCAS], dtype=int)
tdm.write('Orbital electron configuration weighting for CI States \n')
tdm.write('State     Orbital Weights 1...N  \n')
for i in range(NStates):
    for k in range(NCAS):
            if config[i][k] == '0':
                weight = 0
            elif config[i][k] == '1':
                weight = 1
            elif config[i][k] == '2':
                weight = 2
            mo_weight[k] = weight
    tdm.write('  ' + str(i))
    for k in range(NCAS):
        if k == NCAS-1:
            tdm.write('%10d' % (mo_weight[k]) + "\n")
        else:
            tdm.write('%10d' % (mo_weight[k]))
        mo_weight[k] = 0.0 

tdm.write("\n")


#writing the entire dipole moment matrix to 'dipMO_tdcis.txt' - KR

tdm.write('This file contains the dipole moment matrix in MO basis.  \n')
tdm.write('The matrix can be used to compute the total instantaneous \n')
tdm.write('dipole moment of the system. - KR                         \n')
tdm.write('\n')
tdm.write('\n')
tdm.write('\n')

#print 'TDM-x matrix (au): \n', tdip[:,:,0]
#print 'TDM-y matrix (au): \n', tdip[:,:,1]
#print 'TDM-z matrix (au): \n', tdip[:,:,2]

tdm.write('mu_CIx \n \n')
for i in range(NStates):
    for j in range(NStates):
        tdm.write('%14.8f' % tdip[i][j][0])
    tdm.write('\n')

tdm.write('\n \n \n')

tdm.write('mu_CIy \n \n')
for i in range(NStates):
    for j in range(NStates):
        tdm.write('%14.8f' % tdip[i][j][1])
    tdm.write('\n')

tdm.write('\n \n \n')

tdm.write('mu_CIz \n \n \n')
for i in range(NStates):
    for j in range(NStates):
        tdm.write('%14.8f' % tdip[i][j][2])
    tdm.write('\n')

tdm.write('\n')

tdm.close()










import sys
import numpy, math

#
# calculating dipole moment matrix (CI basis) elements for CI states, which will involve
#     expanding the CI basis matrix elements in terms of MO basis matrix elements
#
# NOTE: closed shell cases only
#
# - Karnamohit

sum_dip = 1

n_half = numpy.zeros((NStates), dtype='object')
n_ex = numpy.zeros((NStates,NStates), dtype='object')
e1_indx = numpy.zeros((NStates,NStates,2), dtype='object')
e1_fact = numpy.zeros((NStates), dtype='float')
e1_fact_st = numpy.zeros((NStates,NStates), dtype='object')

# calculating the degree of excitation (WRT other configurations) and the normalization
#     factors (for states with two half-filled orbitals)

for A in range(NStates):
    for B in range(A+1,NStates):
        for i in range(NCAS):
            if (config[A][i] != config[B][i]):
                n_ex[A,B] += abs(float(config[A][i]) - float(config[B][i]))
        n_ex[A,B] /= 2
        n_ex[B,A] = n_ex[A,B]

for A in range(NStates):
    for i in range(NCAS):
        if (float(config[A][i]) == 1):
            n_half[A] += 1
    N = int(n_half[A])
    k = int(n_half[A]/2)
#    if (N == 0):
#        NCk = 1
#    else:
#        NCk = float(numpy.math.factorial(N))/float((numpy.math.factorial(N-k)*numpy.math.factorial(k)))
#    print A,n_half[A],n_ex[0,A]
#    e1_fact[A] = 1./numpy.sqrt(NCk)
    e1_fact[A] = 1./numpy.sqrt(numpy.math.factorial(N))

for A in range(NStates):
    for B in range(A+1,NStates):
        if n_ex[A,B] == 1:
            e1_fact_st[A,B] = e1_fact[A]*e1_fact[B]
            e1_fact_st[B,A] = e1_fact_st[A,B]

print 'config matrix: \n',config
print 'n_ex matrix: \n',n_ex,'\n'
print 'determinant normalization factor matrix: \n',e1_fact,'\n'
print 'cross-term factor matrix: \n',e1_fact_st,'\n'

# keeping a track of mismatched indices (pairs of indices involved in single-excitations)
#
# NOTE: only MO dip. mom. integrals with #(mismatches) <= 1 will survive

for A in range(NStates):
    for B in range(A+1,NStates):
        j = 0
        if (n_ex[A,B] == 1):
            for i in range(NCAS):
                if ((j < 2) and (config[A][i] != config[B][i])):
                    e1_indx[A][B][j] = i
                    e1_indx[B][A][j] = e1_indx[A][B][j]
                    j += 1

#print 'indices for determinants with mismatch of 1 index: \n'
#for A in range(NStates):
#    for B in range(NStates):
#        print 'states ',A,' ',B,':\n'
#        print 'mismatched indices: ',e1_indx[A][B][0],' ',e1_indx[A][B][1],'\n'
#        print 'degree of excitation = ',n_ex[A,B]

# contribution from the frozen orbitals and expectation values of the dipole moment 
#     operator WRT states corresponding to the same configuration (needs to be multiplied 
#     by proper coefficients in the actual sum done to get the CI dip. mom. matrix
#     elements): Slater-Condon rules 

dip_frz = numpy.zeros((NStates,3), dtype='object')

for A in range(NStates):
    if (NFRZ > 0):
        for i in range(NFRZ):
            dip_frz[A][0] += 2*dipXMO[i][i]
            dip_frz[A][1] += 2*dipYMO[i][i]
            dip_frz[A][2] += 2*dipZMO[i][i]
    for i in range(NCAS):
        dip_frz[A,0] += float(config[A][i])*dipXMO[i+NFRZ][i+NFRZ]
        dip_frz[A,1] += float(config[A][i])*dipYMO[i+NFRZ][i+NFRZ]
        dip_frz[A,2] += float(config[A][i])*dipZMO[i+NFRZ][i+NFRZ]

if (sum_dip == 0):
    exit()

# summing over all the MO basis matrix elements/MO integrals to get the CI basis dipole
#     moment matrix elements

print 'CI vector: \n',CI_vect[:,:]
print 'MO basis dipole moment matrix (Z): \n',dipZMO[:,:]

ci_dip = numpy.zeros((NStates,NStates,3), dtype='float')

for A in range(NStates):
    # diagonal elements of CI dip. mom. matrix: ci_dip[A][A]
    for C in range(NStates):
        # contribution from pair of zero mismatched-index determinants
        ci_dip[A][A][:] += (CI_vect[A,C]**2)*dip_frz[C][:]
        if (A == 2):
            print 'determinant: ',C
            print 'CI coefficient: ',CI_vect[A,C]
            print 'dipole moment matrix element (MO basis, diagonal contribution): ',dip_frz[C,2]
            print 'normalization factor: ',e1_fact[C]
        for Ne in range(NELECT):
            for D in range(C+1,NStates):
                # contribution from pair of one mismatched-index determinants
                if (int(n_ex[C,D]) == 1):
                    i = e1_indx[C][D][0] + NFRZ
                    j = e1_indx[C][D][1] + NFRZ
                    if (A == 2):
#                    print 'states: ',A,' ',A
                        print 'determinants: ',C,' ',D,' with mismatched indices: ',i,' ',j
                        print 'degree of excitation: ',n_ex[C,D]
                        print 'CI coefficients: ',CI_vect[A,C],' ',CI_vect[A,D]
                        print 'dipole moment matrix element (MO basis): ',dipZMO[i,j]
                        print 'normalization factor: ',e1_fact[C]*e1_fact[D]*2
                    ci_dip[A][A][0] += 2*CI_vect[A,C]*CI_vect[A,D]*e1_fact_st[C,D]*dipXMO[i][j]
                    ci_dip[A][A][1] += 2*CI_vect[A,C]*CI_vect[A,D]*e1_fact_st[C,D]*dipYMO[i][j]
                    ci_dip[A][A][2] += 2*CI_vect[A,C]*CI_vect[A,D]*e1_fact_st[C,D]*dipZMO[i][j]
    
    for B in range(A+1,NStates):
        # off-diagonal elements of CI dip. mom. matrix: ci_dip[A][B]
        for C in range(NStates):    
            # contribution from pair of zero mismatch-index determinants
            ci_dip[A][B][:] += CI_vect[A,C]*CI_vect[B,C]*dip_frz[C][:]
            for Ne in range(NELECT):
                for D in range(C+1,NStates):
                    # contribution from pair of one mismatched-index determinants
                    if (int(n_ex[C,D]) == 1):
                        i = e1_indx[C][D][0] + NFRZ
                        j = e1_indx[C][D][1] + NFRZ
#                    print 'states: ',A,' ',B
#                    print 'determinants: ',C,' ',D,' with mismatched indices: ',i,' ',j
#                    print 'degree of excitation: ',n_ex[C,D]
#                        ci_dip[A][B][0] += 2*(CI_vect[A,C]*CI_vect[B,D])*e1_fact_st[C,D]*dipXMO[i][j]
#                        ci_dip[A][B][1] += 2*(CI_vect[A,C]*CI_vect[B,D])*e1_fact_st[C,D]*dipYMO[i][j]
#                        ci_dip[A][B][2] += 2*(CI_vect[A,C]*CI_vect[B,D])*e1_fact_st[C,D]*dipZMO[i][j]
                        ci_dip[A][B][0] += (CI_vect[B,C]*CI_vect[A,D]+CI_vect[A,C]*CI_vect[B,D])*e1_fact_st[C,D]*dipXMO[i][j]
                        ci_dip[A][B][1] += (CI_vect[B,C]*CI_vect[A,D]+CI_vect[A,C]*CI_vect[B,D])*e1_fact_st[C,D]*dipYMO[i][j]
                        ci_dip[A][B][2] += (CI_vect[B,C]*CI_vect[A,D]+CI_vect[A,C]*CI_vect[B,D])*e1_fact_st[C,D]*dipZMO[i][j]
#        ci_dip[B,A,:]=ci_dip[A,B,:]
    if (A == 2):
        print 'S2 CI dip. mom.: ',ci_dip[A,A,2]

S0_dip = 2*(CI_vect[0,0]**2)*dipZMO[0,0]+(CI_vect[0,1]**2)*(dipZMO[0,0]+dipZMO[1,1])+2*(CI_vect[0,2]**2)*dipZMO[1,1]
S0_dip += 2*(CI_vect[0,1]*CI_vect[0,0]+CI_vect[0,1]*CI_vect[0,2])*dipZMO[1,0]
print 'S0 state dipole moment = ',S0_dip,' au'

S2_dip = 2*(CI_vect[2,0]**2)*dipZMO[0,0]+(CI_vect[2,1]**2)*(dipZMO[0,0]+dipZMO[1,1])+2*(CI_vect[2,2]**2)*dipZMO[1,1]
S2_dip += 2*numpy.sqrt(2)*(CI_vect[2,1]*CI_vect[2,0]+CI_vect[2,1]*CI_vect[2,2])*dipZMO[1,0]
print 'S2 state dipole moment = ',S2_dip,' au'

print '\n \n i, j = ', r,', ',t
print 'c_01(S2)**2 = ', CI_vect[2,1]**2
print 'normalization factor = ', e1_fact[1]
print 'MO dipole moment matrix element (ij, ji)= ',dipZMO[r,t],', ',dipZMO[t,r]
print '\n \n '

#print 'MO basis dip. mom. matrix (x): \n',dipXMO[:,:],'\n'
#print 'MO basis dip. mom. matrix (y): \n',dipYMO[:,:],'\n'
#print 'MO basis dip. mom. matrix (z): \n',dipZMO[:,:],'\n'

print 'CI basis dip. mom. matrix (x) (Brad): \n',tdip[:,:,0],'\n'
print 'CI basis dip. mom. matrix (y) (Brad): \n',tdip[:,:,1],'\n'
print 'CI basis dip. mom. matrix (z) (Brad): \n',tdip[:,:,2],'\n'

print 'CI basis dip. mom. matrix (x): \n',ci_dip[:,:,0],'\n'
print 'CI basis dip. mom. matrix (y): \n',ci_dip[:,:,1],'\n'
print 'CI basis dip. mom. matrix (z): \n',ci_dip[:,:,2],'\n'

# End of CI dip. mom. matrix calculation

print 'MO overlap matrix: \n',MO.T.dot(MO)
