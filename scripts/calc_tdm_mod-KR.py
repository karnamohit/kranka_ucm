#!/usr/bin/env python
"""
April 17, 2014
@author B.F. Habenicht (bfh)
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

- bfh
NOTE!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

^ This is not true any more. Code has been adjusted to save all the line-numbers with 'FINAL COEFFICIENT
   MATRIX' and pick out the one corresponding to the last occurrence of this string. - Karnamohit (KR).

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
'''
 Important note: if using "Gaussian DV (version i10+)", without specifying the multiplicity of calculated excited states, only states with
                 the same multiplicity as that of the ground state wlll be calculated. With "Gaussian DV (version i14+)", this is not true
                 anymore, and if starting with a singlet state, triplet excited states will also be calculated (at least for H2, HeH+).
                  - KR, 04/26/2019.
'''
cas = '23'
sN = 's5'
ext = 'log'
Logfile = 'casscf'+cas+'_'+sN+'_lih.'+ext   # Gaussian output file using "GDV i10+":
#p CASSCF(2,2,nroot=3)/sto-3g scf(tight) nosymm IOp(3/33=3,3/36=1,4/33=3,5/33=3,6/33=3,9/33=3)

gdv = "default"		# Gaussian DV version
#gdv = "i10+"

NMOs=6			# total number of MOs
NCAS=3			# size of complete active space (CAS)
NFRZ=1			# number of orbitals with frozen occupations
NOCC=4			# NCAS + NFRZ (size of CAS plus orbitals with frozen occupations)
NELECT=2		# number of electrons in CAS
NAtoms=2		# number of atoms in the molecule


'''
 calculating total number of configurations, NStates; bcoz i10+ is used, only N-tuplet configurations (where N is the input
 multiplicity  of the system) are calculated. If this LOG file was created thru i14+, the total number of states for singlet
 H2 or HeH+ would be 4, including one triplet state, instead of 3. - KR, 04/26/2019.
'''

'''
 Weyl's dimension formula for possible no., of determinants from CAS(NMO,Ne,S):

 # dets = ((2S+1)/(NMO+1))*[(NMO+1) choose ((Ne/2)-S)]*[(NMO+1) choose ((Ne/2)+S+1)],

  where S is (Ne_alpha - Ne_beta)/2 (the total spin quantum number of Ne electrons occupying NMO molecular orbitals)

 (refer to C. David Sherrill's notes on Configuration Interaction)
'''
SPIN = 0            	# 0 -> singlet determinants only; 1 -> doublet determinants only; and so on.
N1 = NCAS+1
k1 = NELECT/2 - SPIN
k2 = NELECT/2 + 1 + SPIN
weyl = ((2.*SPIN+1.)/(NCAS+1))*numpy.math.factorial(N1)/(numpy.math.factorial(N1-k1)*numpy.math.factorial(k1))
weyl *= numpy.math.factorial(N1)/(numpy.math.factorial(N1-k2)*numpy.math.factorial(k2))


if (gdv == "i10+"):
    NStates = 3		# for singlet GS HeH+ or H2, CAS(2,2)
elif (gdv == "i14+"):
    NStates = 4		# for singlet GS HeH+ or H2, CAS(2,2)
else:
    NStates = int(weyl)	# default no. of states is determined by Weyl's formula for no. of possible determinants with a given multiplicity.

#print('# configurations = ',NStates)


#nuc_dipz=0.72943428	# nuclear dipole moment along z-axis; nuc. dip. mom. (z) for HeH+ = 0.7294 a. u. - KR, 04/26/2019
nuc_dipx=0.0		# nuclear dipole moment along x-axis of the system in the input
nuc_dipy=0.0		# nuclear dipole moment along y-axis of the system in the input
nuc_dipz=2.89128097	# nuclear dipole moment along z-axis of the system in the input; nuc. dip. mom. (z) for LiH = 2.8913 a. u.



#===============
debug = 0      # my debug flag set = 1 or 2 to print stuff
#===============

## open log file for reading
#infile = open(Logfile, 'r')

outp = 'logfile.tmp'

# Copying the original output file to 'logfile.tmp', without having to do this manually as previously required - KR
shutil.copy(Logfile, outp)

out_file = open(outp, "w")

# Gaussian prints exponentials with the form 1.0D+00, which python doesn't recognize, so we 
#  go through the log file and change it to 1.0E+00 - bfh
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



# find dipoles and fill arrays
# Gaussian prints the dipole matrices as only lower triangle and the matrix is
#  5 columns across and NMOs rows down. If have less than 6 atoms, don't need to
#  worry about sections of matrix being printed later. we'll fix that part later. - bfh

# initialize array for molecular dipole matrices
dipX = numpy.zeros([NMOs,NMOs], numpy.float64)
dipY = numpy.zeros([NMOs,NMOs], numpy.float64)
dipZ = numpy.zeros([NMOs,NMOs], numpy.float64)

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
    last = NMOs % 5
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
                    for j in range(end):
                        elements=lines[n+shift+k*1+i+2].split()
                        s = k*5 + j
                        m = i + k*5
                        dipY[m][s] = float(elements[j+1])
                        if i != j:
                            dipY[s][m] = dipY[m][s]
                shift = shift + irange

#dipole Z
if (NMOs < 6):
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
    last = NMOs % 5
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
                    for j in range(end):
                        elements=lines[n+shift+k*1+i+2].split()
                        s = k*5 + j
                        m = i + k*5
                        dipZ[m][s] = float(elements[j+1])
                        if i != j:
                            dipZ[s][m] = dipZ[m][s]
                shift = shift + irange

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


# initialize MO matrix natoms x natoms
MO = numpy.zeros([NMOs,NMOs], numpy.float64)

# pulling out the line number of the last instance of "FINAL COEFFICIENT MATRIX"
#  in the log file to extract the MO coefficients - KR
dummy1 = []
for (n, line) in enumerate(lines):
    if ('FINAL COEFFICIENT MATRIX' in line):
        dummy1.append(n)
n = dummy1[-1]

# again less than 10 MOs is easy
if (NMOs < 11):
# need to get MO coefs- in Gaussian stored as rows
# this is the part that requires changing the log file,            	<< can be
#  only one instance of 'FINAL COEFFICIENT MATRIX' can exist - bfh	<< ignored - KR
#    for (n, line) in enumerate(lines):
#        if ('FINAL COEFFICIENT MATRIX' in line):
    for i in range(NMOs):
        for j in range(NMOs):
            elements = lines[n+i*2+2].split()                     
            MO[j][i] = float(elements[j])
else:
    loops = NMOs / 10 + 1
    last = NMOs % 10
#    for (n,line) in enumerate(lines):
#        if ('FINAL COEFFICIENT MATRIX' in line):
    for i in range(NMOs):
        for j in range(loops):
            if j == loops-1:
                end = last
            else:
                end = 10
                for k in range(end):
                    s = k+10*j
                    elements=lines[n+i+2+j+loops*i].split()
#                    print "s",s,"i",i,"j",j,"k",k,"end",end,"line",(n+i+j+ii+2)
#                    print elements[k] 
                    MO[s][i] = elements[k]


# the following block of code is unnecessary, hence commented out - KR
'''
# need to truncate MOs so only those in occupied space (NOCC) have values

if (NOCC < NMOs):
    for i in range(NOCC,NMOs):
        for j in range(NOCC,NMOs):
            MO[i][j] = 0.0
            MO[j][i] = 0.0 
            dipX[i][j] = 0.0
            dipX[j][i] = 0.0
            dipY[i][j] = 0.0
            dipY[j][i] = 0.0
            dipZ[i][j] = 0.0
            dipZ[j][i] = 0.0
'''


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

# rotate dipole matrix into MO basis with C'muC - bfh
dipXMO = MO.T.dot(dipX).dot(MO)
dipYMO = MO.T.dot(dipY).dot(MO)
dipZMO = MO.T.dot(dipZ).dot(MO)

if (debug == 2): 
    print 'MO Dipole x,y,z'
    print dipXMO
    print dipYMO
    print dipZMO


#initialize matrix for CI vectors
CI_vect = numpy.zeros([NStates,NStates], numpy.float64)
CI_E = numpy.zeros([NStates], numpy.float64)

# less than 6 MOs as a first go
if (NStates < 6):
    for (n, line) in enumerate(lines):
# note that AND is now spelled ANE due to sed call at beginning - bfh
        if ('FINAL EIGENVALUES ANE EIGENVECTORS' in line):
            for i in range(NStates):
                for j in range(NStates):
                    elements = lines[n+i*2+4].split()
                    k = j+2
                    CI_E[i] = elements[1]
                    CI_vect[i,j] = elements[k]
else:
    loops = NStates / 5 + 1
    last = NStates % 5
    for (n,line) in enumerate(lines):
      if ('FINAL EIGENVALUES ANE EIGENVECTORS' in line):
        for i in range(NStates):
            for j in range(loops):
                if j == 0:
# may need to edit the following depending on different matrix sizes:
#  interchange b/w "...i*(loops+1)..." and "...i*(loops)..." below (including in
#  'else:' block), whenever python encounters errors when reading 'elements'. - KR
                    elements = lines[n+4+i*(loops+1)].split()
                    CI_E[i] = elements[1]
                    for k in range(5):
                         s = k + 2 
                         CI_vect[i,k] = elements[s]
                else: 
                    if j == loops - 1:
                        end = last
                    else:
                        end = 5
                    for k in range(end):
                        elements = lines[n+4+i*(loops+1)+j].split()
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

# initialize arrays to read and analyze MO configuration for each configuration - bfh
tmp = numpy.zeros((NStates), dtype='object')
config = numpy.zeros((NStates,NCAS), dtype='object')
multi = numpy.zeros((NStates,NStates), dtype='int')
alpha = numpy.zeros((NStates), dtype='int')
config_str = numpy.zeros((NStates), dtype='object')

for (n, line) in enumerate(lines):
    if ('BOTTOM WEIGHT' in line):
        for i in range(NStates):
            elements = lines[n + i + 1].split()
            config_str[i] = elements[4]
            tmp = list(elements[4])
            for j in range(NCAS):
               if tmp[j] == '1':
                   tmp[j] = '2'
               if tmp[j] == 'a' or tmp[j] == 'b':
                   tmp[j] = '1'
               config[i][j] = tmp[j]
               if (tmp[j] == 'a'):
                  alpha[i] = alpha[i] + 1

#coupling = numpy.zeros([NStates,NStates], dtype='int')
if (debug == 1):
  print 'config = ', 
  for i in range(NStates):
    for j in range(NCAS):
       print i,j,config[i][j]

#
# NOTE: The following bit of code has been commented out primarily bcoz it is tough to
#        comprehend and modify. A separate, more straight-forward piece of code to
#        to calculate the full CI dipole moment matrix follows, which gives the correct
#        dipole moments for all CI states, including the frozen orbital contribution.
#  - KR
#

'''
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
#  to figure out which CI coefficients should be multiplied by which MOs and
#  what the prefactor needs to be. - bfh
#
# do diagonal terms first, MO#1 = MO#2 - bfh
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


# get off-diagonal MOs and CI coefficients - bfh
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
                                   ifactor = 5  # note that 5 = square root of two later in the code, integer placeholder here - bfh
                               else:
                                   ifactor = 1
                               orb_info[ct][k][sp+2] = ifactor
                               orb_info[k][ct][sp+2] = ifactor
                               ct_orb[k][ct] = sp + 3
                               ct_orb[ct][k] = sp + 3
                               ct_trans = ct_trans + 1

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

mem_trans = num_trans*2*NCAS*NCAS
trans_den = numpy.zeros([mem_trans,NMOs+1,NMOs+1], numpy.float64) 
index = numpy.zeros([4,mem_trans])


# now we calculate the transition dipole moment - bfh
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
#                    if (left == 2 and left == right and ci1 == 1 and ci2 == ci1):
#                        print '\n \n i, j = ', r,', ',t
#                        print 'c_01(S2)**2 = ', tsumij
#                        print 'normalization factor = ', factor
#                        print 'MO dipole moment matrix element (ij, ji)= ',dipZMO[r,t],', ',dipZMO[t,r] 
#                        print '\n \n '
                else:
                    tsumij = CI_vect[left,ci1]*CI_vect[right,ci2] 
                    tsumji = CI_vect[right,ci1]*CI_vect[left,ci2]
                    tdip[n][k][0] = tdip[n][k][0] + (tsumij*dipXMO[t][r] + tsumji*dipXMO[t][r])*factor
                    tdip[n][k][1] = tdip[n][k][1] + (tsumij*dipYMO[t][r] + tsumji*dipYMO[t][r])*factor
                    tdip[n][k][2] = tdip[n][k][2] + (tsumij*dipZMO[t][r] + tsumji*dipZMO[t][r])*factor
#                    print 'tsumij, tsumji, MO dipz, tdipz = ', tsumij, tsumji, dipZMO[r][t], tdip[n][k][2]
#                    if (left == 2 and right == left and ci1 == 1 and ci2 == 2):
#                        print '\n \n i, j = ', r,', ', t
#                        print 'c_01(S2)*c_0011(S2) = ', tsumij
#                        print 'c_0011(S2)*c_01(S2) = ', tsumji
#                        print 'normalization factor = ', factor
#                        print 'MO dipole moment matrix element (ij, ji)= ',dipZMO[r,t],', ',dipZMO[t,r]
#                        print '\n \n '

# symmetrizing the CI dip. mom. matrix. Unnecessary (propagation script reads in only the upper triangular matrix,
#  but done for the sake of completeness. - KR
for i in range(NStates):
    for j in range(i,NStates):
        tdip[j,i,:] = tdip[i,j,:]
'''

#
# NOTE: As promised, the straight-forward piece of code to calculate the full tdip array follows.
#        Hopefully, it is a bit easier to comprehend. - KR
#

# calculating dipole moment matrix (CI basis) elements for CI states, which will involve
#  expanding the CI basis matrix elements in terms of MO basis matrix elements
#
# NOTE: closed shell cases only
#
#  - KR

n_ex = numpy.zeros((NStates,NStates), dtype='int')
e1_indx = numpy.zeros((NStates,NStates,2), dtype='object')
e1_fact = numpy.zeros((NStates), dtype='float')
n_mismatch = numpy.zeros((NStates,NStates), dtype='object')

# calculating the degree of excitation of a configuration (WRT other configurations) - KR

for A in range(NStates):
    for B in range(A,NStates):
        for i in range(NCAS):
            if (config[A][i] != config[B][i]):
                n_ex[A,B] += int(abs(float(config[A][i]) - float(config[B][i])))
        n_ex[A,B] /= 2
        n_ex[B,A] = n_ex[A,B]

# normalization factors, if needed (Gaussian already includes them) - KR

for A in range(NStates):
    e1_fact[A] = 1.#/numpy.sqrt(numpy.math.factorial(n_ex[0,A]))

# keeping track of mismatched indices (pairs of indices involved in index-mismatches b/w relatively 
#  singly-substituted configurations)
#
# NOTE: Only MO dip. mom. integrals with #(mismatches) <= 1 will survive, bcoz of
#        Slater-Condon rules
#  - KR

for A in range(NStates):
    for B in range(A+1,NStates):
        j = 0
#        if (n_ex[A,B] == 1):
        for i in range(NCAS):
            if ((config[A][i] != config[B][i])):
                if ((j < 2)):
                    e1_indx[A][B][j] = i
                    e1_indx[B][A][j] = e1_indx[A][B][j]
                j += 1
        n_mismatch[A,B] = j #/2
        n_mismatch[B,A] = n_mismatch[A,B]

# un-comment the following bit of code to print info about relevant configurations. - KR
'''
print 'matrix of degree of relative substitution:'
print n_mismatch

print 'matrix of degree of relative excitation:'
print n_ex

print 'indices for(',NELECT,'-tuply)-substituted determinants (WRT each other): \n'
for A in range(NStates):
    for B in range(A+1,NStates):
        if (int(n_ex[A,B]) == 1):
            print 'configurations ',A+1,' and ',B+1,':'
            print 'substituted indices: ',e1_indx[A][B][0],',',e1_indx[A][B][1]
            print 'degree of excitation = ',int(n_ex[A,B])
            print '****************************************','\n'
        else:
            print 'configurations ',A+1,' and ',B+1,':'
#            print 'substituted indices: ',e1_indx[A][B][0],' ',e1_indx[A][B][1],'\n'
            print 'degree of excitation = ',int(n_ex[A,B])#,'\n'
            print 'degree of index mis-match = ',int(n_mismatch[A,B]),'\n'

print 'configuration matrix:'
print config_str
'''

# contribution from the frozen orbitals, and expectation values of the dipole moment
#  operator WRT states corresponding to the same configuration (needs to be multiplied
#  by proper CI expansion coefficients in the actual sum done to get the CI dip. mom. 
#  matrix elements), per Slater-Condon rules - KR

dip_frz = numpy.zeros((NStates,3), dtype='object')

for A in range(NStates):
    if (NFRZ > 0):
        for i in range(NFRZ):
            dip_frz[A,0] += 2*dipXMO[i,i]
            dip_frz[A,1] += 2*dipYMO[i,i]
            dip_frz[A,2] += 2*dipZMO[i,i]
    for i in range(NCAS):
        dip_frz[A,0] += float(config[A,i])*dipXMO[i+NFRZ,i+NFRZ]
        dip_frz[A,1] += float(config[A,i])*dipYMO[i+NFRZ,i+NFRZ]
        dip_frz[A,2] += float(config[A,i])*dipZMO[i+NFRZ,i+NFRZ]

# summing over all the MO basis matrix elements/MO integrals to get the CI basis dipole
#     moment matrix elements

tdip = numpy.zeros((NStates,NStates,3), numpy.float64)
# diagonal elements of CI dip. mom. matrix: tdip[A,A]
for A in range(NStates):
    for C in range(NStates):
# contribution from pair of zero mismatch determinants
        tdip[A,A,:] += (CI_vect[A,C]**2)*dip_frz[C,:]*(e1_fact[C])**2
        for D in range(C+1,NStates):
# contribution from pairs of one mismatch determinants
            if (int(n_ex[C,D]) == 1):
                i = e1_indx[C,D,0] + NFRZ
                j = e1_indx[C,D,1] + NFRZ
                tdip[A,A,0] += 2*CI_vect[A,C]*CI_vect[A,D]*e1_fact[C]*e1_fact[D]*dipXMO[i,j]
                tdip[A,A,1] += 2*CI_vect[A,C]*CI_vect[A,D]*e1_fact[C]*e1_fact[D]*dipYMO[i,j]
                tdip[A,A,2] += 2*CI_vect[A,C]*CI_vect[A,D]*e1_fact[C]*e1_fact[D]*dipZMO[i,j]
# off-diagonal elements of CI dip. mom. matrix: tdip[A,B]
    for B in range(A+1,NStates):
        for C in range(NStates):
# contribution from pair of zero mismatch determinants
            tdip[A,B,:] += CI_vect[A,C]*CI_vect[B,C]*dip_frz[C,:]*(e1_fact[C])**2
            for D in range(C+1,NStates):
# contribution from pairs of one mismatch determinants
                if (int(n_ex[C,D]) == 1):
                    i = e1_indx[C,D,0] + NFRZ
                    j = e1_indx[C,D,1] + NFRZ
                    tdip[A,B,0] += (CI_vect[B,C]*CI_vect[A,D]+CI_vect[A,C]*CI_vect[B,D])*e1_fact[C]*e1_fact[D]*dipXMO[i,j]
                    tdip[A,B,1] += (CI_vect[B,C]*CI_vect[A,D]+CI_vect[A,C]*CI_vect[B,D])*e1_fact[C]*e1_fact[D]*dipYMO[i,j]
                    tdip[A,B,2] += (CI_vect[B,C]*CI_vect[A,D]+CI_vect[A,C]*CI_vect[B,D])*e1_fact[C]*e1_fact[D]*dipZMO[i,j]
# symmetrizing the CI dip. mom. matrix. Inconsequential, given the propagation program reads only
#  the upper triangular part of this matrix (upper triangular parts of tdip[:,:,0/1/2] matrices). - KR
        tdip[B,A,:]=tdip[A,B,:]

#
# end of "straight-forward"ly CI dip. mom. matrix calculating piece of code. - KR
#

# converting the dipole moment matrix elements from a. u. to debyes, and factoring in the nuclear dipole moments.
#  the negative sign is bcoz Gaussian calculaes the matrix elements with their signs flipped. - KR
tdip[:,:,0] = (-tdip[:,:,2]+nuc_dipx)/0.393430
tdip[:,:,1] = (-tdip[:,:,2]+nuc_dipy)/0.393430
tdip[:,:,2] = (-tdip[:,:,2]+nuc_dipz)/0.393430

print 'CI basis dip. mom. for ',sN,' (z): \n',tdip[NStates-1,NStates-1,2],' D'



# have now calculated the Transition dipole moments. time to get atom coordinates and dump data
#  to file so TD_CIS.py can read it in. - bfh

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

# input SCF energy- assume energy of first NState - bfh
tdm.write(' \n ')
tdm.write('SCF Energy \n')
tdm.write('SCF Done: E(RHF) = ' + str(CI_E[0]) + "\n")
tdm.write(' \n ')

# Excited State energies - bfh
tdm.write('Excited state energies \n')
tdm.write(' \n')
for i in range(1,NStates):
    prt_E = (CI_E[i] - CI_E[0])*27.211396 
    tdm.write('Excited State   ' + str(i) + '     Singlet-BU     ' + str(prt_E) + ' eV  \n ')

# Ground state to excited state transition dipole moments - bfh
tdm.write(' \n')
tdm.write('\n')
tdm.write('Ground state to excited state transition dipole moments \n')
tdm.write(' \n ')
tdm.write('Ground to excited state transition electric dipole moments \n ')
tdm.write('     state      X       Y     Z     Dip. S.      Osc.     \n ' )
for i in range(1,NStates):
    tdm.write('%-10d %-10.6f %-10.6f %-10.6f %-10.6f %-10.6f' % (i, tdip[0][i][0], tdip[0][i][1], tdip[0][i][2], 0, 0 ) + "\n")

# Excited to excited state transition dipole moments - bfh
tdm.write(' \n')
tdm.write('\n')
tdm.write('Excited state to excited state transition dipole moments \n')
tdm.write(' \n ')
tdm.write('Excited to excited state transition electric dipole moments \n ')
tdm.write('  state I    state J      X       Y     Z     Dip. S.      Osc.     \n ' )
for i in range(2,NStates):
    for j in range(1,i):
        tdm.write('%10d %10d %-10.6f %-10.6f %-10.6f %-10.6f %-10.6f' % (i, j, tdip[j][i][0], tdip[j][i][1], tdip[j][i][2], 0, 0 ) + "\n")
 
tdm.write(' \n ')
tdm.write(' \n ')
 
# Finally, need to write out information on CI configs and MOs occupations - bfh
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
# - bfh

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





