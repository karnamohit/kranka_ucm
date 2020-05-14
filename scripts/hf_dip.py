#!/usr/bin/env python

# Script to extract data from the Gaussian .log AND .fchk files (.fchk for the MO-coeffs.) - KR

import numpy as np
import sys
import os
import shutil
import subprocess

#
# Parameters
#

NMOs = 2  	# no. of MOs

gausslog = 'h2_s0.log'
#tmplog = 'log.tmp'
log = 'log.tmp'
#shutil.copy(gausslog, log)
fchk = 'h2_s0.fchk'

logfile = open(log, 'w')
sub = subprocess.call(['sed', 's/D/E/g', gausslog], stdout=logfile)
logfile.close()

logfile = open(log, 'r')

lines = logfile.readlines()

#
# Reading the matrices in the AO-basis, the MO co-efficient matrix
#

#
# AO basis overlap
#

SAO = np.zeros([NMOs,NMOs], np.float64)

if (NMOs < 6):
    for (n,line) in enumerate(lines):
        if ('*** Overlap ***' in line):
            for i in range (NMOs):
                j = 0
                while (j <= i):
                    elements = lines[n+i+2].split()
                    k = j + 1
                    SAO[i][j] = float(elements[k])
                    if ( i != j):
                        SAO[j][i] = SAO[i][j]
                    j=j+1
else:
    loops = NMOs / 5 + 1
#    print 'loops = ', loops
    last = NMOs % 5
#    print 'last = ', last
    shift = 0
    for (n,line) in enumerate(lines):
        if ('*** Overlap ***' in line):
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
                        SAO[m][s] = float(elements[j+1])
                        if i != j:
                            SAO[s][m] = SAO[m][s]
                shift = shift + irange

print 'Overlap: \n', SAO

#
# Dipole moment matrices
#

dipX = np.zeros([NMOs,NMOs], np.float64)
dipY = np.zeros([NMOs,NMOs], np.float64)
dipZ = np.zeros([NMOs,NMOs], np.float64)

# dipole X

if (NMOs < 6):
    for (n,line) in enumerate(lines):
        if ('Multipole matrices IBuc=  518 IX=    1' in line):
            for i in range (NMOs):
                j = 0
                while (j <= i):
                    elements = lines[n+i+2].split()
                    k = j + 1
                    dipX[i][j] = float(elements[k])
                    if ( i != j):
                        dipX[j][i] = dipX[i][j]
                    j=j+1
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

if (NMOs < 6):
    for (n,line) in enumerate(lines):
        if ('Multipole matrices IBuc=  518 IX=    2' in line):
            for i in range (NMOs):
                j = 0
                while (j <= i):
                    elements = lines[n+i+2].split()
                    k = j + 1
                    dipY[i][j] = float(elements[k])
                    if ( i != j):
                        dipY[j][i] = dipY[i][j]
                    j=j+1
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

# dipole Z

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

print 'dip. mom. matrix (x): \n', dipX
print 'dip. mom. matrix (y): \n', dipY
print 'dip. mom. matrix (z): \n', dipZ


# AO-basis density

densAO = np.zeros([NMOs,NMOs], np.float64)

if (NMOs < 6):
    for (n,line) in enumerate(lines):
        if ('Eensity Matrix:' in line):
            for i in range (NMOs):
                j = 0
                while (j <= i):
                    elements = lines[n+i+2].split()
                    k = j + 4
                    densAO[i][j] = float(elements[k])
                    if ( i != j):
                        densAO[j][i] = densAO[i][j]
                    j=j+1
                    #densAO[i][j] = float(elements[j+1])
                    #if ( i != j):
                        #densAO[j][i] = densAO[i][j]
else:
    loops = NMOs / 5 + 1
#    print 'loops = ', loops
    last = NMOs % 5
#    print 'last = ', last
    shift = 0
    for (n,line) in enumerate(lines):
        if ('Eensity Matrix:' in line):
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
                        densAO[m][s] = float(elements[j+4])
                        if i != j:
                            densAO[s][m] = densAO[m][s]
                shift = shift + irange

print 'AO dens. matrix: \n',densAO

#
# MO-coeffs. from the .fchk file
#

alphaMO = []
MO = np.zeros([NMOs,NMOs], np.float64)

fchk_file = open(fchk,'r')
fchk_lines = fchk_file.readlines()

for (n,line) in enumerate(fchk_lines):
    if ('Alpha MO coefficients' in line):
        loops = NMOs*NMOs / 5 + 1
        x = 0
        for i in range(loops):
            elements = fchk_lines[n+1+x].split()
            alphaMO.extend(elements)
            x += 1
        print 'elements \n', elements, '\n \n', 'alphaMO \n', alphaMO

for i in range(NMOs*NMOs):
    alphaMO[i] = float(alphaMO[i])

for i in range(NMOs):
    for j in range(NMOs):
        MO[i][j] = alphaMO[i*NMOs+j]

MO = MO.T

print 'MO coeffs. \n', MO

densMO = MO.T.dot(SAO.dot(densAO.dot(SAO.dot(MO))))

dipZMO = MO.T.dot(dipZ.dot(MO))

print 'dip. mom. matrix along the z-axis (MO): \n', dipZMO

muZ = -1*np.trace(densAO.dot(dipZ))

muZp = -1*np.trace(densMO.dot(dipZMO))

SMO = MO.T.dot(SAO.dot(MO))

print 'MO overlap: \n', SMO

conv_dip = 0.393430307

print 'muZ_e = ', muZ/conv_dip, ' D'
print 'muZp_e = ', muZp/conv_dip, ' D'


