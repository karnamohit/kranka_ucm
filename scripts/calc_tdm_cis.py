#!/home/kranka/local/miniconda3/bin/python

import numpy as np
import os
import shutil
import subprocess
# import sys
# import re
# import time

# sys.path.append('/path/to/kranka_ucm/scripts/')
from gauss_hf import *

Logfile = './cis-30roots_sys1-h_sto-3g.log'
print('file input:', Logfile)
print('detecting whether file exists...')
if os.path.isfile(Logfile) == True:
    print('file exists.')
elif os.path.isfile(Logfile) == False:
    print('error: file does not exist.')
    quit()

def find_linenum(s,txtfile): # locate matching string's line number(s) using `grep`
    strg = subprocess.Popen(['grep','-n', s, txtfile], stdout=subprocess.PIPE).stdout.read().decode("utf-8")
    lst = strg.split("\n")[:-1]
    lst = [int(lst[i].split(":")[0]) - 1 for i in range(len(lst))]
    return lst

# dipole moment conversion factor: 1 debye = 0.393430307 a.u.
DtoAU = 0.39343
AUtoD = 1./DtoAU
# energy conversion factor: 1 a.u. = 27.211396 eV
AUtoEV = 27.211396
EVtoAU = 1./AUtoEV
rt2 = np.sqrt(2)

terms = Logfile.split('/')[-1].split('_')
ci_str = list(terms[0])
if ci_str[1] == 'i':
    ci = 'cis'
else:
    print('not the right .LOG file for this code.')
    quit()

roots = terms[0].split('-')[1]
basis = terms[-1].split('.')[0]
molecule = terms[-2]
ext = '.'+terms[-1].split('.')[1]

outp = 'logfile.tmp'
shutil.copy(Logfile, outp)
datafile = log_data(outp)
log_lines = datafile.loglines

error = False
try:
    NMOs = datafile.nao        # total number of basis fns/MOs with double occupancy
    NELECT_A = datafile.n_a    # number of alpha electrons in the system
    NELECT_B = datafile.n_b    # number of beta electrons in the system
    NOCC = NELECT_A            # total number of occupied MOs in the reference
    NDOCC = NELECT_B           # number of doubly occupied MOs in the reference
    NSOCC = NELECT_A - NELECT_B    # number of singly occupied MOs in the reference
    NELECT = NELECT_A + NELECT_B
    NOCC = NSOCC + NDOCC
    lst = find_linenum('NAtoms',Logfile)
    for n in lst:
        elements = log_lines[n].split()
        NAtoms = int(float(elements[1]))
    lst = find_linenum('Range of M.O.s used for correlation',Logfile)
    for n in lst:
        elements = log_lines[n].split()
        NMOLow = int(float(elements[-2]))
        NMOHigh = int(float(elements[-1]))
        NFRZ = NMOLow - 1
        NCAS = NMOHigh - NMOLow + 1
    if (error == False):
        linenum = []
        lst = find_linenum('Excited State  ',Logfile)
        for n in lst:
            elements = log_lines[n].split()
            linenum.append(n)
        count = -1
        n = linenum[count]
        elements = log_lines[n].split()[2].split(':')
        NStates = int(float(elements[0])) + 1     # number of CI states/CSFs
    lst = find_linenum('Nuclear    moments (au)',Logfile)
    for n in lst:
        elements = log_lines[n+1].split()
        nuc_dipx = AUtoD*float(elements[1])
        nuc_dipy = AUtoD*float(elements[2])
        nuc_dipz = AUtoD*float(elements[3])
except (ValueError, IndexError, TypeError, NameError) as e:
    print(e)
    pass

dipX = datafile.get_dipole_x_AO()
dipY = datafile.get_dipole_y_AO()
dipZ = datafile.get_dipole_z_AO()

MO = datafile.get_MOcoeffs_AO()

M = NStates
N = int((NOCC - NMOLow + 1)*(NMOHigh - NOCC)) + 1
# forming an array of compound indices. This array will be used as reference to determine index of the singly-substituted determinant within the CI-vector array
cmpd_indx = []
cmpd_indx_dict = {}
k = 1
for i in range((NMOLow-1),NOCC):
    for a in range(NOCC,NMOHigh):
        ia = int((a+1)*((a+1)+1)/2 + (i+1))
        cmpd_indx.append(ia)
        cmpd_indx_dict[str(k)] = [i, a, ia]
        k += 1

CI_vect = np.zeros([M,N], np.float64)
CI_E = np.zeros([M], np.float64)
lst = find_linenum('E(RHF) =',Logfile)
for n in lst: # first element of CI_E will be the HF energy
    elements = log_lines[n].split()
    CI_E[0] = float(elements[4])
CI_vect[0,0] = np.float64(1)    # HF determinant contribution to GS is 1.0
ia = 0
lst = find_linenum('Excited State  ',Logfile)
for n in lst:
    elements = log_lines[n].split()[2].split(':')
    st_ci = int(float(elements[0]))
    CI_E[st_ci] = float(log_lines[n].split()[4])*EVtoAU + CI_E[0]
    for i in range(N):
        try:
            elements = log_lines[n+1+i].split('->')
            mo_i = int(float(elements[0]))
            mo_a = int(float(elements[1].split()[0]))
            c_ia = float(elements[1].split()[1])
            ia = int(mo_a*(mo_a+1)/2 + mo_i)
            CI_vect[st_ci,(cmpd_indx.index(ia)+1)] = c_ia*rt2
        except (ValueError, IndexError, TypeError):
            break

# print('normalization check:\n{}'.format(np.sum(np.diag(np.matmul(CI_vect, CI_vect.T))))); print('shape of CI_vect:\n{}'.format(CI_vect.shape));  quit()

# read MO configurations
config = np.zeros((M,NCAS), np.float64)
for i in range(M):
    for j in range(NCAS):
        if ((j+NFRZ) < NOCC):
            config[i,j] = 2
        else:
            config[i,j] = 0
    if (i > 0):
        for j in cmpd_indx_dict:
            indx = int(j)
            index_i = cmpd_indx_dict[str(j)][0]
            index_a = cmpd_indx_dict[str(j)][1]
            index_ia = int((index_a+1)*((index_a+1)+1)/2 + (index_i+1))
            occ_i = 1*CI_vect[i,indx]*np.conjugate(CI_vect[i,indx])
            occ_a = 1*CI_vect[i,indx]*np.conjugate(CI_vect[i,indx])
            config[i,(index_i-NMOLow+1)] -= occ_i
            config[i,(index_a-NMOLow+1)] += occ_a

print('calculating CI state dipole moments using RhoCI density-matrices...')
dipXCI = np.zeros([M,M], np.float64)
dipYCI = np.zeros([M,M], np.float64)
dipZCI = np.zeros([M,M], np.float64)
densCI_AO = np.zeros((NMOs, NMOs), np.float64)
last = NMOs % 5
if (last == 0):
    loops = int(NMOs / 5)
else:
    loops = int(NMOs / 5) + 1
lst = find_linenum('Alpha density for excited state',Logfile)
for n in lst:   # reading excited CI state density matrices
    st1 = int(float(log_lines[n].split()[-1]))
    shift = 0
    for k in range(loops):
        try:
            irange = NMOs - k*5
            for i in range(irange):
                if (k == (loops - 1)):
                    if (i < last):
                        end = i + 1
                    else:
                        end = last
                else:
                    if (i <= 4):
                        end = i + 1
                    else:
                        end = 5
                elements=log_lines[n+2+k+shift+i].split()
                for j in range(end):
                    s = k*5 + j
                    m = i + k*5
                    densCI_AO[m,s] = float(elements[j+1])
                    if (i != j):
                        densCI_AO[s,m] = densCI_AO[m,s]
            shift += irange
        except (IndexError, ValueError):
            break
    densCI_AO *= 2
    dipXCI[st1,st1] = -1*np.trace(np.matmul(densCI_AO,dipX))
    dipYCI[st1,st1] = -1*np.trace(np.matmul(densCI_AO,dipY))
    dipZCI[st1,st1] = -1*np.trace(np.matmul(densCI_AO,dipZ))
    dens_file = 'dens_ci_'+str(st1)+'_'+str(st1)+'.npz'
    # np.savez(dens_file, densCI_AO)
lst = find_linenum('Final density matrix:',Logfile)
for n in lst:   # reading ground state (HF) density
    st1 = 0
    shift = 0
    for k in range(loops):
        try:
            irange = NMOs - k*5
            for i in range(irange):
                if (k == (loops - 1)):
                    if (i < last):
                        end = i + 1
                    else:
                        end = last
                else:
                    if (i <= 4):
                        end = i + 1
                    else:
                        end = 5
                elements=log_lines[n+2+k+shift+i].split()
                for j in range(end):
                    s = k*5 + j
                    m = i + k*5
                    densCI_AO[m,s] = float(elements[j+1])
                    if (i != j):
                        densCI_AO[s,m] = densCI_AO[m,s]
            shift += irange
        except (IndexError, ValueError):
            break
    densCI_AO *= 2
    dipXCI[st1,st1] = -1*np.trace(np.matmul(densCI_AO,dipX))
    dipYCI[st1,st1] = -1*np.trace(np.matmul(densCI_AO,dipY))
    dipZCI[st1,st1] = -1*np.trace(np.matmul(densCI_AO,dipZ))
    dens_file = 'dens_ci_'+str(st1)+'_'+str(st1)+'.npz'
    # np.savez(dens_file, densCI_AO)
print('reading CI transition dipole moments...')
lst = find_linenum('Excited to excited state transition electric dipole moments',Logfile)
for n in lst:
    for i in range(int((M-2)*(M-1)/2)):
        try:
            elements = log_lines[n+2+i].split()
            st1 = int(float(elements[0]))
            st2 = int(float(elements[1]))
            dipXCI[st1,st2] = float(elements[2])
            dipXCI[st2,st1] = dipXCI[st1,st2]
            dipYCI[st1,st2] = float(elements[3])
            dipYCI[st2,st1] = dipYCI[st1,st2]
            dipZCI[st1,st2] = float(elements[4])
            dipZCI[st2,st1] = dipZCI[st1,st2]
        except (IndexError, ValueError):
            break
lst = find_linenum('Ground to excited state transition electric dipole',Logfile)
for n in lst:
    st1 = 0
    for i in range(1,M):
        try:
            elements = log_lines[n+1+i].split()
            st2 = int(float(elements[0]))
            dipXCI[st1,st2] = float(elements[1])
            dipXCI[st2,st1] = dipXCI[st1,st2]
            dipYCI[st1,st2] = float(elements[2])
            dipYCI[st2,st1] = dipYCI[st1,st2]
            dipZCI[st1,st2] = float(elements[3])
            dipZCI[st2,st1] = dipZCI[st1,st2]
        except (IndexError, ValueError):
            break

'''
densCI_AO = 0 # reading transition density matrices b/w ground and excited states
lst = find_linenum('Alpha transition density to state',Logfile)
for n in lst:
    elements = log_lines[n].split()
    st1 = 0
    st2 = int(float(elements[-1]))
    densCI_AO_a = np.zeros([NMOs,NMOs], np.float64)
    densCI_AO_b = np.zeros([NMOs,NMOs], np.float64)
    for k in range(loops):
        for i in range(NMOs):
            try:
                if (k == (loops - 1)):
                    end = last
                else:
                    end = 5
                dum1 = n+2+k*(1+NMOs)+i
                dum2 = dum1+1+loops*(NMOs+1)
                elements_a = log_lines[dum1].split()
                elements_b = log_lines[dum2].split()
                for j in range(end):
                    s = k*5 + j
                    m = i
                    densCI_AO_a[m,s] = elements_a[j+1]
                    densCI_AO_b[m,s] = elements_b[j+1]
            except (IndexError, ValueError):
                break
    densCI_AO = (densCI_AO_a + densCI_AO_b)/rt2
    dens_file = 'dens_ci_'+str(st1)+'_'+str(st2)+'.npz'
    np.savez(dens_file, densCI_AO)

densCI_AO = 0 # reading transition density matrices b/w excited states
lst = find_linenum('Alpha transition density matrix between excited states',Logfile)
for n in lst:
    elements = log_lines[n].split()
    st1 = int(float(elements[-1]))
    st2 = int(float(elements[-3]))
    densCI_AO_a = np.zeros([NMOs,NMOs], np.float64)
    densCI_AO_b = np.zeros([NMOs,NMOs], np.float64)
    for k in range(loops):
        for i in range(NMOs):
            try:
                if (k == (loops - 1)):
                    end = last
                else:
                    end = 5
                dum1 = n+2+k*(1+NMOs)+i
                dum2 = dum1+2+loops*(NMOs+1)
                elements_a = log_lines[dum1].split()
                elements_b = log_lines[dum2].split()
                for j in range(end):
                    s = k*5 + j
                    m = i
                    densCI_AO_a[m,s] = float(elements_a[j+1])
                    densCI_AO_b[m,s] = float(elements_b[j+1])
            except (IndexError, ValueError):
                break
    densCI_AO = densCI_AO_a + densCI_AO_b
    dens_file = 'dens_ci_'+str(st1)+'_'+str(st2)+'.npz'
    np.savez(dens_file, densCI_AO)
# '''

# writing output
outname = 'tdm_tdci.txt' # name of the output file, to be read by TDCI script
tdm = open(outname, 'w')
tdm.write('===========================================================================================\n')
tdm.write('   This file contains the information one needs to run a time-dependent CIS simulation.    \n')
tdm.write('===========================================================================================\n')
tdm.write('An initial Gaussian calculation was used to calculate the CI eigenvalues and eigenvectors  \n')
tdm.write('as well as the CI state configurations. The program calc_tdm.py was then used to calculate \n')
tdm.write('the transition dipole moments between the CI states as well as the project of the CI states\n')
tdm.write('onto the molecular orbitals. This file was generated by calc_tdm.py as input to TD_CIS.py. \n')
tdm.write(' -BFH 9 May 2014 (modified by KR 8 June 2019)                                              \n')
tdm.write('\n\n')
tdm.write('CI calculation parameters\n\n') # parameters for the TDCI script
tdm.write(' title   = {}\n'.format(molecule+'_'+basis))
tdm.write(' method  = {}\n'.format(ci))
tdm.write(' NMOLow  = {}\n'.format(NMOLow))
tdm.write(' NMOHigh = {}\n'.format(NMOHigh))
tdm.write(' NMOs    = {}\n'.format(NMOs))
tdm.write(' NCAS    = {}\n'.format(NCAS))
tdm.write(' NFRZ    = {}\n'.format(NFRZ))
tdm.write(' NOCC    = {}\n'.format(NOCC))
tdm.write(' NELECT  = {}\n'.format(NELECT))
tdm.write(' NAtoms  = {}\n'.format(NAtoms))
tdm.write(' NStates = {}\n'.format(NStates))
tdm.write('\n\n')
tdm.write('Nuclear dipole moments (a.u.) \n\n')
tdm.write(' x = {: 4.3f} \n'.format(nuc_dipx*DtoAU))
tdm.write(' y = {: 4.3f} \n'.format(nuc_dipy*DtoAU))
tdm.write(' z = {: 4.3f} \n'.format(nuc_dipz*DtoAU))
tdm.write('\n\n')
tdm.write('Setting up the atomic coordinates \n\n')
tdm.write(' Standard orientation: \n')
tdm.write(' ------------------------------------------------------------------- \n')
tdm.write(' Center   Atomic     Atomic            Coordinates (Angstroms)       \n')
tdm.write(' Number   Number      Type            X            Y          Z      \n')
tdm.write(' ------------------------------------------------------------------- \n')
coords = np.zeros([NAtoms,3], np.float64)
atom_info = np.zeros([NAtoms,3], np.float64)
for (n, line) in enumerate(log_lines):
    if (' Standard basis:' in line):
        for i in range(NAtoms):
            elements = log_lines[n + i + 6].split()
            atom_info[i,0] = float(elements[0])
            atom_info[i,1] = float(elements[1])
            atom_info[i,2] = float(elements[2])
            coords[i,0] = float(elements[3])
            coords[i,1] = float(elements[4])
            coords[i,2] = float(elements[5])
            tdm.write('  %-10s %-10s %-10s %-10s %-10s %-10s' % (atom_info[i][0], atom_info[i][1], atom_info[i][2], coords[i][0], coords[i][1], coords[i][2]) + "\n")
tdm.write('\n\n')
tdm.write('SCF Energy \n\n')
tdm.write(' SCF Done: E(RHF) = ' + str(CI_E[0]) + '  Ha \n')
tdm.write('\n\n')
tdm.write('Excited state energies (eV) \n\n')
for i in range(1,NStates):
    prt_E = (CI_E[i] - CI_E[0])*AUtoEV
    tdm.write(' Excited State {:3d} \t {:4.6f}\n'.format(i, prt_E))
tdm.write('\n\n')
tdm.write('CI electric dipole moments (a.u.) \n\n')
tdm.write(' CI state and transition dipole moments\n')
tdm.write(' state I  state J \t      X \t\t      Y \t\t      Z \t\t   Dip. S. \t\t    Osc.\n' )
for i in range(NStates):
    for j in range(NStates):
        wij = abs(CI_E[i] - CI_E[j])
        dij = np.array([dipXCI[i,j], dipYCI[i,j], dipZCI[i,j]])
        dipsij = np.linalg.norm(dij)**2
        fij = (2/3) * wij * dipsij
        tdm.write('   {} \t    {} \t\t  {:3.7f}  \t  {:3.7f}  \t  {:3.7f}  \t  {:3.7f}  \t  {:3.7f}\n'.format(i, j, dij[0], dij[1], dij[2], dipsij, fij))
tdm.write('\n\n')
tdm.write('Output CI Vectors \n\n')
tdm.write(' CI state       CI coeffs 1...N \n')
for i in range(M):
    tdm.write('  %4d \t\t' % i )
    for j in range(N):
        tdm.write('%6.6f\t' % CI_vect[i,j])
    tdm.write('\n')
tdm.write('\n\n')
mo_weight = np.zeros([NCAS], np.float64)
tdm.write('Orbital electron configuration weighting for CI States \n\n')
tdm.write(' CI State \t\t Orbital Weights 1...NMOs  \n')
for i in range(M):
    for k in range(NCAS):
        mo_weight[k] = float(config[i,k])
    tdm.write('  {:3d}\t\t'.format(i))
    for k in range(NCAS):
        if (k == (NCAS - 1)):
            tdm.write('{:4.6f}\n'.format(mo_weight[k]))
        else:
            tdm.write('{:4.6f}\t'.format(mo_weight[k]))
        mo_weight[k] = 0
tdm.write('\n\n')
tdm.close()
