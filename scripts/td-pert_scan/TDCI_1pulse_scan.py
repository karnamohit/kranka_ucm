import numpy as np
import scipy.linalg as sl
import os
# import scipy.integrate as si
# import matplotlib.pyplot as plt
# import sys
# import shutil
# import time

# sys.path.append('/path/to/kranka_ucm/scripts/')
from gauss_hf import *

#====================
# CONVERSION FACTORS
#====================
# dipole moment: 1 debye = 0.393430307 a.u.
DtoAU = 0.39343
AUtoD = 1./DtoAU
# energy: 1 Ha = 27.211396 eV
#              = 219474.6305 /cm
#              = 6.579691506102500e15 Hz
HAtoEV = 27.211396
EVtoHA = 1./HAtoEV
HAtoHZ = 6.5796915061025e15
HZtoHA = 1./HAtoHZ
# time: 1 a.u. = 0.02418884254 fs
AUtoFS = 0.02418884
FStoAU = 1./AUtoFS
#===========
# CONSTANTS
#===========
# atomic {int(number): [str(symbol), float(mass)]} (up to 2nd row)
AtomicMass = {1: ['H', 1.0079], 2: ['He', 4.0026], 3: ['Li', 6.941], \
              4: ['Be', 9.0122], 5: ['B', 10.811], 6: ['C', 12.0107], \
              7: ['N', 14.0067], 8: ['O', 15.9994], 9: ['F', 18.9984], \
              10: ['Ne', 20.1797] \
             }

# this is the output file produced by "calc_tdm_<cis or casscf>_latest"
Logfile = 'tdm_tdci.txt'

# check whether "tdm_tdci.txt" exists
if os.path.isfile(Logfile) == True:
    print('file exists.')
elif os.path.isfile(Logfile) == False:
    print('error: file does not exist.')
    quit()

data_file = log_data(logfile=Logfile, nonlog_error_msg=False)
log_lines = data_file.loglines

for (n,line) in enumerate(log_lines):
    try:
        if ('title' in line):
            # system name for output directory
            molsys = log_lines[n].split('=')[-1].strip()
        if ('method' in line):
            # CI method
            ci = log_lines[n].split('=')[-1].strip()
        if ('NMOLow' in line):
            # Lowest MO in the active space of CIS
            NMOLow = int(float(log_lines[n].split('=')[-1]))
        if ('NMOHigh' in line):
            # Highest MO in the active space of CIS
            NMOHigh = int(float(log_lines[n].split('=')[-1]))
        if ('NMOs    =' in line):
            # Total number of MOs, set to NCAS later
            NMOs = int(float(log_lines[n].split('=')[-1]))
        if ('NCAS' in line):
            # Number of MOs in the CI active space
            NCAS = int(float(log_lines[n].split('=')[-1]))
        if ('NFRZ' in line):
            # Number of MOs w/ frozen occupation
            NFRZ = int(float(log_lines[n].split('=')[-1]))
        if ('NOCC' in line):
            # Number of (doubly) occupied MOs
            NOCC = int(float(log_lines[n].split('=')[-1]))
        if ('NELECT' in line):
            # Number of electrons
            #   (= 2 x # alpha electrons for restricted reference)
            NELECT = int(float(log_lines[n].split('=')[-1]))
        if ('NAtoms' in line):
            # Number of atoms in the molecular system
            NAtoms = int(float(log_lines[n].split('=')[-1]))
        if ('NStates' in line):
            # Number of CI states (including the ground state)
            NStates = int(float(log_lines[n].split('=')[-1]))
    except (IndexError, ValueError) as e:
        print(e)
        print('Error encountered while reading parameters.')
        print('Check the contents of "tdm_tdci.txt"')
        quit()
NMOs = NCAS        # Number of MOs in CAS
NSim = NStates     # Number of states in simulation (redundant at the moment)

if (ci == 'cis'):
    M = NStates
    N = int((NOCC - NMOLow + 1)*(NMOHigh - NOCC)) + 1
else:
    print('CIS parameters not found.')
    quit()

AtNum = np.zeros([NAtoms], np.float64)
AtMass = np.zeros([NAtoms], np.float64)
CoordsX = np.zeros([NAtoms], np.float64)
CoordsY = np.zeros([NAtoms], np.float64)
CoordsZ = np.zeros([NAtoms], np.float64)
for (n, line) in enumerate(log_lines):
    if ('Standard orientation:' in line): 
        for i in range(NAtoms):
            elements = log_lines[n+i+5].split()
            AtNum[i] = int(float(elements[1]))
            AtMass[i] = AtomicMass[AtNum[i]][1]
            CoordsX[i] = float(elements[3])
            CoordsY[i] = float(elements[4])
            CoordsZ[i] = float(elements[5])
for (n, line) in enumerate(log_lines):
    if ('Nuclear dipole moments (a.u.)' in line):
        i = n + 2
        NucDipX = float(log_lines[i].split('=')[-1])
        NucDipY = float(log_lines[i+1].split('=')[-1])
        NucDipZ = float(log_lines[i+2].split('=')[-1])
for (n, line) in enumerate(log_lines): # Get SCF Energy
    if ('SCF Eone:' in line):
        elements = log_lines[n].split()
        ESCF = float(elements[4])
CI_Energies = np.zeros([M], np.float64) # Get CIS energy
for (n,line) in enumerate(log_lines):
    if('Excited State ' in line):
        elements = log_lines[n].split()
        st = int(float(elements[2]))
        CI_Energies[st] = float(elements[3])*EVtoHA
CI_Energies[:] += ESCF

CI_Dipoles = np.zeros((M, M), dtype='object')
linenum = []
for (n,line) in enumerate(log_lines):
    try:
        if ('CI state and transition dipole moments' in line):
            linenum.append(n)
    except (IndexError, ValueError):
        break
count = -1
n = linenum[count]
k = 2
while (k <= int(M**2 + 2)):
    try:
        elements = log_lines[n+k].split()
        st1 = int(float(elements[0]))
        st2 = int(float(elements[1]))
        dipvec = [float(elements[2]), float(elements[3]), float(elements[4])] 
        CI_Dipoles[st1,st2] = dipvec
        k += 1
    except (IndexError, ValueError):
        break

CI_vect = np.zeros([M, N], np.float64)
linenum = []
for (n,line) in enumerate(log_lines):
    try:
        if ('Output CI Vectors' in line):
            linenum.append(n)
    except (IndexError, ValueError):
        break
n = linenum[-1]
k = 3
while (k <= (M + 3)):
    try:
        elements = log_lines[n+k].split()
        i = int(float(elements[0]))
        for j in range(N):
            CI_vect[i,j] = float(elements[j+1])
        k += 1
    except (IndexError, ValueError):
        break

config = np.zeros([M,NMOs], np.float64) # get state-configuration info
linenum = []
for (n,line) in enumerate(log_lines):
    if ('Orbital electron configuration' in line):
        linenum.append(n)
n = linenum[-1]
k = 3
while (k < (M + 3)):
    try:
        elements = log_lines[n+k].split()
        i = int(float(elements[0]))
        for j in range(NMOs):
            try:
                config[i,j] = float(elements[j+1])
            except (IndexError, ValueError):
                break
        k += 1
    except (IndexError, ValueError):
        break

st_ini = 0 # starting and targeted state-indices: chosen from the set of CI states (int)
st_trg = 1
deltaE1 = abs(CI_Energies[st_ini] - CI_Energies[st_trg]) # excitation energy (in au) (float)
delta = 0.08268 # ~0.0020 fs    # step-size (in au) (float)
# delta = 2.067   # ~0.05 fs
steptofs = delta * AUtoFS # conversion factor for time-steps to fs (float)
w, tdperturb = deltaE1, 'ress0s1' # frequency of field (in Ha) (float)
ncyc = REPLACE_CYC    # cycles for which field is switched on (int)
phase = 0.0
period_au = 2.*np.pi / w
period_tsteps = int(period_au/delta)
duration = int(ncyc*period_tsteps) # duration, in time-steps, for which the field is switched on (int)
duration_au = duration*delta
duration_fs = duration_au * AUtoFS
emax = REPLACE_AMP      # electric field amplitude (in au) (float)
ton = 0                 # step number to turn on field (int)
toff = ton + duration   # step number to turn off field (int)
maxstep = 10000 + toff  # number of integration steps (int)
envelope_apply = False  # using field-envelope for the field? (logical, 1 or 0)
envelope = 'gaussian'   # the envelope shape ('trapezoidal', 'gaussian', 'cosine')

if (envelope_apply):
    print('{} envelope applied.'.format(envelope))
    # split "duration" into cycles within the envelope
    cyc = 1                         # can be "ncyc"
    chirp = 0.01                    # chirp (in au)
    cyc += 2                        # in case cyc is set to 1 or 0, to avoid breaking the code
    period_tsteps_scaled = int(duration/cyc)    # scaled period in time-steps
    decay = 0.05                    # exponent that determines sharpness of inflection (smooth trapezoid)
    y= np.zeros((maxstep), np.float64)
    if(envelope == 'trapezoidal'):  # envelope = 'trapezoidal'
        for t in range(ton,toff):
            if(t < toff):
                t_red_up = t-ton-period_tsteps_scaled/2
                t_red_down = t_red_up-(cyc-1)*period_tsteps_scaled
                if (t < ((cyc - 1)*period_tsteps_scaled + ton)):
                    if (t < (period_tsteps_scaled + ton)):
                        y[t] = 1/(1+np.exp(-decay*t_red_up))
                    else:
                        y[t] = 1/(1+np.exp(-decay*(period_tsteps_scaled/2)))
                else:
                    y[t] = 1-1/(1+np.exp(-decay*t_red_down))
            else:
                y[t] = 0
    elif(envelope == 'gaussian'):   # envelope = 'gaussian'
        for t in range(ton,toff):
            if (t < toff):
                y[t] = np.exp(-16*np.log(2)*((toff - t)/duration-1/2)**2)
                y[t] -= 1/16
                y[t] /= 15/16
            else:
                y[t] = 0
    elif(envelope == 'cosine'):     # envelope = 'cosine'
        for t in range(ton,toff):
            if (t < toff):
                y[t] = 0.5 - np.cos(2*(np.pi)*(toff - t)/duration)/2
            else:
                y[t] = 0
efield = np.zeros(maxstep)
if (envelope_apply == True):
    for t in range(ton,toff):
        w_chirp = w + chirp * (t/(ncyc*period_tsteps) - 0.5)    # apply an envelope and a chirp-freq to the sinusoidal field
        efield[t] = emax * y[t] * (np.sin(w_chirp * delta * t + phase))
elif (envelope_apply == False):
    for t in range(ton,toff):
        efield[t] = emax * (np.sin(w * delta * t + phase))      # apply bare sinusoidal field

# direction for the field
x = 1.0     # +ve, -ve or 0.0
y = 0.0     # +ve, -ve or 0.0
z = 0.0     # +ve, -ve or 0.0
direction = [x,y,z]
norm_direc = direction/(np.linalg.norm(direction))   # normalization

# initialize time-dependent arrays
dip = np.zeros((maxstep), np.float64)           # dipole moment
v = np.zeros((maxstep, NSim), np.complex128)    # vector of CI states that form the system's wave-function
v_norm = np.zeros((maxstep), np.float64)
h0 = np.diag(CI_Energies[:] - CI_Energies[0])   # (diagonal) Hamiltonian in the CI state basis (static in time)
d0 = np.zeros((NSim, NSim), np.float64)         # dipole moment matrix in the CI state basis
for i in range(NSim):
    try:
        for j in range(NSim):
            d0[i,j] = np.dot(norm_direc, CI_Dipoles[i,j])
    except (TypeError, ValueError):
        print('Check the "CI_Dipoles" or "norm_direc" matrix!')
        print('(make sure they are filled with correct data.)')
        break
# initial conditions
v[0,st_ini] = 1
v_norm[0] = 1
norm0 = np.linalg.norm(v[0,:])
ctmp1 = v[0,:]/norm0
dip[0] = np.real(np.dot(np.dot(np.conjugate(ctmp1), d0), ctmp1))

# TDCI propagation
for tstep in range(maxstep-1):
    # compute a unitary propagator using the **matrix** exponential
    propagator = sl.expm(-1J*delta*(h0 + -1*efield[tstep]*d0))
    # advance time by "delta" units
    v[tstep+1,:] = propagator @ v[tstep,:]
    v_norm[tstep+1] = np.linalg.norm(v[tstep+1,:])
    # compute the dipole moment
    dip[tstep+1] = np.real(np.dot(np.dot(np.conjugate(v[tstep+1,:]), d0), v[tstep+1,:]))

title = ' REPLACE_TTL '
title = title.strip()

#output time-dependent TDCI coefficients
prop_file = 'time_coeffs.'+title+'.txt'
name_file = prop_file
file = name_file
bcoef = open(file, 'w')
for j in range(2,maxstep):
    t_fs = j*steptofs
    bcoef.write('{:16.8f}'.format(t_fs))
    for i in range(NSim):
        vcoef = v[j,i]
        bcoef.write('\t\t{:16.8f}'.format(vcoef))
    bcoef.write('\n')
bcoef.close()

# output time-dependent electric dipole moment
prop_file = 'time_dipole.'+title+'.txt'
name_file = prop_file
file = name_file
idm = open(file, 'w')
prop_file = 'time_dip2fft.'+title+'.txt'
name_file = prop_file
file = name_file
sdm = open(file, 'w')
# dip *= AUtoD    # Dipole moment in debyes
for i in range(maxstep):
    t_fs = i*steptofs
    t_s = t_fs*1e-15
    idm.write('{:14.8f}\t{: 14.8f}\n'.format(t_fs, dip[i]))
    sdm.write('{:14.8e}\t{: 14.8f}\n'.format(t_s, dip[i]))
idm.close(); sdm.close()

# output time-dependent MO occupations
prop_file = 'time_occup.'+title+'.txt'
name_file = prop_file
file = name_file
tcoef = open(file, 'w')
mo_coef = np.zeros([NMOs], np.float64)
for j in range(maxstep):
    t_fs = j*steptofs
    tcoef.write('{:16.8f}'.format(t_fs))
    for k in range(NMOs):
        mo_coef[k] = 0.0
        for n in range(NSim):
            vcoef = np.conjugate(v[j,n]) * v[j,n]
            mo_coef[k] += vcoef*config[n,k]
        tcoef.write('\t{:14.8f}'.format(mo_coef[k]))
    totN = np.sum(mo_coef)  # total number of electrons at each time-step: a quality-check; should be constant.
    tcoef.write('\t\t{:4.4f}\n'.format(totN))
tcoef.close()

# output time-dependent electric field
prop_file = 'time_field.'+title+'.txt'
name_file = prop_file
file = name_file
fldprt = open(file, 'w')
td_field = np.zeros([3,maxstep])
for i in range(1,maxstep):
    t_fs = i*steptofs
    fldprt.write('{:16.8f}\t\t{: 16.8f}\t\t{: 16.8f}\t\t{: 16.8f}\n'.format(t_fs, direction[0]*efield[i],direction[1]*efield[i],direction[2]*efield[i]))
fldprt.close()

'''
NAOs = NCAS + NFRZ
densTot = np.zeros([maxstep,NAOs,NAOs], np.complex128)
for i in range(NSim):
    for j in range(i,NSim):
        filename = 'dens_ci_'+str(i)+'_'+str(j)+'.npz'
        dens = np.load(filename, allow_pickle=True)
        densCI = dens['arr_0']
        for t in range(maxstep):
            if (i != j):
                # densCI is real, so transpose same as c.c.
                densTot[t,:,:] += v[t,i]*np.conjugate(v[t,j]) * densCI
                densTot[t,:,:] += v[t,j]*np.conjugate(v[t,i]) * densCI.T
            else:
                densTot[t,:,:] += (abs(v[t,i])**2) * densCI
densTot *= 0.5
save_file = saveat + 'td_dens_re+im_'+std_prefix+'_' + tdperturb + '_s0_' + title+'.npz'
print(save_file)
np.savez(save_file, td_dens_im_data=densTot.imag, td_dens_re_data=densTot.real)
# '''
