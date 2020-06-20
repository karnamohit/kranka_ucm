#!/usr/bin/env python
"""Extracts AO matrices from a Gaussian .LOG file. Provides log_data.

log_data reads the .LOG file and stores lines, which its methods can
use to extract the AO matrices. Also evaluates an effective one-elec-
tron matrix for electron-electron interaction operator.
"""

import sys
import os
import re
import numpy as np

__author__ = "Karnamohit Ranka"
__copyright__ = "N/A"
__credits__ = ["Karnamohit Ranka", "Brad Habenicht"]
__license__ = "N/A"
__version__ = "N/A"
__maintainer__ = "Karnamohit Ranka"
__email__ = "kranka@ucmerced.edu"
__status__ = "N/A"

class log_data:
    #
    def __init__(self, logfile=None, print_val=None):
        #
        if (print_val != None):
            if (print_val == False):
                print('set print_val=True to print info about this module\n')
            else:
                print_info(print_val)
        if (logfile == None):
            print('\nNo file specified!')
            print('Exiting...')
            return
        else:
            # check if logfile exists
            print('\nDetecting whether the specified file exists...')
            if os.path.isfile(logfile) == True:
                print('\nFile detected. Specified file: "{}"\n'.format(logfile))
                self.logfile = logfile
            elif os.path.isfile(logfile) == False:
                print('\nError: Specified file does not exist.')
                print('Exiting...')
                return
        #
        # reading "logfile" contents
        log_file = open(self.logfile, 'r')
        print('\nReading data from {}...'.format(self.logfile))
        self.loglines = log_file.readlines()
        log_file.close()
        self.loglines = [re.sub(r'D','E', s) for s in self.loglines]
        #
        error = False
        for (n,line) in enumerate(self.loglines):
            # read # basis fns and electrons
            try: 
                if ('primitive gaussians' in line):
                    elements = self.loglines[n].split()
                    elements2 = self.loglines[n+1].split()
                    #print(elements[0], elements2[0], elements2[3])
                    self.nao = int(float(elements[0]))
                    self.n_a = int(float(elements2[0]))
                    self.n_b = int(float(elements2[3]))
                    #
                    self.nel = self.n_a + self.n_b
                    self.vir_a = self.nao - self.n_a
                    self.vir_b = self.nao - self.n_b
            except (IndexError, ValueError, TypeError):
                print('\n.LOG file not as expected.\nError occurred while trying to read data.')
                print('\nExiting...')
                error = True
        #
        if error == True:
            return
        elif error == False:
            print('\nData read.')
    #
    def get_overlap_AO(self):
        #
        line_num = []
        for (n, line) in enumerate(self.loglines):
            if ('*** Overlap ***' in line):
                line_num.append(n)
        #
        overlap_data = np.zeros((self.nao, self.nao), np.float64)
        #
        count = 0
        n = line_num[count]
        shift = 0
        loops = int(self.nao / 5) + 1
        last = int(self.nao % 5)
        for k in range(loops):
            try:
                irange = self.nao - k*5
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
                    elements=self.loglines[n+2+k+shift+i].split()
                    for j in range(end):
                        s = k*5 + j
                        m = i + k*5
                        overlap_data[m,s] = float(elements[j+1])
                        if (i != j):
                            overlap_data[s,m] = overlap_data[m,s]
                shift += irange
            except (IndexError, ValueError):
                break
        #
        return overlap_data
    #
    def get_kinetic_AO(self):
        #
        line_num = []
        for (n, line) in enumerate(self.loglines):
            if ('*** Kinetic Energy ***' in line):
                line_num.append(n)
        #
        ke_data = np.zeros((self.nao, self.nao), np.float64)
        #
        count = 0
        n = line_num[count]
        shift = 0
        loops = int(self.nao / 5) + 1
        last = int(self.nao % 5)
        for k in range(loops):
            try:
                irange = self.nao - k*5
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
                    elements=self.loglines[n+2+k+shift+i].split()
                    for j in range(end):
                        s = k*5 + j
                        m = i + k*5
                        ke_data[m,s] = float(elements[j+1])
                        if (i != j):
                            ke_data[s,m] = ke_data[m,s]
                shift += irange
            except (IndexError, ValueError):
                break
        #
        return ke_data
    #
    def get_potential_AO(self):
        #
        line_num = []
        for (n, line) in enumerate(self.loglines):
            if ('***** Potential Energy *****' in line):
                line_num.append(n)
        #
        pe_data = np.zeros((self.nao, self.nao), np.float64)
        #
        count = 0
        n = line_num[count]
        shift = 0
        loops = int(self.nao / 5) + 1
        last = int(self.nao % 5)
        for k in range(loops):
            try:
                irange = self.nao - k*5
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
                    elements=self.loglines[n+2+k+shift+i].split()
                    for j in range(end):
                        s = k*5 + j
                        m = i + k*5
                        pe_data[m,s] = float(elements[j+1])
                        if (i != j):
                            pe_data[s,m] = pe_data[m,s]
                shift += irange
            except (IndexError, ValueError):
                break
        #
        return pe_data
    #
    def get_density_AO(self):
        #
        line_num = []
        for (n, line) in enumerate(self.loglines):
            if ('Final density matrix:' in line):
                line_num.append(n)
        #
        dens_data = np.zeros((self.nao, self.nao), np.float64)
        #
        count = 0
        n = line_num[count]
        shift = 0
        loops = int(self.nao / 5) + 1
        last = int(self.nao % 5)
        for k in range(loops):
            try:
                irange = self.nao - k*5
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
                    elements=self.loglines[n+2+k+shift+i].split()
                    for j in range(end):
                        s = k*5 + j
                        m = i + k*5
                        dens_data[m,s] = float(elements[j+1])
                        if (i != j):
                            dens_data[s,m] = dens_data[m,s]
                shift += irange
            except (IndexError, ValueError):
                break
        #
        return dens_data
    #
    def get_ee_twoe_AO(self):
        #
        linenum = []
        read_error = True
        for (n,line) in enumerate(self.loglines):
            if ('*** Eumping Two-Electron integrals ***' in line):
                read_error = False
                linenum.append(n)
        #
        if (read_error == False):
            twoe_AO_4D = np.zeros([self.nao, self.nao, self.nao, self.nao], np.float64)
            count = -1
            n = linenum[count]
            k = 7
            while (k < self.nao**4):
                try:
                    elements = self.loglines[n+k].split()
                    # the two-electron integrals are in Mulliken/chemist's notation: [uv|ls] (cf. Szabo, Ostlund: Table 2.2)
                    u = int(float(elements[1])) - 1 # in the AO basis index notation, \mu
                    v = int(float(elements[3])) - 1 # \nu
                    l = int(float(elements[5])) - 1 # \lambda
                    s = int(float(elements[7])) - 1 # \sigma
                    uvls = float(elements[9])
                    twoe_AO_4D[u,v,l,s] = uvls
                    # 8-fold permutation symmetry of 2-e integrals (w/ real basis fns)
                    twoe_AO_4D[l,s,u,v] = twoe_AO_4D[u,v,l,s]
                    twoe_AO_4D[v,u,s,l] = twoe_AO_4D[u,v,l,s]
                    twoe_AO_4D[s,l,v,u] = twoe_AO_4D[u,v,l,s]
                    twoe_AO_4D[v,u,l,s] = twoe_AO_4D[u,v,l,s]
                    twoe_AO_4D[s,l,u,v] = twoe_AO_4D[u,v,l,s]
                    twoe_AO_4D[u,v,s,l] = twoe_AO_4D[u,v,l,s]
                    twoe_AO_4D[l,s,v,u] = twoe_AO_4D[u,v,l,s]
                    k += 1
                except (IndexError, ValueError, TypeError, NameError):
                    break
            return twoe_AO_4D
        else:
            print('Two-electron integrals not found.')
            return
    #
    def get_ee_onee_AO(self, dens, ee_twoe, exchange=True):
        #
        assert dens.shape[0] == dens.shape[1], 'Density matrix (problem with axes 0 and 1, all axis-dimensions must be the same!)'
        assert ee_twoe.shape[0] == ee_twoe.shape[1], 'ERIs (problem with axes 0 and 1, all axis-dimensions must be the same!)'
        assert ee_twoe.shape[2] == ee_twoe.shape[3], 'ERIs (problem with axes 2 and 3, all axis-dimensions must be the same!)'
        assert ee_twoe.shape[0] == ee_twoe.shape[2], 'ERIs (problem with axes 0 and 2, all axis-dimensions must be the same!)'
        e = True
        if (dens.shape[0] == ee_twoe.shape[0]):
            nbas = dens.shape[0]
            vee_data = np.zeros((nbas, nbas), np.float64)
            e = False
            if (exchange == True):
                for u in range(nbas):
                    for v in range(u,nbas):
                        for l in range(nbas):
                            for s in range(nbas):
                                # coulomb - 0.5*exchange
                                vee_data[u,v] += dens[l,s]*(ee_twoe[u,v,l,s])
                                vee_data[u,v] -= dens[l,s]*(0.5*ee_twoe[u,l,v,s])
                        vee_data[v,u] = np.conjugate(vee_data[u,v])
            elif (exchange == False):
                for u in range(nbas):
                    for v in range(u,nbas):
                        for l in range(nbas):
                            for s in range(nbas):
                                # coulomb
                                vee_data[u,v] += dens[l,s]*(ee_twoe[u,v,l,s])
                        vee_data[v,u] = np.conjugate(vee_data[u,v])
            return vee_data
        elif (e == True):
            print('\nError: Shapes of density and ERI tensors are not compatible.')
            return
    #
    def get_core_AO(self):
        #
        core = self.get_kinetic_AO() + self.get_potential_AO()
        return core

def print_info(logic):
    if (logic == True):
        print('|======================gauss_hf.py=====================|')
        print('|   Class:                                             |')
        print('|******************************************************|')
        print('|   log_data("file.log"):                              |')
        print('|       reads text from "file.log" and extracts data   |')
        print('|       accessible via the following methods.          |')
        print('|       returns:                                       |')
        print('|           nao, # of AO basis fns;                    |')
        print('|           n_a, # of alpha electrons;                 |')
        print('|           n_b, # of beta electrons.                  |')
        print('|******************************************************|')
        print('|   Methods:                                           |')
        print('|******************************************************|')
        print('|   get_overlap_AO():                                  |')
        print('|       returns the overlap matrix, S, in AO basis.    |')
        print('|******************************************************|')
        print('|   get_kinetic_AO():                                  |')
        print('|       returns the kinetic energy operator matrix,    |')
        print('|       in AO basis.                                   |')
        print('|******************************************************|')
        print('|   get_potential_AO():                                |')
        print('|       returns the electron-nuclear potential energy  |')
        print('|       operator matrix, in AO basis.                  |')
        print('|******************************************************|')
        print('|   get_core_AO():                                     |')
        print('|       returns the core Hamiltonian operator matrix,  |')
        print('|       in AO basis.                                   |')
        print('|******************************************************|')
        print('|   get_density_AO():                                  |')
        print('|       returns the density matrix, in AO basis.       |')
        print('|******************************************************|')
        print('|   get_ee_twoe_AO():                                  |')
        print('|       returns the electron-electron repulsion inte-  |')
        print('|       grals, in a rank-4 tensor form, in AO basis.   |')
        print('|******************************************************|')
        print('|   get_ee_onee_AO(dens=dens_bas, ee_twoe=ee_twoe_bas, |')
        print('|                  exchange=True):                     |')
        print('|       returns an effective one-electron matrix for   |')
        print('|       electron-electron interaction, evaluated using |')
        print('|       the ERIs stored in the 4-rank tensor           |')
        print('|       "ee_twoe_bas" and the density matrix "dens_bas"|')
        print('|       (same basis set is assumed for both). If       |')
        print('|       "exchange" set to "True", includes the exchange|')
        print('|       integrals in the evaluation, otherwise Coulomb-|')
        print('|       only.                                          |')
        print('|======================================================|')
    else:
        print('call as print_info(True) to print info about this module')
    return

if (__name__ == '__main__'):
    print('|======================================================|')
    print('|Use "gauss_hf.py" as a python module:                 |')
    print('|>>> import gauss_hf                                   |')
    print('|===================gauss_hf.py========================|')
    print('|     Use "gauss_hf.py" ONLY if interested in the      |')
    print('|     information extracted by this script, from       |')
    print('| a Gaussian .LOG file of a Hartree-Fock calculation.  |')
    print('|======================================================|')
    print('| IMPORTANT NOTE:                                      |')
    print('|======================================================|')
    print('| It is assumed the following keywords and IOp flags   |')
    print('| are specified in the route section of Gaussian       |')
    print('| input:                                               |')
    print('|******************************************************|')
    print('|   #P ... scf(tight,maxcyc=1500,conventional)     &   |')
    print('|     iop(3/33=6) extralinks(l316,l308) noraff     &   |')
    print('|     symm=noint iop(3/33=3) pop(full)             &   |')
    print('|     iop(6/8=1,3/36=1,4/33=3,5/33=3,6/33=3,9/33=3)    |')
    print('|======================================================|')
    print_info(True)