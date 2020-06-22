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
    def __init__(self, logfile=None, nonlog_error_msg=True):
        #
        if (logfile == None):
            print('\nNo file specified!')
            print('Exiting...')
            return
        else:
            # check if logfile exists
            if os.path.isfile(logfile) == True:
                self.logfile = logfile
            elif os.path.isfile(logfile) == False:
                print('\nError: Specified file ({}) does not exist.'.format(logfile))
                print('Exiting...')
                return
        #
        # reading "logfile" contents
        log_file = open(self.logfile, 'r')
        print('\nReading data from {}...\n'.format(self.logfile))
        self.loglines = log_file.readlines()
        log_file.close()
        self.loglines = [re.sub(r'D','E', s) for s in self.loglines]
        #
        gauss_log = True
        for (n,line) in enumerate(self.loglines):
            # read no. of basis fns and electrons
            if ('primitive gaussians' in line):
                elements = self.loglines[n].split()
                elements2 = self.loglines[n+1].split()
                self.nao = int(float(elements[0]))
                self.n_a = int(float(elements2[0]))
                self.n_b = int(float(elements2[3]))
                #
            else:
                gauss_log = False
        #
        if (gauss_log == False):            
            print('\nText file may not as expected (expecting Gaussian .LOG file).\n')
            print('  "n_a", "n_b", "nao" instance variables will NOT be available.')
            print('  Further, only "get_matrix_lowtri_AO()" and "get_ee_onee_AO()"')
            print('  methods may be accessible without errors. Use "help()" method')
            print('  to get more information about this module.\n')
        elif (gauss_log == True):
            print('\nGaussian .LOG data read.\n')
        #
        return
    #
    def help(self):
        #
        print_info(True)
        #
        return
    #
    def get_matrix_lowtri_AO(self, string, nbasis, skip, columns, imaginary=False, Hermitian=True, start_inplace=False, n_0=None):
        #
        if (start_inplace == False):
            n_0 = 0
        elif (start_inplace == True):
            assert n_0 != None, 'feed an integer value for "n_0".'
            n_0 = int(n_0)
        #
        nbasis, skip, columns = int(nbasis), int(skip), int(columns)
        #
        line_num = []
        for (n, line) in enumerate(self.loglines[n_0:]):
            if (string in line):
                line_num.append(n)
        #
        data = np.zeros((nbasis, nbasis), np.float64)
        #
        count = 0
        n = line_num[count]
        shift = 0
        #
        loops = int(nbasis / columns) + 1
        last = int(nbasis % columns)
        for k in range(loops):
            try:
                irange = nbasis - k*columns
                for i in range(irange):
                    if (k == (loops - 1)):
                        if (i < last):
                            end = i + 1
                        else:
                            end = last
                    else:
                        if (i <= (columns - 1)):
                            end = i + 1
                        else:
                            end = columns
                    elements=self.loglines[n+skip+k+shift+i].split()
                    for j in range(end):
                        s = k*columns + j
                        m = i + k*columns
                        data[m,s] = float(elements[j+1])
                        if (i != j):
                            if (imaginary == True and Hermitian == True):
                                data[s,m] = -1*data[m,s]
                            else:
                                data[s,m] = data[m,s]
                shift += irange
            except (IndexError, ValueError):
                break
        #
        if (start_inplace == True):
            return data, n
        else:
            return data
    #
    def get_overlap_AO(self):
        #
        string = '*** Overlap ***'
        overlap_data = self.get_matrix_lowtri_AO(string, self.nao, 2, 5)
        #
        return overlap_data
    #
    def get_kinetic_AO(self):
        #
        string = '*** Kinetic Energy ***'
        ke_data = self.get_matrix_lowtri_AO(string, self.nao, 2, 5)
        #
        return ke_data
    #
    def get_potential_AO(self):
        #
        string = '***** Potential Energy *****'
        pe_data = self.get_matrix_lowtri_AO(string, self.nao, 2, 5)
        #
        return pe_data
    #
    def get_density_AO(self):
        #
        string = 'Final density matrix:'
        dens_data = self.get_matrix_lowtri_AO(string, self.nao, 2, 5)
        #
        return dens_data
    #
    def get_dipole_x_AO(self):
        #
        string = 'Multipole matrices IBuc=  518 IX=    1'
        dipAO_data = self.get_matrix_lowtri_AO(string, self.nao, 2, 5)
        #
        return dipAO_data
    #
    def get_dipole_y_AO(self):
        #
        string = 'Multipole matrices IBuc=  518 IX=    2'
        dipAO_data = self.get_matrix_lowtri_AO(string, self.nao, 2, 5)
        #
        return dipAO_data
    #
    def get_dipole_z_AO(self):
        #
        string = 'Multipole matrices IBuc=  518 IX=    3'
        dipAO_data = self.get_matrix_lowtri_AO(string, self.nao, 2, 5)
        #
        return dipAO_data
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
        assert len(dens.shape) == 2
        assert len(ee_twoe.shape) == 4
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
        string = '****** Core Hamiltonian ******'
        core_data = self.get_matrix_lowtri_AO(string, self.nao, 2, 5)
        #
        return core_data

def print_info(logic):
    if (logic == True):
        print('|=====================gauss_hf.py======================|')
        print('|   Class:                                             |')
        print('|******************************************************|')
        print('|   log_data("/path/to/file.log"):                     |')
        print('|       reads text from "file.log" and extracts data   |')
        print('|       accessible via the methods listed below;       |')
        print('|       expects "file.log" to be a Gaussian .LOG file. |')
        print('|       instance variables:                            |')
        print('|           logfile: path+name of the .LOG file        |')
        print('|           loglines: text, through readlines() method |')
        print('|           nao: # of AO basis fns                     |')
        print('|           n_a: # of alpha electrons                  |')
        print('|           n_b: # of beta electrons                   |')
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
        print('|   get_dipole_x_AO():                                 |')
        print('|       returns the electric dipole moment matrix for  |')
        print('|       Cartesian dipole operator "x" in AO basis.     |')
        print('|******************************************************|')
        print('|   get_dipole_y_AO():                                 |')
        print('|       returns the electric dipole moment matrix for  |')
        print('|       Cartesian dipole operator "y" in AO basis.     |')
        print('|******************************************************|')
        print('|   get_dipole_z_AO():                                 |')
        print('|       returns the electric dipole moment matrix for  |')
        print('|       Cartesian dipole operator "z" in AO basis.     |')
        print('|******************************************************|')
        print('|   get_ee_twoe_AO():                                  |')
        print('|       returns the electron-electron repulsion inte-  |')
        print('|       grals, in a rank-4 tensor form, in AO basis.   |')
        print('|******************************************************|')
        print('|   get_ee_onee_AO(dens_bas, ee_twoe_bas,              |')
        print('|                  exchange=True):                     |')
        print('|       returns an effective one-electron matrix for   |')
        print('|       electron-electron interaction, evaluated using |')
        print('|       the ERIs stored in the 4-rank tensor           |')
        print('|       "ee_twoe_bas" and the density matrix "dens_bas"|')
        print('|       (same basis set is assumed for both). If       |')
        print('|       "exchange" set to "True", includes the exchange|')
        print('|       integrals in the evaluation, otherwise Coulomb-|')
        print('|       only.                                          |')
        print('|******************************************************|')
        print('|   get_matrix_lowertri_AO(string, nbasis, skip,       |')
        print('|                          columns, imaginary=False,   |')
        print('|                          Hermitian=True,             |')
        print('|                          start_inplace=False,        |')
        print('|                          n_0=None):                  |')
        print('|       reads a lower triangular matrix in AO basis.   |')
        print('|       NOTE: mainly for internal use, unless familiar |')
        print('|======================================================|')
    else:
        print('\tCall as print_info(True) to print info about this module.')
    return

if (__name__ == '__main__'):
    print('|======================================================|')
    print('|Use "gauss_hf.py" as a python module:                 |')
    print('|>>> import gauss_hf                                   |')
    print('|>>> from gauss_hf import *                            |')
    print('|====================gauss_hf.py=======================|')
    print('|     Use "gauss_hf.py" ONLY if interested in the      |')
    print('|     information extracted by this script, from       |')
    print('| a Gaussian .LOG file of a Hartree-Fock calculation.  |')
    print('|======================================================|')
    print('| IMPORTANT NOTE:                                      |')
    print('|******************************************************|')
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