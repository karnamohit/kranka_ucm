#!/usr/bin/env python

"""Extracts matrices from a Gaussian .LOG file. Provides the class "log_data".

log_data reads the .LOG file and stores lines, which its methods can use to 
extract the AO basis matrices. Also evaluates an effective one-elecron matrix
for the electron-electron Coulombic interaction operator.
"""

import sys
import os
import subprocess
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

def find_linenum(s,txtfile): # locate matching string's line number(s) using `grep`
    strg = subprocess.Popen(['grep','-n', s, txtfile], stdout=subprocess.PIPE).stdout.read().decode("utf-8")
    lst = strg.split("\n")[:-1]
    lst = [int(lst[i].split(":")[0]) for i in range(len(lst))]
    return lst

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
        print('\nReading data from {}...'.format(self.logfile))
        self.loglines = log_file.readlines()
        log_file.close()
        self.loglines = [re.sub(r'D','E', s) for s in self.loglines]
        #
        gauss_log = True
        dum = 0
        try:
            # read no. of basis fns and electrons
            lst1 = find_linenum('primitive gaussians',logfile)
            for n in lst1[0]:
                elements = self.loglines[n].split()
                elements2 = self.loglines[n+1].split()
                self.nao = int(float(elements[0]))
                self.n_a = int(float(elements2[0]))
                self.n_b = int(float(elements2[3]))
                dum = 1
            # read atomic information
            lst2 = find_linenum('NAtoms',logfile)
            for n in lst2[0]:
                try:
                    elements = self.loglines[n].split()
                    self.NAtoms = int(float(elements[1]))  # total number of atoms in the system
                except ValueError:
                    elements = self.loglines[n].split('=')
                    self.NAtoms = int(float(elements[1]))
        except:
            lst1, lst2 = [], []
            for (n,line) in enumerate(self.loglines):
                if ('primitive gaussians' in line):
                    lst1.append(n)
                    elements = self.loglines[n].split()
                    elements2 = self.loglines[n+1].split()
                    self.nao = int(float(elements[0]))
                    self.n_a = int(float(elements2[0]))
                    self.n_b = int(float(elements2[3]))
                    dum = 1
                if ('NAtoms' in line):
                    lst2.append(n)
                    try:
                        elements = self.loglines[n].split()
                        self.NAtoms = int(float(elements[1]))
                    except ValueError:
                        elements = self.loglines[n].split('=')
                        self.NAtoms = int(float(elements[1]))
        finally:
            if (len(lst1) == 0 or len(lst2) == 0):
                print('''
                \tlooking for basis set, electronic, and atomic info, but not found in file!
                ''')
                pass
        #
        if (dum == 0):
            gauss_log = False
        #
        if (gauss_log == False):
            s='''
            \tText file may not be as expected (expecting Gaussian .LOG file):
            \t"n_a", "n_b", "nao", "NAtoms" might NOT be available.
            \tFurther, only "find_linenum()", "get_ee_onee_AO()" functions and
            \t"get_matrix_lowtri_AO()" method may be accessible without errors. 
            \tUse "help()" method to get more information about this module.\n
            '''
            print(s)
        elif (gauss_log == True):
            print('\tGaussian .LOG data read.\n')
        #
        return
    #
    def help(self):
        #
        print_info(True)
        #
        return
    #
    def get_molecule(self):
        NAtoms = self.NAtoms
        coords = np.zeros([NAtoms,3], np.float64)
        atom_info = np.zeros([NAtoms,3], np.float64)
        try:
            lst = find_linenum(' Standard basis:',self.logfile)
            for n in lst[0]:
                for i in range(NAtoms):
                    elements = self.loglines[n + i + 6].split()
                    atom_info[i,0] = int(float(elements[0]))  # atomic center number
                    atom_info[i,1] = int(float(elements[1]))  # atomic number
                    atom_info[i,2] = float(elements[2])  # atomic type
                    # coordinate values in angstroms
                    coords[i,0] = float(elements[3])  # atomic center along x-axis
                    coords[i,1] = float(elements[4])  # atomic center along y-axis
                    coords[i,2] = float(elements[5])  # atomic center along z-axis
        except:
            lst = []
            for (n, line) in enumerate(self.loglines):
                if (' Standard basis:' in line):
                    lst.append(n)
                    for i in range(NAtoms):
                        elements = self.loglines[n + i + 6].split()
                        atom_info[i,0] = int(float(elements[0]))
                        atom_info[i,1] = int(float(elements[1]))
                        atom_info[i,2] = float(elements[2])
                        coords[i,0] = float(elements[3])
                        coords[i,1] = float(elements[4])
                        coords[i,2] = float(elements[5])
        finally:
            if (len(lst) == 0):
                print('looking for molecular info but not found in file. Exiting...')
                return
        #
        return atom_info, coords
    #
    def get_matrix_lowtri_AO(self, string, nbasis, skip, columns, imaginary=False, Hermitian=True, start_inplace=False, n_0=0):
        #
        if start_inplace:
            if type(n_0) != int:
                print('proper value for "n_0" not provided. Setting n_0=0 (default value)')
                n_0 = 0
            elif n_0 == 0:
                print('warning: detected n_0=0. might want to specify different value of "n_0" if start_inplace=True.')
        #
        nbasis, skip, columns = int(nbasis), int(skip), int(columns)
        #
        try:
            lst = find_linenum(string,self.logfile)
            line_num = [lst[i] + n_0 for i in range(len(lst))]
        except:
            line_num = []
            for (n, line) in enumerate(self.loglines[n_0:]):
                if (string in line):
                    line_num.append(n+n_0)
        finally:
            if (len(line_num) == 0):
                print('looking for "{}" but not found in file. Exiting...'.format(string))
                return
        #
        data = np.zeros((nbasis, nbasis), np.float64)
        #
        count = 0
        n = line_num[count]
        m = n + skip
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
                    m += 1
                shift += irange
            except (IndexError, ValueError):
                break
        #
        if (start_inplace == True):
            return data, m
        else:
            return data
    #
    def get_MOcoeffs_AO(self):
        # read SCF MO coefficients as column-vectors
        NAOs = self.nao
        cMO = np.zeros([NAOs,NAOs], np.float64)
        #
        try:
            line_num = find_linenum('Alpha MOs:',self.logfile)
        except:
            line_num = []
            for (n, line) in enumerate(self.loglines):
                if ('Alpha MOs:' in line):
                    line_num.append(n)
        finally:
            if (len(line_num) == 0):
                print('looking for "{}" but not found in file. Exiting...'.format('Alpha MOs:'))
                return
        #
        count = -1
        nline = line_num[count]
        loops = int(NAOs / 5) + 1
        last = NAOs % 5
        for k in range(loops):
            for i in range(NAOs):
                try:
                    if (k == (loops - 1)):
                        end = last
                    else:
                        end = 5
                    dum1 = nline+3+k*(2+NAOs)+i
                    elements=self.loglines[dum1].split()
                    for j in range(end):
                        s = k*5 + j
                        m = i
                        cMO[m,s] = float(elements[j+1])
                except (IndexError, ValueError):
                    break
        #
        return cMO
    #
    def get_CASSCF_MOcoeffs_AO(self):
        #
        NAOs = self.nao
        cMO_CI = np.zeros([NAOs,NAOs], np.float64)
        #
        try:
            line_num = find_linenum('FINAL COEFFICIENT MATRIX',self.logfile)
        except:
            line_num = []
            for (n, line) in enumerate(self.loglines):
                if ('FINAL COEFFICIENT MATRIX' in line):
                    line_num.append(n)
        finally:
            if (len(line_num) == 0):
                print('looking for "{}" but not found in file. Exiting...'.format('FINAL COEFFICIENT MATRIX'))
                return
        #
        count = -1
        nline = line_num[count]
        loops = int(NAOs / 10) + 1
        last = NAOs % 10
        for i in range(NAOs):
            for k in range(loops):
                try:
                    if (k == (loops - 1)):
                        end = last
                    else:
                        end = 10
                    dum1 = nline+1+i*(1+loops)+(k+1)
                    elements=self.loglines[dum1].split()
                    for j in range(end):
                        s = k*10 + j
                        m = i
                        cMO_CI[s,m] = float(elements[j])
                except (IndexError, ValueError):
                    break
        return cMO_CI
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
        # getting the reference alpha spin-density matrix
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
        try:
            line_num = find_linenum('*** Dumping Two-Electron integrals ***',self.logfile)
        except:
            line_num = []
            for (n,line) in enumerate(self.loglines):
                if ('*** Eumping Two-Electron integrals ***' in line):
                    line_num.append(n)
        finally:
            if (len(line_num) == 0):
                print('Two-electron integrals not found.')
                return
        #
        twoe_AO_4D = np.zeros([self.nao, self.nao, self.nao, self.nao], np.float64)
        count = -1
        n = line_num[count]
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
        #
        return twoe_AO_4D
    #
    def get_core_AO(self):
        #
        string = '****** Core Hamiltonian ******'
        core_data = self.get_matrix_lowtri_AO(string, self.nao, 2, 5)
        #
        return core_data

def get_ee_onee_AO(dens, ee_twoe, exchange=True, rhf=True):
    #
    try:
        assert len(dens.shape) == 2, '\tDensity matrix must be a 2-index object'
        assert len(ee_twoe.shape) == 4, '\tERI tensor must be a 4-index object'
        assert dens.shape[0] == dens.shape[1], '\tDensity matrix: problem with axis-dimensions; they must all be the same!\n'
        for i in range(4):
            for j in range(i+1,4):
                s = '''
                \tERI tensor: problem with axes {} and {} (axis-dimensions must all be the same!)\n
                '''
                assert ee_twoe.shape[i] == ee_twoe.shape[j], s.format(i,j)
    except AssertionError as e:
        print(e)
        return
    # the no-4-FOR-loops way: thanks to Dr. Harish Bhat (https://github.com/hbhat4000)
    nbas = dens.shape[0]
    vee_data = np.zeros((nbas, nbas), np.float64)
    ee_twoe_coul = ee_twoe.reshape((nbas, nbas, nbas**2))
    if (exchange == True):
        ee_twoe_ex = np.swapaxes(ee_twoe, 1, 2).reshape((nbas, nbas, nbas**2))
        for u in range(nbas):
            for v in range(u,nbas):
                # coulomb - 0.5*exchange
                # for l in range(nbas):
                #     for s in range(nbas):
                #         vee_data[u,v] += 2*dens[l,s]*(ee_twoe[u,v,l,s])
                #         vee_data[u,v] -= 2*dens[l,s]*(0.5*ee_twoe[u,l,v,s])
                vee_data[u,v] = 2 * ee_twoe_coul[u,v,:] @ dens.reshape((-1))
                vee_data[u,v] -= 2 * 0.5 * ee_twoe_ex[u,v,:] @ dens.reshape((-1))
                vee_data[v,u] = np.conjugate(vee_data[u,v])
    elif (exchange == False):
        for u in range(nbas):
            for v in range(u,nbas):
                # coulomb
                # for l in range(nbas):
                #     for s in range(nbas):
                #         vee_data[u,v] += 2*dens[l,s]*(ee_twoe[u,v,l,s])
                vee_data[u,v] = 2 * ee_twoe_coul[u,v,:] @ dens.reshape((-1))
                vee_data[v,u] = np.conjugate(vee_data[u,v])
    #
    return vee_data

def print_info(logic):
    if logic:
        s='''
        \t|=====================gauss_hf.py======================|
        \t|   Class:                                             |
        \t|******************************************************|
        \t|   log_data("/path/to/file.log"):                     |
        \t|       reads text from "file.log" and extracts data   |
        \t|       accessible via the methods listed below;       |
        \t|       expects "file.log" to be a Gaussian .LOG file. |
        \t|******************************************************|
        \t|   Instance variables:                                |
        \t|******************************************************|
        \t|   logfile: path+name of the .LOG file                |
        \t|   loglines: text, read through readlines() method    |
        \t|   nao: # of AO basis fns                             |
        \t|   n_a: # of alpha electrons                          |
        \t|   n_b: # of beta electrons                           |
        \t|******************************************************|
        \t|   Methods:                                           |
        \t|******************************************************|
        \t|   get_molecule():                                    |
        \t|       returns two 2-D arrays containing information  |
        \t|       about molecular geometry. First array contains |
        \t|       atomic character info; second array contains   |
        \t|       Cartesian coordinates, in Å, of the            |
        \t|       corresponding atoms.                           |
        \t|******************************************************|
        \t|   get_MOcoeffs_AO():                                 |
        \t|       returns the matrix of alpha MO coefficients    |
        \t|       (as column-vectors) in AO basis.               |
        \t|******************************************************|
        \t|   get_CASSCF_MOcoeffs_AO():                          |
        \t|       returns the matrix of alpha MO coefficients    |
        \t|       (as column-vectors), for a CASSCF calculation, |
        \t|       in AO basis.                                   |
        \t|******************************************************|
        \t|   get_overlap_AO():                                  |
        \t|       returns the overlap matrix, S, in AO basis.    |
        \t|******************************************************|
        \t|   get_kinetic_AO():                                  |
        \t|       returns the kinetic energy operator matrix,    |
        \t|       in AO basis.                                   |
        \t|******************************************************|
        \t|   get_potential_AO():                                |
        \t|       returns the electron-nuclear potential energy  |
        \t|       operator matrix, in AO basis.                  |
        \t|******************************************************|
        \t|   get_core_AO():                                     |
        \t|       returns the core Hamiltonian operator matrix,  |
        \t|       in AO basis.                                   |
        \t|******************************************************|
        \t|   get_density_AO():                                  |
        \t|       returns the alpha spin-density matrix, in AO   |
        \t|       basis. Assumes restricted reference.           |
        \t|******************************************************|
        \t|   get_dipole_x_AO():                                 |
        \t|       returns the electric dipole moment matrix for  |
        \t|       Cartesian dipole operator "x" in AO basis.     |
        \t|******************************************************|
        \t|   get_dipole_y_AO():                                 |
        \t|       returns the electric dipole moment matrix for  |
        \t|       Cartesian dipole operator "y" in AO basis.     |
        \t|******************************************************|
        \t|   get_dipole_z_AO():                                 |
        \t|       returns the electric dipole moment matrix for  |
        \t|       Cartesian dipole operator "z" in AO basis.     |
        \t|******************************************************|
        \t|   get_ee_twoe_AO():                                  |
        \t|       returns the electron-electron repulsion inte-  |
        \t|       grals, in a rank-4 tensor form, in AO basis.   |
        \t|******************************************************|
        \t|   Functions:                                         |
        \t|******************************************************|
        \t|   find_linenum(string, filename):                    |
        \t|       returns a list of line-numbers (0-indexed),    |
        \t|       that have text matching <string>, for          |
        \t|       <filename>. This function uses `grep` so using |
        \t|       this script on a Linux machine is preferred.   |
        \t|******************************************************|
        \t|   get_ee_onee_AO(dens_bas, ee_twoe_bas,              |
        \t|                  exchange=True, , rhf=True):         |
        \t|       returns an effective one-electron matrix for   |
        \t|       electron-electron interaction, evaluated using |
        \t|       the ERIs stored in the 4-rank tensor           |
        \t|       "ee_twoe_bas" and the density matrix "dens_bas"|
        \t|       (same basis set is assumed for both). If       |
        \t|       "exchange" set to "True", includes the exchange|
        \t|       integrals in the evaluation, otherwise Coulomb-|
        \t|       only. Needs alpha spin density, assumes        |
        \t|       restricted reference.                          |
        \t|******************************************************|
        \t|   get_matrix_lowertri_AO(string, nbasis, skip,       |
        \t|                          columns, imaginary=False,   |
        \t|                          Hermitian=True,             |
        \t|                          start_inplace=False,        |
        \t|                          n_0=None):                  |
        \t|       reads a lower triangular matrix in AO basis.   |
        \t|       NOTE: mainly for internal use, unless familiar |
        \t|======================================================|\n
        '''
        print(s)
    else:
        print('''
        \tCall as `print_info(True)` to print information about this
        \tmodule.\n
        ''')
    return

if (__name__ == '__main__'):
    s='''
    \t|======================================================|
    \t| Using "gauss_hf.py" as a python module:              |
    \t| 1. Copy to the same folder as your .py/.ipynb file.  |
    \t|                        OR                            |
    \t|    Add the path to "gauss_hf.py" to your system path:|
    \t|    >>> import sys                                    |
    \t|    >>> sys.path.append('/path/to/file/')             |
    \t|                                                      |
    \t| 2. Try either of these:                              |
    \t|    >>> import gauss_hf                               |
    \t|    >>> from gauss_hf import *                        |
    \t|====================gauss_hf.py=======================|
    \t|     Use "gauss_hf.py" ONLY if interested in the      |
    \t|     information extracted by this script, from       |
    \t| a Gaussian .LOG file of a Hartree-Fock calculation.  |
    \t|======================================================|
    \t| Python libraries required:                           |
    \t|     numpy: may need to install separately.           |
    \t|     re; os; sys, subprocess: usually included with   |
    \t|            standard Python installation.             |
    \t|======================================================|
    \t| IMPORTANT NOTE:                                      |
    \t|******************************************************|
    \t| It is assumed the following keywords and IOp flags   |
    \t| are specified in the route section of Gaussian       |
    \t| input:                                               |
    \t|******************************************************|
    \t|   #P ... scf(tight,maxcyc=1500,conventional)     &   |
    \t|     iop(3/33=6) extralinks(l316,l308) noraff     &   |
    \t|     symm=noint iop(3/33=3) pop(full)             &   |
    \t|     iop(6/8=1,3/36=1,4/33=3,5/33=3,6/33=3,9/33=3)    |
    \t|======================================================|\n
    \t|                                                      |\n
    '''
    print(s)
    print_info(True)