#!/usr/bin/env python

import numpy as np
import os
import time
import sys
# from functools import partial 
from multiprocessing import Pool
sys.path.append('/mnt/c/Users/karna/GitHub/kranka_ucm/scripts/')
sys.path.append('/mnt/c/Users/kranka/GitHub/kranka_ucm/scripts/')
sys.path.append('/mnt/d/GitHub/kranka_ucm/scripts/')
sys.path.append('/home/kranka/proj_isborn/kranka_ucm/scripts/')
from gauss_hf import *
#
# following is where you store the output of the current code
location = './chempaper_data/'; '''
location = './test/trajdata/'
# '''
print('data-files saved at: {}'.format(os.path.dirname(location)))

tdmethod_list = {'hf': 'rt-tdexx', 'h': 'rt-tdh', 'scf': 'tdci'}

pert1 = ['ress0s11cyc','halfress0s11cyc','doubress0s11cyc','ndlaser']
amp1 = ['0.05']
pert_arr = []
for i in pert1:
    for j in amp1:
        pert_arr.append(i+'='+j)

initLog = sys.argv[1]
tdperturb = sys.argv[2]

# print(initLog,tdperturb)

filepath='./logfiles/'
tdfilepath='./td_logfiles/'
ext_tdmethod='.txt'
#
filename = initLog
file = filepath+filename
elements = filename.split('_')
# print(elements); quit()
try:
    if (True in ['field-delta' in j for j in elements]): # detect parameters
        method = elements[0]                            # 'hf', 'h'
        state = elements[2]                             # 's0'
        molecule = elements[3]                          # 'c2h4', 'lih', 'heh+', 'h2'
        basis = elements[4].split('.')[0]               # 'sto-3g', '6-31g'
        perturb = elements[1].split('-')[1]
        tdperturb = perturb
        # print(tdperturb); quit()
    else:
        if ('casscf' in elements[0] or 'hf' in elements[0]):
            method = 'hf'
        else:
            method = elements[0]
        state = elements[1]
        molecule = elements[2]
        basis = elements[3].split('.log')[0]
        perturb = ''
except FileNotFoundError as e:
    print(e)
#
molsys = state+'_'+molecule+'_'+basis
if perturb != '':
    title = perturb+'_'+molsys
elif perturb == '':
    title = molsys
#
# checks whether the given .LOG file exists and prints stuff appropriately
print('file input:\n{}\n'.format(file))
print('\ndetecting whether the specified file exists...\n')
if os.path.isfile(file) == True:
    print('\ngoing through the Gaussian .LOG file:\n{}\n'.format(file))
elif os.path.isfile(file) == False:
    print('\nERROR: the specified file does not exist.')
#
# read the non-TD method .LOG file
data_file = log_data(file,nonlog_error_msg=False)
#
NAOs = data_file.nao
overlap_data = data_file.get_overlap_AO()
ke_data = data_file.get_kinetic_AO()
pe_data = data_file.get_potential_AO()
ee_twoe_data = data_file.get_ee_twoe_AO()
dip_data = [data_file.get_dipole_x_AO(), data_file.get_dipole_y_AO(), data_file.get_dipole_z_AO()]
datafile = location + 'ke+en+overlap+ee_twoe+dip_'+method+'_'+title+'.npz'
#
print('data for T, Ven, S, Vee, electric dipole moment matrices stored in .npz file as:\n',datafile)
np.savez(datafile, ke_data=ke_data, en_data=pe_data, overlap_data=overlap_data, ee_twoe_data=ee_twoe_data, dip_data=dip_data)
#
tdlocation = './td_logfiles/'
td_suffix = tdperturb+'_'+molsys
ext_tdmethod = '.txt'
#
# process TD field and dipole moment matrix data
file_dip = tdfilepath + 'tdip_' + td_suffix + ext_tdmethod
file_fld = tdfilepath + 'efield_' + td_suffix + ext_tdmethod
filetd_1 = log_data(file_dip,nonlog_error_msg=False)
log_lines_1 = filetd_1.loglines
filetd_2 = log_data(file_fld,nonlog_error_msg=False)
log_lines_2 = filetd_2.loglines
#
efield_data = np.array(([[]]), np.float64)
td_dip_data = np.array(([[]]), np.float64)
#
for log_lines in [log_lines_1, log_lines_2]:
    step = 0
    for (n, line) in enumerate(log_lines):
        try:
            if ('Real-Time Eipole Moment' in line):
                elements = log_lines[n].split('|')
                dip = np.array([[float(elements[1]),float(elements[3]),float(elements[5])]])
                if step != 0:
                    td_dip_data = np.append(td_dip_data, dip, axis=0)
                else:
                    td_dip_data = np.append(td_dip_data, dip, axis=1)
                # print(elements,end='\r')
                step += 1
            elif ('Current electric field:' in line):
                elements = log_lines[n+1].split()
                fld = np.array([[float(elements[2]),float(elements[6]),float(elements[10])]])
                if step != 0:
                    efield_data = np.append(efield_data, fld, axis=0)
                else:
                    efield_data = np.append(efield_data, fld, axis=1)
                # print(elements,end='\r')
                step += 1
        except (IndexError, ValueError) as e:
            print(e)
            break
#
filetd = 'td_efield+dipole_'+tdmethod_list[method]+'_'+td_suffix+'.npz'
td_datafile = location + filetd
print('the corresponding .NPZ file has been stored as:\n',td_datafile,'\n')
np.savez(td_datafile, td_efield_data=efield_data, td_dipole_data=td_dip_data)
#
# process TD density matrix data
dens_arr = ['real', 'imag'] # order is IMPORTANT: 'real' must be iterated over before 'imag'
#
td_filename1 = 'dens_real_'+td_suffix
td_filename2 = 'dens_imag_'+td_suffix
td_denreal = tdfilepath+td_filename1+ext_tdmethod
td_denimag = tdfilepath+td_filename2+ext_tdmethod
# checks whether the given .TXT file exists and prints stuff appropriately
print('file input (real):', td_denreal)
print('file input (imag):', td_denimag)
print('detecting whether files exist...')
for td_file in [td_denreal, td_denimag]:
    if os.path.isfile(td_file) == True:
        print('the output of grep (on the Gaussian .LOG file of TD calculation) is being read from:\n',td_file)
    elif os.path.isfile(td_file) == False:
        print('ERROR: file does not exist.')
#
filetd_real = log_data(td_denreal,nonlog_error_msg=False)
log_lines_re = filetd_real.loglines
filetd_imag = log_data(td_denimag,nonlog_error_msg=False)
log_lines_im = filetd_imag.loglines
#
line_num_re = []
for (n, line) in enumerate(log_lines_re):
    if ('Eensity Matrix' in line):
        line_num_re.append(n)
Ntot = len(line_num_re); '''
Ntot = 50
# '''
#
series_data = np.zeros((Ntot, NAOs, NAOs), np.complex128)
# print('series_data set to 0...')
#
# returns corresponding (imaginary and real parts of) density matrices for a given index, from the grepped txt files
def denmat(indx):
    mat_re, _ = filetd_real.get_matrix_lowtri_AO('Density Matrix',NAOs,2,5,imaginary=False,start_inplace=True,n_0=line_num_re[indx])
    mat_im, _ = filetd_imag.get_matrix_lowtri_AO('Density Matrix',NAOs,2,5,imaginary=True,start_inplace=True,n_0=line_num_re[indx])
    return indx, mat_re, mat_im
#
if __name__ == '__main__':
    CPUs = os.cpu_count() # - 1 # "- 1" for personal computers
    pool = Pool(CPUs)
    start_time = time.time()
    with pool:
        try:
            for i in pool.imap(denmat, range(Ntot)):
                series_data[i[0],:,:] = i[1] + 1J*i[2]
        except Exception as e:
            print(e)
        finally:
            pool.close()
            pool.join()
        #
    filetd = 'td_dens_re+im_'+tdmethod_list[method]+'_'+td_suffix+'.npz'; '''
    filetd = 'test2.npz'
    # '''
    td_datafile = location + filetd
    # print(series_data[0,:,:])
    print('the corresponding .NPZ file has been stored as:\n',td_datafile,'\n')
    np.savez(td_datafile, td_dens_im_data=series_data.imag, td_dens_re_data=series_data.real)
    print('time elapsed: {:.2f} s'.format(time.time() - start_time))
