import sys
import os
import shutil
import subprocess
import numpy as np
import re

class el_info:
    
    def __init__(self, logfile, fchkfile=None):
        
        self.fchkname = fchkfile
        self.logname = logfile
        
        print('getting electronic info from the Gaussian files', \
              '"',self.logname,'"','and','"',self.fchkname,'"...')
        
        #print(fchkfile,logfile)
                        
        file1 = open(self.logname,'r')
        self.loglines = file1.readlines()
        file1.close()
        self.loglines = [re.sub(r'D','E', s) for s in self.loglines]
        
        for (n,line) in enumerate(self.loglines):
            if ('primitive gaussians' in line):
                Nbas = int(float(self.loglines[n].split()[0]))
                NELECT = int(float(self.loglines[n+1].split()[0]))
        
        self.nel = NELECT
        self.nmo = Nbas
        
        if (self.fchkname is not None):
            file1 = open(self.fchkname,'r')
            self.fchklines = file1.readlines()
            file1.close()
    
    def overlapAO(self):
        
        SAO = np.zeros([self.nmo,self.nmo], np.float64)
        
        line_num = []
	
        for (n, line) in enumerate(self.loglines):
            if ('*** Overlap ***' in line):
                line_num.append(n)
        
        loops = int(NMOs / 5) + 1
        last = NMOs % 5
        
        count = 0
        n = line_num[count]
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
                    elements=self.loglines[n+2+k+shift+i].split()
                    for j in range(end):
                        s = k*5 + j
                        m = i + k*5
                        SAO[m][s] = float(elements[j+1])
                        if (i != j):
                            SAO[s][m] = SAO[m][s]
                shift += irange
            except (IndexError, ValueError):
                break        
        return SAO
    
    def densityAO(self):
        
        rhoAO = np.zeros([self.nmo,self.nmo], np.float64)
        
        if (self.nmo < 6):
            for (n,line) in enumerate(self.loglines):
                if ('Eensity Matrix:' in line):
                    for i in range(self.nmo):
                        j = 0
                        while (j <= i):
                            elements = self.loglines[n+i+2].split()
                            k = j + 4
                            rhoAO[i,j] = float(elements[k])
                            if (i != j):
                                rhoAO[j,i] = rhoAO[i,j]
                            j += 1
        else:
            loops = self.nmo / 5 + 1
            last = self.nmo % 5
            shift = 0
            for (n,line) in enumerate(self.loglines):
                if ('Eensity Matrix:' in line):
                    for k in range(loops):
                        irange = self.nmo - k * 5
                        for i in range(irange):
                            if (k == loops):
                                if (i < last):
                                    end = i + 1
                                else:
                                    end = last
                            else:
                                if (i < 4):
                                    end = i + 1
                                else:
                                    end = 5
                            for j in range(end):
                                elements = self.loglines[n+shift+k+i+2].split()
                                s = k * 5 + j
                                m = i + k * 5
                                rhoAO[m,s] = float(elements[j+4])
                                if (i != j):
                                    rhoAO[s,m] = rhoAO[m,s]
                        shift += irange
        
        return rhoAO
    
    def coefficientsMOalpha(self):
	    
        # read SCF MO coefficients as column-vectors
        
        NMOs = self.nmo
        cMO = np.zeros([NMOs,NMOs], np.float64)
         
        if (self.fchkname != None):
            alphaMO = []
            for (n,line) in enumerate(self.fchklines):
                if ('Alpha MO coefficients' in line):
                    loops = int(self.nmo**2 / 5 + 1)
                    k = 0
                    for i in range(loops):
                        elements = self.fchklines[n+1+k].split()
                        alphaMO.extend(elements)
                        k += 1
            
            alphaMO = float(alphaMO)
            
            for i in range(self.nmo):
                for j in range(self.nmo):
                    cMO[i,j] = alphaMO[i*self.nmo+j]
                
            cMO = cMO.T
        else:
            line_num = []
            for (n, line) in enumerate(self.loglines):
                if ('Alpha MOs' in line):
                    line_num.append(n)
            count = -1
            nline = line_num[count]
            loops = int(NMOs / 5) + 1
            last = NMOs % 5
            for k in range(loops):
                for i in range(NMOs):
                    try:
                        if (k == (loops - 1)):
                            end = last
                        else:
                            end = 5
                        dum1 = nline+3+k*(2+NMOs)+i
                        elements=self.loglines[dum1].split()
                        for j in range(end):
                            s = k*5 + j
                            m = i
                            cMO[m][s] = float(elements[j+1])
                    except (IndexError, ValueError):
                        break
        
        return cMO
    
    def densityMO_HF(self):
        
        densAO = self.densityAO()
        MO = self.coefficientsMOalpha()[:,:self.nel]
        SAO = self.overlapAO()
        
        rhoMO = MO.T.dot(SAO.dot(densAO.dot(SAO.dot(MO))))
        
        return rhoMO



if (__name__ == '__main__'):
    print('Use this file as a library. Use hf_dip.py if only interested in the information extracted by this script. \n')

