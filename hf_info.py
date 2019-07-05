import sys
import os
import shutil
import subprocess
import numpy as np

class el_info:
    
    def __init__(self,NELECT, NMOs, fchkfile, logfile):
        
        print 'getting electronic info from the Gaussian log file... \n'
        
        self.nel = NELECT
        self.nmo = NMOs
        self.fchkname = fchkfile
        self.logname = logfile
        
        file1 = open(self.fchkname,'r')
        self.fchklines = file1.readlines()
        file1.close()
        
        file1 = self.logname
        log = 'log.tmp'
        file2 = open(log, 'w')
        sub = subprocess.call(['sed', 's/D/E/g', file1], stdout=file2)
        file2.close()
        file2 = open(log,'r')
        self.loglines = file2.readlines()
        file2.close()
    
    def overlapAO(self):
        
        SAO = np.zeros([self.nmo,self.nmo], np.float64)
        
        if (self.nmo < 6):
            for (n,line) in enumerate(self.loglines):
                if ('*** Overlap ***' in line):
                    for i in range(self.nmo):
                        j = 0
                        while (j <= i):
                            elements = self.loglines[n+i+2].split()
                            k = j + 1
                            SAO[i,j] = float(elements[k])
                            if (i != j):
                                SAO[j,i] = SAO[i,j]
                            j += 1
        else:
            loops = self.nmo / 5 + 1
            last = self.nmo % 5
            shift = 0
            for (n,line) in enumerate(self.loglines):
                if ('*** Overlap ***' in line):
                    for k in range(loops):
                        irange = self.nmo - k * 5
                        for i in range(irange):
                            if (k == loops):
                                if (i < last):
                                    end = i + 1
                                else:
                                    end = last
                            else:
                                if (i <= 4):
                                    end = i + 1
                                else:
                                    end = 5
                            for j in range(end):
                                elements = self.loglines[n+shift+k+i+2].split()
                                s = k * 5 + j
                                m = i + k * 5
                                SAO[m,s] = float(elements[j+1])
                                if (i != j):
                                    SAO[s,m] = SAO[m,s]
                        shift += irange
        
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
        
        alphaMO = []
        cMO = np.zeros([self.nmo,self.nmo], np.float64)
        
        for (n,line) in enumerate(self.fchklines):
            if ('Alpha MO coefficients' in line):
                loops = self.nmo**2 / 5 + 1
                k = 0
                for i in range(loops):
                    elements = self.fchklines[n+1+k].split()
                    alphaMO.extend(elements)
                    k += 1
        
        for i in range(len(alphaMO)):
            alphaMO[i] = float(alphaMO[i])
        
        for i in range(self.nmo):
            for j in range(self.nmo):
                cMO[i,j] = alphaMO[i*self.nmo+j]
        
        cMO = cMO.T
        
        return cMO
    
    def densityMO(self):
        
        densAO = self.densityAO()
        MO = self.coefficientsMOalpha()
        SAO = self.overlapAO()
        
        rhoMO = MO.T.dot(SAO.dot(densAO.dot(SAO.dot(MO))))
        
        return rhoMO



if (__name__ == '__main__'):
    print 'Use this file as a library. Use hf_dip.py if only interested in the information extracted by this script. \n'

