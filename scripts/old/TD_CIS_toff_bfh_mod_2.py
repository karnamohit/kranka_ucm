#!/usr/bin/env python
"""
July 2nd, 2013
@author: C. M. Isborn, Noriyuki Tani
"""

"""
Extract data (coordinates, transition dipoles, oscillator strengths, excitation energies) from a CIS computation. 
For now this only gets the molecular coordinates. We will need the other data, too. 
Eventually we want to add the applied electric field and couple it to the transition dipole moments, 
then propagate the electron density in time. 
"""

#=============================================================================================
# IMPORTS
#=============================================================================================

import sys
import numpy # for numerical work
import os    # For interacting with the operating system
#import matplotlib.pyplot as plt #for plotting  
#=============================================================================================
# PARAMETERS
#=============================================================================================

#LogFile = 'tdm_tdcis.txt'      # Gaussian output file with keyword cis=(nstates=50,alltransitiondensities)
LogFile = 'tdm_tdcis.txt'      # output from calc_tdm.py 
NAtoms = 14                   # Number of atoms in the molecule
NStates = 6                  # Number of excited(?) states (what is required by the code is the TOTAL # of states computed - KR)
NMOs = 3           # Number of MOs in CAS
NSim = 6           # Number of states in simulation

stindx = 0         # Initial state from which to begin propagation (ranges from 0, i.e. GS, to (NStates-1), i.e. last excited state) - KR


emax = 0.0500        # Field strength (in au)
maxstep = 10452     # Number of itegration steps
delta = 0.083      # Stepsize (in au)
w = 0.337          # Frequency of field (in au, hartrees)
#
# envelope is turned off at moment- search for ENVELOPE in this document. 
#     but we do still need to define the envelope, it's just not multiplied times the field
#     this is kind of lazy, but it makes it easy to turn the field on if needed later.
#
envelope = "cosine"             # The envelope that could be created (trapezoidal,gaussian or cosine)
chirp = 0.00                   # chirp (in au)
ncyc = 2                        # number of cycles in envelope
phase = 0 #*numpy.pi/2            # phase of electric field
ton = 0
toff = 452 
#ton = 0                         # step number to turn on field
#toff = 560                      # step number to turn off field
#=============================================================================================
# Read file 
#=============================================================================================

# Open file for reading.
infile = open(LogFile, 'r')

# Read all lines.
lines = infile.readlines()

# Close file.
infile.close()

# Set up empty array for coordinates
CoordsX = numpy.zeros([NAtoms], numpy.float64)
CoordsY = numpy.zeros([NAtoms], numpy.float64)
CoordsZ = numpy.zeros([NAtoms], numpy.float64)

# Get coordinates 
for (n, line) in enumerate(lines):
  if ('Standard orientation:' in line):     # Finds line n that has string key 
    for i in range(NAtoms):                 # Loop over the number of atoms
       elements = lines[n + i + 5].split()  # Skip 5 lines before the data starts
       
       CoordsX[i] = float(elements[3])
       CoordsY[i] = float(elements[4])
       CoordsZ[i] = float(elements[5])

# print to screen
for i in range(NAtoms):                     # Loop over the number of atoms
    print "%d   %f    %f    %f" % (i, CoordsX[i], CoordsY[i], CoordsZ[i])

#Read SCF Energy
ESCF = 0
# Looking for the SCF Energy
for (n,line) in enumerate(lines):
    if ('SCF Done' in line):
        elements = lines[n].split()
        ESCF = elements[4]              

print 'ESCF =',ESCF,' Ha'

# Get CIS energy
Energies = numpy.zeros([NStates + 1], numpy.float64)          #initializing Energy

i = 0
for (n,line) in enumerate(lines):
    if('Excited State ' in line):
        # Filling in the values
        elements = lines[n].split()
        Energies[i+1] = float(elements[4])/27.212
        i = i + 1

print 'Energies (eV)'
print Energies

#Read transition dipoles

# Creating the first row of info (info(1))
#Ground_to_e_state = numpy.zeros([NStates])     
Ground_to_e_X = numpy.zeros([NStates])
Ground_to_e_Y = numpy.zeros([NStates])
Ground_to_e_Z = numpy.zeros([NStates])
#Ground_to_e_osc = numpy.zeros([NStates])       oscillation

for (n,line) in enumerate(lines):
    if('Ground to excited state transition electric dipole moments' in line):
        for i in range(NStates-1):
            elements = lines[n+i+2].split()
            # putting values into array
            #Ground_to_e_state[i] = float(elements[0])
            Ground_to_e_X[i] = float(elements[1])
            Ground_to_e_Y[i] = float(elements[2])
            Ground_to_e_Z[i] = float(elements[3])
            #Ground_to_e_osc[i] = float(elements[4])

print 'Ground to excited Z'
print Ground_to_e_Z 

#Creating the second row of info (info(2))
Excited_to_e_state1 = numpy.zeros([NStates*(NStates - 1)/2])        # the first state
Excited_to_e_state2 = numpy.zeros([NStates*(NStates - 1)/2])        # the second state
Excited_to_e_X = numpy.zeros([NStates*(NStates - 1)/2])
Excited_to_e_Y = numpy.zeros([NStates*(NStates - 1)/2])
Excited_to_e_Z = numpy.zeros([NStates*(NStates - 1)/2])
#Excited_to_e_osc = numpy.zeros([NStates*(NStates - 1)/2])          oscillation

trct = 0
for i in range(NStates-1):
    for j in range(i):
        trct = trct + 1

for (n,line) in enumerate(lines):
    if('Excited to excited state transition electric dipole moments' in line):
        for i in range (trct):
            elements = lines[n+i+2].split()
            # putting values into array
            Excited_to_e_state1[i] = (elements[0])
            Excited_to_e_state2[i] = (elements[1])
            Excited_to_e_X[i] = float(elements[2])
            Excited_to_e_Y[i] = float(elements[3])
            Excited_to_e_Z[i] = float(elements[4]) 
            #Excited_to_e_osc[i] = float(elements[7])

# get CI Vectors for later data analysis
ci_vect = numpy.zeros([NStates,NStates], dtype='float')
for (n,line) in enumerate(lines):
    if ('CI State Number' in line):
        for i in range(NStates):
            for j in range(NStates):
                elements = lines[n+1+i].split() 
                ci_vect[i][j] = float(elements[j+1])
print 'CI vectors'
print ci_vect   

# finally, get state configuration info
config = numpy.zeros([NStates,NMOs], dtype='int')
for (n,line) in enumerate(lines):
    if ('Orbital electron configuration' in line):
        for i in range(NStates):
            for j in range(NMOs):
                elements = lines[n+2+i].split()
                print "i",i,"j",j,"line #",(n+2+i)
                print elements
                #print "j",(j),elements[j]
                config[i][j] = int(elements[j+1])

# Initializing TransDipoles
# First create the (NStates + 2) x (NStates + 2) matrix that can take in an array as the element
TransDipoles = numpy.zeros((NStates + 2 ,NStates + 2), dtype='object')
#TransDipoles = numpy.zeros((NStates + 2 ,NStates + 2), dtype='float')       
for i in range(NStates + 2):
    for j in range(NStates + 2):
        TransDipoles[i,j] = numpy.zeros((3,1))      # initializing each element to a 3x1 matrix filled with 0

# Beginning the first for loop
# Assigning X,Y,Z of groand state
for k in range(NStates):
    TransDipoles[k+1,0] = [Ground_to_e_X[k],Ground_to_e_Y[k],Ground_to_e_Z[k]]  
    TransDipoles[0,k+1] = [Ground_to_e_X[k],Ground_to_e_Y[k],Ground_to_e_Z[k]]
    
# Beginning the second for loop
# Assigning X,Y,Z of excited state
for k in range(NStates*(NStates - 1)/2):
    i = Excited_to_e_state1[k]
    j = Excited_to_e_state2[k]
    
    TransDipoles[i,j] = [Excited_to_e_X[k],Excited_to_e_Y[k],Excited_to_e_Z[k]]
    TransDipoles[j,i] = [Excited_to_e_X[k],Excited_to_e_Y[k],Excited_to_e_Z[k]]

#print 'TransDipoles'
#print TransDipoles
#print 'TDMZ' 
#print Excited_to_e_Z 

# CIS Specttrum
Trans_dot = numpy.zeros((NStates+1),dtype = float)
for i in range(NStates):
    Trans_dot[i+1] = numpy.dot(TransDipoles[i+1,0],TransDipoles[i+1,0])

#print 'Trans_dot'
#print Trans_dot

# Plotting the CIS Spectrum
"""
fig = plt.figure(figsize=(6, 4))                    # Creating the figure
vax = fig.add_subplot(121)  
vax.scatter(Energies,Trans_dot)                     # Scatter plot of Energy to Trans_dot
vax.vlines(Energies,[0],Trans_dot)                  # Creating vertical line to scatter plot
plt.grid(True)                                      # adding grid
plt.show()
"""

# Ground state polarizability calculated using the CIS energies and transition dipoles
polarCIS = 0
for i in range (1,NStates-1):
    polarCIS = polarCIS + 2*(numpy.outer(TransDipoles[0,i],TransDipoles[0,i])/Energies[i])

#print 'polarCIS'
#print polarCIS

# Initilization
#NSim = 3           # Number of states in simulation
#emax = 0.05         # Field strength (in au)
#emax = 0.089
#maxstep = 30000     # Number of itegration steps
#delta = 0.05        # Stepsize (in au)

print "Simulation uses ", NSim, " states, ", maxstep, " integration steps with stepsize = ", delta, "au"

#envelope = "cosine"             # The envelope that could be created (trapezoidal,gaussian or cosine)
#w = 0.9767                        # The frequency (in au)
#w = 1.6369
#chirp = -0.01                   # chirp (in au)
#ncyc = 7                        # number of cycles in envelope
#phase = 0*numpy.pi/2 

# Based on envelope, the printed words will change
if (envelope == 'trapzoidal'):
    print ncyc, " cycle trapezoidal pulse (linear ramp on during first cyce, linear ramp off during last cycle)"
elif (envelope == 'gaussian'):
    print ncyc, " cycle gaussian pulse withFWHM = ", (ncyc/2), " cycles"
elif (envelope == 'cosine'):
    print ncyc, " cycle cosine pulse"

if (abs(chirp) < 0.000001):
    print"Maximum field strength ", emax, " au, frequency = ", w, " au" 
else:
    print"Maximum field strength " ,emax, " au, linear chirp: initial frequency ", w - chirp, " au and final frequency ", w + chirp, " au"

# graphing the envelope

# Creating our y value for the graph,which is not stated in mathematica because 
# when plotting in mathematica, the x and y does not need to be initialized
y= numpy.zeros((maxstep),dtype = float)
period = 2*numpy.pi/(w*delta)

#for t in range(0,maxstep):
#
#bfh- changed to toff instead of maxstep
#
for t in range(0,toff):
    # Getting the y values for envelope = 'trapezoidal'
    if(envelope == 'trapezoidal'):
        if(t < (ncyc*period)):
            if (t < (ncyc - 1)*period):
                if (t < period):
                    y[t] = (t/period)
                else:
                    y[t] = 1
            else: 
                y[t] = ((ncyc*period - t)/period)
        else: 
            y[t] = 0
    # Getting the y values for envelope = 'gaussian'
    elif(envelope == 'gaussian'):
        if (t < (ncyc*period)):
            y[t] = ((numpy.exp(-16*numpy.log(2)*(t/(ncyc*period) - 1/2)**2) - 1/16)/(15/16))
        else:
            y[t] = 0
    # Getting the y value for envelope = 'cosine'
    elif(envelope == 'cosine'):          
        if (t < ncyc*period):
            y[t] = 0.5 - numpy.cos(2*(numpy.pi)*t/(ncyc*period))/2
        else:
            y[t] = 0
            
# defining the time steps
t = range(0,maxstep)

# Plotting the envelope
"""
#plotting
plt.plot(t,y)
plt.grid(True)
plt.show()
"""

# Direction for the field
direction = [CoordsX[1]-CoordsX[0],CoordsY[1]-CoordsY[0],CoordsZ[1]-CoordsZ[0]]
norm_direc = direction/(numpy.linalg.norm(direction))   # normalized direction

#print 'direction'
#print norm_direc

# Initialization
dip = numpy.zeros((maxstep),dtype = numpy.float)        # creating an array that is filled with zeros
norm = numpy.zeros((maxstep),dtype = numpy.float)       # creating an array filled with zeros
v = numpy.zeros((maxstep,NSim),dtype = numpy.complex)   # creating the table(matrix) in mathematica code 


norm[0] = 1                         
norm[1] = 1
v[0,stindx] = 1
v[1,stindx] = 1

efield = numpy.zeros(maxstep)                           # creating array for efield
for t in range(1,toff):                              # putting values into efield
#
# ENVELOPE
# BFH- the commented out efield[t] includes the envelope function
#      it has been removed for simplicity. to include it, uncomment
#      the first efield[t] and comment out the second efield[t]
#
#    efield[t] = emax * y[t] * (numpy.sin((w + chirp * (t/(ncyc*period) - 0.5)) * delta * t + phase))
    efield[t] = emax * (numpy.sin(w*delta * t + phase))

h0 = numpy.zeros((NSim,NSim),dtype = numpy.float)       # creating empty matrix(NSim x NSim) called h0
i = m = 0           
for k in range(0,NSim):                                 # putting values into h0
    h0[i,m] = Energies[i]  
    i = i + 1                                           # These two counters are incrimented every loop
    m = m + 1                                           # so that the diagonal of matrix is filled

d0 = numpy.zeros((NSim,NSim),dtype = numpy.float)       # creaing an empty matrix(NSim x NSim) d0
for i in range(0,NSim):
    for j in range(0,NSim):
        d0[i,j] = numpy.dot(norm_direc,TransDipoles[i,j])                


#print 'd0-matrix: \n', d0
#print d0
#print 'h0'
#print h0
# Integration
h = h0 + efield[1] * d0                                 # creating matrix h that is (NSim x NSim)

ctmp1=numpy.zeros((NSim),dtype=complex)                 # initilizing ctmp1  
ctmp = numpy.zeros((NSim), dtype= complex)              # since ctmp and ctmp2 will have complex number
ctmp2 = numpy.zeros((NSim), dtype = complex)            # dtype = complex

ctmp1[stindx] = 1

comp_imag = (numpy.dot(1J*delta*h,ctmp1))               # calculating the complex portion of ctmp

ctmp = ctmp1 + comp_imag                                # putting values into ctmp
ctmp2 = ctmp/numpy.sqrt(abs(numpy.dot(numpy.conjugate(ctmp),ctmp)))

#print "ctmp2:"
#print ctmp2
#print "comp_imag:"
#print comp_imag
#print "ctmp1:"
#print ctmp1

#=============================================================================================
# Propagation: first-order approximation to exp(-1J*h*delta)
#=============================================================================================

# Index for the loop should be correct
for i in range(2,(maxstep)):        
    h = h0 + efield[i] * d0                             # Creating H_ij (Hamiltonian)
    comp_imag = (numpy.dot(2J*delta*h,ctmp2))           # calculating complex portion of ctmp
        
    ctmp  = ctmp1 + comp_imag
    
    ctmp1 = list(ctmp2)                                 # copying/assigning ctmp1 as ctmp2
    norm[i] = numpy.sqrt(abs(numpy.dot(numpy.conjugate(ctmp),ctmp)))
    ctmp2 = ctmp/norm[i]
    v[i,:] = list(ctmp2)                                # cmopying/assigning v[i,:] as ctmp2
    #print "v-matrix:"
    #print v
    dip[i] = numpy.real(numpy.dot(numpy.dot(numpy.conjugate(ctmp2),d0),ctmp2))

#
#NOTE: the instantaneous dipole moment expression above is missing the dipole moment of the
#      states. ONLY transition dipole moments are accounted for. Thus, (in pseudo-code)
#
# dip[t] = TDM terms +
#   sum(conjugate(v[t][i])*v[t][i]*conjugate(ci_vect[i][i])*ci_vect[i][i]
#                                                               for i in range(0,NStates))
#
#      The above is the TOTAL instantaneous dipole moment of the system. - KR
#


#print (max(norm) - 1, min(norm) - 1)

# Adiabatic energies as a function of the field
# two arrays come out for eigenvalues
                                      
#***problem section

nplot = 10
eig = numpy.zeros((21,2),dtype = 'object')                   #creating table
for i in range(1,22):                                      #computing the eigenvalues
    value, vector = numpy.linalg.eig(h0 + (i-1) *(emax/21)*d0)
#   eig[i] = sorted(value[0])
    idx = sorted(value)
#    eig[i-1] = sorted(value)
    eig[i-1,0] = idx[0]
    eig[i-1,1] = idx[1]	

trans_eig = numpy.zeros((21),dtype = 'object')
trans_eig = numpy.transpose(eig)                          #transposing the table

#print 'idx'
#print idx 
#taking care of listplot section

"""
xi = range(0,21)
plt.plot(xi, trans_eig[0])             
plt.plot(xi, trans_eig[1])
for i in range(0,21):
    plt.plot(i,trans_eig[0,i])

plt.axis(0,(emax))

plt.grid(True)
plt.show()
"""


# Plotting section in the Mathematica code 
# Electric field of the pulse
"""
x_e = range(maxstep)
plt.plot(x_e,efield)
plt.grid(True)
plt.show()
"""
# Instantaneous dipole                          
"""
x_d = range(maxstep)                  #setting the domain
plt.plot(x_d,-1*dip)
plt.grid(True)
plt.show()
"""
# Coefficient of the excited state "istate"     
"""
istate = 1
x_i = range(maxstep)                  #setting the domain
plt.plot(x_i,abs(v[:,istate]))
plt.grid(True)
plt.show()
"""
# Coefficients at the end of the simulation     
"""
y_coe_v = numpy.zeros(51)        
for i in range(NSim):                     
    y_coe_v[i] = abs(v[maxstep-1,i])                # putting values into y_coe_v

fig = plt.figure(figsize=(6, 4))                    # Creating the figure    
v_plot = fig.add_subplot(121)                       
v_plot.scatter(Energies[1:NSim],y_coe_v[1:NSim])    # Scatter Plot of Energies to y_coe_v
v_plot.vlines(Energies[1:NSim],[0],y_coe_v[1:NSim]) # Making it so that there is a vertical line to scatter plot
v_plot.grid(True)                                   # Adding grid 
plt.show()
"""
# output instanteous dipole moment
idm = open('tdci_dipole_st'+str(stindx)+'.txt', 'w')
sdm = open('tdci_dipole2fft_st'+str(stindx)+'.txt', 'w')
steptofs = 0.0241884*delta     # conversion factor for au to fs
dip = -dip/0.4      # Dipole moment in debyes
for i in range(2,maxstep):
    t_fs = i*steptofs
    t_s = t_fs*0.000000000000001
    idm.write('%14.8f %14.8f' % (t_fs, dip[i]) + "\n")
    sdm.write('%14.8e %14.8f' % (t_s, dip[i]) + "\n")
idm.close

# need to project state coefficients back onto MOs and print 
test = open('test.txt', 'w')
mo_coef = numpy.zeros([NMOs],dtype=numpy.float64)
tcoef = open('tdci_pops_st'+str(stindx)+'.txt', 'w')
for j in range(2,maxstep):
    t_fs = j*steptofs
    tcoef.write('%16.8f' % t_fs)
    for k in range(NMOs):
        mo_coef[k] = 0.0
    for i in range(NSim):
        tmp = numpy.conjugate(v[j][i])
        tmp2 = v[j][i]
        vcoef = tmp*tmp2
        if i==1:
            test.write('%12.8f' % vcoef) 
        for n in range(NSim):            
            for k in range(NMOs):
                mo_coef[k] = mo_coef[k] + vcoef*config[n][k]*ci_vect[i][n]*ci_vect[i][n]
    for k in range(NMOs):
        tcoef.write('%14.8f' % mo_coef[k])
    tcoef.write(' \n ')
tcoef.close

bcoef = open('tdci_coefs_st'+str(stindx)+'.txt', 'w')
for j in range(2,maxstep):
    t_fs = j*steptofs
    bcoef.write('%16.8f' % t_fs)
    for i in range(NSim):
        tmp = numpy.conjugate(v[j][i])
        tmp2 = v[j][i]
        vcoef = tmp*tmp2
        bcoef.write('%16.8f' % vcoef)
    bcoef.write("\n")
bcoef.close

fldprt = open('tdci_field_st'+str(stindx)+'.txt', 'w')
for i in range(1,maxstep):
    t_fs = i*steptofs
    fldprt.write('%16.8f %16.8f' % (t_fs, efield[i]) + "\n")
fldprt.close




#print 'v = ' 
#print v

#print 'efield = '
#print efield 

#print 'dipole = '
#print dip
