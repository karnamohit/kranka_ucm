import numpy
import sys
from hf_info import el_info

# running this script:
#
# $ python grid.py blah.cube gradient=<true/false> N1 N2 N3
#
#
# N1, N2, N3: # of pts along x-, y- and z-axes
# N1*N2*N3 is the amt of grid-points used for the density and the gradients

NELECT = 2
NMOs = 2
fchkfile = 'h2_s0.fchk'
logfile = 'h2_s0.log'

filename = sys.argv[1]+'.cube'

str1 = sys.argv[2]
str2 = str1.split('=')
densgrad=str2[1]

xpts=int(float(sys.argv[3]))
ypts=int(float(sys.argv[4]))
zpts=int(float(sys.argv[5]))

file = open(filename,'r')
lines = file.readlines()

stepsize = numpy.zeros((3), numpy.float64)

# calculating the volume element

nskip = 8
elements = []
m = 0
l = 0
last_indx = lines.index(lines[-1])
print 'last index ',last_indx
for (n,line) in enumerate(lines):
    while (l <= 2):
        str1 = lines[3+l].split()
        print str1
        str1 = filter(None, str1)
        stepsize[l] = float(str1[l+1])
        l += 1
    while ((nskip+m) <= last_indx):
        str1 = lines[nskip+m]
#        print 'line number ',nskip+m
#        print 'string ',list(lines[nskip+m])
        str2 = str1.split(' ')
        str3 = str2[-1]
        str4 = str3.split('\n')
        str3 = str4[0]
        str2[-1] = str3
        str2 = filter(None, str2)	# gettting rid of the null elements, like ' ' or 0
#        print 'elements ',list(str2)
        elements.extend(str2)
        m += 1

volume = 1

for i in range(3):
    volume *= stepsize[i]

print 'dimensions of the volume element: ',stepsize
#print '1st elements ',lines[nskip+0]
#print 'number of elements ',len(elements)

dens = numpy.zeros((xpts,ypts,zpts), dtype='object')

m = 1
if (densgrad == 'true'):
    
    grad = numpy.zeros((xpts,ypts,zpts,3), dtype='object')
    
    for x in range(xpts):
        for y in range(ypts):
            for z in range(zpts):
                dens[x,y,z] = elements[x*ypts*zpts*4+y*zpts*4+z*4]
                dens[x,y,z] = float(dens[x,y,z])
                for w in range(3):
                    grad[x,y,z,w] = elements[x*ypts*zpts*4+y*zpts*4+z*4+w+1]
    
    print 'last grad element ',grad[(xpts-1),(ypts-1),(zpts-1),2]
    print 'last dens element ',dens[(xpts-1),(ypts-1),(zpts-1)]
    print grad.shape,dens.shape
    
else:
    
    for x in range(xpts):
        for y in range(ypts):
            for z in range(zpts):
                dens[x,y,z] = elements[x*ypts*zpts+y*zpts+z]
                dens[x,y,z] = float(dens[x,y,z])

print ' max. value of density: ',numpy.amax(dens)
elect_info = el_info(NELECT,NMOs,fchkfile,logfile)

SAO = elect_info.overlapAO()
PAO = elect_info.densityAO()
cMO = elect_info.coefficientsMOalpha()
PMO = elect_info.densityMO()

print 'density matrix:\n',PMO
cMOT = cMO.T

grid_el = numpy.sum(dens)*volume

#for x in range(xpts):
#    for y in range(ypts):
#        for z in range(zpts):
#            for i in range(NMOs):
#                griddensMO += cMOT[i,:].dot(SAO.dot(SAO.dot(cMO[:,i])))*float(dens[x,y,z])

densMO = numpy.trace(PMO)

print '# electrons from the grid density = ', grid_el
print 'trace of density matrix = ',densMO

