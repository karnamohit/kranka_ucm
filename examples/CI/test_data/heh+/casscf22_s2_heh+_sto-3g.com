%chk=casscf22_s2.chk
#P CASSCF(2,2,fulldiag,NRoot=4,SaveGEDensities)/sto-3g scf(tight,maxcyc=50) nosymm iop(6/8=1,3/33=3,3/36=1,4/33=3,5/33=3,6/33=3,9/33=3,5/72=-3)

Preliminary HF calculation to determine orbital symmetry

1 1
 H       0.000000    0.000000   -0.386000
 He      0.000000    0.000000    0.386000

