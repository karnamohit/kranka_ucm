%chk=cis_s0_lih_sto-3g.chk
#p hf/sto-3g CIS(NStates=30,Root=1,AllTransitionDensities,Singlets) nosymm IOp(3/33=3,3/36=1,4/33=3,5/33=3,5/72=-3,6/8=1,6/33=3,8/10=90,9/33=3,9/40=3)

CIS calculation

0 1
 H       0.000000    0.000000   -0.765000
 Li      0.000000    0.000000    0.765000

