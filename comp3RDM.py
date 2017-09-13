#!/usr/bin/env python
# Compares the 3RDM from Dice and DMRG for N2 (6e,6o)
import sys
import numpy as np

### Arguments ###
# Script takes four command line arguments: 
#    1) norm threshold
#    2) wd for files organized in the manner shown below
#    3) nelec
#    4) norbs
#
# Ex) python comp3RDM.py 1e-3 n2 10 8

### Settings and command line args ###
passCompTest = True

thold = float(sys.argv[1])
wd = str(sys.argv[2])
nelec = int(sys.argv[3])
norbs = int(sys.argv[4])
norbs2 = norbs * norbs; norbs6 = norbs*norbs*norbs*norbs*norbs*norbs

dmrgFileName = "./%s/dmrg/spatial_threepdm.0.0.txt" % wd
diceFileName = "./%s/spatial3RDM.0.0.txt" % wd
dice2RDMName = "./%s/spatialRDM.0.0.txt" % wd

### Initialize Data Arrays  ###
dmrg = np.zeros( (norbs,)*6 )
dice = np.zeros( (norbs,)*6 )
dice2RDM = np.zeros( (norbs,)*4 )

### Dice 3RDM ###
with open(diceFileName) as f:
    content = f.readlines()

for i in range(1,len(content)):
    C0, C1, C2, D0, D1, D2, val = content[i].split()
    dice[int(C0),int(C1),int(C2),int(D0),int(D1),int(D2)] = float(val)

### Dice 2RDM ###    
with open(dice2RDMName) as f:
    content = f.readlines()

for i in range(1,len(content)):
    c0, c1, d1, d0, val = content[i].split()
    dice2RDM[int(c0),int(c1),int(d0),int(d1)] = float(val)
    dice2RDM[int(c1),int(c0),int(d1),int(d0)] = float(val)

### DMRG 3RDM ###
with open(dmrgFileName) as f:
    content = f.readlines()

for i in range(1,len(content)):
    C0, C1, C2, D0, D1, D2, val = content[i].split()
    dmrg[int(C0),int(C1),int(C2),int(D0),int(D1),int(D2)] = float(val)

### Comparison ###
dm2 = np.einsum('ijkklm->ijlm',dmrg)
dm2 /= (nelec-2)
dm1 = np.einsum('ijjk->ik', dm2)
dm1 /= (nelec-1)

d2 = np.einsum('ijkklm->ijlm',dice)
d2 /= (nelec-2)
d1 = np.einsum('ijjk->ik', d2)
d1 /= (nelec-1)


if ( np.linalg.norm(dmrg-dice) > thold ):
    passCompTest = False
    for i in range(norbs):
        for j in range(norbs):
            for k in range(norbs):
                for l in range(norbs):
                    for m in range(norbs):
                        for n in range(norbs):
                            
                            if ( dmrg[i,j,k,l,m,n] != 0 and \
                                 abs( (dmrg[i,j,k,l,m,n]-dice[i,j,k,l,m,n])/\
                                      dmrg[i,j,k,l,m,n] ) > 1e-3 ):
                                if ( i < 1 and j < 2 and l < 2 ): 
                                    print i,j,k,l,m,n, "   DMRG: ", \
                                        dmrg[i,j,k,l,m,n], "   Dice: ", \
                                        dice[i,j,k,l,m,n]


### Output ###
print "Testing species in %s directory" % wd
if (passCompTest): print "...PASSED\n......Norm(dmrg3RDM-dice3RDM)<%f\n" % thold
else: print "...FAILED\n......Norm(dmrg4RDM-dice4RDM)<%f\n" % thold
print "DMRG Einsum: ", np.einsum('ii',dm1)
print "Dice Einsum: ", np.einsum('ii',d1)

print "Norm between 3RDM: ", np.linalg.norm(dmrg-dice)
print "Norm between 2RDM from Dice: ", np.linalg.norm(dice2RDM-d2)
print "Norb between DRMG2RDM and Dice ", np.linalg.norm(dm2-d2)
print "\n"
