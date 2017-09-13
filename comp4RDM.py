#!/usr/bin/env python
import numpy as np
import sys

### Arguments ###
# Script takes four command line arguments: 
#    1) norm threshold
#    2) wd for files organized in the manner shown below
#    3) nelec
#    4) norbs
#
# Ex) python comp4RDM.py 1e-3 n2 10 8

thold = float(sys.argv[1])
wd = str(sys.argv[2])
nelec = int(sys.argv[3])
norbs = int(sys.argv[4])

outputLvl = "debug" # "debug" or "normal" anything else is equivalent to none
displayMax = 50 # Maximum number of incorrect 4RDM indices to be displayed
passCompTest = True

### Functions ###
def read4RDM(rdm,fileName):
    with open(fileName) as f:
        content = f.readlines()

    for i in range(1,len(content)):
        c0,c1,c2,c3,d0,d1,d2,d3,val = content[i].split()
        rdm[int(c0),int(c1),int(c2),int(c3),
            int(d0),int(d1),int(d2),int(d3)] = float(val)


def read3RDM(rdm,fileName):
    with open(fileName) as f:
        content = f.readlines()

    for i in range(1,len(content)):
        C0, C1, C2, D0, D1, D2, val = content[i].split()
        rdm[int(C0),int(C1),int(C2),int(D0),int(D1),int(D2)] = float(val)
        
            
### Load DMRG 4RDM ###
dmrg = np.zeros( (norbs,)*8 )
import os
os.system("ls -lst ./%s/dmrg/"%wd)
read4RDM(dmrg,"./%s/dmrg/spatial_fourpdm.0.0.txt") % wd
    
### Load Dice 4RDM ###
d4 = np.zeros( (norbs,)*8 )
read4RDM(d4,"./%s/spatial4RDM.0.0.txt") % wd

### Load Dice 3RDM ###
dice3 = np.zeros( (norbs,)*6 )
read3RDM(dice3,"./%s/spatial3RDM.0.0.txt") % wd

### DMRG 3RDM ###
dmrg3 = np.zeros( (norbs,)*6 )
read3RDM(dmrg3,"./%s/dmrg/spatial_threepdm.0.0.txt") % wd

### Comparison ###
if ( np.linalg.norm(d4-dmrg) > thold ):
    passCompTest = False
    counter = 0
    if (outputLvl=="debug"):
        for i in range(norbs):
            for j in range(norbs):
                for k in range(norbs):
                    for l in range(norbs):
                        for m in range(norbs):
                            for n in range(norbs):
                                for o in range(norbs):
                                    for p in range(norbs):
                                        if ( np.isclose(dmrg[i,j,k,l,m,n,o,p],
                                                        d4[i,j,k,l,m,n,o,p]) ): pass
                                        elif ( counter < displayMax ):
                                            print i,j,k,l,m,n,o,p, " ",\
                                                "DMRG: ",dmrg[i,j,k,l,m,n,o,p], \
                                                "  DICE: ",d4[i,j,k,l,m,n,o,p]
                                            counter+=1
                                        else:
                                            break
                                    
                                    if (counter > displayMax ): break
                                if(counter > displayMax): break
                            if(counter > displayMax): break
                        if (counter > displayMax): break    
                    if (counter > displayMax): break    
                if (counter > displayMax): break    
            if (counter > displayMax): break    
    

# Check num electrons correct
dm3 = np.einsum('ijkllmno-> ijkmno',dmrg)/(nelec-3)
dm2 = np.einsum('ijkklm -> ijlm',dm3)/(nelec-2)
dm1 = np.einsum('ijjk-> ik',dm2)/(nelec-1)


d3 = np.einsum('ijkllmno-> ijkmno',d4)/(nelec-3)
d2 = np.einsum('ijkklm -> ijlm',d3)/(nelec-2)
d1 = np.einsum('ijjk-> ik',d2)/(nelec-1)

### Output ###
# ' in print statements indicates that they were generating with einsum
print "Testing species in %s directory" % wd
if (passCompTest): print "...PASSED\n......Norm(dmrg4RDM-dice4RDM)<%f\n" % thold
else: print "...FAILED\n......Norm(dmrg4RDM-dice4RDM)<%f\n" % thold
print "Test DMRG nelec: ", np.einsum('ii',dm1)
print "Test DICE nelec: ", np.einsum('ii',d1)
print "Dice-DMRG Norm. ", np.linalg.norm(d4-dmrg)

print "Dice3'-Dice3 Norm. ", np.linalg.norm(d3-dice3)
print "Dice3'-DMRG3 Norm. ", np.linalg.norm(d3-dmrg3)
print "Maximum absolute difference: ", np.absolute(d4-dmrg).max()
print "\n"
