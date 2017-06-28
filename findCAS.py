#!/usr/bin/env python
#
# Author: James Smith <orbames.smith9113@gmail.com>
'''
Find MOs with certain AO character, e.g. C pz, Fe dz2.
'''

import numpy, sys

def findCASfromOut( filename, aos_tofind, num_best_orbs=20, **kwargs ):
    '''Read output file after mf.analyze() used to aid in identifying orbitals
        intended for CASSCF calculations. The user must call mf.analyze(ncol=1).
        The indices are printed in decreasing order of how well they match the
        AO.

    Args:
        filename: str
            Name of output from PySCF calculation where mf.analyze() results
            are stored.
        aos_tofind: list of strings
            The list of AOs for find, e.g.   ['Fe 3dxz', 'C 2pz'].

    Kwargs:
        num_best_orbs: int
            The number of "best" MOs to print/return. Default is 20.

    Returns:
        mo_idx: list of int
            Returns the index of the single best MO for the desired AO.

    Examples:
        >>> mf.kernel()
        >>> mf.analyze(ncol=1)
        >>> num_best_orbs = 3
        >>> aos_tofind = ['Fe 3dxz']
        >>> findCAS( filename, num_atoms, num_orbs, mo_start, mo_end,
        .. aos_tofind )
        Fe   3dxz
        ============

        81.0
        68.0
        91.0
    '''
    ncasao = len(aos_tofind)
    mo_idx = [0]*ncasao
    mo_coeff = [0]*ncasao
    this_mo_coeff = [0]*ncasao


    # These are like to become unstable in the future if PySCF changes output
    # format.
    start_seq = " ** MO coefficients (expansion on meta-Lowdin AOs) **"
    start_seq = start_seq.split()
    end_seq = " ** Mulliken pop on meta-lowdin orthogonal AOs  **"
    end_seq = end_seq.split()
    gto_seq = "number of NR cGTOs ="; gto_seq = gto_seq.split()
    natom_seq = "[INPUT] num atoms ="; natom_seq = natom_seq.split()

    for i in range(len(aos_tofind)):
        aos_tofind[i] = aos_tofind[i].split()

    with open( filename, 'r' ) as f:
        data = f.readlines()

    # Automatically find start and end of mf.analyze() MOs
    # TODO: Tighten this up and work it into the main loop.
    for i in range(len(data)):
        d_split = data[i].split()
        d_split2 = d_split[0:len(d_split)-1]

        if ( d_split2 == natom_seq ): num_atoms = int(d_split[-1])
        if ( d_split2 == gto_seq ): num_orbs = int(data[i].split()[-1])
        if ( d_split == start_seq ): mo_start = i + 1
        if ( d_split == end_seq  ): mo_end = i - 1


    # Data structures dependent on num_orbs, mo_start, and mo_end
    pos_mos = numpy.zeros((num_orbs,2*ncasao))
    data = data[mo_start:mo_end] # Only keep necessary part of output file.
    mo_counter = 0

    for i in range(len(data)):
        data[i] = data[i].split()
        ncols = len(data[i]) - 3
        if ( data[i][0][0] == '#' ):
            if ( i > 0 ):
                for orb in range(len(aos_tofind)):
                    if ( this_mo_coeff[orb] > mo_coeff[orb] ):
                        # Track the best MOs
                        mo_idx[orb] = mo_counter
                        mo_coeff[orb] = this_mo_coeff[orb]

                    # Keep track of all MOs
                    pos_mos[mo_counter][orb*2] = mo_counter
                    pos_mos[mo_counter][orb*2 + 1] = this_mo_coeff[orb]
                    this_mo_coeff[orb] = 0 # Resets for next MO

            mo_counter += 1
            continue

        # Check to see if AO is in aos_tofind.
        for orb in range(len(aos_tofind)):
            if (data[i][1] == aos_tofind[orb][0] and
                data[i][2] == aos_tofind[orb][1]):
                this_mo_coeff[orb] += abs(float(data[i][3]))


    # print mo_idx

    for orb in range(ncasao):
        print( aos_tofind[orb][0] + " " + aos_tofind[orb][1] + "\n============" )

        best_mos = pos_mos[:,orb*2:orb*2 +2]
        best_mos = best_mos[best_mos[:,1].argsort()]
        best_mos = best_mos[len(best_mos)-num_best_orbs:len(best_mos)]
        print_mos = []
        for i in range(num_best_orbs):
            print_mos.append( int(best_mos[num_best_orbs-1-i,0]) )
            # print( "%i" % best_mos[num_best_orbs-1-i,0] )

        print( print_mos )
        print( "\n" )

    return mo_idx

def prindIrrepData( irreps, symData ):
    for i in range(len(irreps)):
        print( irreps[i], ": ", symData[i])

def getOccAndSym( casOrbs, pointGroup, filename, printMOs ):
    # Data objects
    orderedCAS = sorted( casOrbs )
    irreps = getPointGroupIrreps( pointGroup )
    casOcc = [0]*len(irreps)
    if ( pointGroup != "c1" ):
        casSym = [0]*len(irreps)
        coreSym = [0]*len(irreps)



    # Get Data
    gto_seq = "number of NR cGTOs ="; gto_seq = gto_seq.split()

    with open( filename, 'r' ) as f:
        data = f.readlines()
    for i in range(len(data)):
        if ( data[i] == "**** MO energy ****\n" ):
            startLine = i + 2
            continue
        if ( data[i] == " ** MO coefficients (expansion on meta-Lowdin AOs) **\n" ):
            endLine = i
            continue

        d_split = data[i].split()
        d_split2 = d_split[0:len(d_split)-1]
        if ( d_split2 == gto_seq ):
            num_orbs = int(data[i].split()[-1])
            continue

    # Sift through data and collect the relevant parts
    orbI = 0 # Current orb we're looking at
    naselec = 0 # Number of active space orbitals

    if ( pointGroup != "c1" ):
        for line in range(startLine,endLine):
            orbI += 1
            lineSplit = data[line].split()

            # Numerical representation of irrep
            thisIrrep = lineSplit[2].strip("(")

            for i in range(len(irreps)):
                if (irreps[i] == thisIrrep ):
                    thisIrrep = i

            thisOcc = int(lineSplit[-1])
            # thisString = "#" + str( orbI )

            # CAS Data
            if ( orbI in casOrbs ):
                naselec += int(lineSplit[-1])
                casSym[thisIrrep] += 1
                if (int(lineSplit[-1]) > 0):
                    casOcc[thisIrrep] += int(lineSplit[-1])
                if printMOs: print (lineSplit)

            # Core Data
            else:
                if ( thisOcc > 0 ):
                    coreSym[thisIrrep] += thisOcc
    else:
        for line in range(startLine,endLine):
            orbI += 1
            lineSplit = data[line].split()

            thisOcc = int(lineSplit[-1])

            # CAS Data
            if ( orbI in casOrbs ):
                naselec += int(lineSplit[-1])
                if printMOs: print (lineSplit)


    if ( printMOs ):
        print ("CAS Data\n=============")
        print ("Number of electrons in active space: ", naselec)
        print ("Number of orbitals in active space: ", len(casOrbs))

        if ( pointGroup != "c1" ):
            print ("Symmetry of active space orbitals: ")
            prindIrrepData(irreps, casSym )
            print ("Electrons in active space orbitals: ")
            prindIrrepData(irreps, casOcc )
            print ("Ordered CAS list: ", ordered_CAS)
            print ("\nCore Data\n=============")
            print ("Electrons in core orbitals by symmetry: ")
            prindIrrepData(irreps, coreSym )

    if ( pointGroup != "c1" ):
        return  naselec, casSym, casOcc, coreSym
    else:
        return naselec

def getPointGroupIrreps( pointGroup ):
    symLists = {
        "d2h": ["Ag", "B1g", "B2g", "B3g", "Au", "B1u", "B2u", "B3u"],
        "c1": ["A"]
    }

    try:
        return symLists[pointGroup]

    except ValueError:
        print( "Symmetry group %s not implemented at this time" % pointGroup)



if __name__ == '__main__':

    # Make sure the user calls mf.analyze(ncol=1)
    rawArgs = sys.argv
    filename = str(rawArgs[1])
    num_best_orbs = 25

    # aos_tofind = ['C 2pz', 'Mn 3dx2-y2', 'Mn 3dz2', 'Mn 3dxy','O 2px', 'O 2py', 'O 2pz', 'N 2px', 'N 2py']

    # aos_tofind = ['Mn 3dz2', 'Mn 3dx2-y2', 'Mn 3dxy', 'Mn 3dxz', 'Mn 3dyz']
    # aos_tofind = ['C 2pz','N 2pz','O 2pz']
    # aos_tofind = ['N 2px','O 2px','N 2py','O 2py']
    aos_tofind = ['Cl 3pz','Cl 3px','Cl 3py']

    findCASfromOut( filename, aos_tofind, num_best_orbs=num_best_orbs )

    casOrbs = [69,75,70,79,50,71,46,44,64,60,23,10,3,30,94,185,110,84,38,44,11,17,27] # See page 85 of LNb

    getOccAndSym( casOrbs, "c1", filename, True )

    # aos_tofind = ['Fe 3dxz', 'Fe 3dxy', 'Fe 3dx2-y2', 'Fe 3dz^2', 'Fe 3dyz',
    #     'Fe 4dxz', 'Fe 4dxy', 'Fe 4dx2-y2', 'Fe 4dz^2', 'Fe 4dyz']
    #
    # findCASfromOut( filename, num_atoms, num_orbs, aos_tofind,
    #     num_best_orbs=num_best_orbs )
    #
    # aos_tofind = ['N 2px', 'N 2py', 'Fe 4px', 'Fe 4py']
    #
    # findCASfromOut( filename, num_atoms, num_orbs, aos_tofind,
    #     num_best_orbs=num_best_orbs )
