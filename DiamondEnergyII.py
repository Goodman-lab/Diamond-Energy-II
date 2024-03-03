#################################################################
## Diamond Energy version 2.0
##
## Input of just an InChI will do a conformation search which
## lists all low-energy conformations and an estimate of their
## energies, approximately in kJ/mol.
##
## Input of an InChI and a conformation number will report the
## energy and generate an .sdf file of the structure of that
## conformation
##
#################################################################

import sys
import array
import time
from time import process_time
# from IPython.core.display import display
from rdkit import Chem
# from rdkit.Chem.Draw import IPythonConsole
# IPythonConsole.ipython_useSVG=True
from rdkit.Chem import AllChem, rdMolDescriptors
from rdkit.Chem import Draw
from rdkit.Chem import rdmolfiles
from rdkit.Chem import Descriptors
# from IPython.display import display, Markdown, HTML
import re

start_time = process_time()


def takeSecond(elem):
    return elem[1]


def atom_counts(molecule, atom_symbol):
    """Return the count of a specific atom type in a molecule."""
    return sum(1 for atom in molecule.GetAtoms() if atom.GetSymbol() == atom_symbol)


def get_alcohol_oxygens(molecule):
    """Return a list of oxygen atom indices that belong to alcohol groups."""
    alcohols = []
    for atom in molecule.GetAtoms():
        if atom.GetSymbol() == 'O' and atom.GetDegree() == 2:
            neighbors = [n.GetSymbol() for n in atom.GetNeighbors()]
            if 'C' in neighbors and 'H' in neighbors:
                alcohols.append(atom.GetIdx())
    return alcohols


def extract_and_convert(s):
    """Extracts the numeric part of a string and converts it to an integer."""
    match = re.search(r'(\d+)', s)
    if match:
        return int(match.group(1))
    return None


def printatoms():
    print("                I#,At,Di,Tw,Nd, x, y, z, Co ")
    for i in range(number_of_skeleton):
        print("%3d " % (i), atom[i * 9:i * 9 + 9])


def printxyzatoms():
    # this is an XYZ file
    print(number_of_skeleton)
    print(sys.argv[1])
    scalefactor = 0.85
    for i in range(number_of_skeleton):
        if list_of_connectivity_int[i] in oxygen_atom:
            print("O   %8.4f  %8.4f  %8.4f" % (float(atom[i * 9 + 5]) * 0.885, float(atom[i * 9 + 6]) * 0.885,
                                               float(change_enantiomer * atom[i * 9 + 7]) * 0.885))
        else:
            print("C   %8.4f  %8.4f  %8.4f" % (float(atom[i * 9 + 5]) * 0.885, float(atom[i * 9 + 6]) * 0.885,
                                               float(change_enantiomer * atom[i * 9 + 7]) * 0.885))


def printmol():
    # NB must be in .mol order, not InChI order, for a molfile
    # this is .sdf file single molecule .mol format
    with open("Mol_conf"+str(conf_num)+".sdf", "a+") as f:
        f.write(sys.argv[1])
        f.write("\n")
        f.write("\n")
        f.write("\n")
        ## Here the bond number info need to add the number of bonds connecting the first and final ring atom as well
        f.write("%3d%3d  0  0  0  0  0  0  0  0999 V2000" % (
        number_of_skeleton, number_of_skeleton - 1 + len(first_atom_InChI)))
        f.write("\n")
        # scalefactor = 0.85s
        # write coordinates information
        i = 0
        while i < number_of_skeleton:
            if list_of_connectivity_int[i] in oxygen_atom:
                f.write("%10.4f%10.4f%10.4f O   0  0  0  0  0  0  0  0  0  0  0  0" % (
                float(atom[i * 9 + 5]) * 0.885, float(atom[i * 9 + 6]) * 0.885, float(atom[i * 9 + 7]) * 0.885))
                f.write("\n")
            else:
                f.write("%10.4f%10.4f%10.4f C   0  0  0  0  0  0  0  0  0  0  0  0" % (
                float(atom[i * 9 + 5]) * 0.885, float(atom[i * 9 + 6]) * 0.885, float(atom[i * 9 + 7]) * 0.885))
                f.write("\n")
            i += 1
        # write bond connectivity information
        j = 0
        bond = []
        while j < number_of_skeleton:
            if (atominfo(inchi2wn[int(list_of_connectivity_int[j])], 1) == 0 and list_of_connectivity_int[
                j] not in first_atom_InChI):
                pass
            else:
                if list_of_connectivity_int[j] in first_atom_InChI:
                    current_ring_index = first_atom_InChI.index(list_of_connectivity_int[j])
                    final_atom = final_atom_InChI[current_ring_index]
                    bond_value = []
                    bond_value.append(list_of_connectivity_int[j])
                    bond_value.append(atominfo(inchi2wn[list_of_connectivity_int[j]], 1))
                    bond_value.sort()
                    bond.append(bond_value)
                    bond_value_1 = []
                    bond_value_1.append(list_of_connectivity_int[j])
                    bond_value_1.append(final_atom)
                    bond_value_1.sort()
                    bond.append(bond_value_1)
                else:
                    bond_value_2 = []
                    bond_value_2.append(list_of_connectivity_int[j])
                    bond_value_2.append(atominfo(inchi2wn[list_of_connectivity_int[j]], 1))
                    bond_value_2.sort()
                    bond.append(bond_value_2)
            j += 1
        # bond.sort()
        bond_num = 0
        while bond_num < len(bond):
            if (inchi2wn[bond[bond_num][0]] + 1) != (inchi2wn[bond[bond_num][1]] + 1):
                f.write("%3d%3d  1  0" % (inchi2wn[bond[bond_num][0]] + 1, inchi2wn[bond[bond_num][1]] + 1))
                f.write("\n")
            else:
                pass
            bond_num += 1
        f.write("M  END")
        f.write("\n")
        f.write("$$$$")
        f.write("\n")
    f.close()


def atominfo(atomnumber, information):
    # Return data for atom 'atomnumber'
    # 0: InchiNumber
    # 1: Previous Attached Atom
    # 2: Was Direction; now energy of atom
    # 3: Twist
    # 4: New Direction
    # 5,6,7: x,y,z coordinates
    # 8: conformation - number of the relevant torsion angle
    # print("atominfo working",atomnumber,information,atomnumber*8+information)
    return atom[atomnumber * 9 + information]


def setatominfo(atomnumber, information, datapoint):
    atom[atomnumber * 9 + information] = datapoint


def GetRingSystems(mol, includeSpiro=False):
    ri = mol.GetRingInfo()
    systems = []
    for ring in ri.AtomRings():
        ringAts = set(ring)
        nSystems = []
        for system in systems:
            nInCommon = len(ringAts.intersection(system))
            if nInCommon and (includeSpiro or nInCommon > 1):
                ringAts = ringAts.union(system)
            else:
                nSystems.append(system)
        nSystems.append(ringAts)
        systems = nSystems
    return systems


def mol_with_atom_index(mol):
    for atom in mol.GetAtoms():
        print(atom.GetIdx())
        atom.SetAtomMapNum(atom.GetIdx())
    return mol


#################################################################
## From a conformation number, recalculate newdirection and
## calculate the x, y, z coordinates of the conformation
## The current conformation is define by list conformation[]
## NB: direction (atominfo(#, 2) is currently unused
## Rest of system is unchanged
## old_conformation[] has the list for the energies in the starting structure
## except first time around
def find_all_conf_number():
    # energy:
    # count number of overlaps, adjacent atoms, 1,5 diaxial, gauche
    num_overlap = 0
    num_adjacent = 0
    num_15diaxial = 0
    num_gauche = 0
    num_axeq = 0
    num_13connect = 0
    energy_total = 0
    ##count the number of C_C interraction
    num_15diaxial_C_C = 0
    num_gauche_C_C = 0
    num_axeq_C_C = 0
    ##count the number of C_O interraction
    num_15diaxial_C_O = 0
    num_gauche_C_O = 0
    num_axeq_C_O = 0
    ##count the number of O_O interraction
    num_15diaxial_O_O = 0
    num_gauche_O_O = 0
    num_axeq_O_O = 0
    ##count the number of C_OH interraction
    num_15diaxial_C_OH = 0
    num_gauche_C_OH = 0
    num_axeq_C_OH = 0
    ##count the number of O_OH interraction
    num_15diaxial_O_OH = 0
    num_gauche_O_OH = 0
    num_axeq_O_OH = 0
    ##count the number of OH_OH interraction
    num_15diaxial_OH_OH = 0
    num_gauche_OH_OH = 0
    num_axeq_OH_OH = 0
    # Original Energy components
    energy_axeq = -2
    energy_gauche = 2
    energy_15diaxial = 20
    energy_adjacent = 200
    energy_overlap = energy_adjacent
    # New Energy components
    energy_axeq_C_O = -1
    energy_gauche_C_O = -3
    energy_15diaxial_C_O = 4
    energy_axeq_O_O = -6
    energy_gauche_O_O = -4
    energy_15diaxial_O_O = 11
    energy_axeq_C_OH = -1
    energy_gauche_C_OH = -2
    energy_15diaxial_C_OH = 2
    energy_axeq_O_OH = -1
    energy_gauche_O_OH = -4
    energy_15diaxial_O_OH = -11
    energy_axeq_OH_OH = -9
    energy_gauche_OH_OH = -3
    energy_15diaxial_OH_OH = -16
    energy_adjacent_O_OH = -26
    energy_adjacent_OH_OH = -29
    energy_13connect_O_OH = -21
    energy_13connect_OH_OH = -23



    current_atom_energy = 0
    current_atom_marker = 0
    num_repeat_torsions = -1
    nt = 0
    while nt < number_rotatable_bonds:
        if conformation[nt] == old_conformation[nt]:
            num_repeat_torsions += 1
        else:
            nt = number_rotatable_bonds
        nt += 1
    # print("Number of repeated torsions",num_repeat_torsions+1,conformation,old_conformation)
    # print("Conformation",conformation)

    # do the cs process
    ia = 0
    while ia < number_of_skeleton:
        current_atom_marker = ia
        # print("current_atom_marker", current_atom_marker)
        current_torsion = atominfo(ia, 8)
        # print("TEST!!!current_torsion:", current_torsion)
        # print("TEST!!!atom_info:")
        # printatoms()
        ##if current_torsion > -1:
        if current_torsion > -1 and current_torsion >= num_repeat_torsions:
            # print("TEST!!!current_torsion:", current_torsion)
            if atominfo(inchi2wn[atominfo(inchi2wn[atominfo(ia, 1)], 1)], 1) == 0:
                # print("TEST!!!current_torsion:", current_torsion)
                total_twist = 0
            ## Expect the first and the second atom on the ring,
            ## the rest atoms will maintain untwist status while rotating mol
            ## Consequently, the relative twist value of these atoms will not be changed
            ## Judge whether the current atom is a fourth ring atom onward
            elif atominfo(ia, 0) in ringatoms:
                # print("TEST!!!current_torsion:", current_torsion)
                ring_index = -1
                for ring in inchi_ring_order:
                    ring_index += 1
                    if atominfo(ia, 0) in ring:
                        atom_index = ring.index(atominfo(ia, 0))
                        break
                #print("ring_index:", ring_index)
                #print("atom_index:", atom_index)
                # The first two ring atoms, are treated the same as other non-ring atoms
                if atom_index < 2:
                    total_twist = (atominfo(ia, 3) + conformation[current_torsion]) % 3
                else:
                    if ring_twist_statue[ring_index] == 0:
                        if conformation[current_torsion] == 0:
                            total_twist = atominfo(ia, 3)
                        if conformation[current_torsion] == 1:
                            total_twist = atominfo(ia, 3)
                        # ring flex
                        if conformation[current_torsion] == 2:
                            # since ring flex will affect the branch before the third ring atom
                            if atom_index == 2:
                                total_twist = conformation[current_torsion]
                            if atom_index == 3:
                                total_twist = (conformation[current_torsion] + ring_flip_list[0][0]) % 3
                            if atom_index == 4:
                                total_twist = (conformation[current_torsion] + ring_flip_list[0][1]) % 3
                            if atom_index == 5:
                                total_twist = (conformation[current_torsion] + ring_flip_list[0][2]) % 3
                    if ring_twist_statue[ring_index] == 1:
                        if conformation[current_torsion] == 0:
                            total_twist = atominfo(ia, 3)
                        # ring flex
                        if conformation[current_torsion] == 1:
                            # since ring flex will affect the branch before the third ring atom
                            if atom_index == 2:
                                total_twist = conformation[current_torsion]
                            if atom_index == 3:
                                total_twist = (conformation[current_torsion] + ring_flip_list[1][0]) % 3
                            if atom_index == 4:
                                total_twist = (conformation[current_torsion] + ring_flip_list[1][1]) % 3
                            if atom_index == 5:
                                total_twist = (conformation[current_torsion] + ring_flip_list[1][2]) % 3
                        if conformation[current_torsion] == 2:
                            total_twist = atominfo(ia, 3)

            # If the attached atom on ring atom
            # In this version code, do not consider their movement
            elif (atominfo(ia, 0) not in ringatoms) and (atominfo(ia, 1) in ringatoms):
                ring_index = -1
                for ring in inchi_ring_order:
                    ring_index += 1
                    if atominfo(ia, 1) in ring:
                        attached_atom_index = ring.index(atominfo(ia, 1))
                        break
                #print("ring_index:", ring_index)
                #print("attached_atom_index:", attached_atom_index)
                #print("ring_torsion:", ring_torsion)
                current_ring=ring_torsion[ring_index]
                if ring_twist_statue[ring_index] == 0:
                    if conformation[current_ring] == 0:
                        total_twist = atominfo(ia, 3)
                    if conformation[current_ring] == 1:
                        total_twist = atominfo(ia, 3)
                    if conformation[current_ring] == 2:
                        # test atom_index=3
                        if attached_atom_index==1:
                            total_twist = (atominfo(ia, 3)+2)%3
                        if attached_atom_index==2:
                            total_twist = (atominfo(ia, 3)+1)%3
                        if attached_atom_index==3:
                            total_twist = (atominfo(ia, 3)+2)%3
                        if attached_atom_index==4:
                            total_twist = (atominfo(ia, 3)+1)%3
                        if attached_atom_index==5:
                            total_twist = (atominfo(ia, 3)+2)%3
                        #print("current_ring:", current_ring)
                        #print("ring_twist_statue[ring_index]:", ring_twist_statue[ring_index])
                        #print("current_ring_chekcing_atominfo(ia, 3):", atominfo(ia, 3))
                        #print("current_ring_chekcing_total_twist:", total_twist)
                        #print("current_ring_conformation[current_ring]:", conformation[current_ring])
                if ring_twist_statue[ring_index] == 1:
                    if conformation[current_ring] == 0:
                        total_twist = atominfo(ia, 3)
                    if conformation[current_ring] == 1:
                        # test atom_index=3
                        if attached_atom_index == 1:
                            total_twist = (atominfo(ia, 3)+1)%3
                        if attached_atom_index == 2:
                            total_twist = (atominfo(ia, 3)+2)%3
                        if attached_atom_index == 3:
                            total_twist = (atominfo(ia, 3)+1)%3
                        if attached_atom_index == 4:
                            total_twist = (atominfo(ia, 3)+2)%3
                        if attached_atom_index==5:
                            total_twist = (atominfo(ia, 3)+1)%3
                        #print("current_ring:", current_ring)
                        #print("ring_twist_statue[ring_index]:", ring_twist_statue[ring_index])
                        #print("current_ring_chekcing_atominfo(ia, 3):", atominfo(ia, 3))
                        #print("current_ring_chekcing_total_twist:", total_twist)
                        #print("current_ring_conformation[current_ring]:", conformation[current_ring])
                    if conformation[current_ring] == 2:
                        total_twist = atominfo(ia, 3)

            else:
                total_twist = (atominfo(ia, 3) + conformation[current_torsion]) % 3
                # print("TEST!!!total_twist:", total_twist)
            # print("### total_twist, line 119,total_twist, atominfo(ia,3), conformation[current_torsion], ia, current_torsion")
            # print("### total_twist, line 119",total_twist, atominfo(ia,3), conformation[current_torsion], ia, current_torsion)
            # find direction of penultimate atom
            direction = atominfo(inchi2wn[atominfo(inchi2wn[atominfo(ia, 1)], 1)], 4)
            # print("direction:", direction)
            # find direction of attached atom
            old_direction = atominfo(inchi2wn[atominfo(ia, 1)], 4)
            # print("old_direction:", old_direction)
            if twisting[total_twist * 64 + old_direction * 8 + direction] > 8:
                # print("conf_num, conformation",conf_num,conformation)
                # print("current_torsion",current_torsion,"num_repeat_torsions",num_repeat_torsions)
                # print("atom, attach_atom",ia,atominfo(ia,0),atominfo(ia,3),inchi2wn[atominfo(ia,1)])
                # print("t-1, t, t+1",twisting[total_twist*64+old_direction*8+direction-1],twisting[total_twist*64+old_direction*8+direction],twisting[total_twist*64+old_direction*8+direction+1])
                # print("twisting",twisting[total_twist*64+old_direction*8+direction],"tt,od,d",total_twist,old_direction,direction,"combine",total_twist*64+old_direction*8+direction)
                # printatoms()
                # printxyzatoms()
                print("****** twisting too big ****** ", "current_atom_marker:", current_atom_marker, "old_direction:",
                      old_direction, "direction:", direction, "total_twist:", total_twist)
                print(twisting[total_twist * 64 + old_direction * 8 + direction])
            newdirection = twisting[total_twist * 64 + old_direction * 8 + direction]
            # print("newdirection:", newdirection)
            iattach_atom = inchi2wn[atominfo(ia, 1)]
            # print("starting on atom",ia)
            # print(atominfo(ia,0),atominfo(ia,1),atominfo(ia,2),atominfo(ia,3),atominfo(ia,4),atominfo(ia,5),atominfo(ia,6),atominfo(ia,7),atominfo(ia,8))
            # print("dir, old, new, iattach_atom",direction,old_direction,newdirection,iattach_atom,"tt",total_twist)
            xcoord = atominfo(iattach_atom, 5) + directions[newdirection * 3]
            ycoord = atominfo(iattach_atom, 6) + directions[newdirection * 3 + 1]
            zcoord = atominfo(iattach_atom, 7) + directions[newdirection * 3 + 2]
            # print("old x,y,z",atominfo(iattach_atom,5),atominfo(iattach_atom,6),atominfo(iattach_atom,7),"directions",newdirection,";",directions[newdirection*3],directions[newdirection*3+1],directions[newdirection*3+2])
            # print("new x,y,z",xcoord,ycoord,zcoord)
            setatominfo(ia, 4, newdirection)
            setatominfo(ia, 5, xcoord)
            setatominfo(ia, 6, ycoord)
            setatominfo(ia, 7, zcoord)
            if ia > 2:
                # No energy contribution from first three atoms, which are joined to each other
                # Now calculate the contributions to the energy
                # compare distance with all previous atoms
                # However, can we avoid Pythagoras, because distances so close?
                # overlap - idential coordinates
                # adjacent - one bond length (root 3) away, and not bonded: dx+dy+dz=3, but so does 0,0,3
                # 1,5 - two bond lengths: 2root2, and not bonded: dx+dy+dz=4
                # gauche - three bonds (does bonding matter? probably not) (0,0,0 via 2,2,0 to 1,3,-1) root11
                # make a list of atoms sorted by x, y, z?
                # Should butane and methylpropane have the same energy? Give an energy for methyl, methylene, methine, C ?
                # Compile a list of connections (1,3; 1,4) for each atom?
                # favourable non-bonded interactions?
                # print("worrying about energy")
                current_atom_energy = 0
                ea = 0
                while ea < ia:
                    # for ea in range(0, ia):
                    # if ea>=ia:
                    #  print("ea>=ia",ea,ia)
                    xea = atominfo(ea, 5) - xcoord
                    yea = atominfo(ea, 6) - ycoord
                    zea = atominfo(ea, 7) - zcoord
                    distance2 = xea * xea + yea * yea + zea * zea
                    # print("distance2",ia,ea,"d2",distance2,"x,y,z",xea,yea,zea)
                    if distance2 == 0:
                        num_overlap = num_overlap + 1
                        # no point in going on - energy too high
                        # if it is overlap, then the total energy will be marked lower than -1000
                        energy_total = -1000 - ia
                        test_atom = current_atom_marker
                        while test_atom < number_of_skeleton:
                            if atominfo(test_atom, 8) < current_torsion:
                                ia = current_atom_marker
                                break
                            else:
                                ia = number_of_skeleton
                                ea = ia
                            test_atom += 1
                        current_atom_energy += energy_overlap

                    if distance2 == 3 and not inchi2wn[atominfo(ia, 1)] == ea:
                        # if this adjacent situation happens between the first and the final ring atom on the same ring,
                        # this situation will not be regarded as adjacent
                        # then need to skip this energy calculation round
                        if atominfo(ia, 0) in final_atom_InChI and atominfo(ea, 0) in first_atom_InChI:
                            ia_index = final_atom_InChI.index(atominfo(ia, 0))
                            ia_index = first_atom_InChI.index(atominfo(ea, 0))
                            if ia_index == ia_index:
                                ea += 1
                                continue
                        else:
                            # pass
                            # if it is an alcohol oxygen interact with the other oxygen
                            if (ia in oxygen_atom) and (ea in oxygen_atom):
                                if (wn2inchi[ia] in alcohol_oxygen_atom) and (wn2inchi[ea] in alcohol_oxygen_atom):
                                    current_atom_energy += energy_adjacent_OH_OH
                                else:
                                    current_atom_energy += energy_adjacent_O_OH
                            # else it is adjacent
                            else:
                                # print("adjacent",ia,ea,atominfo(ia,1),inchi2wn[atominfo(ia,1)])
                                num_adjacent = num_adjacent + 1
                                # if it is adjacent
                                # no point in going on - energy too high
                                # if it is adjacent, then the total energy will be marked lower than -1000
                                energy_total = -1000 - ia
                                test_atom = current_atom_marker
                                while test_atom < number_of_skeleton:
                                    if atominfo(test_atom, 8) < current_torsion:
                                        ia = current_atom_marker
                                        break
                                    else:
                                        ia = number_of_skeleton
                                        ea = ia
                                    test_atom += 1
                                current_atom_energy += energy_adjacent

                    if distance2 == 8:
                        # two bonds away, or 1,5 diaxial?
                        # if two bonds away, atominfo(ia,1) must be attached to ea
                        if not (atominfo(ia, 1) == atominfo(ea, 1) or inchi2wn[
                            atominfo(inchi2wn[atominfo(ia, 1)], 1)] == ea):
                            # if they were to be bound to each other
                            # print("1,3 connection",ia,ea,";",atominfo(ia,1),atominfo(ea,1),";",inchi2wn[atominfo(inchi2wn[atominfo(ia,1)],1)],ea)
                            # Good moment to look for gauche interactions???
                            # What is attached to ea, other than atominfo(ia,1) ?
                            # Might be helpful to make a list in advance
                            num_15diaxial = num_15diaxial + 1
                            if (ia in oxygen_atom) or (ea in oxygen_atom):
                                if (ia in oxygen_atom) and (ea in oxygen_atom):
                                    # num_15diaxial_O_O = num_15diaxial_O_O + 1
                                    if (wn2inchi[ia] in alcohol_oxygen_atom) and (wn2inchi[ea] in alcohol_oxygen_atom):
                                        current_atom_energy += energy_15diaxial_OH_OH
                                    elif (wn2inchi[ia] in alcohol_oxygen_atom) or (wn2inchi[ea] in alcohol_oxygen_atom):
                                        current_atom_energy += energy_15diaxial_O_OH
                                    else:
                                        current_atom_energy += energy_15diaxial_O_O
                                else:
                                    # num_15diaxial_C_O = num_15diaxial_C_O + 1
                                    if (wn2inchi[ia] in alcohol_oxygen_atom) and (wn2inchi[ea] in alcohol_oxygen_atom):
                                        current_atom_energy += energy_15diaxial_C_OH
                                    else:
                                        current_atom_energy += energy_15diaxial_C_O
                            else:
                                num_15diaxial_C_C = num_15diaxial_C_C + 1
                                current_atom_energy += energy_15diaxial
                        # in 1,3 connection case, if they are alcohol oxygen interact with the other osygen, then need to add new motif
                        else:
                            num_13connect += 1
                            if (ia in oxygen_atom) and (ea in oxygen_atom):
                                # if both are alcohol oxygen
                                if (wn2inchi[ia] in alcohol_oxygen_atom) and (wn2inchi[ea] in alcohol_oxygen_atom):
                                    current_atom_energy += energy_13connect_OH_OH
                                # else one of atoms is ether oxygen
                                else:
                                    current_atom_energy += energy_13connect_O_OH

                    if distance2 == 11:
                        # possible gauche interaction; only if connected??? Rather unusual if not connected; Risk it!
                        num_gauche = num_gauche + 1
                        if (ia in oxygen_atom) or (ea in oxygen_atom):
                            if (ia in oxygen_atom) and (ea in oxygen_atom):
                                # num_gauche_O_O = num_gauche_O_O + 1
                                if (wn2inchi[ia] in alcohol_oxygen_atom) and (wn2inchi[ea] in alcohol_oxygen_atom):
                                    current_atom_energy += energy_gauche_OH_OH
                                elif (wn2inchi[ia] in alcohol_oxygen_atom) or (wn2inchi[ea] in alcohol_oxygen_atom):
                                    current_atom_energy += energy_gauche_O_OH
                                else:
                                    current_atom_energy += energy_gauche_O_O
                            else:
                                # num_gauche_C_O = num_gauche_C_O + 1
                                if (wn2inchi[ia] in alcohol_oxygen_atom) or (wn2inchi[ea] in alcohol_oxygen_atom):
                                    current_atom_energy += energy_gauche_C_OH
                                else:
                                    current_atom_energy += energy_gauche_C_O
                        else:
                            num_gauche_C_C = num_gauche_C_C + 1
                            current_atom_energy += energy_gauche

                    if distance2 == 16:
                        num_axeq = num_axeq + 1
                        if (ia in oxygen_atom) or (ea in oxygen_atom):
                            if (ia in oxygen_atom) and (ea in oxygen_atom):
                                # num_axeq_O_O = num_axeq_O_O + 1
                                if (wn2inchi[ia] in alcohol_oxygen_atom) and (wn2inchi[ea] in alcohol_oxygen_atom):
                                    current_atom_energy += energy_axeq_OH_OH
                                elif (wn2inchi[ia] in alcohol_oxygen_atom) or (wn2inchi[ea] in alcohol_oxygen_atom):
                                    current_atom_energy += energy_axeq_O_OH
                                else:
                                    current_atom_energy += energy_axeq_O_O
                            else:
                                # num_axeq_C_O = num_axeq_C_O + 1
                                if (wn2inchi[ia] in alcohol_oxygen_atom) or (wn2inchi[ea] in alcohol_oxygen_atom):
                                    current_atom_energy += energy_axeq_C_OH
                                else:
                                    current_atom_energy += energy_axeq_C_O
                        else:
                            num_axeq_C_C = num_axeq_C_C + 1
                            current_atom_energy += energy_axeq
                    ea += 1
                # if the ring is broken after the searing process, then the structure need to be throw away
                # if it is not overlap or adjacent, then the process will continue
            if energy_total > -1000:
                energy_total += current_atom_energy
                # print("current_atom_marker,ia,current_atom_energy",current_atom_marker,ia,current_atom_energy,energy_total)
            else:
                return energy_total
            setatominfo(current_atom_marker, 2, current_atom_energy)
        else:
            energy_total += atominfo(current_atom_marker, 2)
            # print("Adding energy from previous structure",ia,current_atom_marker,atominfo(current_atom_marker,2),energy_total)
            # printatoms()
            # print("===")
        ia += 1
        # if energy_total > -1:
        # gauche interaction is about 2 kJ/mol (two butane conformations)
        # 1,5 diaxial is at least 20 kJ/mol - this is a relaxed diaxial cyclohexane, so should be more
        # overlap and adjacent are very large - no more detail needed
        # print("Energy_total",energy_total)
        # energy_total=energy_gauche*num_gauche+energy_15diaxial*num_15diaxial+energy_adjacent*(num_overlap+num_adjacent)
        # print("Energy of",conf_num,":",num_overlap,num_adjacent,num_15diaxial,num_gauche,num_axeq,"overall",energy_total,conformation)
        # f1.write("\n")
        # f1.write("Energy of conformer"+" "+str(conf_num)+":"+" "+str(num_overlap)+" "+str(num_adjacent)+" "+str(num_axeq_C_C)+" "+str(num_gauche_C_C)+" "+str(num_15diaxial_C_C)+" "+str(num_axeq_C_O)+" "+str(num_gauche_C_O)+" "+str(num_15diaxial_C_O)+" "+str(num_axeq_O_O)+" "+str(num_gauche_O_O)+" "+str(num_15diaxial_O_O)+" "+"overall_energy:"+" "+str(energy_total)+" "+str(conformation))
        # f1.write("\n")
        # printatoms()
        # else:
        # print("*** Energy overload ***",conf_num,energy_total,-1-energy_total,conformation)
    print()
    return energy_total


#################################################################
## Start and calculate preliminary quantities:
## Number of atoms, stereogenic centres, rotatable bonds, etc
#################################################################

print("###########################################")
print("############ Diamond Energy II ############")
print("###########################################")
print(sys.argv[1])


m = Chem.MolFromInchi(sys.argv[1])
number_of_skeleton = m.GetNumAtoms()
number_of_carbons = atom_counts(m, 'C')
number_of_oxygens = atom_counts(m, 'O')
print("Number of skeleton: ", number_of_skeleton)
print("Number of carbons:", number_of_carbons)
print("Number of oxygens:", number_of_oxygens)



ringatoms = [atom + 1 for ring in Chem.GetSymmSSSR(m) for atom in ring]
print("ringatoms:", ringatoms)
oxygen_atom = [atom.GetIdx() + 1 for atom in m.GetAtoms() if atom.GetSymbol() == 'O']
print("oxygen_atom:", oxygen_atom)
alcohol_oxygen_atom = get_alcohol_oxygens(m)
print("alcohol_oxygen_atom:", alcohol_oxygen_atom)




change_enantiomer = 1
len_list_numbers_stereogenic_centres = 0
if sys.argv[1].find("/t") > -1:
    stereogenic_centres = sys.argv[1].split("/")[4].split("t")[1]
    number_stereogenic_centres = rdMolDescriptors.CalcNumAtomStereoCenters(m)
    print("number_stereogenic_centres:", number_stereogenic_centres)
    list_stereogenic_centres = stereogenic_centres.split(",")
    print("list_stereogenic_centres:", list_stereogenic_centres)
    initial_stereoinfo=list_stereogenic_centres.copy()
    list_numbers_stereogenic_centres = array.array("i")
    for i in range(0, len(list_stereogenic_centres)):
        list_numbers_stereogenic_centres.append(int(list_stereogenic_centres[i].replace("+", "").replace("-", "")))
    len_list_numbers_stereogenic_centres = len(list_numbers_stereogenic_centres)
    print(number_stereogenic_centres, "stereogenic_centres_initial:", stereogenic_centres)
    print("list_numbers_stereogenic_centres_initial", list_numbers_stereogenic_centres)
    print("len_list_numbers_stereogenic_centres", len(list_numbers_stereogenic_centres))
    if sys.argv[1].find("/m1") > -1:
        change_enantiomer = -1
else:
    stereogenic_centres = ''
    number_stereogenic_centres = 0



connectivity = sys.argv[1].split("/")[2].split("c")[1]
print("connectivity ", connectivity)
list_of_connectivity = connectivity.replace(")", ")-").replace("(", "-(").replace(",", "-,").split("-")
print("list_of_connectivity ", list_of_connectivity)
# remove duplicates since introduce ring and get all atoms InChI order
temp_list_of_connectivity = connectivity.replace(")", ")-").replace("(", "-(").replace(",", "-,").split("-")
# print("temp_list_of_connectivity ",temp_list_of_connectivity)
temp_all_atom_order = connectivity.replace(")", "-").replace("(", "-").replace(",", "-").split("-")
# print("temp_all_atom_order:", temp_all_atom_order)
all_atom_order = []
for atom in temp_all_atom_order:
    if atom not in all_atom_order:
        all_atom_order.append(atom)
print("all_atom_InChI_order:", all_atom_order)
list_of_connectivity_int = [extract_and_convert(x) for x in all_atom_order if extract_and_convert(x) is not None]
print("list_of_connectivity_int", list_of_connectivity_int)



if sys.argv[1].find("/t") > -1:
    # sort the stereocenter order based on list of connectivity int order
    sorted_indices = sorted(range(len(list_numbers_stereogenic_centres)), key=lambda k: list_of_connectivity_int.index(list_numbers_stereogenic_centres[k]))
    list_numbers_stereogenic_centres = [list_numbers_stereogenic_centres[i] for i in sorted_indices]
    list_stereogenic_centres = [list_stereogenic_centres[i] for i in sorted_indices]
    print("list_stereogenic_centres_sort:", list_stereogenic_centres)
    print("list_numbers_stereogenic_centres", list_numbers_stereogenic_centres)
else:
    pass



# Remove the repeated first atom from list_of_connectivity
first_atom_InChI = []
time = 0
for atom_info in all_atom_order:
    if int(atom_info) in ringatoms:
        time += 1
        if time % 6 == 1:
            # print(time)
            first_atom_InChI.append(int(atom_info))
print("first_atom_InChI:", first_atom_InChI)
# re-sort the order of the temp_inchi_ring_order based on the first_atom_InChI
inchi_ring_order = []
i = 0
num = 0
for atom_info in all_atom_order:
    if int(atom_info) in ringatoms:
        if num % 6 == 0:
            inchi_ring_order.append([int(atom_info)])
            num += 1
        else:
            inchi_ring_order[i].append(int(atom_info))
            num += 1
        if num % 6 == 0:
            i += 1
    # print("num:", num)
    # print("i:", i)
print("inchi_ring_order:", inchi_ring_order)




# inchi2workingnumber translates the inchi atom number to the position in the array
# InChI numbers are 1 - n; working numbers are 0 to n-1
# The atoms are numbered as InChI numbers, and as the position in the atom list
# For example, if an InChI has c1-6(2)7(3,4)5
# Then: inchi2wn[6]=1
# and:  wn2inchi[1]=6
# Then: inchi2wn[7]=3
# and:  wn2inchi[3]=7
inchi2wn = array.array('i')
wn2inchi = array.array('i')
tmp_inchi2workingnumber = all_atom_order
for workingposition in range(0, number_of_skeleton):
    if int(tmp_inchi2workingnumber[workingposition]) not in wn2inchi:
        wn2inchi.append(int(tmp_inchi2workingnumber[workingposition]))
        inchi2wn.append(0)
    else:
        pass
inchi2wn.append(0)
for workingposition in range(0, number_of_skeleton):
    inchi2wn[wn2inchi[workingposition]] = workingposition
tmp_inchi2workingnumber.clear()
# print("wn2inchi", wn2inchi)
# print("inchi2wn", inchi2wn[1: number_of_carbons+1])



# Number of methyl groups; may not need this, but never know...
# number_methyl=int(sys.argv[1].split("/")[3].split("h")[1].split(",").pop().split("-").pop().split("H")[0])
# print(number_methyl,"methyls")

## using NumRotatableBonds function in rdkit to count rotatable bonds, which ignore bonds on ring automatically
## let the ring bond regarbe ded as a rotatable bond to do ring flex in this version
number_rotatable_bonds = Descriptors.NumRotatableBonds(m) + len(inchi_ring_order)
print("number_rotatable_bonds:", number_rotatable_bonds)



## number_rotatable_bonds=len(connectivity.replace(")","-").split("-"))-3(for alkanes)
number_conformations = 3 ** number_rotatable_bonds
print(number_rotatable_bonds, "rotatable bonds;", number_conformations, "conformations")
conformation = array.array("i")
old_conformation = array.array("i")
conf_number = 0
# print("len(sys.argv)",len(sys.argv),sys.argv)
if len(sys.argv) > 4:
    conf_number = int(sys.argv[4])
tmp_conf_number = conf_number
if number_rotatable_bonds == 0:
    conformation.append(0)
else:
    for i in range(0, number_rotatable_bonds):
        conformation.append(0)
old_conformation = conformation[:]
old_conformation[0] = -1
for i in range(0, number_rotatable_bonds):
    # print("conf_sort_out_setup",number_rotatable_bonds,i,tmp_conf_number,conformation)
    conformation[number_rotatable_bonds - i - 1] = tmp_conf_number % 3
    tmp_conf_number = int(tmp_conf_number / 3)
# print("Conformation",conformation)




# The directions array gives the vector for each step between atoms
directions = array.array('i',
                         [ 1,  1,  1,
                          -1, -1,  1,
                          -1,  1, -1,
                           1, -1, -1,
                          -1, -1, -1,
                           1,  1, -1,
                           1, -1,  1,
                          -1,  1,  1])
# print(directions)




# The twistings array gives the twisted vectors
# The 9 should never occur
# With penultimate direction  pendir
# and attached atom direction attdir
# the new direction will be pendir (untwisted)
# twisting[attdir*8+pendir] to twist up
# twisting[64+attdir*8+pendir] to twist down
# newdirection=twisting[twist*64+old_direction*8+direction]
twisting = array.array('i',
                       [9, 9, 9, 9, 9, 5, 6, 7,
                        9, 9, 9, 9, 4, 9, 6, 7,
                        9, 9, 9, 9, 4, 5, 9, 7,
                        9, 9, 9, 9, 4, 5, 6, 9,
                        9, 1, 2, 3, 9, 9, 9, 9,
                        0, 9, 2, 3, 9, 9, 9, 9,
                        0, 1, 9, 3, 9, 9, 9, 9,
                        0, 1, 2, 9, 9, 9, 9, 9,

                        9, 9, 9, 9, 9, 7, 5, 6,
                        9, 9, 9, 9, 6, 9, 7, 4,
                        9, 9, 9, 9, 7, 4, 9, 5,
                        9, 9, 9, 9, 5, 6, 4, 9,
                        9, 2, 3, 1, 9, 9, 9, 9,
                        3, 9, 0, 2, 9, 9, 9, 9,
                        1, 3, 9, 0, 9, 9, 9, 9,
                        2, 0, 1, 9, 9, 9, 9, 9,

                        9, 9, 9, 9, 9, 6, 7, 5,
                        9, 9, 9, 9, 7, 9, 4, 6,
                        9, 9, 9, 9, 5, 7, 9, 4,
                        9, 9, 9, 9, 6, 4, 5, 9,
                        9, 3, 1, 2, 9, 9, 9, 9,
                        2, 9, 3, 0, 9, 9, 9, 9,
                        3, 0, 9, 1, 9, 9, 9, 9,
                        1, 2, 0, 9, 9, 9, 9, 9])
# print(twisting)


# The ring atoms twist value should be store in lists in advance for building the ring structure
# Once the first three ring 'atoms confirm
# The ring atoms also confirm closing the ring
# ring_twist_list_1=[1, 2, 1] # beta
# ring_twist_list_2=[2, 1, 2] # alpha
ring_twist_list = [[1, 2, 1], [2, 1, 2]]


# Correspondingly, there will be two different ring flip relative twist values lists
# for ring_twist_list_1 and ring_twist_list_2 in their situation
ring_flip_list = [[0, 2, 0], [0, 1, 0]]

# Since ring flex will affect the branch before the third ring atom
# need to jump to rebuild the branch attached to the second ring atom
global affected_branch
affected_branch=[]




#################################################################
# Assemble molecule
# First three atoms are all connected as propane
# but need to distinguish ring atoms
# Each atom has the InChI atom number,
# the InChI# of the atom to which it is attached,
# the bond direction (from list of eight possibilities),
# twist (0,1,2) for multiple atoms on same last bond
# The overall new direction
# x,y,z coordinates (signed integers)
# NB: first atom is one, for easy array indexing
# I#,At,Di,Tw,Nd, x, y, z
#
# Need to index flexible torsion angles
# Each is uniquely defined by terminal atom, but need a list
#################################################################
atom1 = int(list_of_connectivity[0])
atom2 = int(list_of_connectivity[1])

# create a list named ring_twist_statue to record how the ring is built, alpha or beta
# The default statue is ring_twist_list_1, ring_twist_statue value should be 0
ring_twist_statue = []
ring_twist_statue.append(0)

firstatom_dir = 0
# create a list to store each ring bond's twist value,
# Record and store where the rotatable bond is ring flex indicator using ring_torsion list
global ring_torsion
ring_torsion = []
global ring_twist
ring_twist=[]
final_atom_InChI = []



# Record the first ring atom have been looked at
# avoid looking at it second time as there will be two first ring number shown in InChI string for ring system
processed_atoms = set()
attach = [atom1, atom2]
atom = array.array('i',
                   [int(list_of_connectivity[0]), 0, 0, 0, 5, 0, 0, 0, -1,
                    int(list_of_connectivity[1]), wn2inchi[0], 0, 0, 0, 1, 1, 1, -1])
# If the first two atoms are not in a ring, so could just re-use the previous method to build the propane
if (atom1 and atom2) not in inchi_ring_order[0]:
    #print("Checking!")
    if not list_of_connectivity[2].isnumeric():
        if list_of_connectivity[3].find(",") > -1:
            atom.extend([int(list_of_connectivity[2].replace("(", "").replace(")", "")), wn2inchi[1], 0, 1, 7, 0, 2, 2, -1])
            atom.extend([int(list_of_connectivity[3].replace("(", "").replace(")", "").replace(",", "")), wn2inchi[1], 0, 2, 6, 2, 0, 2, -1])
            if list_of_connectivity[3].find(")") == -1:
                attach.append(int(list_of_connectivity[3].replace("(", "").replace(")", "")))
            else:
                # wait to test, atom 4 might be on torsion angle 0??
                atom.extend([int(list_of_connectivity[4]), wn2inchi[1], 0, 0, 5, 2, 2, 0, -1])
                attach.append(int(list_of_connectivity[4]))

        else:
            atom.extend([int(list_of_connectivity[2].replace("(", "").replace(")", "")), wn2inchi[1], 0, 1, 7, 0, 2, 2, -1])
            if list_of_connectivity[2].find(")") > -1:
                atom.extend([int(list_of_connectivity[3]), wn2inchi[1], 0, 0, 5, 2, 2, 0, -1])
                attach.append(int(list_of_connectivity[3]))
                if int(list_of_connectivity[3]) in ringatoms:
                    firstatom_dir = 5
                    twist = 0
                    current_twist = twist
                    ring_twist.append([current_twist])
                    # Mark the atom as processed
                    processed_atoms.add(int(list_of_connectivity[3]))
            else:
                attach.append("(")
                attach.append(int(list_of_connectivity[2].replace("(", "").replace(")", "")))
    else:
        atom.extend([int(list_of_connectivity[2]), wn2inchi[1], 0, 0, 5, 2, 2, 0, -1])
        attach.append(int(list_of_connectivity[2]))
        # Once the first ring atom is recognized
        # need to record the info of the current ring when the first ring atom is found
        if int(list_of_connectivity[2]) in ringatoms:
            firstatom_dir = 5
            twist = 0
            current_twist = twist
            ring_twist.append([current_twist])
            # Mark the atom as processed
            processed_atoms.add(int(list_of_connectivity[2]))


# If the first three atoms involve a ring system, need to add ring information at the beginning
else:
    if atom1 in ringatoms:
        ## If atom1 is the first atom of the ring, it represents this is a cyclohexane
        ## So the first three ring atoms are built as propane
        ## Here add the first atom's information, including direction/twist value/append ring_atom list
        firstatom_dir = 5
        twist = 0
        current_twist = twist
        # add first ring atom twist value
        ring_twist.append([current_twist])
        # Mark the atom as processed
        processed_atoms.add(atom1)
        ## Here add the the second atom's information, including direction/twist value/append ring_atom list
        secondatom_dir = 0
        ring_twist[0].append(0)
        ## Here add the third atom's information, including direction/twist value/append ring_atom list/atom_list/attach_list
        atom.extend([int(list_of_connectivity[2]), wn2inchi[1], 0, 0, 5, 2, 2, 0, -1])
        attach.append(int(list_of_connectivity[2]))
        ring_twist[0].append(0)
    else:
        ## If atom2 is the first atom of the ring, it represents this is methyl-cyclohexane
        ## need to add its ring information
        firstatom_dir = 0
        twist = 0
        current_twist = twist
        ring_twist.append([current_twist])
        # Mark the atom as processed
        processed_atoms.add(atom2)
        if not list_of_connectivity[2].isnumeric():
            ## means there are two methyl at the first ring atom
            ## So atom4 is the second ring atom, need to add its info
            atom.extend([int(list_of_connectivity[2].replace("(", "").replace(")", "")), wn2inchi[1], 0, 0, 5, 2, 2, 0, -1])
            if not list_of_connectivity[3].isnumeric():
                pass
            else:
                atom.extend([int(list_of_connectivity[3]), wn2inchi[1], 0, 1, 7, 0, 2, 2, -1])
                attach.append(int(list_of_connectivity[3]))
                secondatom_dir = 5
                ring_twist[0].append(1)
        else:
            atom.extend([int(list_of_connectivity[2]), wn2inchi[1], 0, 0, 5, 2, 2, 0, -1])
            attach.append(int(list_of_connectivity[2]))
            secondatom_dir = 5
            ring_twist[0].append(0)

## current_position is the end of the fixed atoms and the start of the flexible molecule
## NB - 9 is the number of integers per atom
current_position = int(len(atom) / 9)
print("current_position ", current_position)
#print("atom_info:", atom)






#############################################################################
# Now add rest of structure.
# Assume, for the moment, structure zero, so all torsion angles are 180
# Probably hold conformation as base three string:
#############################################################################
# print("attach",attach)

# To sort out each atom, need to find:
# old_direction = actual direction of attached atom: atominfo(inchi2wn[attach_atom],4)
# direction = actual direction of penultimate atom: atominfo(inchi2wn[atominfo(inchi2wn[attach_atom],1)],4)
# new direction = direction + affect of twist, so actual direction = newdirection=twisting[twist*64+old_direction*8+direction]

highest_torsion_angle = -1
# list_of_torsion_angles=[]
# list_of_torsion_angles.append(current_torsion_angle)
# print("list_of_torsion_angles",list_of_torsion_angles)
# The array conformation defines whether the angle should be 180, 60 or -60
# conformation[list_of_torsion_angles[current_torsion_angle]] will be zero, one or two


for workingposition in range(current_position, len(list_of_connectivity)):
    # print("workingposition", wn2inchi[workingposition], list_of_connectivity[workingposition],"  ",attach)
    # print("workingposition", list_of_connectivity[workingposition])
    # print("TEST", type(list_of_connectivity[workingposition]))
    current_atom = int(list_of_connectivity[workingposition].replace("(", "").replace(")", "").replace(",", ""))
    #print("current_atom", current_atom)
    # if it is a first ring atom, and has been built before
    # need to skip this round
    if current_atom in processed_atoms:
        continue

    # if there is ring atom
    # get ring index and atom index for further used
    ring_lookup = {atom: idx for idx, ring in enumerate(inchi_ring_order) for atom in ring}
    if current_atom in ring_lookup:
        ring_index = ring_lookup[current_atom]
        atom_index = inchi_ring_order[ring_index].index(current_atom)
    else:
        pass


    # if list_of_connectivity[workingposition]=="10)":
    # print("CHECKING")
    # Choose arbitrary conformation for base structure; number zero? As extended as possible
    # print("direction",direction,"x,y,z",directions[direction*3],directions[direction*3+1],directions[direction*3+2])
    # What atom (InChI#) is this atom attached to? Put in attach_atom
    if list_of_connectivity[workingposition].isnumeric():
        attach_atom = int(attach[-1])
        attach.append(current_atom)
        old_direction = atominfo(inchi2wn[attach_atom], 4)
        penultimate_atom = atominfo(inchi2wn[attach_atom], 1)
        direction = atominfo(inchi2wn[penultimate_atom], 4)
        twist = 0
        newdirection = direction
        # print("This branch twist:", twist)
        # First to identify whether the current atom is on a ring
        if current_atom in ringatoms:
            # Then to identify this ring atom's index on the ring
            # Creating a dictionary for faster lookups
            # print("atom_index:", atom_index)
            # The first atom found on the ring,
            # treat it the same as other non-ring atoms
            if atom_index == 0 and (current_atom not in processed_atoms):
                if (ring_index > 0) and (current_atom not in processed_atoms):
                    ring_twist_statue.append(0)
                firstatom_dir = newdirection
                # Follow its previous twist value, that is 0 in this condition
                current_twist = twist
                ring_twist.append([current_twist])
                # Mark the atom as processed
                processed_atoms.add(current_atom)
            # The first three ring-atom will maintain previous method
            elif atom_index == 1:
                secondatom_dir = newdirection
                current_twist = twist
                #print("ring_twist:", ring_twist)
                ring_twist[ring_index].append(current_twist)
            # The first three ring-atom will maintain previous method
            elif atom_index == 2:
                #print("current_atom:", current_atom)
                current_twist = twist
                ring_twist[ring_index].append(current_twist)
            # The forth ring atom onward will need special twist to close the ring
            elif atom_index == 3:
                #print("ring_twist", ring_twist)
                #print("ring_index:", ring_index)
                #print("ring_twist[ring_index]", ring_twist[ring_index])
                current_twist = (ring_twist[ring_index][2]+ring_twist_list[0][0]) % 3
                twist = current_twist
                ring_twist[ring_index].append(current_twist)
                newdirection = twisting[twist * 64 + old_direction * 8 + direction]
            # The fifth ring atom will need special twist to close the ring
            elif atom_index == 4:
                current_twist = (ring_twist[ring_index][2]+ring_twist_list[0][1]) % 3
                twist = current_twist
                ring_twist[ring_index].append(current_twist)
                newdirection = twisting[twist * 64 + old_direction * 8 + direction]
            elif atom_index == 5:
                current_twist = (ring_twist[ring_index][2]+ring_twist_list[0][2]) % 3
                twist = current_twist
                ring_twist[ring_index].append(current_twist)
                newdirection = twisting[twist * 64 + old_direction * 8 + direction]
                # Record the last ring atom, and further used to justify whether it connect to the first ring
                final_atom_InChI.append(current_atom)
        else:
            #  Here consider the branch atoms that link to ring atoms
            #  If current atom is attached to an atom on the ring,
            #  this related information also neesd to record
            #  for justify and make cosrrection for the stereocenter on the ring
            if attach_atom in ringatoms:
                # Then to identify this attached atom's index on the ring
                # Creating a dictionary for faster lookups
                ring_lookup = {atom: idx for idx, ring in enumerate(inchi_ring_order) for atom in ring}
                ring_index = ring_lookup[attach_atom]
                print("ring_index:", ring_index)
                attach_atom_index = inchi_ring_order[ring_index].index(attach_atom)
                print("attach_atom_index:", attach_atom_index)

                # if the atom is attached to the first aring atom，this means there are two branches on the first atom ring
                # let it keeped to be considered later
                # now only consider there is only one branch on a certain ring atom
                if attach_atom_index == 0:
                    twist = 1
                if attach_atom_index == 1:
                    twist = 1
                    affected_branch.append(current_atom)
                if attach_atom_index == 2:
                    twist = 1
                if attach_atom_index == 3:
                    twist = 0
                if attach_atom_index == 4:
                    twist = 0
                if attach_atom_index == 5:
                    real_attach = int(list_of_connectivity[workingposition - 1].replace("(", "").replace(")", "").replace(",", ""))
                    if real_attach != attach_atom:
                        if "(" in attach:  # Check if "(" is in the list
                            top_item = attach.pop()
                            while not top_item == "(":
                                top_item = attach.pop()
                            #print("TEST_attach:", attach)
                            attach_atom = int(attach[-1])
                            #print("TEST_attach_atom:", attach_atom)
                            attach.append(current_atom)
                            #print("TEST_attach_later:", attach)
                            twist = 1
                        else:
                            twist = 1
                    else:
                        twist = 0
                newdirection = twisting[twist * 64 + old_direction * 8 + direction]
            else:
                pass

    else:
        # print("not numeric", workingposition, current_atom)
        # This could be open bracket - so a branch
        # Close bracket, so end of a branch
        # Comma, so alternative branch
        # Or a combination of these
        if list_of_connectivity[workingposition].find("(") > -1:
            # print("Found '('")
            #print("current_atom:", current_atom)
            attach_atom = int(attach[-1])
            attach.append("(")
            attach.append(current_atom)
            # print("current_atom", current_atom)
            # print("attach", attach)
            penultimate_atom = atominfo(inchi2wn[attach_atom], 1)
            direction = atominfo(inchi2wn[penultimate_atom], 4)
            old_direction = atominfo(inchi2wn[attach_atom], 4)

            # First to identify whether the current atom is on a ring
            if current_atom in ringatoms:
                # Then to identify this ring atom's index on the ring
                # Creating a dictionary for faster lookups
                # print("atom_index:", atom_index)
                # The first atom found on the ring,
                # treat it the same as other non-ring atoms
                if atom_index == 0 and (current_atom not in processed_atoms):
                    if (ring_index > 0) and (current_atom not in processed_atoms):
                        ring_twist_statue.append(0)
                    firstatom_dir = newdirection
                    # Follow its previous twist value, that is 0 in this condition
                    twist = 1
                    current_twist = twist
                    ring_twist.append([current_twist])
                    # Mark the atom as processed
                    processed_atoms.add(current_atom)
                # The first three ring-atom will maintain previous method
                if atom_index == 1:
                    secondatom_dir = newdirection
                    twist = 1
                    current_twist = twist
                    ring_twist[ring_index].append(current_twist)
                # The first three ring-atom will maintain previous method
                if atom_index == 2:
                    twist = 0
                    current_twist = twist
                    ring_twist[ring_index].append(current_twist)
                # The forth ring atom onward will need special twist to close the ring
                if atom_index == 3:
                    #print("Test_ringtwist:", ring_twist)
                    current_twist = (ring_twist[ring_index][2] + ring_twist_list[0][0]) % 3
                    twist = current_twist
                    ring_twist[ring_index].append(current_twist)
                    newdirection = twisting[twist * 64 + old_direction * 8 + direction]
                # The fifth ring atom will need special twist to close the ring
                if atom_index == 4:
                    #print("Test_ringtwist:", ring_twist)
                    current_twist = (ring_twist[ring_index][2] + ring_twist_list[0][1]) % 3
                    twist = current_twist
                    ring_twist[ring_index].append(current_twist)
                    newdirection = twisting[twist * 64 + old_direction * 8 + direction]
                if atom_index == 5:
                    #print("Test_ringtwist:", ring_twist)
                    current_twist = (ring_twist[ring_index][2] + ring_twist_list[0][2]) % 3
                    twist = current_twist
                    ring_twist[ring_index].append(current_twist)
                    newdirection = twisting[twist * 64 + old_direction * 8 + direction]
                    # Record the last ring atom, and further used to justify whether it connect to the first ring
                    final_atom_InChI.append(current_atom)


            ## if current atom is attached to an atom on ring,
            ## then the overall twist of current atom it is fixed
            ## as there are only limited certian direction for them to go due to stereo heric
            if (current_atom not in ringatoms) and (attach_atom in ringatoms):
                # identify this atom's index of the ring
                attach_atom_index = inchi_ring_order[ring_index].index(attach_atom)
                # if the atom is attached to the first aring atom，this means there are two branches on the first atom ring
                # let it keeped to be considered later
                # now only consider there is only one branch on a certain ring atom
                if attach_atom_index == 0:
                    twist = 1
                if attach_atom_index == 1:
                    twist = 1
                    affected_branch.append(current_atom)
                if attach_atom_index == 2:
                    twist = 0
                # test
                if attach_atom_index == 3:
                    twist = 0
                if attach_atom_index == 4:
                    twist = 0
                if attach_atom_index == 5:
                    twist = 0
                newdirection = twisting[twist * 64 + old_direction * 8 + direction]
            if (current_atom not in ringatoms) and (attach_atom not in ringatoms):
                twist=1


            # NB this arbitrarily assigns stereochemistry; sort out after assigning whole structure
            # print("workingposition; direction, old; twist",workingposition,"; ",direction, old_direction,"; ",twist)
            if list_of_connectivity[workingposition].find(")") > -1:
                # print(workingposition," attach ",attach)
                attach.pop()
                attach.pop()

        else:
            # print("'(' not found; try ',' amd ')'")
            # print("attach", attach)
            if list_of_connectivity[workingposition].find(",") > -1 or list_of_connectivity[workingposition].find(
                    ")") > -1:
                # Go back to last "("
                #
                # Three situations:
                #  ) - go back to last branch
                #  , - attach to last branch but maintain current
                # ,) - attach to last branch and go back to last branch
                #
                # For example, if attach is [ 1, 13, 17, (, 14, (, 4 ] and current atom is 5
                #
                # If the latest character is 5) then should be attached to 4 and close branch
                #   so attach becomes [ 1, 13, 17, (, 14 ] ; attach_atom is 4
                #
                # If the latest character is ,5 then should be attached to 14 and maintain branch
                #   so attach becomes [ 1, 13, 17, (, 14, (, 5 ] ; attach_atom is 14
                #
                # If the latest character is ,5) then should be attached to 14 and close branch
                #   so attach becomes [ 1, 13, 17, (, 14 ] ; attach_atom is 14
                #
                if int(list_of_connectivity[workingposition].split(")")[0]) in first_atom_InChI:
                    # print("CHECKING")
                    attach_atom = int(attach[-1])
                    # print("attach_atom", attach_atom)
                    top_item = attach.pop()
                    # print("top_item", top_item)
                    while not top_item == "(":
                        top_item = attach.pop()
                        # print("workingposition, attach", list_of_connectivity[workingposition], attach)
                    continue
                else:
                    attach_atom = int(attach[-1])
                    # print("attach_atom", attach_atom)
                    top_item = attach.pop()
                    while not top_item == "(":
                        top_item = attach.pop()
                        # print("list_of_connectivity[workingposition],attach",list_of_connectivity[workingposition],attach)
                if list_of_connectivity[workingposition].find(",") > -1:
                    attach_atom = int(attach[-1])
                    twist = 2
                    if list_of_connectivity[workingposition].find(")") == -1:
                        # comma and not )
                        attach.append("(")
                        attach.append(current_atom)
                else:
                    twist = 0
                # attach should now be sorted out
                penultimate_atom = atominfo(inchi2wn[attach_atom], 1)
                direction = atominfo(inchi2wn[penultimate_atom], 4)
                old_direction = atominfo(inchi2wn[attach_atom], 4)
            else:
                print("There is no 'else:'; should never be here")

    # Find current torsion angle
    if highest_torsion_angle < 0:
        if current_atom not in ringatoms:
            current_torsion_angle = 0
            highest_torsion_angle = 0
        else:
            if atom_index == 0 or atom_index == 1 or atom_index == 2:
                current_torsion_angle = 0
                highest_torsion_angle = 0
                if atom_index == 2:
                    ring_torsion.append(current_torsion_angle)
            else:
                current_torsion_angle = -1
                highest_torsion_angle = -1
    else:
        # New torsion, or re-use old one? Has attached atom been used before?
        # print("torsions",current_torsion_angle,highest_torsion_angle,attach_atom,"current_position,workingposition",current_position,workingposition)
        test_current_torsion_angle = -1
        for check_at in range(current_position, int(len(atom) / 9)):
            # print("check_at", check_at)
            # print("len_atom", len(atom))
            if atominfo(check_at, 1) == attach_atom:
                # reuse torsion
                test_current_torsion_angle = atominfo(check_at, 8)
                current_torsion_angle = atominfo(check_at, 8)
                if current_atom in ringatoms:
                    if atom_index == 2:
                        ring_torsion.append(current_torsion_angle)

        if test_current_torsion_angle < 0:
            if atominfo(inchi2wn[attach_atom], 8) == -1:
                current_torsion_angle = 0
            else:
                # Need new torsion angle
                ## if atoms are on a ring,
                ## then (the second ring atom - the fianl ring atom) shoudl be on the same torsional angle
                ## so first need to identify whether current atom is on ring
                ## if so, then get the atom's index on the ring
                ## this can help us to decide which torsion angle it could be classified into
                if current_atom in ringatoms:
                    ## if current atom is the first two atom on ring, need new torsion angle
                    # if atom_index==0 or atom_index==1:
                    # if atom_index == 0:
                    ###### TEST for the first three atoms, all added new torsion
                    if atom_index < 3:
                        highest_torsion_angle += 1
                        current_torsion_angle = highest_torsion_angle
                        if atom_index == 2:
                            ring_torsion.append(current_torsion_angle)
                    ## The atom on ring reuse the second ring atom's torsion angle value
                    else:
                        # print("ring_index:", ring_index)
                        # print("inchi_ring_order[ring_index][1]:", inchi_ring_order[ring_index][1])
                        # print("test", inchi2wn[inchi_ring_order[ring_index][1]])
                        current_torsion_angle = atominfo(inchi2wn[attach_atom], 8)
                ## need to identify whether this atom is attached to an atom on ring
                ## if the attached atom is on ring, then current atom would not add a new value
                elif (current_atom not in ringatoms) and (attach_atom in ringatoms):
                    ring_lookup = {atom: idx for idx, ring in enumerate(inchi_ring_order) for atom in ring}
                    ring_index = ring_lookup[attach_atom]
                    # print("ring_index:", ring_index)
                    attach_atom_index = inchi_ring_order[ring_index].index(attach_atom)
                    if attach_atom_index == 1:
                        highest_torsion_angle += 1
                        current_torsion_angle = highest_torsion_angle
                    else:
                        current_torsion_angle = atominfo(inchi2wn[attach_atom], 8)
                # elif (current_atom not in ringatoms) and (attach_atom in ringatoms):
                # current_torsion_angle=atominfo(inchi2wn[attach_atom],8)
                # elif (current_atom in ringatoms) and (attach_atom in ringatoms):
                # current_torsion_angle=atominfo(inchi2wn[attach_atom],8)
                else:
                    highest_torsion_angle += 1
                    current_torsion_angle = highest_torsion_angle
                # print("highest_torsion_angle",highest_torsion_angle,"current_torsion_angle",current_torsion_angle)

    # coordinates are those of attached atom, plus direction
    # printatoms()
    # print("twist,old_direction,direction",twist,old_direction,direction)
    if (current_atom not in ringatoms) and (attach_atom not in ringatoms):
        newdirection = twisting[twist * 64 + old_direction * 8 + direction]
    else:
        pass
    # if len(inchi_ring_order) > 1:
    # first_atom_num=inchi_ring_order[0]
    # if current_atom==first_atom_num:
    # first_atom_time+=1
    # if first_atom_time > 1:
    # pass
    # else:
    # print("xyz coords",attach_atom,inchi2wn[attach_atom],"old, direction, new, twist",old_direction,direction,newdirection,twist)
    # print("test")
    # print(current_atom)
    # print("newdirection:", newdirection)
    xcoord = atominfo(inchi2wn[attach_atom], 5) + directions[newdirection * 3]
    ycoord = atominfo(inchi2wn[attach_atom], 6) + directions[newdirection * 3 + 1]
    zcoord = atominfo(inchi2wn[attach_atom], 7) + directions[newdirection * 3 + 2]
    # we have a sidechain; find post-sidechain atom
    # atom.extend([current_atom, direction, twist, attach_atom, xcoord, ycoord, zcoord], conformation)
    # if not attach_atom in list_of_torsion_angles:
    #  list_of_torsion_angles.append(attach_atom)
    # current_torsion_angle = list_of_torsion_angles.index(attach_atom)
    # print("current_torsion_angle",current_torsion_angle,list_of_torsion_angles)
    # print("workingposition",workingposition,"; atom extend",current_atom, direction, attach_atom, xcoord, ycoord, zcoord)
    # print(current_atom, attach_atom, direction, twist, newdirection, xcoord, ycoord, zcoord, current_torsion_angle)
    # atom.extend([current_atom, attach_atom, direction, twist, newdirection, xcoord, ycoord, zcoord, current_torsion_angle])
    # print("current_atom:", current_atom)
    # print([current_atom, attach_atom, 0, twist, newdirection, xcoord, ycoord, zcoord, current_torsion_angle])
    atom.extend([current_atom, attach_atom, 0, twist, newdirection, xcoord, ycoord, zcoord, current_torsion_angle])
    #print("atom_info:",atom)
# print("ring_twist:", ring_twist)
# printatoms()
# printxyzatoms()
# print("ring_torsion:",  ring_torsion)
#print("final_atom_InChI:", final_atom_InChI)
global conf_num
#conf_num = 10000000
print("ring_twist_statue:", ring_twist_statue)
#printmol()







#############################################################################
# Check Chiral Centres
# print(list_stereogenic_centres)
# print("list_numbers_stereogenic_centres",list_numbers_stereogenic_centres)
#############################################################################


## need to treat ring stereo centre and alknae stereocenter differently

# if number_stereogenic_centres>0:
if len_list_numbers_stereogenic_centres > 0:
    print("Stereocentres to check")
    print("list_numbers_stereogenic_centres:", list_numbers_stereogenic_centres)
    # printatoms()
    sca = array.array("i")
    for s_centre in range(0, len(list_numbers_stereogenic_centres)):
        del sca
        sca = array.array("i")
        ## if the stereocenter is an alkane atom, all the function will be the same as before version
        if list_numbers_stereogenic_centres[s_centre] not in ringatoms:
            # print("stereocentre: ",list_numbers_stereogenic_centres[s_centre], s_centre)
            # print("inchi_ring_order: ", inchi_ring_order)
            # print("list_numbers_stereogenic_centres",list_numbers_stereogenic_centres[s_centre], list_stereogenic_centres[s_centre])
            # printatoms()
            sca.append(atominfo(inchi2wn[list_numbers_stereogenic_centres[s_centre]], 1))
            smallest = sca[-1]
            largest = sca[-1]
            # for workingposition in range(0, number_of_carbons):
            for workingposition in range(1, number_of_carbons + 1):
                # print("growing sca",list_numbers_stereogenic_centres[s_centre],workingposition,atominfo(inchi2wn[workingposition],1),sca)
                if atominfo(inchi2wn[workingposition], 1) == list_numbers_stereogenic_centres[s_centre]:
                    sca.append(workingposition)
                    if sca[-1] < smallest:
                        smallest = sca[-1]
                    if sca[-1] > largest:
                        largest = sca[-1]
                # print("largest",workingposition,largest,sca,"lnsc",list_numbers_stereogenic_centres)
            second_smallest = largest
            second_largest = smallest
            for workingposition in range(0, len(sca)):
                if sca[workingposition] < second_smallest and not sca[workingposition] == smallest:
                    second_smallest = sca[workingposition]
                if sca[workingposition] > second_largest and not sca[workingposition] == largest:
                    second_largest = sca[workingposition]
            # print("sca, smallest, second_smallest, second_largest, largest",sca, smallest, second_smallest, second_largest, largest)
            if len(sca) > 3:
                smallest = second_smallest

            # print("stereocentre",atominfo(inchi2wn[list_numbers_stereogenic_centres[s_centre]],5),atominfo(inchi2wn[list_numbers_stereogenic_centres[s_centre]],6),atominfo(inchi2wn[list_numbers_stereogenic_centres[s_centre]],7))
            # print("smallest",atominfo(inchi2wn[smallest],5),atominfo(inchi2wn[smallest],6),atominfo(inchi2wn[smallest],7))
            # print("smallest",atominfo(inchi2wn[second_largest],5),atominfo(inchi2wn[second_largest],6),atominfo(inchi2wn[second_largest],7))
            # print("largest",atominfo(inchi2wn[largest],5),atominfo(inchi2wn[largest],6),atominfo(inchi2wn[largest],7))

            # Transfer three largest coordinates to matrix: smallest, second largest, largest
            m = array.array("i")
            m.append(atominfo(inchi2wn[largest], 5) - atominfo(inchi2wn[list_numbers_stereogenic_centres[s_centre]], 5))
            m.append(atominfo(inchi2wn[largest], 6) - atominfo(inchi2wn[list_numbers_stereogenic_centres[s_centre]], 6))
            m.append(atominfo(inchi2wn[largest], 7) - atominfo(inchi2wn[list_numbers_stereogenic_centres[s_centre]], 7))

            m.append(atominfo(inchi2wn[second_largest], 5) - atominfo(inchi2wn[list_numbers_stereogenic_centres[s_centre]], 5))
            m.append(atominfo(inchi2wn[second_largest], 6) - atominfo(inchi2wn[list_numbers_stereogenic_centres[s_centre]], 6))
            m.append(atominfo(inchi2wn[second_largest], 7) - atominfo(inchi2wn[list_numbers_stereogenic_centres[s_centre]], 7))

            m.append(atominfo(inchi2wn[smallest], 5) - atominfo(inchi2wn[list_numbers_stereogenic_centres[s_centre]], 5))
            m.append(atominfo(inchi2wn[smallest], 6) - atominfo(inchi2wn[list_numbers_stereogenic_centres[s_centre]], 6))
            m.append(atominfo(inchi2wn[smallest], 7) - atominfo(inchi2wn[list_numbers_stereogenic_centres[s_centre]], 7))

            # print("m",m)

            # We now have three points, translated so stereogenic centre is at (0,0,0)
            # These are the points with the larges InChI#
            # Construct a vector from the origin, to the average position
            # Construct a vector from the highest priority to the second highest priority
            # (priority increases with InChI#)
            # Cross these vectors together, than dot with vector to smallest point
            # Sign of the outcome should correlate with InChI +/- stereochemical indicator

            m2 = array.array("i")
            m2.append(m[6])
            m2.append(m[7])
            m2.append(m[8])

            m2.append(m[0] + m[3] + m[6])
            m2.append(m[1] + m[4] + m[7])
            m2.append(m[2] + m[5] + m[8])

            m2.append(m[3] - m[0])
            m2.append(m[4] - m[1])
            m2.append(m[5] - m[2])

            sign_stereo = m2[0] * (m2[4] * m2[8] - m2[5] * m2[7]) - m2[1] * (m2[3] * m2[8] - m2[5] * m2[6]) + m2[2] * (m2[3] * m2[7] - m2[4] * m2[6])

            # print("centre: ",list_numbers_stereogenic_centres[s_centre],sign_stereo,list_stereogenic_centres[s_centre])
            # printatoms()

            inchi_sign = 1
            if list_stereogenic_centres[s_centre].find("+") == -1:
                inchi_sign = -1
            if inchi_sign * sign_stereo > 0:
                # need to invert center
                # print("Inverting centre",list_stereogenic_centres[s_centre],sca)
                # If the stereocenter is a methine, may need to alter just one twist
                for workingposition in range(0, len(sca)):
                    i = (3 - atominfo(inchi2wn[sca[workingposition]], 3)) % 3
                    if i > 0:
                        setatominfo(inchi2wn[sca[workingposition]], 3, i)
                        # print("check setatominfo",i,inchi2wn[sca[workingposition]],(3-atom[inchi2wn[sca[workingposition]]*8+3])%3,atominfo(inchi2wn[sca[workingposition]],3))
                        # atom[inchi2wn[sca[workingposition]]*9+3]=(3-atom[inchi2wn[sca[workingposition]]*9+3])%3
                        twist = atominfo(inchi2wn[sca[workingposition]], 3)
                        # why use direction, not newdirection?
                        # could this be replaced by atominfo(inchi2wn[atominfo(inchi2wn[atominfo(inchi2wn[sca[workingposition]],1)],1)],4) ??
                        # direction=atominfo(inchi2wn[sca[workingposition]],2)
                        # if direction == atominfo(inchi2wn[atominfo(inchi2wn[atominfo(inchi2wn[sca[workingposition]],1)],1)],4):
                        #  print("OK!", workingposition, direction)
                        # else:
                        #  print("Oh dear!", workingposition, direction, atominfo(inchi2wn[atominfo(inchi2wn[atominfo(inchi2wn[sca[workingposition]],1)],1)],4))
                        # direction (penultimate direction)
                        # old_direction (direction of attached atom)
                        direction = atominfo(inchi2wn[atominfo(inchi2wn[atominfo(inchi2wn[sca[workingposition]], 1)], 1)], 4)
                        old_direction = atominfo(inchi2wn[atominfo(inchi2wn[sca[workingposition]], 1)], 4)
                        newdirection = twisting[twist * 64 + old_direction * 8 + direction]
                        # print("about to change newdirection",sca[workingposition],inchi2wn[sca[workingposition]],"o,d,n",old_direction,direction,newdirection)
                        setatominfo(inchi2wn[sca[workingposition]], 4, newdirection)
            # printatoms()
            # Coordinates are now inconsistent with twists; and with conformation[]
            # This should set it right - but perhaps not needed at this stage
            for workingposition in range(current_position, number_of_carbons):
                attach_atom = atominfo(workingposition, 1)
                newdirection = atominfo(workingposition, 4)
                # print("prep ",current_position,workingposition,attach_atom,newdirection)
                xcoord = atominfo(inchi2wn[attach_atom], 5) + directions[newdirection * 3]
                ycoord = atominfo(inchi2wn[attach_atom], 6) + directions[newdirection * 3 + 1]
                zcoord = atominfo(inchi2wn[attach_atom], 7) + directions[newdirection * 3 + 2]
                setatominfo(workingposition, 5, xcoord)
                setatominfo(workingposition, 6, ycoord)
                setatominfo(workingposition, 7, zcoord)
                # print("reset",workingposition,";",xcoord,ycoord,zcoord,";",attach_atom,newdirection)


            # If the stereocenter is in a ring, and this center is wrongly built,
            # then need to justify whether should be re-built
            # to correct the first_ring_atom based stereocenter
        else:
            if list_numbers_stereogenic_centres[s_centre] in first_atom_InChI:
                print("need to check whether to invert first ring-atom based centre")
                print("stereocentre: ", list_numbers_stereogenic_centres[s_centre], s_centre)
                # print("inchi_ring_order: ", inchi_ring_order)
                current_workingatom = list_numbers_stereogenic_centres[s_centre]
                print("current_workingatom:", current_workingatom)
                # Find the ring index where the atom is located
                ring_index = next((i for i, ring in enumerate(inchi_ring_order) if list_numbers_stereogenic_centres[s_centre] in ring), None)
                if ring_index is None:
                    raise ValueError("Atom not found in any ring")
                print(f"ring_index: {ring_index}")

                # Get the surrounding heavy atoms around the ring atom to form a stereocenter
                center_sca = inchi_ring_order[ring_index][0]
                temp_largest = inchi_ring_order[ring_index][5]
                temp_second_largest = atominfo(inchi2wn[list_numbers_stereogenic_centres[s_centre]], 1)
                temp_smallest = inchi_ring_order[ring_index][1]

                stereo_list = [temp_largest, temp_second_largest, temp_smallest]
                stereo_list.sort()
                print(f"stereo_list: {stereo_list}")

                # Get the largest, second largest, and smallest values from the stereo_list
                largest, second_largest, smallest = stereo_list[2], stereo_list[1], stereo_list[0]

                # Transfer three largest coordinates to matrix: smallest, second largest, largest
                m = array.array("i")
                m.append(atominfo(inchi2wn[largest], 5) - atominfo(inchi2wn[list_numbers_stereogenic_centres[s_centre]], 5))
                m.append(atominfo(inchi2wn[largest], 6) - atominfo(inchi2wn[list_numbers_stereogenic_centres[s_centre]], 6))
                m.append(atominfo(inchi2wn[largest], 7) - atominfo(inchi2wn[list_numbers_stereogenic_centres[s_centre]], 7))

                m.append(atominfo(inchi2wn[second_largest], 5) - atominfo(inchi2wn[list_numbers_stereogenic_centres[s_centre]], 5))
                m.append(atominfo(inchi2wn[second_largest], 6) - atominfo(inchi2wn[list_numbers_stereogenic_centres[s_centre]], 6))
                m.append(atominfo(inchi2wn[second_largest], 7) - atominfo(inchi2wn[list_numbers_stereogenic_centres[s_centre]], 7))

                m.append(atominfo(inchi2wn[smallest], 5) - atominfo(inchi2wn[list_numbers_stereogenic_centres[s_centre]], 5))
                m.append(atominfo(inchi2wn[smallest], 6) - atominfo(inchi2wn[list_numbers_stereogenic_centres[s_centre]], 6))
                m.append(atominfo(inchi2wn[smallest], 7) - atominfo(inchi2wn[list_numbers_stereogenic_centres[s_centre]], 7))

                # print("m",m)

                # We now have three points, translated so stereogenic centre is at (0,0,0)
                # These are the points with the larges InChI#
                # Construct a vector from the origin, to the average position
                # Construct a vector from the highest priority to the second highest priority
                # (priority increases with InChI#)
                # Cross these vectors together, than dot with vector to smallest point
                # Sign of the outcome should correlate with InChI +/- stereochemical indicator

                m2 = array.array("i")
                m2.append(m[6])
                m2.append(m[7])
                m2.append(m[8])

                m2.append(m[0] + m[3] + m[6])
                m2.append(m[1] + m[4] + m[7])
                m2.append(m[2] + m[5] + m[8])

                m2.append(m[3] - m[0])
                m2.append(m[4] - m[1])
                m2.append(m[5] - m[2])

                sign_stereo = m2[0] * (m2[4] * m2[8] - m2[5] * m2[7]) - m2[1] * (m2[3] * m2[8] - m2[5] * m2[6]) + m2[2] * (m2[3] * m2[7] - m2[4] * m2[6])

                # print("centre: ",list_numbers_stereogenic_centres[s_centre],sign_stereo,list_stereogenic_centres[s_centre])
                # printatoms()
                # if sign_stereo > 0, so by similar to CIP rule it is clockwise
                print("sign_stereo:", sign_stereo)

                inchi_sign = 1
                # if InChI did not find + and it is m0, then by InChI rule it is counterwise
                # need to change inchi_sign to -1
                if list_stereogenic_centres[s_centre].find("+") == -1 and sys.argv[1].find("/m0") > -1:
                    inchi_sign = -1
                # test example InChI=1S/C8H16O/c1-3-8-7(2)5-4-6-9-8/h7-8H,3-6H2,1-2H3/t7-,8+/m1/s1
                if list_stereogenic_centres[s_centre].find("+") > -1 and sys.argv[1].find("/m1") > -1:
                    inchi_sign = -1
                print("inchi_sign:", inchi_sign)

                if inchi_sign * sign_stereo > 0:
                    # InChI rule is against the CIP rule
                    # so they cannot be the same sign
                    print("NEED TO REBUILD RING")
                    current_ring_index = next((i for i, ring in enumerate(inchi_ring_order) if
                                               current_workingatom in ring), None)
                    if ring_index is None:
                        raise ValueError("Atom not found in any ring")
                    print(f"current_ring_index: {current_ring_index}")
                    # need to invert center of the wrongly built ring, the first ring atom
                    ring_twist_chose = 1  # change to statue to rebuild ring
                    print("ring_twist_statue_before:", ring_twist_statue)
                    ring_twist_statue[current_ring_index] = ring_twist_chose
                    print("ring_twist_statue_final:", ring_twist_statue)
                    print("ring_twist_list[1]", ring_twist_list[1])

                    print("current_workingatom:", current_workingatom)
                    starting=inchi2wn[current_workingatom]
                    print("starting:", starting)
                    for workingposition in range(starting, number_of_skeleton):
                        attach_atom = atominfo(workingposition, 1)
                        # direction (penultimate direction)
                        # old_direction (direction of attached atom)
                        old_direction = atominfo(inchi2wn[attach_atom], 4)
                        penultimate_atom = atominfo(inchi2wn[attach_atom], 1)
                        direction = atominfo(inchi2wn[penultimate_atom], 4)

                        ## if the atom is non ring atoms, also needs to re-build the structure
                        if atominfo(workingposition, 0) not in ringatoms:
                            ## read the previous twist value
                            ## if not correct then need to re-set the new value
                            twist = atominfo(workingposition, 3)
                            # direction=atominfo(inchi2wn[atominfo(workingposition, 1)], 4)
                            # old_direction=atominfo(inchi2wn[atominfo(workingposition, 1)], 4)
                            newdirection = twisting[twist * 64 + old_direction * 8 + direction]
                            # print("about to change newdirection",sca[workingposition],inchi2wn[sca[workingposition]],"o,d,n",old_direction,direction,newdirection)
                            setatominfo(workingposition, 4, newdirection)

                        ## Need to re-twist to close the ring
                        ## to make the first ring atom in the right stereo-center
                        else:
                            # Find the current ring index where the atom is located
                            current_ring_index = next((i for i, ring in enumerate(inchi_ring_order) if
                                               atominfo(workingposition, 0) in ring), None)
                            if ring_index is None:
                                raise ValueError("Atom not found in any ring")
                            print(f"current_ring_index: {current_ring_index}")
                            current_atom_index=inchi_ring_order[current_ring_index].index(atominfo(workingposition, 0))
                            print('current_atom_index', current_atom_index)
                            previous_twist = ring_twist[current_ring_index][2]
                            # The first three atom will maintain original twist value
                            # just need to read and build the ring structure
                            # if it is the first ring, then it involves first three atoms, so do not need to re-calculate
                            # just use original structure data
                            if current_atom_index < 3:
                                if current_ring_index==0:
                                    newdirection=atominfo(workingposition, 4)
                                else:
                                    twist = atominfo(workingposition, 3)
                                    newdirection = twisting[twist * 64 + old_direction * 8 + direction]
                                    setatominfo(workingposition, 4, newdirection)
                            if ring_twist_statue[current_ring_index] == 0:
                                if current_atom_index == 3:
                                    current_twist = (previous_twist + ring_twist_list[0][0]) % 3
                                    twist = current_twist
                                    newdirection = twisting[twist * 64 + old_direction * 8 + direction]
                                    setatominfo(inchi2wn[atominfo(workingposition, 0)], 3, twist)
                                    setatominfo(inchi2wn[atominfo(workingposition, 0)], 4, newdirection)
                                if current_atom_index == 4:
                                    current_twist = (previous_twist + ring_twist_list[0][1]) % 3
                                    twist = current_twist
                                    newdirection = twisting[twist * 64 + old_direction * 8 + direction]
                                    setatominfo(inchi2wn[atominfo(workingposition, 0)], 3, twist)
                                    setatominfo(inchi2wn[atominfo(workingposition, 0)], 4, newdirection)
                                if current_atom_index == 5:
                                    current_twist = (previous_twist + ring_twist_list[0][2]) % 3
                                    twist = current_twist
                                    newdirection = twisting[twist * 64 + old_direction * 8 + direction]
                                    setatominfo(inchi2wn[atominfo(workingposition, 0)], 3, twist)
                                    setatominfo(inchi2wn[atominfo(workingposition, 0)], 4, newdirection)
                            if ring_twist_statue[current_ring_index] == 1:
                                if current_atom_index == 3:
                                    current_twist = (previous_twist + ring_twist_list[1][0]) % 3
                                    twist = current_twist
                                    newdirection = twisting[twist * 64 + old_direction * 8 + direction]
                                    setatominfo(inchi2wn[atominfo(workingposition, 0)], 3, twist)
                                    setatominfo(inchi2wn[atominfo(workingposition, 0)], 4, newdirection)
                                if current_atom_index == 4:
                                    current_twist = (previous_twist + ring_twist_list[1][1]) % 3
                                    twist = current_twist
                                    newdirection = twisting[twist * 64 + old_direction * 8 + direction]
                                    setatominfo(inchi2wn[atominfo(workingposition, 0)], 3, twist)
                                    setatominfo(inchi2wn[atominfo(workingposition, 0)], 4, newdirection)
                                if current_atom_index == 5:
                                    current_twist = (previous_twist + ring_twist_list[1][2]) % 3
                                    twist = current_twist
                                    newdirection = twisting[twist * 64 + old_direction * 8 + direction]
                                    setatominfo(inchi2wn[atominfo(workingposition, 0)], 3, twist)
                                    setatominfo(inchi2wn[atominfo(workingposition, 0)], 4, newdirection)


                        print("attach_atom, newdirection", attach_atom, newdirection)
                        xcoord = atominfo(inchi2wn[attach_atom], 5) + directions[newdirection * 3]
                        ycoord = atominfo(inchi2wn[attach_atom], 6) + directions[newdirection * 3 + 1]
                        zcoord = atominfo(inchi2wn[attach_atom], 7) + directions[newdirection * 3 + 2]
                        setatominfo(workingposition, 5, xcoord)
                        setatominfo(workingposition, 6, ycoord)
                        setatominfo(workingposition, 7, zcoord)






            ## or the stereocenter is in other position of the ring
            ## also need to check whether the ring bond is in correct orientation
            # printatoms()
            else:
                # Justify whether the first ring atom is correct in stereocenter##
                print("need to check whether to invert centre")
                print("stereocentre: ", list_numbers_stereogenic_centres[s_centre], s_centre)
                print("inchi_ring_order: ", inchi_ring_order)

                # To get the attached ring atom index for further used
                ring_lookup = {atom: idx for idx, ring in enumerate(inchi_ring_order) for atom in ring}
                ring_index = ring_lookup[list_numbers_stereogenic_centres[s_centre]]
                # print("ring_index:", ring_index)
                atom_index = inchi_ring_order[ring_index].index(list_numbers_stereogenic_centres[s_centre])
                print("ring_twist_statue[ring_index]:", ring_twist_statue[ring_index])

                ## Build a list to store these surrounding atoms
                ## and then we could use this list to justify which of these atoms is the non-ring atoms
                ## In the last step, we can then twist and build a new branch to make this stereocenter works correctly
                surrounding_atoms = []


                for workingposition in range(1, number_of_skeleton):
                    if atominfo(workingposition, 1) == list_numbers_stereogenic_centres[s_centre]:
                        print("wn2inchi[workingposition]:", wn2inchi[workingposition])
                        surrounding_atoms.append(wn2inchi[workingposition])
                surrounding_atoms.append(atominfo(inchi2wn[list_numbers_stereogenic_centres[s_centre]], 1))
                print("final_atom_InChI:", final_atom_InChI)
                if list_numbers_stereogenic_centres[s_centre] in final_atom_InChI:
                    surrounding_atoms.append(first_atom_InChI[ring_index])
                surrounding_atoms.sort()
                print(f"surrounding_atoms: {surrounding_atoms}")
                printatoms()

                # Get the largest, second largest, and smallest values from the stereo_list
                largest, second_largest, smallest = surrounding_atoms[2], surrounding_atoms[1], surrounding_atoms[0]

                print("surrounding_atoms:", largest, second_largest, smallest)
                print("HERE is FINE")
                # printatoms()

                # Transfer three largest coordinates to matrix: smallest, second largest, largest
                m = array.array("i")
                m.append(atominfo(inchi2wn[largest], 5) - atominfo(inchi2wn[list_numbers_stereogenic_centres[s_centre]], 5))
                m.append(atominfo(inchi2wn[largest], 6) - atominfo(inchi2wn[list_numbers_stereogenic_centres[s_centre]], 6))
                m.append(atominfo(inchi2wn[largest], 7) - atominfo(inchi2wn[list_numbers_stereogenic_centres[s_centre]], 7))

                m.append(atominfo(inchi2wn[second_largest], 5) - atominfo(inchi2wn[list_numbers_stereogenic_centres[s_centre]], 5))
                m.append(atominfo(inchi2wn[second_largest], 6) - atominfo(inchi2wn[list_numbers_stereogenic_centres[s_centre]], 6))
                m.append(atominfo(inchi2wn[second_largest], 7) - atominfo(inchi2wn[list_numbers_stereogenic_centres[s_centre]], 7))

                m.append(atominfo(inchi2wn[smallest], 5) - atominfo(inchi2wn[list_numbers_stereogenic_centres[s_centre]], 5))
                m.append(atominfo(inchi2wn[smallest], 6) - atominfo(inchi2wn[list_numbers_stereogenic_centres[s_centre]], 6))
                m.append(atominfo(inchi2wn[smallest], 7) - atominfo(inchi2wn[list_numbers_stereogenic_centres[s_centre]], 7))

                # print("m",m)

                # We now have three points, translated so stereogenic centre is at (0,0,0)
                # These are the points with the larges InChI#
                # Construct a vector from the origin, to the average position
                # Construct a vector from the highest priority to the second highest priority
                # (priority increases with InChI#)
                # Cross these vectors together, than dot with vector to smallest point
                # Sign of the outcome should correlate with InChI +/- stereochemical indicator

                m2 = array.array("i")
                m2.append(m[6])
                m2.append(m[7])
                m2.append(m[8])

                m2.append(m[0] + m[3] + m[6])
                m2.append(m[1] + m[4] + m[7])
                m2.append(m[2] + m[5] + m[8])

                m2.append(m[3] - m[0])
                m2.append(m[4] - m[1])
                m2.append(m[5] - m[2])

                sign_stereo = m2[0] * (m2[4] * m2[8] - m2[5] * m2[7]) - m2[1] * (m2[3] * m2[8] - m2[5] * m2[6]) + m2[2] * (m2[3] * m2[7] - m2[4] * m2[6])

                # print("centre: ",list_numbers_stereogenic_centres[s_centre],sign_stereo,list_stereogenic_centres[s_centre])
                # printatoms()
                print("sign_stereo", sign_stereo)

                inchi_sign = 1
                if list_stereogenic_centres[s_centre].find("+") == -1 and sys.argv[1].find("/m0") > -1:
                    inchi_sign = -1
                # test for this example InChI=1S/C5H10O3/c1-4-2-8-5(6)3-7-4/h4-6H,2-3H2,1H3/t4-,5+/m1/s1
                if list_stereogenic_centres[s_centre].find("+") > -1 and sys.argv[1].find("/m1") > -1:
                    inchi_sign = -1
                print("inchi_sign:", inchi_sign)

                if inchi_sign * sign_stereo >= 0:
                    print("NEED TO RETWIST BRANCH!!!!!!!!!")
                    ## Here the strategy is to twist the branch atom that attached to this ring atom stereocenter
                    ## Find the non-ring atom first
                    for sur_atom in surrounding_atoms:
                        if sur_atom not in ringatoms:
                            atom_need_change = sur_atom
                    ## get the working number of the atom need to be changed(the branch atom that attached to current ring stereocenter)
                    atom_branch_index = inchi2wn[atom_need_change]
                    print("workingposition", workingposition)
                    #printatoms()

                    for workingposition in range(current_position, number_of_skeleton):
                        attach_atom = atominfo(workingposition, 1)
                        # direction (penultimate direction)
                        # old_direction (direction of attached atom)
                        old_direction = atominfo(inchi2wn[attach_atom], 4)
                        penultimate_atom = atominfo(inchi2wn[attach_atom], 1)
                        direction = atominfo(inchi2wn[penultimate_atom], 4)
                        # If it is exactly the first atom of the branch that needs to be changed
                        if workingposition == atom_branch_index:
                            print("HERE DO THE BRANCH CHNAGE")

                            # Find the current ring index where the attached ring atom is located
                            current_ring_index = next((i for i, ring in enumerate(inchi_ring_order) if
                                               attach_atom in ring), None)
                            if current_ring_index is None:
                                raise ValueError("Atom not found in any ring")
                            print(f"current_ring_index: {current_ring_index}")
                            current_atom_index=inchi_ring_order[current_ring_index].index(attach_atom)
                            print('current_atom_index', current_atom_index)
                            previous_twist = ring_twist[current_ring_index][2]

                            twist_old = atominfo(workingposition, 3)
                            print("twist_old:", twist_old)
                            # Test -- statue ==0
                            if ring_twist_statue[current_ring_index]==0:
                                if current_atom_index==1:
                                    twist_new = (atominfo(workingposition, 3) + 1)%3
                                if current_atom_index==2:
                                    twist_new = (atominfo(workingposition, 3) + 2)%3

                                if current_atom_index==3:
                                    twist_new = (atominfo(workingposition, 3) + 2)%3

                                if current_atom_index==4:
                                    twist_new = (atominfo(workingposition, 3) + 2)%3
                                if current_atom_index==5:
                                    twist_new = (atominfo(workingposition, 3) + 2)%3
                                # test index =2
                                #else:
                                #    twist_new = (atominfo(workingposition, 3) + 2)%3
                            # Test -- statue ==1
                            if ring_twist_statue[current_ring_index]==1:
                                if current_atom_index==1:
                                    twist_new = (atominfo(workingposition, 3) + 1)%3
                                # Test
                                if current_atom_index==2:
                                    twist_new = (atominfo(workingposition, 3) + 1)%3

                                # Test
                                if current_atom_index==3:
                                    twist_new = (atominfo(workingposition, 3) + 1)%3

                                if current_atom_index==4:
                                    twist_new = (atominfo(workingposition, 3) + 1)%3
                                if current_atom_index==5:
                                    twist_new = (atominfo(workingposition, 3) + 1)%3
                                # test index =2
                                #else:
                                    #twist_new = (atominfo(workingposition, 3) + 1)%3
                            print("ring_twist_statue[current_ring_index]:", ring_twist_statue[current_ring_index])
                            print("twist_new:", twist_new)



                            # twist_new=atominfo(atominfo(inchi2wn[atom_need_change],1), 3)+2
                            newdirection = twisting[twist_new * 64 + old_direction * 8 + direction]
                            # print("about to change newdirection",sca[workingposition],inchi2wn[sca[workingposition]],"o,d,n",old_direction,direction,newdirection)
                            setatominfo(workingposition, 3, twist_new)
                            setatominfo(workingposition, 4, newdirection)
                        else:
                            twist = atominfo(workingposition, 3)
                            newdirection = twisting[twist * 64 + old_direction * 8 + direction]
                            setatominfo(workingposition, 4, newdirection)
                        xcoord = atominfo(inchi2wn[attach_atom], 5) + directions[newdirection * 3]
                        ycoord = atominfo(inchi2wn[attach_atom], 6) + directions[newdirection * 3 + 1]
                        zcoord = atominfo(inchi2wn[attach_atom], 7) + directions[newdirection * 3 + 2]
                        setatominfo(workingposition, 5, xcoord)
                        setatominfo(workingposition, 6, ycoord)
                        setatominfo(workingposition, 7, zcoord)

#conf_num = 10000001
#printatoms()
print("ring_twist_statue:", ring_twist_statue)
#printmol()




#################################################################
# Everything should now be OK for connectivity, etc
# work out energy and geometry for different conformations
# conf_number=0
#################################################################

energy_list=[]
lowest_energy=999999
lowest_energy_conformation=0


if len(sys.argv) > 100:
  conf_num=conf_number
  energy=find_all_conf_number()
  #if energy < 0:
    #print("Energy of structure",sys.argv[2],"is > 200")
  #else:
    #print("Energy of structure",sys.argv[2],"is",energy)
else:
  # need to do conformation search
  conf_num=0
  #energy_list.append([conf_num, 0])
  while conf_num < number_conformations:
    tmp_conf_number=conf_num
    #num_rot=0
    for i in range(0,number_rotatable_bonds):
      #print("conformation",number_rotatable_bonds,i,number_rotatable_bonds-i,conformation)
      conformation[number_rotatable_bonds-i-1]=tmp_conf_number%3
      tmp_conf_number=int(tmp_conf_number/3)
    #print("Conformation",conf_num,"Energy: ",find_all_conf_number(),conformation)
    #print("conf_num:", conf_num)
    energy=find_all_conf_number()
    if energy < -1000:
        #print("current_energy:", energy)
        stop_atom=-1000-energy
        # Should now be able to skip all structures with same motif
        # need to skip by torsion number not atom number...
        # What is relevant torsion angle from atom number?
        # atominfo(stop_atom,8)
        # comment out following lines to stop skipping:
        stop_torsion=atominfo(stop_atom, 8)
        #print("stop_torsion:", stop_torsion)
        skip_number=int(3**(number_rotatable_bonds-atominfo(stop_atom, 8)-1))-1
        #print("skip_number:", skip_number)
        conf_num+=skip_number
        #print("conf_num:", conf_num)
    else:
        #print("current_energy:", energy)
        if energy < lowest_energy:
            lowest_energy=energy
            lowest_energy_conformation=conf_num
        if energy < 200:
            energy_list.append([conf_num, energy])
            printmol()
            # print("Conformation",conf_num,"Energy: ",energy,conformation)
    old_conformation = conformation[:]
    conf_num+=1

  #printatoms()
  print()
  print("Dihedral angles, using InChI numbering")
  current_dihedral=0
  for i in range(number_of_skeleton):
    #print("i",i,atom[i*9])
    a1=atom[i*9]
    a2=0
    a3=0
    a4=0
    if atom[inchi2wn[a1]*9+1]>0:
      a2=atom[inchi2wn[a1]*9+1]
      if atom[inchi2wn[a2]*9+1]>0:
        a3=atom[inchi2wn[a2]*9+1]
        if atom[inchi2wn[a3]*9+1]>0:
          a4=atom[inchi2wn[a3]*9+1]
          if atom[inchi2wn[a1]*9+8]==current_dihedral:
            current_dihedral +=1
            print("Dihedral",current_dihedral,"atoms:",a1,a2,a3,a4)
    #print(i,"Dihedral tmp",dihedral_number,"atoms:",a1,a2,a3,a4)
  print()
  print("Number of conformations with accessible energies (< 200):",len(energy_list))
  for i in range(0,len(energy_list)):
    tmp_conf_number=energy_list[i][0]
    #print("tmp_conf_number",energy_list[i],energy_list[i][0],tmp_conf_number)
    float_conformation=[]
    for j in range(0,number_rotatable_bonds):
      #print("conformation",number_rotatable_bonds,i,number_rotatable_bonds-i,conformation)
      #print("tmp_conf_number latest",tmp_conf_number,tmp_conf_number%3,120.0*(tmp_conf_number%3))
      float_conformation.append(180.0-120.0*(tmp_conf_number%3))
      #print("float_conformation",float_conformation)
      tmp_conf_number=int(tmp_conf_number/3)
    print("Conformation, Energy",energy_list[i],float_conformation)
  conf_num=lowest_energy_conformation
  tmp_conf_number=conf_num
  for i in range(0,number_rotatable_bonds):
    conformation[number_rotatable_bonds-i-1]=tmp_conf_number%3
    tmp_conf_number=int(tmp_conf_number/3)
  print("conformation:",conformation)
  energy=find_all_conf_number()
  print()
  print("lowest energy",lowest_energy," at conformation",lowest_energy_conformation)


print()
#print("atominfo:", atominfo(0,0))
#print(len(energy_list))
print("###################################################################################")
print("Lowest_energy_conformer:")
printatoms()
printxyzatoms()
#printmol()
#print()
#printmolatoms()
print("Number of conformations with accessible energies (< 200):",len(energy_list))
print()
end_time=process_time()
running_time=end_time-start_time
print("Diamond_running_time:", str(end_time-start_time))
#table.cell(int(sys.argv[2])+1, 7).value=running_time
#data.save('diamond_runtime.xlsx')
print("energy_list:", energy_list)
#energy_list.sort(key=takeSecond)
#if len(energy_list) >= 25:
  #final_list=energy_list[:25]
#else:
  #final_list=energy_list


#with open(sys.argv[2]+"_"+sys.argv[3]+"_energy_list.txt", "a+") as f_compare:
#  for confor_num in energy_list:
#    f_compare.write(str(confor_num[0]))
#    f_compare.write("\n")
#f_compare.close()
