"""
calculate the number of bonds for each node, and the functionality, avoiding
counting linkers with both ends on the same node. User passes in the file
to calculate functionality for, the reference file, the number of nodes, and
the number of linkers.
"""

import os
import numpy as np
import sys

def main():
    output_file = 'functionality.txt'
    l_atoms = ('Rh_x', 'N_x')

    with open(output_file, 'w') as f:
        f.write("functionality   stdev   single   inter   double\n")

    in_files = []
    all_dirs = os.listdir()
    dirs = []
    for d in all_dirs:
        if d[:4] == 'step' and not d[-3:-1] == 'md':
            dirs.append(d)
    dirs = sorted(dirs, key=lambda a: int(a[-3:]))
    for d in dirs:
        in_files.append(d + '/min.lmps')

    ref_file = sys.argv[1]
    nodes = int(sys.argv[2])
    linkers = int(sys.argv[3])

    for func_file in in_files:
        print('Analyzing', func_file)

        vals = calc_func(func_file, ref_file, nodes, linkers, l_atoms)
        func, stdev, single, inter, same_MOP = vals

        with open(output_file, 'a') as f:
            f.write(str(func) + '   ')
            f.write(str(round(stdev, 2)) + '   ')
            f.write(str(single) + '   ')
            f.write(str(inter) + '   ')
            f.write(str(same_MOP) + '\n')

def calc_func(func_file, ref_file, nodes, linkers, l_atoms):
    # Important params
    link_ind = []
    link_set = set()
    linker_array = np.zeros([nodes + linkers + 1, 2])
    new_bond = 9
    #node_array = np.zeros([nodes + linkers + 1, 12])

    with open(ref_file, 'r') as f:
        lines = f.readlines()

    header = lines[:30]
    for line in header:
        # Get total number of atoms
        if line.strip()[-5:] == 'atoms':
            num_atoms = int(line.split()[0])
        # Get indices of unlinked linker atoms
        if l_atoms[0] in line or l_atoms[1] in line:
            link_ind.append(int(line.split()[0]))

    # Keep track of atom number and the associated molecule number
    atoms_array = np.zeros([num_atoms + 1, 1])
    atoms_flag = False
    for line in lines:
        if 'Atoms' in line:
            atoms_flag = True
        if 'Velocity' in line or 'Bonds' in line:
            atoms_flag = False
        if atoms_flag and line != '\n' and 'Atoms' not in line:
            # check if the atom is one of the link atoms
            if int(line.split()[2]) == link_ind[0] or int(line.split()[2]) == link_ind[1]:
                idx = int(line.split()[0])
                link_set.add(idx)
                # Save molecule index of linker atoms
                atoms_array[idx, 0] = line.split()[1]
                if int(line.split()[2]) == link_ind[1]:
                    #print(atoms_array[idx])
                    pass

    with open(func_file, 'r') as f:
        lines = f.readlines()

    bonds_flag = False
    for line in lines:
        if 'Bonds' in line:
            bonds_flag = True
        if 'Angles' in line:
            bonds_flag = False
        if bonds_flag and line != '\n' and 'Bonds' not in line:
            # check if the atom is one of the link atoms
            if int(line.split()[2]) in link_set and int(line.split()[3]) in link_set:
            #if int(line.split()[1] == new_bond):
                # get linker node pair
                pair = (atoms_array[int(line.split()[2]), 0], atoms_array[int(line.split()[3]), 0])
                linker = int(max(pair[0], pair[1]))
                node = int(min(pair[0], pair[1]))
                #print(pair)
                if linker_array[linker, 0] == 0:
                    linker_array[linker, 0] = node
                elif linker_array[linker, 1] == 0:
                    linker_array[linker, 1] = node
                else:
                    print('WARNING: too many bonds for one linker')

    #print(linker_array)

    same_MOP = 0
    inter = 0
    single = 0
    nodes_array = np.zeros([nodes+1,])
    for i in range(nodes+1, nodes+linkers+1):
        #print(linker_array[i])
        if linker_array[i, 0] == 0 and linker_array[i, 1] == 0:
            #print('Either not a linker, or free linker')
            pass
        elif linker_array[i, 0] == 0 or linker_array[i, 1] == 0:
            single += 1
        elif linker_array[i, 0] == linker_array[i, 1]:
            same_MOP += 1
        else:
            inter += 1
            nodes_array[int(linker_array[i, 0])] += 1
            nodes_array[int(linker_array[i, 1])] += 1

    return np.average(nodes_array[1:]), np.std(nodes_array[1:]), single, inter, same_MOP



if __name__ == '__main__':
    main()
