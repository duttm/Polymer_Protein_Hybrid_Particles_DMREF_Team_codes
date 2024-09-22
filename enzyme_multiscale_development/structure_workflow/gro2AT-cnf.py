#!/usr/bin/env python3
# copyright 2023 Mason Hooten mason.simulation@gmail.com MIT license

# Read a gro trajectory file and convert it to a bunch of xyz files suitable
# for use with the A+T pair distribution function routine.

import numpy as np
import sys

def read_frame(f):
    """Read a single frame from a GROMACS trajectory file, returning coordinate data."""
    
    # Input data is in a fixed width format.
    # ; resnumresname atomname atomnum x y z vx vy vz
    #     2CAB     CA    1   2.157   2.021   1.901 [no velocities here] 
    # Note in practice resnum and resname are concatenated as one string.
    #
    # Output will be a list of lists
    # [ ['2CAB', 'HA', 2, 2.0, 3.0, 2.0], ... ]
    
    title_line = f.readline()
    if not title_line:
        return False,False,False,False,False,False,False
    title_data = title_line.split()
    
    #TODO: return time and frame in trajectories
    #if traj:
    #    time  = float(title_data[-3])
    #    frame = int(title_data[-1])
    #    print("frame: ",frame)
    
    atom_count_line = f.readline()
    atom_count_data = atom_count_line.split()
    atom_count      = int(atom_count_data[0])
    
    atom_data = []
    
    for i in range(atom_count):
        data_line   = f.readline()
        f1 = int(data_line[0:5].strip())        # resnum
        f2 = data_line[5:10].strip()            # restype
        f3 = data_line[10:15].strip()           # atomname
        f4 = int(data_line[15:20].strip())      # atomnum
        f5 = float(data_line[20:28].strip())    # x
        f6 = float(data_line[28:36].strip())    # y
        f7 = float(data_line[36:44].strip())    # z
        one_atom = [f1,f2,f3,f4,f5,f6,f7]
        atom_data.append(one_atom)
        
    box_line    = f.readline()
    box_vector  = [float(i) for i in box_line.split()]
    
    max_resnum=-1
    atom_type_list = []
    for atom in atom_data:
        # The modulo terms are to fix the box.
        # Atoms which traveled outside the box and have not yet been
        # reassigned to their SIMULATION NEIGHBORHOOD (by Gromacs) should
        # be counted in the right ANALYSIS NEIGHBORHOOD (by this program).
        # TODO decide if this is the right place for this step.
        atom[3] = float(atom[3]) % box_vector[0]
        atom[4] = float(atom[4]) % box_vector[1]
        atom[5] = float(atom[5]) % box_vector[2]
        
        if atom[0] > max_resnum:
            max_resnum = atom[0]

        if atom[2] not in atom_type_list:
            atom_type_list.append(atom[2])    

    res_sequence = ['' for i in range(max_resnum)]
    for atom in atom_data:
        if res_sequence[atom[0]-1] == '':
            res_sequence[atom[0]-1] = atom[1]

    #atom_types = [atom[2] for atom in atom_data]
    #atom_types = [*set(atom_types)]
    #atom_types.sort()
    
    #atom_type_count = len([*set(atom_types)])
    #for attype in atom_types:
        
    #atom_type_list.sort()

    headers = ['resnum','restype', 'atomname', 'atomnum', 'x', 'y', 'z']
    
    return(atom_count, atom_data, box_vector, res_sequence, max_resnum, atom_type_list, headers)

def main(structure_file) -> None:
    """Read a trajectory."""
    with open(structure_file) as f:
        frame_count = -1
        while True:
            (atom_count,atom_data,box_vector,res_seq,max_resnum,attypes,head) = read_frame(f)
            if not atom_count: break
            frame_count += 1
            with open("cnf."+str(frame_count).zfill(5),'w') as g:
                g.write(str(atom_count)+'\n')
#                g.write(str(box_vector[0]*10)+'\n')
                g.write(str(box_vector[0])+'\n')
                for atom in atom_data:
                    g.write(' '.join([str(i) for i in atom[4:7]])+'\n')
#                    g.write(' '.join([str(i*10) for i in atom[4:7]])+'\n')

TRAJ_FILE = sys.argv[1]

if __name__ == '__main__':
    main(TRAJ_FILE)
