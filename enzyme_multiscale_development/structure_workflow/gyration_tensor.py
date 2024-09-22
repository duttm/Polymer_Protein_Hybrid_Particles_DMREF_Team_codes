#!/usr/bin/env python3

import numpy as np
import sys
np.set_printoptions(precision=3)
import re
import os

VERBOSE=0
VERYBOSE=0

ROUNDING=1
ROUNDVAL=5

def read_frame(f):
    """Read a single frame from a GROMACS trajectory file, returning coordinate data."""
#    np.set_printoptions(precision=3)
    # Input data is in a fixed width format.
    # ; resnumresname atomname atomnum x y z vx vy vz
    #     2CAB     CA    1   2.157   2.021   1.901 [no velocities here] 
    # Note in practice resnum and resname are concatenated as one string.
    #
    # Output will be a list of lists
    # [ ['2CAB', 'HA', 2, 2.0, 3.0, 2.0], ... ]
    
    title_line = f.readline()
    if not title_line:
        return (False, False,False)
    title_data = title_line.split()
    
    #TODO: return time and frame in trajectories
    #if traj:
    time  = str(title_data[-3])
    frame = str(title_data[-1])
#    print("frame: ",frame)
    timeframe=(time,frame)
    
    atom_count_line = f.readline()
    atom_count_data = atom_count_line.split()
    atom_count      = int(atom_count_data[0])
    
#    atom_data = []
    coords = np.zeros(shape=(atom_count,3))
    atom_names = np.chararray((atom_count,1),3)
    
    for i in range(atom_count):
        data_line   = f.readline()
        f1 = int(data_line[0:5].strip())        # resnum
        f2 = data_line[5:10].strip()            # restype
        f3 = data_line[10:15].strip()           # atomname
        f4 = int(data_line[15:20].strip())      # atomnum
        f5 = float(data_line[20:28].strip())    # x
        f6 = float(data_line[28:36].strip())    # y
        f7 = float(data_line[36:44].strip())    # z
#        one_atom = np.array(f1,f2,f3,f4,f5,f6,f7)
#        atom_data.append(one_atom)
        atom_coords = np.array([f5,f6,f7])
        atom_coords = np.round(atom_coords,3)
        coords[i] = atom_coords
        atom_names[i] = f3
#        atom_names[i] = np.array(f3)
        
    box_line    = f.readline()
    box_vector  = [float(i) for i in box_line.split()]
    
    # Apply periodic boundaries
#    box_np      = np.array(box_vector)   
#    coords      = (coords + box_np) % box_np
#    print('coords',coords)
#    print('atom_names\n',atom_names)
#    print('sim time and step\n',timeframe)

    return (coords, atom_names, timeframe)
    
def get_mass(trajectory_type):
    """Retrieve masses from somewhere. Should be topology eventually.
    For now, atoms_file is a txt file like
        ...
        SQ5n
        SP2a
        SC3
        ...
    And topology_file is a txt file containing just the atom_type
    
    """
    masses=[]
    atoms = np.loadtxt('atoms.txt',dtype='<U32')
    if trajectory_type == 'cg': # martini3 cg
        for atom in atoms:
            if atom[0] == 'S':
                mass = 54
            elif atom[0] == 'T':
                mass = 36
            else: # 'R' regular beads are passively identified in Martini
                mass = 72
            masses.append(mass)
    elif trajectory_type == 'aa': # charmm or amber should work. untested
        for atom in atoms:
            if  str(atom)[2] == 'C':
                mass = 12.011
            elif str(atom)[2] == 'O':
                mass = 15.999
            elif str(atom)[2] == 'N':
                mass = 14.007
            elif str(atom)[2] == 'S':
                mass = 32.06
            elif str(atom)[2] == 'H':
                mass = 1.008
            else: # no default, need to make things look wrong
                mass = 9999
            masses.append(mass)
    elif trajectory_type == 'fffcg1':
        for atom in atoms:
            if   atom == 'AMD':
                mass = 43.028
            elif atom == 'CAB':
                mass = 27.044
            elif atom == 'COO':
                mass = 44.01
            elif atom == 'NH3':
                mass = 17.034
            elif atom == 'PHE':
                mass = 77.1
            else: # no default, need to make things look wrong
                mass = 9999
            masses.append(mass)
    elif trajectory_type == 'fffcg2':
        for atom in atoms:
            if   atom == 'AMD1':
                mass = 43.028
            elif atom == 'CA1':
                mass = 13.018
            elif atom == 'COO':
                mass = 44.01
            elif atom == 'NH3':
                mass = 17.034
            elif atom == 'PHA':
                mass = 39.054
            elif atom == 'PHB':
                mass = 26.036
            elif atom == 'PHC':
                mass = 26.036
            else: # no default, need to make things look wrong
                mass = 9999
            masses.append(mass)
    mass_array = np.array(masses)
    mass_array.resize(len(mass_array),1)
    return mass_array
    
def write_header(parametername:str, input_filename:str, output_filename='output.txt', usertext='#'):
    inputfile_absolute_path = os.path.abspath(input_filename)
    with open(output_filename, 'w') as f:
        f.write('# gyration_tensor.py -- '+parametername+'\n')
#        f.write('# full path '++'\n')
        f.write('# calculated from '+inputfile_absolute_path+'\n')
        f.write('# '+usertext+'\n')

def main(trajectory_file,trajectory_type) -> None:
    print("START")
    """Read a trajectory."""
    np.set_printoptions(precision=3)
    frame_count = 0
    with open(trajectory_file) as f:
        while True:
            (coords,atom_names,(time,frame)) = read_frame(f)
            if not coords.all(): break
#            if not coords: break

            if VERYBOSE==1: print(coords)
            if VERYBOSE==1: print(atom_names)

            if frame_count == 0:
                write_header('Radius of Gyration (Rg)', trajectory_file, 'Rg.txt', 'index time[ps] md_step Rg[nm] Rgx Rgz Rgz')
                write_header('Relative Shape Anisotropy kappa^2', trajectory_file, 'k2.txt', 'index time[ps] md_step k')
                write_header('Asphericity b', trajectory_file, 'b.txt', 'index time[ps] md_step b')
                if trajectory_type == 'aa':
                    np.savetxt('atoms.txt', atom_names, '%s4')
            frame_num = str(frame_count)
            frame_count += 1

            masses = get_mass(trajectory_type)
            if ROUNDING==1: masses = masses.round(ROUNDVAL)
#            print(masses)

            cm = sum(coords*masses) / sum(masses)
            if ROUNDING==1: cm = cm.round(ROUNDVAL)
            dev = coords-cm
            if ROUNDING==1: dev = dev.round(ROUNDVAL)

            if VERYBOSE==1: print("cordmass =",cordmass)
            if VERYBOSE==1: print(sum(cordmass))    
            if VERYBOSE==1: print("cm =",cm)
            if VERYBOSE==1: print("dev =",dev)

            # Gyration Tensor
            S = np.zeros((3,3))
            for i in range(3):
                for j in range(3):
                    S[i,j] = sum(dev[:,i]*dev[:,j])
            S = S/len(coords)
            if ROUNDING==1: S = S.round(ROUNDVAL)

            if VERYBOSE==1: print("S =",S)
            
            w,v = np.linalg.eig(S)
            if VERYBOSE==1: print("w =",w) # eigenvals
            if VERYBOSE==1: print("sqrt(sum(w)) =",np.sqrt(sum(w)))

            with open('princ.txt','a') as princfile:
                princfile.write("unsorted "+str(frame_count-1)+" "+str(w)+" "+str(v)+'\n')
            
            # sort the eigen vals in numerical order, with corresponding vectors
            w, v = zip(*sorted(zip(w, v)))
            lam3,lam2,lam1 = w
            if VERYBOSE==1: print(w)
            if VERYBOSE==1: print(lam1,lam2,lam3)
            
            # Radius of Gyration
            Rg = str(np.sqrt(lam1+lam2+lam3).round(5))
            if VERBOSE==1: print("Rg =",Rg)
            Rgx = str(np.sqrt(lam2 + lam3).round(5))
            Rgy = str(np.sqrt(lam1 + lam3).round(5))
            Rgz = str(np.sqrt(lam1 + lam2).round(5))
            
            # Relative Shape Anisotropy
            K2 = str(1 - 3 * (lam1*lam2 + lam2*lam3 + lam3*lam1) / (lam1+lam2+lam3)**2)
            if VERBOSE==1: print("K2 =",K2)

            # Asphericity    
            b = str(lam1 - 1/2 * (lam2+lam3))
            if VERBOSE==1: print("b =",b)
            
            Rg_line = ' '.join([frame_num, time, frame, Rg, Rgx, Rgy, Rgz])
            K2_line = ' '.join([frame_num, time, frame, K2])
            b_line = ' '.join([frame_num, time, frame, b])
            
            with open('Rg.txt','a') as Rgfile:
                Rgfile.write(Rg_line+'\n')
            with open('k2.txt','a') as K2file:
                K2file.write(K2_line+'\n')
            with open('b.txt','a') as bfile:
                bfile.write(b_line+'\n')
            print("wrote frame "+str(frame_count))
##            with open('princ.txt','a') as bfile:
##                bfile.write(str(frame_count-1)+" "+str(w)+" "+str(v)+'\n')

TRAJ_FILE = sys.argv[1]
TRAJ_TYPE = str(sys.argv[2])
if __name__ == '__main__':
    main(TRAJ_FILE,TRAJ_TYPE)
