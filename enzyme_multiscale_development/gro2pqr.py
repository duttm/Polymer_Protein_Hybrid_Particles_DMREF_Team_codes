#!/usr/bin/env python3
# make a pqr file out of a gro file
# plus an atoms file and a charges file. eventually those should come from
# say a topology

import sys

TRAJ_FILE=sys.argv[1]
ATOMS_FILE=sys.argv[2]
CHARGES_FILE=sys.argv[3]
TRAJ_FILE_NOEXT=TRAJ_FILE.split('.')[0]

def main(trajectory_file, modnum=10):
    with open(ATOMS_FILE) as atmhandl:
        atoms = atmhandl.readlines()
    with open(CHARGES_FILE) as chrghandl:
        charges = chrghandl.readlines()
    with open(trajectory_file) as f:
        titleline = f.readline()
        atom_count = int(f.readline().strip())
        print(atom_count)
        for i in range(atom_count):
            data_line   = f.readline()
            
            # read a gro file
            # gro atom C format is "%5d%-5s%5s%5d%8.3f%8.3f%8.3f%8.4f%8.4f%8.4f"
            # resnum
            resnum = int(data_line[0:5].strip())        # resnum
            restype = data_line[5:10].strip()            # restype
            atomname = data_line[10:15].strip()           # atomname
            atomnum = int(data_line[15:20].strip())      # atomnum
            xcoordnm = float(data_line[20:28].strip())    # x
            ycoordnm = float(data_line[28:36].strip())    # y
            zcoordnm = float(data_line[36:44].strip())    # z
            
            charge = float(charges[i])
            
            atomtype = atoms[i][0]
            if atomtype == "S":
                radius = 4.1
            elif atomtype == "T":
                radius = 3.4
            else: radius = 4.7                
            
            # write a pqr file
            # pqr format is whitespace delimited,
            # however one of the programs i'm using (can't remember which)
            # looks for a field starting at some column number, and that can get
            # fouled up by the wrong whitespace.
            # using mostly the pdb format works.
            ATOM = "ATOM"
            printline = f'{ATOM:<6}{atomnum:>5} {atomname:<4} {restype:<3}  {resnum:>4}    {xcoordnm*10:>8.3f}{ycoordnm*10:>8.3f}{zcoordnm*10:>8.3f}{charge:>6.2f}{radius:>6.2f}'
            
            with open(TRAJ_FILE_NOEXT+'.pqr', 'a+') as outhandl:
                outhandl.write(printline+'\n')
#            print(printline)
        with open(TRAJ_FILE_NOEXT+'.pqr','a+') as outhandl:
            outhandl.write('\n'.join(["TER","END"]))
#        print('\n'.join(["TER","END"]))

if __name__ == '__main__':
    main(TRAJ_FILE)
