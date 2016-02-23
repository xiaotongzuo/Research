#!/usr/bin/env python

# Import Rosetta, sys, random
import random
import sys
from rosetta import *

# Initialize Rosetta.
init()

# aa_list is the single letter code list of the 20 normal occuring amino acids.
aa_list="ACDEFGHIKLMNPQRSTVWY"

# this function is to create the alpha-helix
def idealize_helix(n_residues):
    aa_seq=[]
    phi=-60    # the phi torsion angle for alpha helix is -60
    psi=-45    # the psi torsion angle for alpha helix is -45
    for i in range(n_residues):
        random_num=random.randrange(0,len(aa_list))  # create a random number from 0 to 20
        aa_seq.append(aa_list[random_num])    # add amino acid to the sequence randomly
    aa_seq_str=''.join(aa_seq)    # convert the sequence list to a string
    pose=pose_from_sequence(aa_seq_str,"fa_standard")    # create the pose
    for resi_num in range(1,n_residues+1):    # set torsion angles to create an alpha helix
        pose.set_phi(resi_num,phi)
        pose.set_psi(resi_num,psi)
    return pose

def main(args):    # check the input
    if len(args)==2:
        n_residues=int(args[1])
        helix=idealize_helix(n_residues)
        helix.dump_pdb("helix.pdb")
    else:
        print "Wrong input"
        sys.exit(1)

# Check if we're importing this file as a module
if __name__=='__main__':
    main(sys.argv)
else:
    # If we are not importing, run this
    helix = idealize_helix(20)
    # Output result
    helix.dump_pdb("helix.pdb")
