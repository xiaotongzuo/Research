from rosetta import *
init()
import glob
from toolbox import cleanATOM


# step 1, locate the right directory, and cleanATOM
filenames = glob.glob('*.pdb')
for i, filename in enumerate(filenames):
    cleanATOM(filename)

# step 2, calculate the bond length and write into a .txt file
path = '/Users/XT/Downloads/top8000_chains_70/random_10_pdb/'
filenames = glob.glob(path+'*.clean.pdb')
f=open("NClength.txt","w+")
for i, filename in enumerate(filenames):
    #f.write(str(i+1)+'\t'+str(filename)+'\n')
    pose=pose_from_pdb(filename)
    for resi_num in range(1,pose.total_residue()+1):
        N_xyz = pose.residue(resi_num).xyz("N")
        CA_xyz = pose.residue(resi_num).xyz("CA")
        N_CA_vector = CA_xyz-N_xyz
        f.write(str(N_CA_vector.norm)+"\n")
f.close()