#!/usr/bin/env python
'''This python file is for Homework 2, Problem 4. I didn't use the template. I created a torsion angle file in my code but did not use the torsion angles to calculate the L/H/E propensities.
    --by Xiaotong Zuo, Feb. 2016
    '''

# import
from rosetta import *
init()
from toolbox import get_secstruct
from toolbox import cleanATOM
import sys

# use 1m40.pdb as template, first, cleanATOM
cleanATOM("1m40.pdb")

# load pose
pose=pose_from_pdb("1m40.clean.pdb")


### I did not use torsion angles to calculate the propensities!
# create torsion.dat: phi and psi
f=open("/Users/XT/Downloads/1yy8.pdb","r")
g=open("torsion.dat","w+")
for line in f.readlines():
    a=line.split()
    if "ATOM"==a[0] and "CA"==a[2]:
        g.write(str(a[3])+"\t"+str(a[5])+"\t"+str(pose.phi(int(a[5])))+"\t"+str(pose.psi(int(a[5])))+"\t"+str(pose.psi(int(a[5]))+pose.phi(int(a[5])))+'\n')
        print a
f.close()
g.close()

# read torsion.dat, return residue name and torsion angles
def read_torsions(torsion_fn='torsions.dat'):
    "Return the number of residues and corresponding torsion angles."
    f=open("torsion.dat","r")
    n_residues=[]
    torsion_angles=[]
    for line in f.readlines():
        a=line.split()
        n_residues.append(a[1])
        torsion_angles.append(a[2:4])
    f.close()
    return n_residues, torsion_angles
### I did not use torsion angles to calculate the propensities!


# read secondary structure file
output=sys.stdout
outputfile=open("secstructure.dat","w+")
sys.stdout=outputfile
# write the secstruct output to a file.
get_secstruct(pose, space=3000, page=3000)
outputfile.close()
sys.stdout=output


# secstruct file_delete water
f=open("secstruct.dat","r+")
secstruct=[]
for line in f.readlines():
    secstruct.append(line)
f.close()
g=open("secstruct_delwater.dat","w+")
for i in range(len(secstruct[1])):
    if secstruct[1][i]!="w":
        g.write(str(secstruct[1][i])+"\t"+str(secstruct[2][i])+"\n")
g.close()


# Problem 4: ALA propensity
f=open("secstruct_delwater.dat","r+")
aa_list=[]
N_alpha_tot=0
N_beta_tot=0
N_loop_tot=0
N_alpha_A=0
N_beta_A=0
N_loop_A=0
for line in f.readlines():
    if len(line.split())==2:
        aa_list.append(line.split())
f.close()
for i in range(len(aa_list)):
    if aa_list[i][1]=="H":
        N_alpha_tot+=1
    elif aa_list[i][1]=="E":
        N_beta_tot+=1
    else:
        N_loop_tot+=1
for i in range(len(aa_list)):
    if aa_list[i][0]=="A":
        if aa_list[i][1]=="H":
            N_alpha_A+=1
        elif aa_list[i][1]=="E":
            N_beta_A+=1
        else:
            N_loop_A+=1

P_alpha_ALA=float(N_alpha_A)/float(N_alpha_tot)
P_beta_ALA=float(N_beta_A)/float(N_beta_tot)
P_loop_ALA=float(N_loop_A)/float(N_loop_tot)
print P_alpha_ALA, P_beta_ALA, P_loop_ALA
g=open("ALA_PROPENSITY.dat","w+")
g.write(" P_alpha_ALA,  P_beta_ALA,  P_loop_ALA"+'\n')
g.write(str(P_alpha_ALA)+'\t'+ str(P_beta_ALA)+'\t'+str(P_loop_ALA))
g.close()

# bonus 1: all_type_aa_propensities.dat
# L/E/H count for each residue.
N_alpha=0
N_beta=0
N_loop=0
# all Amino Acid types
aa_type="ACDEFGHIKLMNPQRSTVWY"
# write in a file
g=open("all_type_aa_propensities.dat","w+")
g.write("AA"+"      P_alpha "+"P_beta  "+"P_loop"+"\n\n")
# for loop: all AA types
for aa_index in range(len(aa_type)):
    # for loop: all AA in pdb file
    for i in range(len(aa_list)):
        # the name of AA from pdb == the name of AA from type list
        if aa_list[i][0]==aa_type[aa_index]:
            if aa_list[i][1]=="H":
                N_alpha+=1
            elif aa_list[i][1]=="E":
                N_beta+=1
            else:
                N_loop+=1
    # calculate the propensities
    P_alpha=float(N_alpha)/float(N_alpha_tot)
    P_beta=float(N_beta)/float(N_beta_tot)
    P_loop=float(N_loop)/float(N_loop_tot)
    # write in a file
    g.write(str(aa_type[aa_index])+"\t"+str(round(P_alpha,4))+"\t"+str(round(P_beta,4))+"\t"+str(round(P_loop,4))+"\n")
g.close()



