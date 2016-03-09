#!/usr/bin/env python
''' This file is for CompBio Homework 4, Problem 4.
    Please input your:
    1. pose name(like 1B72), or the number of aa and aa name(like 10, A), 
    2. the number of decoys.
    Examples:
    python HW4_4.py 10  A  100
    python HW4_4.py 1B72 100
    
    it will output a text file called 'torsion_angles.txt'
    
    --Xiaotong Zuo, Feb. 24'''
# import
import time
start=time.clock()
from rosetta import *
init(notebook=True)
import random
from math import exp
from toolbox import pose_from_rcsb
import sys



scorefxn=ScoreFunction()   # I copied the weights from full-atom scorefxn
scorefxn.set_weight(fa_atr,0.8)
scorefxn.set_weight(fa_rep,0.44)
scorefxn.set_weight(hbond_sr_bb,1.17)
scorefxn.set_weight(hbond_lr_bb,1.17)
scorefxn.set_weight(hbond_bb_sc,1.17)
scorefxn.set_weight(hbond_sc,1.1)


def torsion_mover(pose):
    '''randomly choose a residue and change its torsion angle randomly.'''
    global score
    score=scorefxn(pose)
    global old_score
    old_score=score         # store score of last accepted pose
    # randomly choose a residue
    global tot_resi
    tot_resi=pose.total_residue()
    resi_num=random.randint(1,tot_resi)  # inclusive
    # randomly choose a torsion angle
    torsion=random.choice(['phi','psi'])
    if torsion=='phi':
        current_angle=pose.phi(resi_num)
        new_angle=random.gauss(current_angle,25)
        pose.set_phi(resi_num,new_angle)   # set a new angle
    else:
        current_angle=pose.psi(resi_num)
        new_angle=random.gauss(current_angle,25)
        pose.set_psi(resi_num,new_angle)
    return pose

def score_pose(pose):
    '''after a move is made, score the pose and do Monte Carlo'''
    # Scoring
    score=scorefxn(pose)
    delta_E=score-old_score
    if delta_E>=0:
        P=exp(-delta_E)*100
        if P<random.randint(1,100):
            pose.assign(save_pose)   # reject pose
            score=scorefxn(pose)
    return pose

def check_input(args):
    '''check the inputs, and create initial pose.'''
    #checking input
    print args
    if len(args)==3:
        number_of_decoys=args[2]
        pose_name=args[1]
        initial_pose=pose_from_rcsb(pose_name)
        g=open(pose_name+'.clean.pdb','r')
        new_pdb=open(pose_name+'.aa_only.pdb','w')
        # delete lines which are not amino acid
        for line in g.readlines():
            if len(line.split()[3])==3:
                new_pdb.write(line)
        g.close()
        new_pdb.close()
        initial_pose=pose_from_file(pose_name+'.aa_only.pdb')
    elif len(args)==4:
        number_of_decoys=args[3]
        aa_num=args[1]
        aa_name=args[2]
        initial_pose=pose_from_sequence(int(aa_num)*str(aa_name))
    return number_of_decoys, initial_pose


def main(args):
    '''main function, generate decoys and output a torsion angle txt file.'''
    #pose_name, number_of_decoys, num_of_iterations=200
    pose=Pose()
    global save_pose
    save_pose=Pose()
    list=check_input(args)
    number_of_decoys=list[0]
    initial_pose=list[1]
    # write torsion angles into a .txt file
    f=open('torsion_angles.txt','w+')
    for arg in args:
        f.write(str(arg)+'\t')
    f.write('\n'+'decoy_num'+'\t'+'phi'+'\t\t'+'psi'+'\n')

    # decoys loop
    for decoy_num in range(1,int(number_of_decoys)+1):
        pose.assign(initial_pose)   # reset pose to the start point everytime
        for iteration in range(100):    # iterations
            save_pose.assign(pose)   # store last accepted pose
            torsion_mover(pose)
            score_pose(pose)
        # pose.dump...
    # phi psi write to a file
        for resi_num in range(1,tot_resi+1):
            phi=pose.phi(resi_num)
            psi=pose.psi(resi_num)
            f.write(str(decoy_num)+'\t'+str(phi)+'\t'+str(psi)+'\n')
    f.close()


if __name__=='__main__':
    main(sys.argv)
end=time.clock()
print end-start