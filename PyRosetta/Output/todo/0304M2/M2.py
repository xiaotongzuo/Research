#!/usr/bin/env python

#for M2 dimer

import random
from rosetta import *
init(extra_options='-include_sugars -override_rsd_type_limit -read_pdb_link_records -write_pdb_link_records')

##### Things can be changed
#initial_pose=pose_from_file('/Users/XT/Dropbox/M2dimer_repacking.pdb')   # change file dir
initial_pose=pose_from_pdb('/home/Xiaotong/jazz_home/M2dimer_repacking.pdb')

scorefxn=get_fa_scorefxn()
scorefxn.set_weight(fa_intra_rep, 0.440)

# Job Distributor
jd=PyJobDistributor("0304_M2_low",1000,scorefxn)

# define movemap
tail_start=140+118-50  # from chain B resi 128, cut head 25*2, 128 or 122
tail_end=140+136-50  # to Chain B sugar ends
sugar_start=140+136-50+1
sugar_end=140+140-50
#####


#####
pose=Pose()
pose.assign(initial_pose)

# MoveMap
mm=MoveMap()
mm.set_bb_true_range(tail_start,sugar_end)   # inclusive
mm.set_branches(True)                 # sugar torsion angles on
mm.set_chi_true_range(tail_start,sugar_end)  ###### add sc chi
mm.set_nu(False)                   # for Ring Mover, turn off
# Mover
shearmover=ShearMover(mm,1,5) # add to 5 times
smallmover=SmallMover(mm,1,5)
minmover=MinMover()
minmover.movemap(mm)
minmover.score_function(scorefxn)
minmover.min_type('dfpmin')
#Monte Carlo
mc_outer_core=MonteCarlo(pose,scorefxn,1)
mc_mid_core=MonteCarlo(pose,scorefxn,1)
mc_inner_core=MonteCarlo(pose,scorefxn,1)
mc_outer_refine=MonteCarlo(pose,scorefxn,1)
mc_mid_refine=MonteCarlo(pose,scorefxn,1)
mc_inner_refine=MonteCarlo(pose,scorefxn,1)
# task_pack
task_pack=standard_packer_task(pose)
task_pack.restrict_to_repacking()
task_pack.temporarily_fix_everything()
for resi in range (tail_start,sugar_end+1):
    task_pack.temporarily_set_pack_residue(resi,True)
#####

###########
def pack_mover_min_renew_taskpack(pose):   # meaningless for sugars, recalculate
    for resi in range(tail_start,tail_end+1):# repacking    3%    ###### add 1 residue, 3 resi now
        pack_mover=PackRotamersMover(scorefxn,task_pack)
        pack_mover.apply(pose)
        global mc_outer_refine
        mc_outer_refine.boltzmann(pose)

def rotamer_trials_renew_taskpack(pose):
    for resi in range (tail_start,sugar_end+1):  # turn on 128~140
        smallmover.apply(pose)
        shearmover.apply(pose)
        current_angle=pose.omega(sugar_start)
        new_angle=random.gauss(current_angle,25)
        pose.set_omega(sugar_start,new_angle)   # set a new angle
        global mc_inner_refine
        mc_inner_refine.boltzmann(pose)
        # Rotamer Trials
        rotamer_trials=RotamerTrialsMover(scorefxn,task_pack)
        rotamer_trials.apply(pose)
        mc_inner_refine.boltzmann(pose)
###########

count=1
while not jd.job_complete:
    pose=Pose()
    pose.assign(initial_pose)
    # Floppytail Algorithm
    ## 1. max angle for shear mover: 180 degree
    max_angle=180
    mc_outer_core.reset(pose)
    for i in range(5):
        shearmover.angle_max('L',max_angle)    # loop for sugar
        smallmover.angle_max('L',max_angle)
        shearmover.angle_max('H',max_angle)    # loop for sugar
        smallmover.angle_max('H',max_angle)
        shearmover.angle_max('S',max_angle)    # loop for sugar
        smallmover.angle_max('S',max_angle)
        mc_mid_core.reset(pose)
        for i in range(150):                   # min mover,50, 5%, 150
            mc_inner_core.reset(pose)
            for ii in range(19):               # shear+ring,95%,19
                smallmover.apply(pose)
                shearmover.apply(pose)
                mc_inner_core.boltzmann(pose)
                new_omega=random.uniform(0,360)
                pose.set_omega(sugar_start, new_omega)  # set the omega of first sugar resi
                mc_inner_core.boltzmann(pose)
            minmover.apply(pose)
            mc_mid_core.boltzmann(pose)
        mc_outer_core.boltzmann(pose)
        max_angle=max_angle/2
        ## mc.show...
        
    '''## 2. max angle for shear mover: 8 degree
    max_angle=4
    shearmover.angle_max('L',max_angle)    # loop for sugar
    smallmover.angle_max('L',max_angle)
    shearmover.angle_max('H',max_angle)    # loop for sugar
    smallmover.angle_max('H',max_angle)
    shearmover.angle_max('S',max_angle)    # loop for sugar
    smallmover.angle_max('S',max_angle)
    mc_outer_refine.reset(pose)
    for i in range(10):   # 10
        mc_mid_refine.reset(pose)
        for i in range(6):  # 6
            mc_inner_refine.reset(pose)
            for i in range(2):  # 2
                rotamer_trials_renew_taskpack(pose)
            minmover.apply(pose)
            mc_mid_refine.boltzmann(pose)
        pack_mover_min_renew_taskpack(pose)
        minmover.apply(pose)
        mc_outer_refine.boltzmann(pose)'''
    jd.output_decoy(pose)
print "done"




