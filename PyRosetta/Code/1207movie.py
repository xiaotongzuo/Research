#!/usr/bin/env python

# 3000+1000  --ok
# repacking & rotamer 2 each time
# add small mover
# rmsd
# time too long

# scorefxn
# sugar turn off something
# ringmover off
import random
from rosetta import *
init(extra_options='-include_sugars -override_rsd_type_limit -read_pdb_link_records -write_pdb_link_records')

initial_pose=pose_from_pdb('/Users/XT/Dropbox/Project/pilus_monomer_repacking.pdb')   # change file
#initial_pose=pose_from_pdb('/home/Xiaotong/Github/Project/pilus_monomer_repacking.pdb')
pose=Pose()

tail_start=132-16
tail_mid=136-16   # this aa Cys is always False
tail_end=139-16
sugar_start=140-16
sugar_end=144-16

# Pymol Mover
pmm=PyMOLMover()
pmm.keep_history(True)
scorefxn=get_score_function()
# MoveMap
mm=MoveMap()
mm.set_bb_true_range(tail_start,tail_end)   # inclusive
mm.set_branches(True)                 # sugar
mm.set_nu(False)                   # for Ring Mover   # turn off!!!

# Mover
shearmover=ShearMover(mm,1,5)
smallmover=SmallMover(mm,1,5)
minmover=MinMover()
minmover.movemap(mm)
minmover.score_function(scorefxn)
minmover.min_type('dfpmin')
#ringmover=RingConformationMover()
#ringmover.movemap(mm)



# Job Distributor
jd=PyJobDistributor("1207output",1,scorefxn)
jd.native_pose =initial_pose


num_accepts=0
###########
def pack_mover_min_renew_taskpack(pose):
    for resi in range(tail_start,sugar_end+1-1):# repacking    3%
        if resi != tail_mid and resi+1 !=tail_mid:   # 10 cycles, 116-127, 116,117,118,121,122,123,124,125,126,127, gap 119,120
            task_pack.restrict_to_repacking()
            task_pack.temporarily_fix_everything()
            task_pack.temporarily_set_pack_residue(resi,True)
            task_pack.temporarily_set_pack_residue(resi+1,True)
            pack_mover=PackRotamersMover(scorefxn,task_pack)
            pack_mover.apply(pose)
            pmm.apply(pose)
            minmover.apply(pose)
            pmm.apply(pose)
            if mc.boltzmann(pose):
                global num_accepts
                num_accepts+=1
                pmm.apply(pose)



def rotamer_trials_renew_taskpack(pose):
    for resi in range (tail_start,sugar_end+1-1):  # turn on 121,128+1  sugar_end+1-1
        if resi != tail_mid and resi+1 !=tail_mid:
            smallmover.apply(pose)
            pmm.apply(pose)
            shearmover.apply(pose)
            pmm.apply(pose)
            #ringmover.apply(pose)
            #pmm.apply(pose)
            task_pack.restrict_to_repacking()
            task_pack.temporarily_fix_everything()
            task_pack.temporarily_set_pack_residue(resi,True)
            task_pack.temporarily_set_pack_residue(resi+1,True)
            # Rotamer Trials
            rotamer_trials=RotamerTrialsMover(scorefxn,task_pack)
            rotamer_trials.apply(pose)
            pmm.apply(pose)
            if mc.boltzmann(pose):
                global num_accepts
                num_accepts+=1
                pmm.apply(pose)
###########


count=0
while not jd.job_complete:
    pose=Pose()
    pose.assign(initial_pose)
    #Monte Carlo
    mc=MonteCarlo(pose,scorefxn,1)
    # task pack
    task_pack=standard_packer_task(pose)
    for resi in range(tail_start, tail_mid):
        task_pack.temporarily_set_pack_residue(resi,True)
    for resi in range(tail_mid+1, sugar_end):
        task_pack.temporarily_set_pack_residue(resi,True)
    # 120 Cys skipped
    count+=1       # cycles
    num_accepts=0
    # floppy tail algorithm
    ## 1. max angle for shear mover: 180 degree
    shearmover.angle_max('L',180)    # loop for sugar   ???
    smallmover.angle_max('L', 180)
    for i in range(2):                   # min mover,50, 5%, 150
        for ii in range(19):               # shear+ring,19, 95%
            smallmover.apply(pose)
            pmm.apply(pose)
            shearmover.apply(pose)
            pmm.apply(pose)
            new_omega=random.uniform(0,360)
            pose.set_omega(sugar_start, new_omega)  # set the omega of first sugar resi
            pmm.apply(pose)
            #mc.boltzmann(pose)
            if mc.boltzmann(pose):
                num_accepts+=1       # acceptance rate ?
                pmm.apply(pose)
            #ringmover.apply(pose)
            #pmm.apply(pose)
            #mc.boltzmann(pose)
            #if mc.boltzmann(pose):
            #   num_accepts+=1
        pmm.apply(pose)
        minmover.apply(pose)
        pmm.apply(pose)
        #mc.boltzmann(pose)
        if mc.boltzmann(pose):
            num_accepts+=1

## 2. max angle for shear mover: 4 degree
    shearmover.angle_max('L',4)
    smallmover.angle_max('L',4)
    for i in range(1):   # 25
        for i in range(2):  # 2
            for i in range(2):  # 16>>2
                rotamer_trials_renew_taskpack(pose)
            minmover.apply(pose)
            pmm.apply(pose)
            if mc.boltzmann(pose):
                num_accepts+=1
                pmm.apply(pose)
        pack_mover_min_renew_taskpack(pose)
    jd.output_decoy(pose)
    
'''# Acceptance Rate
# log_acceptance_rate
f=open('log.txt','a')
acceptance_rate=float(num_accepts)/float(mc.total_trials())
acceptance_rate=float('%.4f'%acceptance_rate)
f.write('count \t num_accepts \t acceptance_rate \n')
f.write(str(count)+'\t'+str( num_accepts)+'\t\t' +str(acceptance_rate)+'\n')
f.close()'''
