#!/usr/bin/env python

#for M2 dimer

import random
from rosetta import *
init(extra_options='-include_sugars -override_rsd_type_limit -read_pdb_link_records -write_pdb_link_records')

def init(initial_pose, job_name, job_number):
    initial_pose=pose_from_pdb(initial_pose)   # change file dir
    #initial_pose=pose_from_pdb('/home/Xiaotong/jazz_home/M2dimer_repacking.pdb')

    scorefxn=get_fa_scorefxn()
    scorefxn.set_weight(fa_intra_rep, 0.440)

    # Job Distributor
    jd=PyJobDistributor("0229_low_M2",2000,scorefxn)
    #jd.native_pose =initial_pose  # rmsd


    pose=Pose()
    pose.assign(initial_pose)

    tail_start_2=140+118-50  # from chain B resi 128, cut head 25*2, 128 or 122
    tail_end_2=140+136-50  # to Chain B sugar ends
    sugar_start=140+136-50+1
    sugar_end=140+140-50

    # MoveMap
    mm=MoveMap()
    mm.set_bb_true_range(tail_start_2,sugar_end)   # inclusive
    mm.set_branches(True)                 # sugar torsion angles on
    mm.set_nu(False)                   # for Ring Mover, turn off
    mm.set_chi_true_range(tail_start_2,sugar_end)  ###### add sc chi
    # Mover
    shearmover=ShearMover(mm,1,5) # add to 5 times
    smallmover=SmallMover(mm,1,5)
    minmover=MinMover()
    minmover.movemap(mm)
    minmover.score_function(scorefxn)
    minmover.min_type('dfpmin')

    #Monte Carlo
    mc_outer_core=MonteCarlo(pose,scorefxn,1)
    mc_inner_core=MonteCarlo(pose,scorefxn,1.1)    ### add acceptance rate
    mc_outer_refine=MonteCarlo(pose,scorefxn,1)
    mc_mid_refine=MonteCarlo(pose,scorefxn,1)
    mc_inner_refine=MonteCarlo(pose,scorefxn,1)


###########
def pack_mover_min_renew_taskpack(pose):   # meaningless for sugars, recalculate
    for resi in range(tail_start_2,tail_end_2-1):# repacking    3%    ###### add 1 residue, 3 resi now
        task_pack=standard_packer_task(pose)
        task_pack.restrict_to_repacking()
        task_pack.temporarily_fix_everything()
        task_pack.temporarily_set_pack_residue(resi,True)
        task_pack.temporarily_set_pack_residue(resi+1,True)
        task_pack.temporarily_set_pack_residue(resi+2,True)
        pack_mover=PackRotamersMover(scorefxn,task_pack)
        pack_mover.apply(pose)
        global mc_outer_refine
        mc_outer_refine.boltzmann(pose)

def rotamer_trials_renew_taskpack(pose):
    for resi in range (tail_start_2,sugar_end-1):  # turn on 128~140
        smallmover.apply(pose)
        shearmover.apply(pose)
        #ringmover.apply(pose)
        task_pack=standard_packer_task(pose)
        task_pack.restrict_to_repacking()
        task_pack.temporarily_fix_everything()
        task_pack.temporarily_set_pack_residue(resi,True)
        task_pack.temporarily_set_pack_residue(resi+1,True)
        task_pack.temporarily_set_pack_residue(resi+2,True)
	# Rotamer Trials
        rotamer_trials=RotamerTrialsMover(scorefxn,task_pack)
        rotamer_trials.apply(pose)
        global mc_inner_refine
        mc_inner_refine.boltzmann(pose)
###########

def main():
    count=1
    init('/Users/XT/Dropbox/Project/pilus_monomer_repacking_delH.pdb',"0229_low_M2",2)
    while not jd.job_complete:
        pose=Pose()
        pose.assign(initial_pose)
        # Floppytail Algorithm
        ## 1. max angle for shear mover: 180 degree
        shearmover.angle_max('L',180)    # loop for sugar
        smallmover.angle_max('L', 180)
        mc_outer_core.reset(pose)
        for i in range(150):                   # min mover,50, 5%, 150
            mc_inner_core.reset(pose)
            for ii in range(20):               # shear+ring,95%,19
                smallmover.apply(pose)
                shearmover.apply(pose)
                new_omega=random.uniform(0,360)
                pose.set_omega(sugar_start, new_omega)  # set the omega of first sugar resi
                mc_inner_core.boltzmann(pose)
            minmover.apply(pose)
            mc_outer_core.boltzmann(pose)
            ## mc.show...
        ## 2. max angle for shear mover: 8 degree
        shearmover.angle_max('L',8)
        smallmover.angle_max('L',8)
        mc_outer_refine.reset(pose)
        for i in range(100):   # 10
            mc_mid_refine.reset(pose)
            for i in range(6):  # 6
                mc_inner_refine.reset(pose)
                for i in range(5):  # 2
                    smallmover.apply(pose)
                    shearmover.apply(pose)
                    #new_omega=random.uniform(0,360)
                    #pose.set_omega(sugar_start, new_omega)  # set the omega of first sugar resi
                    #rotamer_trials_renew_taskpack(pose)
                    mc_inner_refine.boltzmann(pose)
                minmover.apply(pose)
                mc_mid_refine.boltzmann(pose)
            #pack_mover_min_renew_taskpack(pose)
            minmover.apply(pose)
            mc_outer_refine.boltzmann(pose)
        '''sugar_rmsd=non_peptide_heavy_atom_RMSD(pose,initial_pose)
        f=open('log.txt','a')
        #acceptance_rate=float(num_accepts)/float(mc.total_trials())
        #acceptance_rate=float('%.4f'%acceptance_rate)
        f.write(str(count)+'\t'+str(jd.current_num)+'\t'+str(jd.current_name)+'\t'+str(sugar_rmsd)+'\n')
        count+=1
        f.close()'''
        jd.output_decoy(pose)
    print "done"




