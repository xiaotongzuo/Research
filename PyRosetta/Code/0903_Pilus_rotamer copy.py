#!/usr/bin/env python
from rosetta import *
init(extra_options='-include_sugars -override_rsd_type_limit -read_pdb_link_records -write_pdb_link_records')

initial_pose=pose_from_pdb('/Users/XT/Dropbox/Project/pilus_repacking.pdb')
#pose=pose_from_pdb('/Users/XT/Dropbox/Project/pilus_repacking.pdb')
pose=Pose()

# Pymol Mover
pmm=PyMOLMover()

#loop_domain=range(132-16,139+1-16)
scorefxn=get_score_function()
# MoveMap
mm=MoveMap()
mm.set_bb_true_range(132-16,139-16)   # inclusive
mm.set_branches(True)                 # sugar
mm.set_nu(True)                       # for Ring Mover

# Mover
smallmover=SmallMover(mm,1,1)
shearmover=ShearMover(mm,1,1)
minmover=MinMover()
minmover.movemap(mm)   
minmover.score_function(scorefxn)
minmover.min_type('dfpmin')
ringmover=RingConformationMover()
ringmover.movemap(mm)

# Job Distributor
jd=PyJobDistributor("0903decoys",100,scorefxn)
# flappy tail algorithm
f=open('log.txt','a')
f.write('count \t num_accepts \t acceptance_rate \n')
f.close()
## 1. max angle for shear mover: 180 degree
count=0
while not jd.job_complete:
    pose=Pose()
    pose.assign(initial_pose)
    #Monte Carlo
    mc=MonteCarlo(pose,scorefxn,1)         # change kT, half is ideal
    shearmover.angle_max('L',180)
    # task pack
    task_pack=standard_packer_task(pose)
    # 120 Cys skipped
    count+=1
    num_accepts=0
    num_rejects=0
    for i in range(250):                   # min mover,50, 5%      ?set tolerance & type df    250
        for ii in range(19):               # shear+ring,19, 95%
            shearmover.apply(pose)
            #mc.boltzmann(pose)
            if mc.boltzmann(pose):
                num_accepts+=1
            ringmover.apply(pose)
            #mc.boltzmann(pose)
            if mc.boltzmann(pose):
                num_accepts+=1
            #pmm.apply(pose)
        minmover.apply(pose)       
        #pmm.apply(pose)
        #mc.boltzmann(pose)
        if mc.boltzmann(pose):
            num_accepts+=1

## 2. max angle for shear mover: 4 degree
    shearmover.angle_max('L',4)
    for i in range(11):                    # repacking    3%    90/8
        for ii in range(121,128+1):        # repacking    3%    90
            for iii in range(2):           # 2,2 min mover    3%    90
                for j in range(2):
                    for k in range (121,128+1):  # 121,128+1 shear+ring+rotamer   95% 16/8
                        shearmover.apply(pose)
                        ringmover.apply(pose)
                        #mc.boltzmann(pose)
                        if mc.boltzmann(pose):
                            num_accepts+=1
                        #print k
                        task_pack.restrict_to_repacking()
                        task_pack.temporarily_fix_everything()
                        for kk in range(116,119+1):
                            task_pack.temporarily_set_pack_residue(kk,True)
                        task_pack.temporarily_set_pack_residue(k,True)
                        # Rotamer Trials
                        rotamer_trials=RotamerTrialsMover(scorefxn,task_pack)
                        # Rotamer Trials
                        rotamer_trials.apply(pose)     # should add more loops, *8
                        #mc.boltzmann(pose)
                        if mc.boltzmann(pose):
                            num_accepts+=1
                        #pmm.apply(pose)
                minmover.apply(pose)
                #pmm.apply(pose)
                #mc.boltzmann(pose)
                if mc.boltzmann(pose):
                    num_accepts+=1
            #repack
            task_pack.restrict_to_repacking()
            task_pack.temporarily_fix_everything()
            for kk in range(116,119+1):
                task_pack.temporarily_set_pack_residue(kk,True)
            task_pack.temporarily_set_pack_residue(ii,True)
            pack_mover=PackRotamersMover(scorefxn,task_pack)
            pack_mover.apply(pose)                     # should add more loops, *8
            minmover.apply(pose)
            #mc.boltzmann(pose)
            if mc.boltzmann(pose):
                num_accepts+=1
            #pmm.apply(pose)

    jd.output_decoy(pose)


    # Acceptance Rate
    print 'Total num of Monte Carlo trials is ', mc.total_trials()
    print 'Num_Accepts=', num_accepts
    acceptance_rate=float(num_accepts)/float(mc.total_trials())
    float('%.4f'%acceptance_rate)
    print 'Acceptance Rate for count',count,' is ', acceptance_rate
    f=open('log.txt','a')
    f.write(str(count)+'\t'+str( num_accepts)+' \t\t' +str(acceptance_rate)+' \n')
    f.close()

