#!/usr/bin/env python
''' This file is for CompBio Homework 5.
    Ab initio folding algorithm.
    
    --Xiaotong Zuo, Mar.06'''
# import
import time
start=time.clock()
from rosetta import *
init()
# sequence
seq='NFYGELVDLGVKEKLIEKAGAWYSYKGEKIGQGKANATAWLKDNPETAKEIEKKVRELLL'
# load into pose
initial_pose=pose_from_sequence(seq)
pose=Pose()
pose.assign(initial_pose)
# set up movemap, fragment set
kT=1
mm = MoveMap()
mm.set_bb(True)

fragset9 = ConstantLengthFragSet(9)
fragset9.read_fragment_file("aat000_09_05.200_v1_3")
mover_9mer = ClassicFragmentMover(fragset9, mm)
fragset3 = ConstantLengthFragSet(3)
fragset3.read_fragment_file("aat000_03_05.200_v1_3")
mover_3mer = ClassicFragmentMover(fragset3, mm)

for decoy_num in range(10):
    # step 1, low resolution fragment insertion
    # centroid mode
    switch = SwitchResidueTypeSetMover("centroid")
    switch.apply(pose)
    # define scorefxn and Monte Carlo object.
    centroid_scorefxn=get_cen_scorefxn()
    mc_low=MonteCarlo(pose,centroid_scorefxn,kT)
    # 9mer, 3mer
    for iteration in range(1000):
        mover_9mer.apply(pose)
        mc_low.boltzmann(pose)
    for iteration in range(2000):
        mover_3mer.apply(pose)
        mc_low.boltzmann(pose)
    # dump a low resolution pdb file.
    pose.dump_pdb('Low_Res_decoy.pdb')


    # step 2, high resolution torsion mover and min mover
    # full atom mode
    switch = SwitchResidueTypeSetMover("fa_standard")
    switch.apply(pose)
    # define scorefxn, movers and Monte Carlo object.
    scorefxn=get_fa_scorefxn()
    mc=MonteCarlo(pose,scorefxn,kT)
    smallmover=SmallMover(mm,kT,5)
    shearmover=ShearMover(mm,kT,5)
    minmover=MinMover()
    minmover.movemap(mm)
    minmover.score_function(scorefxn)
    minmover.min_type('dfpmin')
    # 100 cycles of small+shear movements, followed by a minmover.
    for minmove in range(10):
        for torsion_mover in range(100):
            smallmover.apply(pose)
            shearmover.apply(pose)
            mc.boltzmann(pose)
        minmover.apply(pose)
        mc.boltzmann(pose)
    # dump a new pdb file.
    pose.dump_pdb('High_Res_decoy.pdb')
end=time.clock()
print end-start