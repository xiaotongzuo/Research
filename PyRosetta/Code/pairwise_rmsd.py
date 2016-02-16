"""
This program writes a file containing a vector with all pairwise backbone 
RMSDs of all or the first n decoys in a directory.

Input: Accepts the directory contaning the decoys as the first argument
and the number of decoys to consider as the second (optional) argument.

For example: python pairwise_rmsd.py out_pdb 100 , or
			 python pairwise_rmsd.py out_pdb

Output: A file names "all_vs_all_rmsd.dat" or "all_vs_all_rmsd_<n>.dat" containing a
vector with all pairwise RMSDs in the format (2,1)(3,1)...(n,1)(3,2)...(n,n-1)
"""

import os
import sys
import time
#import psutil # for memory profiling
import numpy as np
from pyRMSD.RMSDCalculator import RMSDCalculator
#from pyRMSD.condensedMatrix import CondensedMatrix
from prody import *

def parse_decoys(decoy_dir, number_of_decoys = None):  # do not need a decoy_dir
	"""
	Extracts backbone coordinates from first n pdbs in the given directory and 
	returns their coordinates as a numpy array.	If not n is specified, 
	it considers every pdb file.
	decoy_dir = directory with decoys
	number_of_decoys = n (optional)
	"""
    decoy_list = []   # delete all lines which read all the pdb files
	for file in os.listdir(decoy_dir):
		if file.endswith(".pdb"):
			decoy_list.append(os.path.join(decoy_dir,file))
	decoy_list = sorted(decoy_list)

	#selecting first n decoys, if n is specified
	if number_of_decoys != None:
		decoys_considered = decoy_list[:number_of_decoys]
	else:	
		decoys_considered = decoy_list

	#using ProDy to parse PDBs
	pdb = parsePDB(decoys_considered[0])

	#extracting C-alpha coordinates as numpy array
    coords = np.array((pdb.select("name CA C1").getCoordsets()))  # CA C N
	for i in range(1,len(decoys_considered)):
		pdb = parsePDB(decoys_considered[i])
		coords_temp = np.array((pdb.select("name CA C1").getCoordsets()))   # CA C N
		coords = np.vstack((coords, coords_temp))
	number_of_conformations = coords.shape[0]
	number_of_atoms = coords.shape[1]
	print str(number_of_conformations), " decoys loaded, each containing ", str(number_of_atoms), " atoms."

	return coords

def pairwise_rmsd(coords):
	"""
	Returns vector with all pairwise backbone RMSDs from given coordinates
	coords = numpy array of coordinates generated using ProDy
	"""
	calculator = RMSDCalculator("QCP_SERIAL_CALCULATOR", fittingCoordsets  = coords)
	t1 = time.time()
	rmsds = calculator.pairwiseRMSDMatrix()
	t2 = time.time()
	print "With QCP Serial calculator it took: ", t2-t1 ,"s."
	sys.stdout.flush()
	#rmsd_matrix = CondensedMatrix(rmsds)

	return rmsds

def main(args):
	#checking input
	if len(args) == 3: #check if number of decoys is given
		number_of_decoys = int(args[2])
	elif len(args) == 2:
		number_of_decoys = None #just initializing the variable for later
	else:
		print "Error: Incorrect number of arguments given."
		sys.exit(1)

	if os.path.isdir(args[1]):
		decoy_dir = args[1]
	else:
		print "Error: Directory does not exist."
		sys.exit(1)

	#parsing PDBs
	coords = parse_decoys(decoy_dir, number_of_decoys)
	#calculating the backbone rmsd of one decoy vs another
	print "Calculating the backbone rmsd of one decoy vs another."
	sys.stdout.flush()
	rmsd_matrix = pairwise_rmsd(coords)
	print "No. of entries in pairwise RMSD matrix: ",len(rmsd_matrix)

	#saving to file
	t1 = time.time()
	if number_of_decoys != None:
		rmsd_out=open("all_vs_all_rmsd_"+str(number_of_decoys)+".dat",'w')
	else:
		rmsd_out=open("all_vs_all_rmsd.dat",'w')
	for i in range(len(rmsd_matrix)):
		rmsd_out.write(str(rmsd_matrix[i])+"\t")
	rmsd_out.close()
	t2 = time.time()
	print "Vector with all pairwise backbone RMSDs written to file in ", t2-t1 ,"s."

if __name__ == "__main__":
	main(sys.argv)