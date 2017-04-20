#!/share/apps/python-2.7.2/bin/python
#
# Either include pyrosetta in your path or run something like 
# source /home/ateufel/Pyrosetta/PyRosetta.monolith.ubuntu.release-80/SetPyRosettaEnvironment.sh
# so that pyrosetta is in the path
# run this program like: python score_mutant_AIT.py tiny_gene
# where tiny_gene is the name of the pdb file you want to use without the ".pdb" part
# note that your pdb files must start a residue 1, this may mean you have to renumber your pbd file 

# Imports
import sys
import argparse
import random, math, os 
from rosetta import init, pose_from_pdb, get_fa_scorefxn, \
                    standard_packer_task, change_cys_state, \
                    Pose, MoveMap, RotamerTrialsMover, MinMover

from toolbox import mutate_residue
from toolbox import cleanATOM
from time import time


def main():
	#takes name of pdb file without the extention
	args =  sys.argv	
	pdb_file = args[1]
	#set up timer to figure out how long the code took to run
	t0=time()

	# Initialize Rosetta.
	init(extra_options='-mute basic -mute core')

	# Constants
	PACK_RADIUS = 10.0
	#Amino acids, notice there is no C
	AAs = ("A","D","E","F","G","H","I","K","L","M","N","P","Q","R","S","T","V","W","Y")
	#Number of mutations to accept
	max_accept_mut = 5000
	#Population size
	N = 100
	#Beta (temp term)
	beta = 1

	#Prepare data headers
	data = ['Variant,Rosetta Score,"delta-delta-G",Probability,Generation\n']

	#Load and clean up pdb file  
	name=pdb_file+".pdb"
	cleanATOM(name)
	clean_name=pdb_file+".clean.pdb"
	initial_pose = pose_from_pdb(clean_name)

	#Set up ScoreFunction
	sf = get_fa_scorefxn()

	#Set up MoveMap.
	mm = MoveMap()
	mm.set_bb(True)
	mm.set_chi(True)

	#Pack and minimize initial pose to remove clashes.
	pre_pre_packing_score = sf(initial_pose)

	task = standard_packer_task(initial_pose)
	task.restrict_to_repacking()
	task.or_include_current(True)
	pack_rotamers_mover = RotamerTrialsMover(sf, task)
	#pack_rotamers_mover.apply(initial_pose)

	min_mover = MinMover()
	min_mover.movemap(mm)
	min_mover.score_function(sf)
	min_mover.min_type('dfpmin_armijo_nonmonotone')
	#min_mover.apply(initial_pose)

	post_pre_packing_score = sf(initial_pose)

	#Set threshold for selection 
	#threshold = post_pre_packing_score/2
	threshold = args[2]
	print 'threshold:', threshold

	data.append('WT,' + str(post_pre_packing_score) + ',0.0,0.0,0\n')

	#number of residues to select from
	n_res = initial_pose.total_residue()

	#start evolution
	i=0
	gen=0
	while i < max_accept_mut:

		#update the number of generations that have pased
		gen+=1

		#print 'accepts:', i 

		#pick a place to mutate
		mut_location = random.randint(1, n_res)

		#get the amino acid at that position
		res = initial_pose.residue(mut_location)
		
		#don't mess with C, just choose again
		#while(res == 'C'):
		while(res.name1() == 'C'):
			mut_location = random.randint(1, n_res)
			#get the amino acid at that position
			res = initial_pose.residue(mut_location)

		#choose the amino acid to mutate to
		new_mut_key = random.randint(0,len(AAs)-1)
		proposed_res = AAs[new_mut_key]

		#make the mutation
		mutant_pose = mutate_residue(initial_pose, mut_location, proposed_res, PACK_RADIUS, sf)

		#score mutant
		variant_score = sf(mutant_pose)

		#get the probability that the mutation will be accepted
		probability = calc_prob_mh(variant_score, post_pre_packing_score, N, beta, threshold)

		#test to see if mutation is accepted
		if random.random() < probability:
	
			#create a name for the mutant if its going to be kept 
			variant_name = res.name1() + str(initial_pose.pdb_info().number(mut_location)) + str(proposed_res)

			#save name and energy change
			data.append(variant_name + "," + str(variant_score) + "," + str(variant_score - post_pre_packing_score) + "," + str(probability) + "," + str(gen) + "\n")

			#pdb_name=str(i)+".pdb"
			#mutant_pose.dump_pdb(pdb_name)
			if i == (max_accept_mut - 1):
				pdb_name=str(i)+".pdb"	
				mutant_pose.dump_pdb(pdb_name)

			#update the wildtype 
			initial_pose = mutant_pose
			post_pre_packing_score = variant_score

			#update number of accepts
			i+=1

	print '\nMutations and scoring complete.'
	t1 = time()
	# Output results.
	data_filename = pdb_file + '_variant_scores_mh_5000.csv'
	with open(data_filename, "w") as f:
		f.writelines(data)

	print 'Data written to:', data_filename
	print 'program takes %f' %(t1-t0)


###assorted functions that have to do with scoring and prob of acceptance ####


#score functions for met-hastings selection
def calc_prob_mh(stab_mut, stab_org, N, beta, thresholds):

  xi = calc_x(stab_org, beta, thresholds)
  xj = calc_x(stab_mut, beta, thresholds)

  if xj > xi:
    return((1.0))
  else:
    exponent = -2 * float(N) * (xi - xj)
    return(safe_calc(exponent))



#score functions for met-hastings selection
def calc_prob_fix(stab_mut, stab_org, N, beta, thresholds):

  xi = calc_x_fix(stab_org, beta, thresholds)
  xj = calc_x_fix(stab_mut, beta, thresholds)

  if xi == xj:
    return(1/float(N))

  if(xj==0.0):
    return(0.0)

  try:
    p =((1-pow( (xi/xj),2)) /(1-pow( (xi/xj), (2 * float(N) )) ) )
  except OverflowError as e:
    p = 0.0
  return (p)



def calc_x_fix(data, beta, threshold):
  total = 0
  exponent = float(beta) * (float(data) - float(threshold))
  total = 1/(safe_calc(exponent) + 1)
  return(total)

def calc_x(data, beta, threshold):
  total = 0
  exponent = float(beta) * (float(data) - float(threshold))
  total += -math.log(safe_calc(exponent) + 1)
  return(total)


def safe_calc(exponent):
  if exponent > 700:
    #print("system maxed")
    return(sys.float_info.max)
  else:
    return(math.exp(exponent))


#Run main program
if __name__ == '__main__':
   main()
