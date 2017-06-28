from Bio import SeqIO
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.cm as cm
from pylab import *

def cal_divergence(filePath, naturalPath, withGaps):
	file_sequences = SeqIO.parse(open(filePath),'fasta')
	natural_sequences = SeqIO.parse(open(naturalPath),'fasta')
	
	natural_seq = natural_sequences.next().seq
	
	if not withGaps:
		natural_seq = natural_seq.ungap("-")
	
	col_count = len(natural_seq)
		
	result_list = []
	for file_sequence in file_sequences:
		file_seq = file_sequence.seq
		result = 0
		for i in range(0, col_count):
			if file_seq[i] != natural_seq[i]:
				result += 1
		result_list.append(result/float(col_count))
	
	return result_list

PDBS = ["1ci0A", "1g58B", "1hujA", "1ibsA", "1jlwA", "1kzlA", "1m3uA", "1mozA", "1pv1A", "1qmvA", "1riiA", "1v9sB", "1w7wB", "1x1oB", "1ypiA", "1znnA", "2a84A", "2bcgY", "2br9A", "2cjmC", "2esfA", "2fliA", "2gv5D", "1b4tA", "1efvB", "1gv3A", "1hurA", "1ky2A", "1okcA", "1r6mA", "1xtdA", "1ysbA", "1zwkA", "2aiuA", "2cfeA", "2cnvA", "2eu8A", "2g0nB"]
#PDBS = ["1b4tA"]
FASTA_PATH = "/Users/qian/Desktop/evol_sim_vs_rosetta/sequences/"
NATURAL_FILE_PATH = "aligned_sequences/"
NATURAL_FILE_POSFIX = "Aligned_Sequences.fasta"
ROSETTA_FILE_PATH = "designed_sequences_fasta/"
ROSETTA_FILE_POSFIX = "designed_seqs.fasta"
EVOLVED_FILE_PATH = "designed_sequences_fasta/"
EVOLVED_FILE_POSFIX = "evolved_seqs.fasta"

rosetta_array = []
evolved_array = []
natural_array = []

pdb_names = []
chain_names = []

for protein in PDBS:
	pdb_id = protein[0:4].upper()
	chain_id = protein[4:5].upper()
	pdb_names.append(pdb_id)
	chain_names.append(chain_id)
		
	natural_file = FASTA_PATH + NATURAL_FILE_PATH + pdb_id + "_" + chain_id + "_" + NATURAL_FILE_POSFIX
	rosetta_file = FASTA_PATH + ROSETTA_FILE_PATH + pdb_id + "_" + chain_id + "_" + ROSETTA_FILE_POSFIX
	evolved_file = FASTA_PATH + EVOLVED_FILE_PATH + pdb_id + "_" + chain_id + "_" + EVOLVED_FILE_POSFIX
	
	rosetta = cal_divergence(rosetta_file, natural_file, False)
	evolved = cal_divergence(evolved_file, natural_file, False)
	natural = cal_divergence(natural_file, natural_file, True)
	
	rosetta_array.append(mean(rosetta))
	evolved_array.append(mean(evolved))
	natural_array.append(mean(natural))
	
all_data = []
all_data.append(rosetta_array)
all_data.append(evolved_array)
all_data.append(natural_array)

#open file name
divergence_filename = "mean_seq_divergence.csv"
divergence_file = open(divergence_filename, "w")
divergence_file.write("PDB\tchain\trosetta\tevolved\tnatural\n")
divergence_length = len(pdb_names)
length_counter = 0
while(length_counter<divergence_length):
    rosetta = rosetta_array[length_counter]
    evolved = evolved_array[length_counter]
    natural = natural_array[length_counter]
    divergence_filestring = pdb_names[length_counter] +"\t" + chain_names[length_counter] + "\t" + str(rosetta) + "\t" + str(evolved) + "\t" + str(natural) + "\n"
    divergence_file.write(divergence_filestring)
    length_counter = length_counter + 1
divergence_file.close()

fig = plt.figure(dpi = 500)
rcParams['figure.figsize'] = [8,6]
rcParams['font.size'] = 20     
rcParams['lines.linewidth'] = 2


ax = axes([0.15, 0.15, 0.8, 0.8])
b1 = ax.boxplot(all_data, sym = 'ko')
setp(b1['whiskers'], color = 'black', linestyle = '-')
setp(b1['boxes'], color =  'black')
setp(b1['caps'], color = 'black')
setp(b1['medians'], color = 'black')
setp(b1['fliers'], color  = 'black')
ylabel("Mean sequence divergence")
ax.spines['top'].set_visible(False)
ax.spines['right'].set_visible(False)
ax.get_xaxis().tick_bottom()
ax.get_yaxis().tick_left()
plt.ylim(0, 1)
plt.xticks([1, 2, 3], ["designed", "evolved", "natural"])
gca().set_yticklabels(['{:.0f}%'.format(x*100) for x in gca().get_yticks()]) 

save_fig_title = "mean_divergence" + ".eps"
savefig(save_fig_title, format = None)

#plt.show()
