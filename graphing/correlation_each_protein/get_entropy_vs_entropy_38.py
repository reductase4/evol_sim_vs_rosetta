import sys, os, math, string, re, gzip, urllib, shutil, Bio
import subprocess
import cStringIO
from numpy import *
from pylab import *
import matplotlib.pyplot as plt
import analysis_functions
import matplotlib.cm as cm

def make_array(list_of_values):
    new_value_list = []
    for value in list_of_values:
        new_value = value.rstrip()
        #print new_value
        new_value = float(new_value)
        new_value_list.append(new_value)
    value_array = array(new_value_list)
    return value_array
    
def plt_plot_color1(x, y, ax, xlabel, ylabel):
	x = np.array(x, dtype='|S32')
	y = np.array(y, dtype='|S32')
	x = x.astype(np.float)
	y = y.astype(np.float)

	d1, = ax.plot(x, y, 'o', color="cyan", markersize=1, markeredgecolor="cyan")
	ax.spines['top'].set_visible(False)
	ax.spines['right'].set_visible(False)
	ax.spines['bottom'].set_linewidth(0.5)
	ax.spines['left'].set_linewidth(0.5)
	ax.get_xaxis().tick_bottom()
	ax.get_yaxis().tick_left()
	ax.tick_params(labelsize=6, size=1)
	#ax.set_xlabel(xlabel)
	#ax.set_ylabel(ylabel)
	ax.set_yticks([0.5, 1.0, 2.0, 3.0, 4.0, 5.0])
	ax.set_xticks([0.5, 1.0, 2.0, 3.0, 4.0, 5.0])
	#ax.set_xticklabels(["0", "0.25", "0.5", "0.75", "1"])
	ax.set_xticklabels([])
	ax.set_yticklabels([])
	ax.set_xlim(0.5,5.0)
	ax.set_ylim(0.5,5.0)
	return d1
	
def plt_plot_color2(x, y, ax, xlabel, ylabel):
	x = np.array(x, dtype='|S32')
	y = np.array(y, dtype='|S32')
	x = x.astype(np.float)
	y = y.astype(np.float)

	#d2, = ax.plot(x, y, 'o', color="wheat", markersize=1, markeredgecolor="wheat")
	d2, = ax.plot(x, y, 'o', color="magenta", markersize=1, markeredgecolor="magenta")
	ax.spines['top'].set_visible(False)
	ax.spines['right'].set_visible(False)
	ax.spines['bottom'].set_linewidth(0.5)
	ax.spines['left'].set_linewidth(0.5)
	ax.get_xaxis().tick_bottom()
	ax.get_yaxis().tick_left()
	ax.tick_params(labelsize=6, size=1)
	#ax.set_xlabel(xlabel)
	#ax.set_ylabel(ylabel)
	ax.set_yticks([0.5, 1.0, 2.0, 3.0, 4.0, 5.0])
	ax.set_xticks([0.5, 1.0, 2.0, 3.0, 4.0, 5.0])
	#ax.set_xticklabels(["0", "0.25", "0.5", "0.75", "1"])
	ax.set_xticklabels([])
	ax.set_yticklabels([])
	ax.set_xlim(0.5,5.0)
	ax.set_ylim(0.5,5.0)
	return d2


PDBS = ["1ci0A", "1g58B", "1hujA", "1ibsA", "1jlwA", "1kzlA", "1m3uA", "1mozA", "1pv1A", "1qmvA", "1riiA", "1v9sB", "1w7wB", "1x1oB", "1ypiA", "1znnA", "2a84A", "2bcgY", "2br9A", "2cjmC", "2esfA", "2fliA", "2gv5D", "1b4tA", "1efvB", "1gv3A", "1hurA", "1ky2A", "1okcA", "1r6mA", "1xtdA", "1ysbA", "1zwkA", "2aiuA", "2cfeA", "2cnvA", "2eu8A", "2g0nB"]
#PDBS = ["1ci0A"]
#PDBS = ["1ci0A", "1g58B"]
count = 1 #Count used to plot figures

A = rand(8,8)    
w, h = figaspect(A)
fig1 = plt.figure(dpi = 400, figsize=(w, h))
#rcParams['font.size'] = 10
#rcParams['lines.linewidth'] = 1


for protein in PDBS:
	pdb_part = protein[0:4].upper()
	chain_part = protein[4]
	file = "graph_data_" + pdb_part + "_" + chain_part + "_natural.csv"
	natural_proteins = file #These lines get the chain, pdb_id filename, 
	fileparts = re.split("_", file)
	pdb_id = fileparts[2].upper()
	chain_id = fileparts[3]

	natural_file = open(file)
	natural_file_data = natural_file.readlines()
	natural_data = []
	for line in natural_file_data:
		data = re.split("\t", line)
		natural_data.append(data)

	natural_entropy = natural_data[0]
	natural_entropy.pop(0)
	natural_entropy_array = analysis_functions.make_array(natural_entropy)  

	designed_filename_rosetta = "graph_data_" + pdb_id + "_" + chain_id + "_rosetta.csv"
	[designed_entropy_array_rosetta, designed_RSA_array_rosetta, designed_KL_array_rosetta] = analysis_functions.get_designed_graph_data(designed_filename_rosetta)
	designed_filename_evolved = "graph_data_" + pdb_id + "_" + chain_id + "_evolved.csv"
	[designed_entropy_array_evolved, designed_RSA_array_evolved, designed_KL_array_evolved] = analysis_functions.get_designed_graph_data(designed_filename_evolved)
	

	p1 = plt.subplot(8,10,(1+(count-1)*2))
	p2 = plt.subplot(8,10,(2+(count-1)*2))


	subplots_adjust(left=0.08, bottom=0.1, right=0.96, top=0.92, wspace=0.2,hspace=0.7)
	
	d1 = plt_plot_color1(natural_entropy_array, designed_entropy_array_rosetta, p1, "", "")
	#p1.text(0, 3, "FB", fontweight = 'bold', ha = 'center', va = 'center', fontsize = 18)
	p1.set_title(pdb_id + '_' + chain_id, fontsize=6, fontweight = 'bold')

	d2 = plt_plot_color2(natural_entropy_array, designed_entropy_array_evolved, p2, "", "")
	#p2.text(0, 3, "ES", fontweight = 'bold', ha = 'center', va = 'center', fontsize = 18)
	
	
	if count%5 == 1:
		p1.set_ylabel('$\mathregular{e^H}$',fontsize=6, fontweight = 'bold')
		p1.set_yticklabels([" ", "1", "", "3", "", "5"], fontsize=5)
	
	if count == 34 or count == 35 or count == 36 or count == 37 or count == 38:
		p1.set_xlabel('$\mathregular{e^H}$(NS)',fontsize=6, fontweight = 'bold')
		p1.set_xticklabels([" ", "1", "", "3", "", "5"], fontsize=5)
		p2.set_xlabel('$\mathregular{e^H}$(NS)',fontsize=6, fontweight = 'bold')
		p2.set_xticklabels([" ", "1", "", "3", "", "5"], fontsize=5)
		
	print count
	count = count + 1

led = fig1.legend((d1, d2), ('FB', 'ES'), 'lower center', fontsize=6, ncol=3)
led.get_frame().set_linewidth(0.5)	

save_fig_title = "entropy_vs_entropy_38" + ".pdf"
savefig(save_fig_title, format = None)

#plt.show()	
	