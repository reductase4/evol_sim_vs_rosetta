import sys, os, math, string, re, gzip, urllib, shutil, Bio
import subprocess
import cStringIO
from numpy import *
from pylab import *
import matplotlib.pyplot as plt
import analysis_functions as af

#Date Last Updated: Jun 20 2017
#Description: This is one of the scripts that graphs the data for the paper.

def get_plot_format(yaxis_label, ax):
    ax.set_ylabel(yaxis_label)
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    ax.get_xaxis().tick_bottom()
    ax.get_yaxis().tick_left()
    ax.set_xticks([1, 2, 3])
    ax.set_xticklabels(["designed", "evolved \n from design", "natural"])

# mean KL data
mean_KL_method_file = open("graph_mean_KL_all_method_data.csv", "r")
mean_KL_method_data = mean_KL_method_file.readlines()
mean_KL_method_file.close()
header = mean_KL_method_data.pop(0)

all_method_data = []
all_method_mean_KL_data_array = [] 
for line in mean_KL_method_data:
    data = re.split("\t", line)
    data.pop(0)
    data_array = af.make_array(data)
    all_method_data.append(data_array)
    mean_KL_method_values_array = af.make_array(data)
all_method_mean_KL_data_array = array(all_method_data)

# effective number data
all_temp_entropy_values = []

protein_file_name = "graph_mean_data_natural.csv"
[natural_mean_RSA_values, natural_mean_entropy_values, natural_cor_entropy_RSA_values, natural_mean_split_KL_values, natural_cor_entropy_icn_values, natural_cor_entropy_iwcn_values] = af.get_mean_designed_data(protein_file_name)
   
protein_file_name = "graph_mean_data_rosetta.csv"
[designed_mean_RSA_values_rosetta, designed_mean_entropy_values_rosetta, designed_cor_entropy_RSA_values_rosetta, designed_mean_KL_values_rosetta, designed_cor_entropy_icn_values_rosetta, designed_cor_entropy_iwcn_values_rosetta] = af.get_mean_designed_data(protein_file_name)

  
protein_file_name = "graph_mean_data_evolved.csv"
[designed_mean_RSA_values_evolved, designed_mean_entropy_values_evolved, designed_cor_entropy_RSA_values_evolved, designed_mean_KL_values_evolved, designed_cor_entropy_icn_values_evolved, designed_cor_entropy_iwcn_values_evolved] = af.get_mean_designed_data(protein_file_name)


all_temp_entropy_values.append(designed_mean_entropy_values_rosetta)
all_temp_entropy_values.append(designed_mean_entropy_values_evolved)
all_temp_entropy_values.append(natural_mean_entropy_values)

all_temp_entropy_values_array = []


for element in all_temp_entropy_values:
    new_array = []
    for num in element:
        value = num
        value = value.rstrip()
        new_value = float(value)
        new_array.append(new_value)
    all_temp_entropy_values_array.append(new_array)

# RSA - Effective # correlation data
all_method_cor_entropy_RSA_values = []

all_method_cor_entropy_RSA_values.append(designed_cor_entropy_RSA_values_rosetta)
all_method_cor_entropy_RSA_values.append(designed_cor_entropy_RSA_values_evolved)
all_method_cor_entropy_RSA_values.append(natural_cor_entropy_RSA_values)


all_method_cor_entropy_RSA_values_length = len(all_method_cor_entropy_RSA_values[0])
all_method_cor_entropy_RSA_values_array = []


for element in all_method_cor_entropy_RSA_values:
    new_array = []
    for num in element:
        value = num
        value = value.rstrip()
        new_value = float(value)
        new_array.append(new_value)
    all_method_cor_entropy_RSA_values_array.append(new_array)

# effective  number correlation data
protein_file_name = "graph_entropy_corr.csv"
[natural_rosetta_corr_values, natural_evolved_corr_values, rosetta_evolved_corr_values] = af.get_entropy_corr_data(protein_file_name)



# plot
A = rand(8,8)    
w, h = figaspect(A)
fig2 = plt.figure(dpi = 500, figsize=(w, h))

rcParams['font.size'] = 8      
rcParams['lines.linewidth'] = 1

p1 = plt.subplot(221)
p2 = plt.subplot(222)
p3 = plt.subplot(223)
p4 = plt.subplot(224)
subplots_adjust(left=0.12, bottom=0.1, right=0.96, top=0.96, wspace=0.35, hspace=0.35)

#Make boxplot for method vs Mean KL with unordered boxplot
get_plot_format("Mean KL Divergence", p1)
b1 = p1.boxplot(all_method_mean_KL_data_array, sym = "ko")
p1.text(0, 5.0, "A", fontweight = 'bold', ha = 'center', va = 'center', fontsize = 8)
setp(b1['whiskers'], color = 'black', linestyle = '-')
setp(b1['boxes'], color =  'black')
setp(b1['caps'], color = 'black')
setp(b1['medians'], color = 'black')
setp(b1['fliers'], color  = 'black', markersize=4)
p1.set_yticks([0.0, 1.0, 2.0, 3.0, 4.0, 5.0]) 
p1.set_yticklabels(["0.0", "1.0", "2.0", "3.0", "4.0", "5.0"])

#Make boxplot for effective number correlation
get_plot_format("Mean effective # of amino acids", p2)
b2 = p2.boxplot(all_temp_entropy_values_array, sym = 'ko')
p2.text(0, 6.0, "B", fontweight = 'bold', ha = 'center', va = 'center', fontsize = 8)
setp(b2['whiskers'], color = 'black', linestyle = '-')
setp(b2['boxes'], color =  'black')
setp(b2['caps'], color = 'black')
setp(b2['medians'], color = 'black')
setp(b2['fliers'], color  = 'black', markersize=4)
p2.set_ylim(0, 6)

# effective  number correlation
def plt_plot_diagonal_single(x, y, ax, xlabel, ylabel):
	x = np.array(x, dtype='|S32')
	y = np.array(y, dtype='|S32')
	x = x.astype(np.float)
	y = y.astype(np.float)

	ax.plot(x, y, 'ko', markersize=4)
	ax.plot([-0.2, 0.6], [-0.2, 0.6], '--k')
	ax.spines['top'].set_visible(False)
	ax.spines['right'].set_visible(False)
	ax.get_xaxis().tick_bottom()
	ax.get_yaxis().tick_left()
	ax.set_xlim(-0.4, 0.8)
	ax.set_ylim(-0.4, 0.8)
	ax.set_xlabel(xlabel)
	ax.set_ylabel(ylabel)
	ax.set_yticks([-0.4, -0.2, 0.0, 0.2, 0.4, 0.6, 0.8])
	ax.set_xticks([-0.2, 0.0, 0.2, 0.4, 0.6, 0.8])
	ax.set_yticklabels(["-0.4", "-0.2", "0.0", "0.2", "0.4", "0.6", "0.8"])
	ax.set_xticklabels(["-0.2", "0.0", "0.2", "0.4", "0.6", "0.8"])
	
plt_plot_diagonal_single(natural_rosetta_corr_values, natural_evolved_corr_values, p3, 'Correlation \n' + r'$\mathregular{n_{eff}}$ (natural) vs. $\mathregular{n_{eff}}$ (designed)', 'Correlation ' + r'$\mathregular{n_{eff}}$ (natural) vs.' + '\n' + r'$\mathregular{n_{eff}}$ (evolved from designed)')
p3.text(-0.6, 0.8, "C", fontweight = 'bold', ha = 'center', va = 'center', fontsize = 8)

# RSA - Effective # correlation


get_plot_format("Correlation RSA vs. Effective # ", p4)
b4 = p4.boxplot(all_method_cor_entropy_RSA_values_array, sym = 'ko')
p4.text(0, 0.8, "D", fontweight = 'bold', ha = 'center', va = 'center', fontsize = 8)
setp(b4['whiskers'], color = 'black', linestyle = '-')
setp(b4['boxes'], color =  'black')
setp(b4['caps'], color = 'black')
setp(b4['medians'], color = 'black')
setp(b4['fliers'], color  = 'black', markersize=4)
p4.set_ylim(-0.4, 0.8)
p4.set_yticks([-0.4, -0.2, 0.0, 0.2, 0.4, 0.6, 0.8])

save_fig_title = "combo_resimulated" + ".eps"
savefig(save_fig_title, format = None)


