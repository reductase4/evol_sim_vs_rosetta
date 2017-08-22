import sys, os, math, string, re, gzip, urllib, shutil, Bio
import subprocess
import cStringIO
from numpy import *
from pylab import *
import matplotlib.pyplot as plt
import analysis_functions as af

#Date Last Updated: Aug 18 2017
#Description: This is one of the scripts that graphs the data for the paper.

#matplotlib.rcParams['backend'] = "Qt4Agg"


def get_plot_format_cor_plot(xaxis_label, yaxis_label, ax):
 
    ax.set_xlabel(xaxis_label)
    ax.set_ylabel(yaxis_label)
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)

    ax.get_xaxis().tick_bottom()
    ax.get_yaxis().tick_left()
    ax.set_xlim(-0.3, 0.4)
    ax.set_ylim(-0.4, 0.8)
    ax.set_yticks([-0.4, -0.2, 0.0, 0.2, 0.4, 0.6, 0.8])

# proteins with matched stability
all_method_cor_entropy_RSA_values_matched = []

protein_file_name = "graph_mean_data_rosetta_matched.csv"
[designed_mean_RSA_values_rosetta, designed_mean_entropy_values_rosetta, designed_cor_entropy_RSA_values_rosetta_matched, designed_mean_KL_values_rosetta, designed_cor_entropy_iwcn_values_rosetta] = af.get_mean_designed_data(protein_file_name)

designed_cor_entropy_RSA_values_array_rosetta_matched = array(designed_cor_entropy_RSA_values_rosetta_matched)  

protein_file_name = "graph_mean_data_evolved_matched.csv"
[designed_mean_RSA_values_evolved, designed_mean_entropy_values_evolved, designed_cor_entropy_RSA_values_evolved_matched, designed_mean_KL_values_evolved, designed_cor_entropy_iwcn_values_evolved] = af.get_mean_designed_data(protein_file_name)

designed_cor_entropy_RSA_values_array_evolved_matched = array(designed_cor_entropy_RSA_values_evolved_matched)  


all_method_cor_entropy_RSA_values_matched.append(designed_cor_entropy_RSA_values_array_rosetta_matched)
all_method_cor_entropy_RSA_values_matched.append(designed_cor_entropy_RSA_values_array_evolved_matched)

all_method_cor_entropy_RSA_values_array_matched = []


for element in all_method_cor_entropy_RSA_values_matched:
    new_array = []
    for num in element:
        value = num
        value = value.rstrip()
        new_value = float(value)
        new_array.append(new_value)
    all_method_cor_entropy_RSA_values_array_matched.append(new_array)

# evolved proteins with mean and max scores
all_method_cor_entropy_RSA_values_max = []

protein_file_name = "graph_mean_data_rosetta.csv"
[designed_mean_RSA_values_rosetta, designed_mean_entropy_values_rosetta, designed_cor_entropy_RSA_values_rosetta, designed_mean_KL_values_rosetta, designed_cor_entropy_iwcn_values_rosetta] = af.get_mean_designed_data(protein_file_name)

designed_cor_entropy_RSA_values_array_rosetta = array(designed_cor_entropy_RSA_values_rosetta)  

protein_file_name = "graph_mean_data_evolved_mean.csv"
[designed_mean_RSA_values_evolved, designed_mean_entropy_values_evolved, designed_cor_entropy_RSA_values_evolved_mean, designed_mean_KL_values_evolved, designed_cor_entropy_iwcn_values_evolved] = af.get_mean_designed_data(protein_file_name)

designed_cor_entropy_RSA_values_array_evolved_mean = array(designed_cor_entropy_RSA_values_evolved_mean)  

protein_file_name = "graph_mean_data_evolved_max.csv"
[designed_mean_RSA_values_evolved, designed_mean_entropy_values_evolved, designed_cor_entropy_RSA_values_evolved_max, designed_mean_KL_values_evolved, designed_cor_entropy_iwcn_values_evolved] = af.get_mean_designed_data(protein_file_name)

designed_cor_entropy_RSA_values_array_evolved_max = array(designed_cor_entropy_RSA_values_evolved_max)  

all_method_cor_entropy_RSA_values_max.append(designed_cor_entropy_RSA_values_array_rosetta)
all_method_cor_entropy_RSA_values_max.append(designed_cor_entropy_RSA_values_array_evolved_mean)
all_method_cor_entropy_RSA_values_max.append(designed_cor_entropy_RSA_values_array_evolved_max)

all_method_cor_entropy_RSA_values_array_max = []


for element in all_method_cor_entropy_RSA_values_max:
    new_array = []
    for num in element:
        value = num
        value = value.rstrip()
        new_value = float(value)
        new_array.append(new_value)
    all_method_cor_entropy_RSA_values_array_max.append(new_array)

# plot
A = rand(5,10)    
w, h = figaspect(A)
fig2 = plt.figure(dpi = 500, figsize=(w, h))

rcParams['font.size'] = 16      
rcParams['lines.linewidth'] = 1

p1 = plt.subplot(121)
p2 = plt.subplot(122)

subplots_adjust(left=0.1, bottom=0.1, right=0.96, top=0.96, wspace=0.3, hspace=0.35)


#Make boxplot for RSA correlation of proteins with matched stability
get_plot_format_cor_plot("", "Correlation RSA vs. Effective #", p1)
b1 = p1.boxplot(all_method_cor_entropy_RSA_values_array_matched, sym = 'ko')
p1.text(0.15, 0.8, "A", fontweight = 'bold', ha = 'center', va = 'center', fontsize = 16)
setp(b1['whiskers'], color = 'black', linestyle = '-')
setp(b1['boxes'], color =  'black')
setp(b1['caps'], color = 'black')
setp(b1['medians'], color = 'black')
setp(b1['fliers'], color  = 'black')
p1.set_xticks([1, 2])
p1.set_xticklabels(["designed", "evolved"])

#Make boxplot for RSA correlation of proteins with mean and max scores
get_plot_format_cor_plot("", "", p2)
b2 = p2.boxplot(all_method_cor_entropy_RSA_values_array_max, sym = 'ko')
p2.text(0, 0.8, "B", fontweight = 'bold', ha = 'center', va = 'center', fontsize = 16)
setp(b2['whiskers'], color = 'black', linestyle = '-')
setp(b2['boxes'], color =  'black')
setp(b2['caps'], color = 'black')
setp(b2['medians'], color = 'black')
setp(b2['fliers'], color  = 'black')
p2.set_xticks([1, 2, 3])
p2.set_xticklabels(["designed", "evolved, \n mean score", "evolved, \n max score"])


save_fig_title = "Cor_Mean_entropy_RSA_matched_max_combo" + ".eps"
savefig(save_fig_title, format = None)

