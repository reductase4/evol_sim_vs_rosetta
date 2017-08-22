import sys, os, math, string, re, gzip, urllib, shutil, Bio
import subprocess
import cStringIO
from numpy import *
from pylab import *
import matplotlib.pyplot as plt
import analysis_functions as af

#Date Last Updated: Apr 15 2017
#Description: This is one of the scripts that graphs the data for the paper.

#matplotlib.rcParams['backend'] = "Qt4Agg"

def get_plot_color(natural_cor_entropy_RSA_value):
    if(natural_cor_entropy_RSA_value < 0.1): 
        color = cm.hot_r(natural_cor_entropy_RSA_value + 0.2)
        #color = cm.PuRd(natural_cor_entropy_RSA_value + 0.2)
    else: 
        color = cm.hot_r(natural_cor_entropy_RSA_value)
        #color = cm.PuRd(natural_cor_entropy_RSA_value)
    return color

def get_plot_format_cor_plot(xaxis_label, yaxis_label, ax):
 
    #xlabel(xaxis_label)
    #ylabel(yaxis_label)
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)

    ax.get_xaxis().tick_bottom()
    ax.get_yaxis().tick_left()
    plt.xlim(-0.3, 0.4)
    plt.ylim(-0.4, 0.8)
    #plt.xticks([-0.2, rosetta, 0.2, 0.4, 0.6, 0.8, 1.0, 1.2, 1.4], ["FB", "rosetta","0.2", "0.4", "0.6", "0.8", "1.0", "1.2", "NS"])
    plt.xticks([-0.2, 0.05, 0.3], ["designed \n sequences", "evolved \n sequences", "natural \n sequences"])
    plt.yticks([-0.4, -0.2, 0.0, 0.2, 0.4, 0.6, 0.8])

count = 1

pdb_names = []
natural_data = []
natural_cor_entropy_RSA_values = []
all_method_cor_entropy_RSA_values = []

modified_methods = [-0.2, 0.05, 0.3]
modified_method_array = array(modified_methods)

protein_file = open("graph_mean_data_natural.csv", "r")
natural_protein_data = protein_file.readlines()
protein_file.close()

header = natural_protein_data.pop(0)
for line in natural_protein_data:
    data = re.split("\t", line)
    natural_data.append(data)

for data in natural_data:
    pdb_names.append(data[0])
    natural_cor_entropy_RSA_values.append(data[4])

natural_cor_entropy_RSA_values_array = af.make_array(natural_cor_entropy_RSA_values)

protein_file_name = "graph_mean_data_rosetta.csv"
[designed_mean_RSA_values_rosetta, designed_mean_entropy_values_rosetta, designed_cor_entropy_RSA_values_rosetta, designed_mean_KL_values_rosetta, designed_cor_entropy_iwcn_values_rosetta] = af.get_mean_designed_data(protein_file_name)

designed_cor_entropy_RSA_values_array_rosetta = array(designed_cor_entropy_RSA_values_rosetta)  

protein_file_name = "graph_mean_data_evolved.csv"
[designed_mean_RSA_values_evolved, designed_mean_entropy_values_evolved, designed_cor_entropy_RSA_values_evolved, designed_mean_KL_values_evolved, designed_cor_entropy_iwcn_values_evolved] = af.get_mean_designed_data(protein_file_name)

designed_cor_entropy_RSA_values_array_evolved = array(designed_cor_entropy_RSA_values_evolved)  


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

all_method_cor_entropy_values_transpose = transpose(all_method_cor_entropy_RSA_values_array)

#This makes the combo boxplot and line plot for the correlation between RSA and Entropy
(m,n) = all_method_cor_entropy_values_transpose.shape

#Make boxplot for RSA correlation
fig = plt.figure(dpi = 500)
rcParams['figure.figsize'] = [8,6]
rcParams['font.size'] = 20     
rcParams['lines.linewidth'] = 2
ax = axes([0.15, 0.12, 0.8, 0.8])

get_plot_format_cor_plot("methods", "RSA - Effective # Correlation", ax)
b3 = boxplot(all_method_cor_entropy_RSA_values_array, sym = 'ko')
#text(0, 0.8, "A", fontweight = 'bold', ha = 'center', va = 'center', fontsize = 20)
setp(b3['whiskers'], color = 'black', linestyle = '-')
setp(b3['boxes'], color =  'black')
setp(b3['caps'], color = 'black')
setp(b3['medians'], color = 'black')
setp(b3['fliers'], color  = 'black')
ylabel("Correlation RSA vs. Effective #")
ax.spines['top'].set_visible(False)
ax.spines['right'].set_visible(False)
ax.get_xaxis().tick_bottom()
ax.get_yaxis().tick_left()
plt.xticks([1, 2, 3], ["designed", "evolved", "natural"])
save_fig_title = "Cor_Mean_entropy_RSA" + ".eps"
savefig(save_fig_title, format = None)

'''
A = rand(6,14)    
w, h = figaspect(A)
fig1 = plt.figure(dpi = 500, figsize=(w, h))

rcParams['font.size'] = 20      
rcParams['figure.figsize'] = [14,6]
rcParams['lines.linewidth'] = 2

i = 0
ax = axes([0.1, 0.15, 0.38, 0.75])
get_plot_format_cor_plot("methods", "RSA - Effective # Correlation", ax)
b3 = boxplot(all_method_cor_entropy_RSA_values_array, sym = 'ko')
text(0, 0.8, "A", fontweight = 'bold', ha = 'center', va = 'center', fontsize = 20)
setp(b3['whiskers'], color = 'black', linestyle = '-')
setp(b3['boxes'], color =  'black')
setp(b3['caps'], color = 'black')
setp(b3['medians'], color = 'black')
setp(b3['fliers'], color  = 'black')
ylabel("Correlation")
ax.spines['top'].set_visible(False)
ax.spines['right'].set_visible(False)
ax.get_xaxis().tick_bottom()
ax.get_yaxis().tick_left()
plt.xticks([1, 2, 3], ["designed \n sequences", "evolved \n sequences", "natural \n sequences"])


while (i < m): 
    #print i
    PDB = pdb_names[i]
    cor_entropy_RSA_array = all_method_cor_entropy_values_transpose[i,:]
    natural_cor_entropy_RSA = cor_entropy_RSA_array[2]
    
    ax2 = axes([0.6, 0.15, 0.38, 0.75])   
    get_plot_format_cor_plot("methods", "RSA - Effective # Correlation", ax2)
    p1 = plot(modified_method_array, cor_entropy_RSA_array, color = get_plot_color(natural_cor_entropy_RSA), linestyle = "-", marker = "o")
    i = i + 1
text(-0.41, 0.8, "B", fontweight = 'bold', ha = 'center', va = 'center', fontsize = 20)
save_fig_title = "Cor_Mean_entropy_RSA" + ".eps"
#save_fig_title = "Cor_Mean_entropy_RSA" + ".pdf"
savefig(save_fig_title, format = None)
count = count + 1
#plt.show()
'''