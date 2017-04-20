import sys, os, math, string, re, gzip, urllib, shutil, Bio
import subprocess
import cStringIO
from numpy import *
from pylab import *
import matplotlib.pyplot as plt
import analysis_functions as af
import matplotlib.cm as cm

all_temp_entropy_values = []
all_temp_entropy_values = []
all_temp_entropy_values = []


natural_mean_entropy_values = []
designed_mean_entropy_values_rosetta = []
designed_mean_entropy_values_evolved = []

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


fig = plt.figure(dpi = 500)
rcParams['figure.figsize'] = [8,6]
rcParams['font.size'] = 20     
rcParams['lines.linewidth'] = 2


ax = axes([0.12, 0.15, 0.8, 0.8])
b1 = ax.boxplot(all_temp_entropy_values_array, sym = 'ko')
setp(b1['whiskers'], color = 'black', linestyle = '-')
setp(b1['boxes'], color =  'black')
setp(b1['caps'], color = 'black')
setp(b1['medians'], color = 'black')
setp(b1['fliers'], color  = 'black')
xlabel("Methods")
#ylabel("Mean entropy # of amino acids")
ylabel(r'$\mathbf{e^H}$')
ax.spines['top'].set_visible(False)
ax.spines['right'].set_visible(False)
ax.get_xaxis().tick_bottom()
ax.get_yaxis().tick_left()
plt.ylim(0, 6)
plt.xticks([1, 2, 3], ["FB", "ES", "NS"])

save_fig_title = "mean_entropy" + ".pdf"
savefig(save_fig_title, format = None)

#plt.show()
