import sys, os, math, string, re, gzip, urllib, shutil, Bio
import subprocess
import cStringIO
from numpy import *
from pylab import *
import matplotlib.pyplot as plt
import analysis_functions

#Date Last Updated: Apr 15 2017
#Description: This is one of the scripts that graphs the data for the paper.


mean_KL_method_file = open("graph_mean_KL_all_method_data.csv", "r")
mean_KL_method_data = mean_KL_method_file.readlines()
mean_KL_method_file.close()
header = mean_KL_method_data.pop(0)

mean_KL_method_ordered_file = open("graph_mean_KL_all_method_data_ordered.csv", "r")
mean_KL_method_ordered_data = mean_KL_method_ordered_file.readlines()
mean_KL_method_ordered_file.close()
ordered_header = mean_KL_method_ordered_data.pop(0)

all_method_data = []
all_method_mean_KL_data_array = [] 
for line in mean_KL_method_data:
    data = re.split("\t", line)
    data.pop(0)
    data_array = analysis_functions.make_array(data)
    all_method_data.append(data_array)
    mean_KL_method_values_array = analysis_functions.make_array(data)
all_method_mean_KL_data_array = array(all_method_data)

all_method_ordered_data = []
all_method_ordered_mean_KL_data_array = [] 
for line in mean_KL_method_ordered_data:
    data = re.split("\t", line)
    data.pop(0)
    data_array = analysis_functions.make_array(data)
    all_method_ordered_data.append(data_array)
    mean_KL_ordered_method_values_array = analysis_functions.make_array(data)
all_method_ordered_mean_KL_data_array = array(all_method_ordered_data)

#Make boxplot for method vs Mean KL with unordered boxplot
fig = plt.figure(dpi = 500)
rcParams['figure.figsize'] = [8,6]
rcParams['font.size'] = 20     
rcParams['lines.linewidth'] = 2
#ax1 = axes([0.12, 0.15, 0.8, 0.8])
ax1 = axes([0.15, 0.12, 0.8, 0.8])
b1 = boxplot(all_method_mean_KL_data_array, sym = "ko")
#text(0.15, 5.0, "A", fontweight = 'bold', ha = 'center', va = 'center', fontsize = 20)
setp(b1['whiskers'], color = 'black', linestyle = '-')
setp(b1['boxes'], color =  'black')
setp(b1['caps'], color = 'black')
setp(b1['medians'], color = 'black')
setp(b1['fliers'], color  = 'black')
#xlabel("methods")
ylabel("Mean KL Divergence")
ax1.spines['top'].set_visible(False)
ax1.spines['right'].set_visible(False)
ax1.get_xaxis().tick_bottom()
ax1.get_yaxis().tick_left()
plt.yticks([0.0, 1.0, 2.0, 3.0, 4.0, 5.0], ["0.0", "1.0", "2.0", "3.0", "4.0", "5.0"])
plt.xticks([1, 2, 3], ["designed", "evolved", "natural"])

save_fig_title = "KL" + ".eps"
savefig(save_fig_title, format = None)

'''
#This creates the combined boxplot figure for the Mean KL
rcParams['font.size'] = 20     
rcParams['figure.figsize'] = [14,6]
rcParams['lines.linewidth'] = 2

fig2 = plt.figure(dpi = 500)
ax1 = axes([0.07, 0.118, 0.43, 0.849])
grid()
b1 = boxplot(all_method_mean_KL_data_array, sym = "ko")
text(0.15, 5.0, "A", fontweight = 'bold', ha = 'center', va = 'center', fontsize = 20)
setp(b1['whiskers'], color = 'black', linestyle = '-')
setp(b1['boxes'], color =  'black')
setp(b1['caps'], color = 'black')
setp(b1['medians'], color = 'black')
setp(b1['fliers'], color  = 'black')
#xlabel("methods")
ylabel("Mean KL Divergence")
ax1.spines['top'].set_visible(False)
ax1.spines['right'].set_visible(False)
ax1.get_xaxis().tick_bottom()
ax1.get_yaxis().tick_left()
plt.yticks([0.0, 1.0, 2.0, 3.0, 4.0, 5.0], ["0.0", "1.0", "2.0", "3.0", "4.0", "5.0"])
plt.xticks([1, 2, 3], ["designed \n sequences", "evolved \n sequences", "natural \n sequences"])


#Make boxplot for method vs Mean KL with ordered boxplot  
grid()
ax2 = axes([0.56, 0.118, 0.43, 0.849])
b2 = boxplot(all_method_ordered_mean_KL_data_array, sym = "ko")
text(0.15, 5.0, "B", fontweight = 'bold', ha = 'center', va = 'center', fontsize = 20)
setp(b2['whiskers'], color = 'black', linestyle = '-')
setp(b2['boxes'], color =  'black')
setp(b2['caps'], color = 'black')
setp(b2['medians'], color = 'black')
setp(b2['fliers'], color  = 'black')
#xlabel("methods")
ax2.get_xaxis().tick_bottom()
ax2.get_yaxis().tick_left()
ax2.spines['top'].set_visible(False)
ax2.spines['right'].set_visible(False)
plt.yticks([0.0, 1.0, 2.0, 3.0, 4.0, 5.0], ["0.0", "1.0", "2.0", "3.0","4.0", "5.0"])
plt.xticks([1, 2, 3], ["designed \n sequences", "evolved \n sequences", "natural \n sequences"])
save_fig_title = "KL" + ".eps"
savefig(save_fig_title, format = None)
'''


