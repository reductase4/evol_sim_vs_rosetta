import sys, os, math, string, re, gzip, urllib, shutil, Bio
import subprocess
import cStringIO
from numpy import *
from pylab import *
import matplotlib.pyplot as plt
import analysis_functions
import matplotlib.cm as cm


def plt_plot_diagonal(x, y, ax, xlabel, ylabel):
	x = np.array(x, dtype='|S32')
	y = np.array(y, dtype='|S32')
	x = x.astype(np.float)
	y = y.astype(np.float)

	ax.plot(x, y, 'ko', markersize=3)
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

def plt_plot_diagonal_single(x, y, ax, xlabel, ylabel):
	x = np.array(x, dtype='|S32')
	y = np.array(y, dtype='|S32')
	x = x.astype(np.float)
	y = y.astype(np.float)

	ax.plot(x, y, 'ko', markersize=5)
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


protein_file_name = "graph_entropy_corr.csv"
[natural_rosetta_corr_values, natural_evolved_corr_values, rosetta_evolved_corr_values] = analysis_functions.get_entropy_corr_data(protein_file_name)

protein_file_name = "graph_mean_data_natural.csv"
[natural_mean_RSA_values, natural_mean_entropy_values, natural_cor_entropy_RSA_values, natural_mean_split_KL_values, natural_cor_entropy_icn_values, natural_cor_entropy_iwcn_values] = analysis_functions.get_mean_designed_data(protein_file_name)
   
protein_file_name = "graph_mean_data_rosetta.csv"
[designed_mean_RSA_values_rosetta, designed_mean_entropy_values_rosetta, designed_cor_entropy_RSA_values_rosetta, designed_mean_KL_values_rosetta, designed_cor_entropy_icn_values_rosetta, designed_cor_entropy_iwcn_values_rosetta] = analysis_functions.get_mean_designed_data(protein_file_name)

  
protein_file_name = "graph_mean_data_evolved.csv"
[designed_mean_RSA_values_evolved, designed_mean_entropy_values_evolved, designed_cor_entropy_RSA_values_evolved, designed_mean_KL_values_evolved, designed_cor_entropy_icn_values_evolved, designed_cor_entropy_iwcn_values_evolved] = analysis_functions.get_mean_designed_data(protein_file_name)


A = rand(8,8)    
w, h = figaspect(A)
fig1 = plt.figure(dpi = 400, figsize=(w, h))

rcParams['font.size'] = 20      
rcParams['lines.linewidth'] = 2

ax = axes([0.2, 0.2, 0.7, 0.7])

#cor entropy-entropy
plt_plot_diagonal_single(natural_rosetta_corr_values, natural_evolved_corr_values, ax, r'$\mathregular{e^H}$' + "(NS - FB)", r'$\mathregular{e^H}$' + "(NS - ES)")
ax.text(-0.95, 0.8, "A", fontweight = 'bold', ha = 'center', va = 'center', fontsize = 8)

save_fig_title = "Cor_entropy" + ".pdf"
savefig(save_fig_title, format = None)

A = rand(8,8)    
w, h = figaspect(A)
fig2 = plt.figure(dpi = 400, figsize=(w, h))

#fig1 = plt.figure(dpi = 400)
rcParams['font.size'] = 8      
#rcParams['figure.figsize'] = [8,6]
rcParams['lines.linewidth'] = 0.5

p1 = plt.subplot(331)
p2 = plt.subplot(332)
p3 = plt.subplot(333)
p4 = plt.subplot(334)
p5 = plt.subplot(335)
p6 = plt.subplot(336)
p7 = plt.subplot(337)
p8 = plt.subplot(338)
p9 = plt.subplot(339)

subplots_adjust(left=0.12, bottom=0.118, right=0.96, top=0.92, wspace=0.5, hspace=0.56)

#cor entropy-RSA
plt_plot_diagonal(designed_cor_entropy_RSA_values_rosetta, designed_cor_entropy_RSA_values_evolved, p1, r'$\mathregular{e^H}$'+"- RSA (FB)", r'$\mathregular{e^H}$'+"- RSA (ES)")
p1.text(-0.77, 0.9, "A", fontweight = 'bold', ha = 'center', va = 'center', fontsize = 8)

plt_plot_diagonal(natural_cor_entropy_RSA_values, designed_cor_entropy_RSA_values_evolved, p2, r'$\mathregular{e^H}$'+"- RSA (NS)", r'$\mathregular{e^H}$'+"- RSA (ES)")
p2.text(-0.77, 0.9, "B", fontweight = 'bold', ha = 'center', va = 'center', fontsize = 8)

plt_plot_diagonal(designed_cor_entropy_RSA_values_rosetta, natural_cor_entropy_RSA_values, p3, r'$\mathregular{e^H}$'+"- RSA (FB)", r'$\mathregular{e^H}$'+"- RSA (NS)")
p3.text(-0.77, 0.9, "C", fontweight = 'bold', ha = 'center', va = 'center', fontsize = 8)

#cor entropy-icn
plt_plot_diagonal(designed_cor_entropy_icn_values_rosetta, designed_cor_entropy_icn_values_evolved, p4, r'$\mathregular{e^H}$'+"- iCN (FB)", r'$\mathregular{e^H}$'+"- iCN (ES)")
p4.text(-0.77, 0.9, "D", fontweight = 'bold', ha = 'center', va = 'center', fontsize = 8)

plt_plot_diagonal(natural_cor_entropy_icn_values, designed_cor_entropy_icn_values_evolved, p5, r'$\mathregular{e^H}$'+"- iCN (NS)", r'$\mathregular{e^H}$'+"- iCN (ES)")
p5.text(-0.77, 0.9, "E", fontweight = 'bold', ha = 'center', va = 'center', fontsize = 8)

plt_plot_diagonal(designed_cor_entropy_icn_values_rosetta, natural_cor_entropy_icn_values, p6, r'$\mathregular{e^H}$'+"- iCN (FB)", r'$\mathregular{e^H}$'+"- iCN (NS)")
p6.text(-0.77, 0.9, "F", fontweight = 'bold', ha = 'center', va = 'center', fontsize = 8)

#cor entropy-iwcn
plt_plot_diagonal(designed_cor_entropy_iwcn_values_rosetta, designed_cor_entropy_iwcn_values_evolved, p7, r'$\mathregular{e^H}$'+"- iWCN (FB)", r'$\mathregular{e^H}$'+"- iWCN (ES)")
p7.text(-0.77, 0.9, "G", fontweight = 'bold', ha = 'center', va = 'center', fontsize = 8)

plt_plot_diagonal(natural_cor_entropy_iwcn_values, designed_cor_entropy_iwcn_values_evolved, p8, r'$\mathregular{e^H}$'+"- iWCN (NS)", r'$\mathregular{e^H}$'+"- iWCN (ES)")
p8.text(-0.77, 0.9, "H", fontweight = 'bold', ha = 'center', va = 'center', fontsize = 8)

plt_plot_diagonal(designed_cor_entropy_iwcn_values_rosetta, natural_cor_entropy_iwcn_values, p9, r'$\mathregular{e^H}$'+"- iWCN (FB)", r'$\mathregular{e^H}$'+"- iWCN (NS)")
p9.text(-0.77, 0.9, "I", fontweight = 'bold', ha = 'center', va = 'center', fontsize = 8)

save_fig_title = "Cor_combo" + ".pdf"
savefig(save_fig_title, format = None)
