Written/Last Updated on Apr 20, 2017 by Qian Jiang
 
Directories:

analysis/
	Contains all of the scripts used to analyze the yeast-proteins data set
  unordered_analysis/
		Contains the scripts that are used to analyze the yeast-proteins data set (This is for the unranked data)
	ordered_analysis/
		Contains the scripts that are used to analyze the yeast-proteins data set (This is for the ranked data)
	analysis_functions.py
		This a helper file that contains a series of functions used in the analysis				
	calculate_distribution_files.py
		This is a script that calculates all of the amino acid count data for the designed yeast proteins.
	calculate_natural_distribution_files.py
		This is a script that calculates all of the amino acid count data for the natural yeast proteins.
	generate_graph_data.py
		This is a script that calculates all of the entropy, KL-Divergence, Entropy-Entropy, Entropy-RSA, Entropy-iCN and Entropy-iWCN correlation analysis
	

graphing/
	A series of python scripts that analysis the data and produce the figures for the paper.
	analysis_functions.py
		This a helper file that contains a series of functions used in the analysis				
	get_mean_effective_num_plot.py
    This script is the script that is used to create Figure 1.
  get_cor_combo.py
    This script is the script that is used to create Figures 2, and S1.
  get_cor_entropy_RSA_comb.py
    This script is the script that is used to create Figure 3.
  get_combo_freq_plots.py
    This script is the script that is used to create Figure 4.
  get_KL_boxplots.py
    This script is the script that is used to create Figure 5.
  get_entropy_vs_rsa_38.py
    This script is the script that is used to create Figure S2.
  get_entropy_vs_entropy_38.py
    This script is the script that is used to create Figure S3.
    
t_test/
  Contains all of the R scripts used in the analysis
	t_test.R 
		This is a script that performs the statistics in the manuscript. 

rosetta_scripts/
	Contains the two shell scripts that were used to run RosettaDesign and Backrub
	ucsf_yeash.sh
		This is a shell script that was used to create the flexible backbon protein designs
		For arguments you must give it the PDB name and the temperature used for Backrub
	ucsf_yeash_fixed.sh
		This the shell script that was used to create the fixed backbone protein design
		For arguments you must give it the PDB name and the temperature used for Backrub.
		Since it is a fixed design we just set the temperature to 0.0 purely for identification purposes.

sequences/
	Contains all of the sequences (both designed and natural) for the project
	aligned_sequences/
		These are the natural alignments
	designed_sequences/
		These are the designed sequences extracted from the yeast-protein designed proteins 
	

structures/
	This contains the PDB Files that were used as templates in the yeast-proteins design process
	
