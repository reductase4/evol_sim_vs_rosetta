Written/Last Updated on Aug 22, 2017 by Qian Jiang
 
Directories:

analysis/
    Contains all of the scripts used to analyze the yeast-proteins data set

    analysis_functions.py
        This a helper file that contains a series of functions used in the analysis				
    calculate_distribution_files.py
        This is a script that calculates the amino acid count data for the designed proteins.
    calculate_natural_distribution_files.py
        This is a script that calculates the amino acid count data for the natural proteins.
    generate_graph_data.py
        This is a script that calculates the entropy, KL-Divergence, Entropy-Entropy, Entropy-RSA, Entropy-iCN and Entropy-iWCN correlation.
	

graphing/
    A series of python scripts that analysis the data and produce the figures for the paper.

    analysis_functions.py
        This a helper file that contains a series of functions used in the analysis
    get_seq_divergence.py
        This script is the script that is used to create Figure 1.
    get_KL_boxplots.py
        This script is the script that is used to create Figure 2.
    get_mean_effective_num_plot.py
        This script is the script that is used to create Figure 3.
    get_cor_combo.py
        This script is the script that is used to create Figures 4, 7 and S2.
    get_cor_entropy_RSA.py
        This script is the script that is used to create Figure 5.
    get_re_simulated_combo.py
        This script is the script that is used to create Figure 6.
    get_combo_freq_plots.py
        This script is the script that is used to create Figure S1.
    density.R
        This script is the script that is used to create Figure S3.
    get_cor_entropy_RSA_comb.py
        This script is the script that is used to create Figure S5.
    get_seq_divergence_ten.py
        This script is the script that is used to create Figure S7.

matched_stability/
    find_matched_seq.py
        This is a script that finds score matched sequences
    density.R
        This script is the script that is used to create Figure S4.

max_score_evol_sim/
    density.R
        This script is the script that is used to create Figure S6.

t_test/
    Contains all of the R scripts used in the analysis

    t_test.R 
        This is a script that performs the statistics in the manuscript. 

methods_scripts/
    Contains scripts that were used to run minimization, RosettaDesign and evolutionary simulation

    score_mutant_AIT_v1_1.py
        This is a script that was used to minimize native structures
        For arguments you must give it the PDB name.
    design_protein_fixed.sh
        This the shell script that was used to create the fixed backbone protein design
        For arguments you must give it the PDB name and the temperature used for Backrub.
        For a fixed design, we just set the temperature to 0.0 purely for identification purposes.
    score_mutant_AIT_v1_1.py
        This is a script that was used to minimize native structures
        For arguments you must give it the PDB name and score extracted from minimized structures.

sequences/
    Contains all of the sequences (both designed and natural) for the project
    
    aligned_sequences/
        These are the natural alignments
    designed_sequences/
        These are the designed sequences extracted from the yeast-protein designed proteins 	

sequences_ten/
    Contains all of the sequences (both designed and natural) for evolved_from_design analysis
    
    aligned_sequences/
        These are the natural alignments
    designed_sequences/
        These are the designed sequences extracted from the yeast-protein designed proteins 	

structures/
    This contains the PDB Files that were used as templates in the yeast-proteins design process
