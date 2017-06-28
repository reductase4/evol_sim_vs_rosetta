#!/usr/bin/python

import sys, os, math, string, re, gzip, urllib, shutil, Bio
import subprocess
import cStringIO
from numpy import *
from scipy.stats import spearmanr
import analysis_functions 
import calc_wcn

#Last Date Updated: Apr 15, 2017
#Description: This is a script that takes all of the *.dat files with the generated AA count and RSA values and then calculates the necessary values for comparision. Ex. Entropy, KL Divergence, iCN, iWCN, effective # of amino acid

#Get files 
#This searches for all of the unaligned sequences 
searchStr = "align_natural_data_array_" + "[a-zA-Z0-9_\.\-]*" +  ".dat"
#find all csv files that match the search string
method = "rosetta"
method_array = ["rosetta", "evolved"]
data = []
count = 1 
#These array that are used to store various types of data
natural_mean_RSA_values = []
natural_mean_entropy_values = []
natural_mean_effective_values = []
natural_cor_entropy_RSA_values = []
natural_cor_entropy_icn_values = []
natural_cor_entropy_iwcn_values = []
natural_cor_effective_RSA_values = []

#Mean RSA
designed_mean_RSA_values_rosetta = []
designed_mean_RSA_values_evolved = []

natural_mean_RSA_values = []
natural_mean_entropy_values = []

mean_iCN13_values = []
mean_iWCN_values = []

#KL-Divergence
KL_array_rosetta = []
KL_array_evolved = []

#Mean KL-Divergence
designed_mean_KL_values_rosetta = []
designed_mean_KL_values_evolved = []

#Mean Entropy
designed_mean_entropy_values_rosetta = []
designed_mean_entropy_values_evolved = []


#The correlation between RSA and Entropy at sites
designed_cor_entropy_RSA_values_rosetta = [] 
designed_cor_entropy_RSA_values_evolved = []
designed_cor_entropy_icn_values_rosetta = [] 
designed_cor_entropy_icn_values_evolved = []
designed_cor_entropy_iwcn_values_rosetta = [] 
designed_cor_entropy_iwcn_values_evolved = []

split_KL_array = []
natural_mean_split_KL_values = []

#Mean KL- Divergence for Buried Sites
buried_mean_KL_values_rosetta = []
buried_mean_KL_values_evolved = []

#Mean KL-Divergence for Partially Buried Sites
intermediate_mean_KL_values_rosetta = []
intermediate_mean_KL_values_evolved = []

#Mean KL-Divergence for Surface Sites
surface_mean_KL_values_rosetta = []
surface_mean_KL_values_evolved = []

#Mean Entropy for Buried Sites
buried_mean_entropy_values_rosetta = []
buried_mean_entropy_values_evolved = []

#Mean Entropy for Partially Buried Sites
intermediate_mean_entropy_values_rosetta = []
intermediate_mean_entropy_values_evolved = []

#Mean Entropy for Surface Sites
surface_mean_entropy_values_rosetta = []
surface_mean_entropy_values_evolved = []

buried_natural_mean_KL_values = []
intermediate_natural_mean_KL_values = []
surface_natural_mean_KL_values = []

buried_natural_mean_entropy_values = []
intermediate_natural_mean_entropy_values = []
surface_natural_mean_entropy_values = []

#Entropy correlation
natural_rosetta_entropy_corr_values = []
natural_evolved_entropy_corr_values = []
rosetta_evolved_entropy_corr_values = []

pdb_names = []
chain_names = []
for path, names, filename in os.walk('.',False): #Searchs for all the *.dat files for the #natural protein 
    for file in filename:
        if(re.search(searchStr, file)!=None):
            fileparts = re.split("_",file)
            pdb_id = fileparts[4].upper() #Gets the PDB Name and the chain_id
            chain_id = fileparts[5]
            chain_id = chain_id[0]
            print "Processsing file: " + file	
            pdb_names.append(pdb_id)
            chain_names.append(chain_id)

            natural_proteins = file #Open the files with results
            designed_proteins_rosetta = "align_data_array_" + pdb_id + "_" + chain_id + "_" + str("rosetta")  + ".dat"
            designed_proteins_evolved = "align_data_array_" + pdb_id + "_" + chain_id + "_" + str("evolved")  + ".dat"

            split_natural_1 = "align_natural_sample1_data_array_" + pdb_id + "_" + chain_id + ".dat"
            split_natural_2 = "align_natural_sample2_data_array_" + pdb_id + "_" + chain_id + ".dat"

            #Calculates all of the data for comparison (ex. entropy)
            natural_distribution = analysis_functions.get_AA_distribution_KL(natural_proteins)     
            natural_entropy = analysis_functions.get_native_entropy(natural_proteins)
            natural_entropy_array = analysis_functions.make_array(natural_entropy)
            natural_RSA = analysis_functions.get_RSA_Values(natural_proteins)
            natural_RSA_array = analysis_functions.make_array(natural_RSA)
            natural_mean_RSA_values.append(mean(natural_RSA_array)) 
            natural_mean_entropy_values.append(mean(natural_entropy_array))
            
            
            #Calculates cn & wcn
#             cn13_data = analysis_functions.get_cn13_values(pdb_id, chain_id)
#             iCN13 = cn13_data[0]
#             iCN13_array = analysis_functions.make_array(cn13_data)
#             mean_iCN13_values.append(mean(iCN13_array))
            
            iwcn_data = calc_wcn.get_iwcn_values(pdb_id, chain_id)
            iWCN_array = analysis_functions.make_array(iwcn_data)
            mean_iWCN_values.append(mean(iWCN_array))

            designed_distribution_rosetta = analysis_functions.get_AA_distribution_KL(designed_proteins_rosetta)        
            designed_entropy_rosetta = analysis_functions.get_native_entropy(designed_proteins_rosetta)
            designed_entropy_array_rosetta = analysis_functions.make_array(designed_entropy_rosetta)
            designed_RSA_rosetta = analysis_functions.get_RSA_Values(designed_proteins_rosetta)
            designed_RSA_array_rosetta = analysis_functions.make_array(designed_RSA_rosetta)
            designed_mean_RSA_values_rosetta.append(mean(designed_RSA_array_rosetta)) 
            designed_mean_entropy_values_rosetta.append(mean(designed_entropy_array_rosetta))


            designed_distribution_evolved = analysis_functions.get_AA_distribution_KL(designed_proteins_evolved)
            designed_entropy_evolved = analysis_functions.get_native_entropy(designed_proteins_evolved)
            designed_entropy_array_evolved = analysis_functions.make_array(designed_entropy_evolved)
            designed_RSA_evolved = analysis_functions.get_RSA_Values(designed_proteins_evolved)
            designed_RSA_array_evolved = analysis_functions.make_array(designed_RSA_evolved)
            designed_mean_RSA_values_evolved.append(mean(designed_RSA_array_evolved)) 
            designed_mean_entropy_values_evolved.append(mean(designed_entropy_array_evolved))

			#Calculates the RSA & entropy correlation
            [natural_RSA_entropy_corr,p_value] = spearmanr(natural_RSA_array, natural_entropy_array) 
            natural_RSA_entropy_corr = float(natural_RSA_entropy_corr)
            natural_cor_entropy_RSA_values.append(natural_RSA_entropy_corr)        

            [designed_RSA_entropy_corr_evolved,p_value] = spearmanr(designed_RSA_array_evolved, designed_entropy_array_evolved) 
            designed_RSA_entropy_corr_evolved = float(designed_RSA_entropy_corr_evolved)
            designed_cor_entropy_RSA_values_evolved.append(designed_RSA_entropy_corr_evolved)            

            [designed_RSA_entropy_corr_rosetta,p_value] = spearmanr(designed_RSA_array_rosetta, designed_entropy_array_rosetta) 
            designed_RSA_entropy_corr_rosetta = float(designed_RSA_entropy_corr_rosetta)
            designed_cor_entropy_RSA_values_rosetta.append(designed_RSA_entropy_corr_rosetta)            
            
            
            #Calculates the iWCN & entropy correlation
            [natural_iwcn_entropy_corr,p_value] = spearmanr(iWCN_array, natural_entropy_array) 
            natural_iwcn_entropy_corr = float(natural_iwcn_entropy_corr)
            natural_cor_entropy_iwcn_values.append(natural_iwcn_entropy_corr)
            

            [designed_iwcn_entropy_corr_evolved,p_value] = spearmanr(iWCN_array, designed_entropy_array_evolved) 
            designed_iwcn_entropy_corr_evolved = float(designed_iwcn_entropy_corr_evolved)
            designed_cor_entropy_iwcn_values_evolved.append(designed_iwcn_entropy_corr_evolved)

            [designed_iwcn_entropy_corr_rosetta,p_value] = spearmanr(iWCN_array, designed_entropy_array_rosetta) 
            designed_iwcn_entropy_corr_rosetta = float(designed_iwcn_entropy_corr_rosetta)
            designed_cor_entropy_iwcn_values_rosetta.append(designed_iwcn_entropy_corr_rosetta)

			#Calculates the entropy correlation
            [natural_rosetta_entropy_corr,p_value] = spearmanr(natural_entropy_array, designed_entropy_array_rosetta)
            natural_rosetta_entropy_corr = float(natural_rosetta_entropy_corr)
            natural_rosetta_entropy_corr_values.append(natural_rosetta_entropy_corr)
			
            [natural_evolved_entropy_corr,p_value] = spearmanr(natural_entropy_array, designed_entropy_array_evolved)
            natural_evolved_entropy_corr = float(natural_evolved_entropy_corr)
            natural_evolved_entropy_corr_values.append(natural_evolved_entropy_corr)
			
            [rosetta_evolved_entropy_corr,p_value] = spearmanr(designed_entropy_array_rosetta, designed_entropy_array_evolved)
            rosetta_evolved_entropy_corr = float(rosetta_evolved_entropy_corr)
            rosetta_evolved_entropy_corr_values.append(rosetta_evolved_entropy_corr)
			
            #Calculates the KL Divergence Values 
            KL_rosetta = analysis_functions.get_Kullback_Leibler(natural_distribution, designed_distribution_rosetta)
            KL_evolved = analysis_functions.get_Kullback_Leibler(natural_distribution, designed_distribution_evolved)

            KL_array_rosetta = array(KL_rosetta)
            KL_array_evolved = array(KL_evolved)
        
            designed_mean_KL_values_rosetta.append(mean(KL_array_rosetta))
            designed_mean_KL_values_evolved.append(mean(KL_array_evolved))
            
            split_1_distribution = analysis_functions.get_AA_distribution_KL(split_natural_1)
            split_2_distribution = analysis_functions.get_AA_distribution_KL(split_natural_2)
            split_KL = analysis_functions.get_Kullback_Leibler(split_1_distribution, split_2_distribution)
            split_KL_array = array(split_KL)
            natural_mean_split_KL_values.append(mean(split_KL_array))

            #These lines seperate the data into three classes: buried, partially buried, and surface sites
            [buried_entropy_values_rosetta, buried_KL_values_rosetta, intermediate_entropy_values_rosetta, intermediate_KL_values_rosetta, surface_entropy_values_rosetta, surface_KL_values_rosetta] = analysis_functions.get_position_dependent_data(designed_RSA_rosetta, designed_entropy_rosetta, KL_rosetta)

            [buried_entropy_values_evolved, buried_KL_values_evolved, intermediate_entropy_values_evolved, intermediate_KL_values_evolved, surface_entropy_values_evolved, surface_KL_values_evolved] = analysis_functions.get_position_dependent_data(designed_RSA_evolved, designed_entropy_evolved, KL_evolved)

            [natural_buried_entropy_values, natural_buried_KL_values, natural_intermediate_entropy_values, natural_intermediate_KL_values, natural_surface_entropy_values, natural_surface_KL_values] = analysis_functions.get_position_dependent_data(natural_RSA, natural_entropy, split_KL)
             
            #Store the data that was just calculated 
            buried_KL_array_rosetta = array(buried_KL_values_rosetta)
            buried_KL_array_evolved = array(buried_KL_values_evolved)
             
            buried_entropy_array_rosetta = array(buried_entropy_values_rosetta)
            buried_entropy_array_evolved = array(buried_entropy_values_evolved)
            
            buried_mean_entropy_values_rosetta.append(mean(buried_entropy_array_rosetta))
            buried_mean_entropy_values_evolved.append(mean(buried_entropy_array_evolved))

            buried_mean_KL_values_rosetta.append(mean(array(buried_KL_array_rosetta)))
            buried_mean_KL_values_evolved.append(mean(array(buried_KL_array_evolved)))

            intermediate_KL_array_rosetta = array(intermediate_KL_values_rosetta)
            intermediate_KL_array_evolved = array(intermediate_KL_values_evolved)
            
            intermediate_entropy_array_rosetta = array(intermediate_entropy_values_rosetta)
            intermediate_entropy_array_evolved = array(intermediate_entropy_values_evolved)

            intermediate_mean_entropy_values_rosetta.append(mean(intermediate_entropy_array_rosetta))
            intermediate_mean_entropy_values_evolved.append(mean(intermediate_entropy_array_evolved))
            
            intermediate_mean_KL_values_rosetta.append(mean(array(intermediate_KL_array_rosetta)))  
            intermediate_mean_KL_values_evolved.append(mean(array(intermediate_KL_array_evolved))) 

            surface_KL_array_rosetta = array(surface_KL_values_rosetta)
            surface_KL_array_evolved = array(surface_KL_values_evolved)
            
            surface_entropy_array_rosetta = array(surface_entropy_values_rosetta)
            surface_entropy_array_evolved = array(surface_entropy_values_evolved)

            surface_mean_KL_values_rosetta.append(mean(array(surface_KL_array_rosetta)))
            surface_mean_KL_values_evolved.append(mean(array(surface_KL_array_evolved)))

            surface_mean_entropy_values_rosetta.append(mean(surface_entropy_array_rosetta))
            surface_mean_entropy_values_evolved.append(mean(surface_entropy_array_evolved))

            buried_natural_mean_KL_values.append(mean(array(natural_buried_KL_values)))
            intermediate_natural_mean_KL_values.append(mean(array(natural_intermediate_KL_values)))
            surface_natural_mean_KL_values.append(mean(array(natural_surface_KL_values)))

            buried_natural_mean_entropy_values.append(mean(array(natural_buried_entropy_values)))
            intermediate_natural_mean_entropy_values.append(mean(array(natural_intermediate_entropy_values)))
            surface_natural_mean_entropy_values.append(mean(array(natural_surface_entropy_values)))

            #open file name
            natural_graphing_filename = "graph_data_" + pdb_id + "_" + chain_id + "_natural" + ".csv"
            natural_file = open(natural_graphing_filename, "w")
            natural_file.write("entropy\t" + analysis_functions.dump_csv_line3(natural_entropy_array))
            natural_file.write("RSA\t" + analysis_functions.dump_csv_line3(natural_RSA_array))
#             natural_file.write("iwcn_ca\t" + analysis_functions.dump_csv_line3(iCN13_array))
            natural_file.write("iWCN\t" + analysis_functions.dump_csv_line3(iWCN_array))
            natural_file.close()
        
            method_array = ["rosetta", "evolved"] 
            for method in method_array: #These lines write the calculated site data to to the graph data files for each methoderature treatment
                graphing_filename = "graph_data_" + pdb_id + "_" + chain_id + "_" + str(method) + ".csv"
                graph_file = open(graphing_filename, "w")
                if method == "rosetta":      
                    graph_file.write("entropy\t" + analysis_functions.dump_csv_line3(designed_entropy_array_rosetta))
                    graph_file.write("RSA\t" + analysis_functions.dump_csv_line3(designed_RSA_array_rosetta))
                    graph_file.write("KL\t" + analysis_functions.dump_csv_line3(KL_array_rosetta))
        
                elif method == "evolved": #
                    graph_file.write("entropy\t" + analysis_functions.dump_csv_line3(designed_entropy_array_evolved))
                    graph_file.write("RSA\t" + analysis_functions.dump_csv_line3(designed_RSA_array_evolved))
                    graph_file.write("KL\t" + analysis_functions.dump_csv_line3(KL_array_evolved))
                graph_file.close()

#These files print the MEAN DATA FOR EACH PDB (Seperated by tabs)
natural_mean_RSA_values_array = array(natural_mean_RSA_values)
natural_mean_entropy_values_array = array(natural_mean_entropy_values)
natural_cor_entropy_RSA_values_array = array(natural_cor_entropy_RSA_values)       
natural_mean_split_KL_array = array(natural_mean_split_KL_values)
natural_mean_effective_values_array = array(natural_mean_effective_values)

natural_mean_buried_KL_array = array(buried_natural_mean_KL_values)
natural_mean_intermediate_KL_array = array(intermediate_natural_mean_KL_values)
natural_mean_surface_KL_array = array(surface_natural_mean_KL_values)

natural_mean_buried_entropy_array = array(buried_natural_mean_entropy_values)
natural_mean_intermediate_entropy_array = array(intermediate_natural_mean_entropy_values)
natural_mean_surface_entropy_array = array(surface_natural_mean_entropy_values)

#open file name
entropy_corr_filename = "graph_entropy_corr.csv"
entropy_corr_file = open(entropy_corr_filename, "w")
entropy_corr_file.write("PDB\tchain\tnatural_rosetta_corr\tnatural_evolved_corr\trosetta_evolved_corr\n")
entropy_corr_length = len(pdb_names)
length_counter = 0
while(length_counter<entropy_corr_length):
    natural_rosetta_corr = natural_rosetta_entropy_corr_values[length_counter]
    natural_evolved_corr = natural_evolved_entropy_corr_values[length_counter]
    rosetta_evolved_corr = rosetta_evolved_entropy_corr_values[length_counter]
    entropy_corr_filestring = pdb_names[length_counter] +"\t" + chain_names[length_counter] + "\t" + str(natural_rosetta_corr) + "\t" + str(natural_evolved_corr) + "\t" + str(rosetta_evolved_corr) + "\n"
    entropy_corr_file.write(entropy_corr_filestring)
    length_counter = length_counter + 1
entropy_corr_file.close()

#open file name
natural_mean_graphing_filename = "graph_mean_data_natural.csv"
natural_mean_file = open(natural_mean_graphing_filename, "w")
natural_mean_file.write("PDB\tchain\tmean_RSA\tmean_entropy\tcor_entropy_RSA\tsplit_mean_KL\tcor_entropy_iwcn\n")
natural_name_length = len(pdb_names)
length_counter = 0
while(length_counter<natural_name_length):
    natural_mean_RSA = natural_mean_RSA_values[length_counter]
    natural_mean_entropy = natural_mean_entropy_values[length_counter]
    natural_cor_entropy_RSA = natural_cor_entropy_RSA_values[length_counter] 
    natural_mean_split_KL = natural_mean_split_KL_array[length_counter]
    #natural_cor_entropy_icn = natural_cor_entropy_icn_values[length_counter]
    natural_cor_entropy_iwcn = natural_cor_entropy_iwcn_values[length_counter]
    natural_filestring = pdb_names[length_counter] +"\t" + chain_names[length_counter] + "\t" + str(natural_mean_RSA) + "\t" + str(natural_mean_entropy) + "\t" + str(natural_cor_entropy_RSA) + "\t" + str(natural_mean_split_KL) + "\t" + str(natural_cor_entropy_iwcn) + "\n"
    natural_mean_file.write(natural_filestring)
    length_counter = length_counter + 1
natural_mean_file.close()

#open file name
natural_mean_graphing_filename = "graph_mean_buried_data_natural.csv"
natural_mean_file = open(natural_mean_graphing_filename, "w")
natural_mean_file.write("PDB\tchain\tmean_entropy\tsplit_mean_KL\n")
natural_name_length = len(pdb_names)
length_counter = 0
while(length_counter<natural_name_length):
    buried_natural_mean_entropy = natural_mean_buried_entropy_array[length_counter]
    buried_natural_mean_split_KL = natural_mean_buried_KL_array[length_counter]
    natural_filestring = pdb_names[length_counter] +"\t" + chain_names[length_counter] + "\t" + str(buried_natural_mean_entropy) + "\t" + str(natural_mean_split_KL) + "\n"
    natural_mean_file.write(natural_filestring)
    length_counter = length_counter + 1
natural_mean_file.close()


designed_mean_RSA_values_array_rosetta = array(designed_mean_RSA_values_rosetta)
designed_mean_entropy_values_array_rosetta = array(designed_mean_entropy_values_rosetta)
designed_cor_entropy_RSA_values_array_rosetta = array(designed_cor_entropy_RSA_values_rosetta)  
designed_mean_KL_values_array_rosetta = array(designed_mean_KL_values_rosetta)

buried_mean_KL_values_array_rosetta = array(buried_mean_KL_values_rosetta)
intermediate_mean_KL_values_array_rosetta = array(intermediate_mean_KL_values_rosetta)
surface_mean_KL_values_array_rosetta = array(surface_mean_KL_values_rosetta)

buried_mean_entropy_values_array_rosetta = array(buried_mean_entropy_values_rosetta)
intermediate_mean_entropy_values_array_rosetta = array(intermediate_mean_entropy_values_rosetta)
surface_mean_entropy_values_array_rosetta = array(surface_mean_entropy_values_rosetta)

#open file name
length_counter = 0
designed_mean_graphing_filename = "graph_mean_data_rosetta.csv"
designed_mean_file = open(designed_mean_graphing_filename, "w")
designed_mean_file.write("PDB\tchain\tmean_RSA\tmean_entropy\tcor_entropy_RSA\tmean_KL\tcor_entropy_iwcn\n")
designed_name_length = len(pdb_names)
while(length_counter<designed_name_length):
    designed_mean_RSA = designed_mean_RSA_values_rosetta[length_counter]
    designed_mean_entropy = designed_mean_entropy_values_rosetta[length_counter]
    designed_cor_entropy_RSA = designed_cor_entropy_RSA_values_rosetta[length_counter] 
    designed_mean_KL = designed_mean_KL_values_rosetta[length_counter]
    designed_cor_entropy_iwcn = designed_cor_entropy_iwcn_values_rosetta[length_counter]
    designed_filestring = pdb_names[length_counter] + "\t" + chain_names[length_counter] + "\t" + str(designed_mean_RSA) + "\t" + str(designed_mean_entropy) + "\t" + str(designed_cor_entropy_RSA) + "\t" + str(designed_mean_KL) + "\t" + str(designed_cor_entropy_iwcn) + "\n"
    designed_mean_file.write(designed_filestring)
    length_counter = length_counter + 1
designed_mean_file.close()

#open file name
length_counter = 0
designed_mean_graphing_filename = "graph_mean_buried_data_rosetta.csv"
designed_mean_file = open(designed_mean_graphing_filename, "w")
designed_mean_file.write("PDB\tchain\tmean_entropy\tmean_KL\n")
designed_name_length = len(pdb_names)
while(length_counter<designed_name_length):
    designed_mean_entropy = buried_mean_entropy_values_rosetta[length_counter]
    designed_mean_KL = buried_mean_KL_values_rosetta[length_counter]
    designed_filestring = pdb_names[length_counter] +"\t" + chain_names[length_counter] + "\t" + str(designed_mean_entropy)+ "\t" + str(designed_mean_KL) + "\n"
    designed_mean_file.write(designed_filestring)
    length_counter = length_counter + 1
designed_mean_file.close()


designed_mean_RSA_values_array_evolved = array(designed_mean_RSA_values_evolved)
designed_mean_entropy_values_array_evolved = array(designed_mean_entropy_values_evolved)
designed_cor_entropy_RSA_values_array_evolved = array(designed_cor_entropy_RSA_values_evolved)
designed_mean_KL_values_array_evolved = array(designed_mean_KL_values_evolved) 

buried_mean_KL_values_array_evolved = array(buried_mean_KL_values_evolved)
intermediate_mean_KL_values_array_evolved = array(intermediate_mean_KL_values_evolved)
surface_mean_KL_values_array_evolved = array(surface_mean_KL_values_evolved)

buried_mean_entropy_values_array_evolved = array(buried_mean_entropy_values_evolved)
intermediate_mean_entropy_values_array_evolved = array(intermediate_mean_entropy_values_evolved)
surface_mean_entropy_values_array_evolved = array(surface_mean_entropy_values_evolved)
#print buried_mean_entropy_values_array_evolved
#open file name
length_counter = 0
designed_mean_graphing_filename = "graph_mean_data_evolved.csv"
designed_mean_file = open(designed_mean_graphing_filename, "w")
designed_mean_file.write("PDB\tchain\tmean_RSA\tmean_entropy\tcor_entropy_RSA\tmean_KL\tcor_entropy_iwcn\n")
designed_name_length = len(pdb_names)
while(length_counter<designed_name_length):
    designed_mean_RSA = designed_mean_RSA_values_evolved[length_counter]
    designed_mean_entropy = designed_mean_entropy_values_evolved[length_counter]
    designed_cor_entropy_RSA = designed_cor_entropy_RSA_values_evolved[length_counter] 
    designed_mean_KL = designed_mean_KL_values_evolved[length_counter]
    designed_cor_entropy_iwcn = designed_cor_entropy_iwcn_values_evolved[length_counter]
    designed_filestring = pdb_names[length_counter] +"\t" + chain_names[length_counter] + "\t" + str(designed_mean_RSA) + "\t" + str(designed_mean_entropy) + "\t" + str(designed_cor_entropy_RSA) + "\t" + str(designed_mean_KL) + "\t" + str(designed_cor_entropy_iwcn) + "\n"
    designed_mean_file.write(designed_filestring)
    length_counter = length_counter + 1
designed_mean_file.close()

#open file name
length_counter = 0
designed_mean_graphing_filename = "graph_mean_buried_data_evolved.csv"
designed_mean_file = open(designed_mean_graphing_filename, "w")
designed_mean_file.write("PDB\tchain\tmean_entropy\tmean_KL\n")
designed_name_length = len(pdb_names)
while(length_counter<designed_name_length):
    designed_mean_entropy = buried_mean_entropy_values_array_evolved[length_counter]
    designed_mean_KL = buried_mean_KL_values_array_evolved[length_counter]
    designed_filestring = pdb_names[length_counter] +"\t" + chain_names[length_counter] + "\t" + str(designed_mean_entropy)+ "\t" + str(designed_mean_KL) + "\n"
    designed_mean_file.write(designed_filestring)
    length_counter = length_counter + 1
designed_mean_file.close()

mean_KL_method_file = open("graph_mean_KL_all_method_data.csv", "w")
mean_KL_method_file.write("PDB\tmean_KL_method_rosetta\tmean_KL_method_evolved\tmean_KL_split_natural\n")
length_counter = 0
while(length_counter<designed_name_length):
    string = pdb_names[length_counter] + "\t" + str(designed_mean_KL_values_rosetta[length_counter]) + "\t" + str(designed_mean_KL_values_evolved[length_counter]) + "\t"+ str(natural_mean_split_KL_values[length_counter]) + "\n"
    mean_KL_method_file.write(string)
    length_counter = length_counter + 1
mean_KL_method_file.close() 

mean_KL_method_file = open("graph_mean_KL_buried_method_data.csv", "w")
mean_KL_method_file.write("PDB\tmean_KL_method_rosetta\tmean_KL_method_evolved\tmean_KL_split_natural\n")
length_counter = 0
while(length_counter<designed_name_length):
    string = pdb_names[length_counter] + "\t" + str(buried_mean_KL_values_rosetta[length_counter]) + "\t" + str(buried_mean_KL_values_evolved[length_counter]) + "\t"+ str(buried_natural_mean_KL_values[length_counter]) + "\n"
    mean_KL_method_file.write(string)
    length_counter = length_counter + 1
mean_KL_method_file.close()  

mean_KL_method_file = open("graph_mean_KL_intermediate_method_data.csv", "w")
mean_KL_method_file.write("PDB\tmean_KL_method_rosetta\tmean_KL_method_evolved\tmean_KL_split_natural\n")
length_counter = 0
while(length_counter<designed_name_length):
    string = pdb_names[length_counter] + "\t" + str(intermediate_mean_KL_values_rosetta[length_counter]) + "\t" + str(intermediate_mean_KL_values_evolved[length_counter]) + "\t"+ str(intermediate_natural_mean_KL_values[length_counter]) + "\n"
    mean_KL_method_file.write(string)
    length_counter = length_counter + 1
mean_KL_method_file.close()  

mean_KL_method_file = open("graph_mean_KL_surface_method_data.csv", "w")
mean_KL_method_file.write("PDB\tmean_KL_method_rosetta\tmean_KL_method_evolved\tmean_KL_split_natural\n")
length_counter = 0
while(length_counter<designed_name_length):
    string = pdb_names[length_counter] + "\t" + str(surface_mean_KL_values_rosetta[length_counter]) + "\t" + str(surface_mean_KL_values_evolved[length_counter]) + "\t"+ str(surface_natural_mean_KL_values[length_counter]) + "\n"
    mean_KL_method_file.write(string)
    length_counter = length_counter + 1
mean_KL_method_file.close()  

mean_entropy_method_file = open("graph_mean_entropy_buried_method_data.csv", "w")
mean_entropy_method_file.write("PDB\tmean_entropy_method_rosetta\tmean_entropy_method_evolved\tmean_entropy_split_natural\n")
length_counter = 0
while(length_counter<designed_name_length):
    string = pdb_names[length_counter] + "\t" + str(buried_mean_entropy_values_rosetta[length_counter]) + "\t" + str(buried_mean_entropy_values_evolved[length_counter]) + "\t"+ str(buried_natural_mean_entropy_values[length_counter]) + "\n"
    mean_entropy_method_file.write(string)
    length_counter = length_counter + 1
mean_entropy_method_file.close() 

mean_entropy_method_file = open("graph_mean_entropy_intermediate_method_data.csv", "w")
mean_entropy_method_file.write("PDB\tmean_entropy_method_rosetta\tmean_entropy_method_evolved\tmean_entropy_split_natural\n")
length_counter = 0
while(length_counter<designed_name_length):
    string = pdb_names[length_counter] + "\t" + str(intermediate_mean_entropy_values_rosetta[length_counter]) + "\t" + str(intermediate_mean_entropy_values_evolved[length_counter]) + "\t"+ str(intermediate_natural_mean_entropy_values[length_counter]) + "\n"
    mean_entropy_method_file.write(string)
    length_counter = length_counter + 1
mean_entropy_method_file.close() 

mean_entropy_method_file = open("graph_mean_entropy_surface_method_data.csv", "w")
mean_entropy_method_file.write("PDB\tmean_entropy_method_rosetta\tmean_entropy_method_evolved\tmean_entropy_split_natural\n")
length_counter = 0
while(length_counter<designed_name_length):
    string = pdb_names[length_counter] + "\t" + str(surface_mean_entropy_values_rosetta[length_counter]) + "\t" + str(surface_mean_entropy_values_evolved[length_counter]) + "\t"+ str(surface_natural_mean_entropy_values[length_counter]) + "\n"
    mean_entropy_method_file.write(string)
    length_counter = length_counter + 1
mean_entropy_method_file.close()  

stats_file = open("stats_data.csv", "w")
stats_file.write("PDB\tChain\tmean_entropy_method_rosetta\tmean_entropy_method_evolved\tnatural_cor_entropy_RSA_value\n")
length_counter = 0
while(length_counter<designed_name_length):
    #print str(natural_cor_entropy_RSA_values[length_counter])
    designed_filestring = pdb_names[length_counter] +"\t" + chain_names[length_counter] + "\t" + str(designed_mean_entropy_values_rosetta[length_counter]) +  "\t" + str(designed_mean_entropy_values_evolved[length_counter]) + "\t"+ str(natural_cor_entropy_RSA_values[length_counter]) + "\n"
    stats_file.write(designed_filestring)
    length_counter = length_counter + 1
stats_file.close()




