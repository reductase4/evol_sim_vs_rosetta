import sys, os, math, string, re, gzip, urllib, shutil, Bio
import subprocess
import cStringIO
from numpy import *
from scipy.stats import pearsonr
import analysis_functions

#Date last updated: August 7, 2evolved3
#Description: This is a script that takes all of the *.dat files with the generated AA count and RSA values and then calculates the necessary values for comparision. Ex. Entropy, KL Divergence
#Get files 
#This searches for all of the frequency data files
searchStr = "align_natural_data_array_ordered_" + "[a-zA-Z0-9_\.\-]*" +  ".dat"
method = "rosetta"
method_array = ["rosetta", "evolved"]
data = []
count = 1 

#This is a series of arrays that are used to store data
natural_mean_RSA_values = []
natural_mean_entropy_values = []
natural_cor_entropy_RSA_values = []

designed_mean_RSA_values_rosetta = []
designed_mean_RSA_values_evolved = []

natural_mean_RSA_values = []
natural_mean_entropy_values = []

KL_array_rosetta = []
KL_array_evolved = []

designed_mean_KL_values_rosetta = []
designed_mean_KL_values_evolved = []


designed_mean_entropy_values_rosetta = []
designed_mean_entropy_values_evolved = []

designed_cor_entropy_RSA_values_rosetta = []
designed_cor_entropy_RSA_values_evolved = []

split_KL_array = []
natural_mean_split_KL_values = []

buried_mean_KL_values_rosetta = []
buried_mean_KL_values_evolved = []

intermediate_mean_KL_values_rosetta = []
intermediate_mean_KL_values_evolved = []

surface_mean_KL_values_rosetta = []
surface_mean_KL_values_evolved = []

buried_mean_entropy_values_rosetta = []
buried_mean_entropy_values_evolved = []

intermediate_mean_entropy_values_rosetta = []
intermediate_mean_entropy_values_evolved = []

surface_mean_entropy_values_rosetta = []
surface_mean_entropy_values_evolved = []

buried_natural_mean_KL_values = []
intermediate_natural_mean_KL_values = []
surface_natural_mean_KL_values = []

buried_natural_mean_entropy_values = []
intermediate_natural_mean_entropy_values = []
surface_natural_mean_entropy_values = []


pdb_names = []
chain_names = []
for path, names, filename in os.walk('.',False): #Searchs for all the *.dat files for the #natural protein 
    for file in filename:
        if(re.search(searchStr, file)!=None):
            fileparts = re.split("_",file)
            pdb_id = fileparts[5].upper() #Gets the PDB Name and the chain_id
            chain_id = fileparts[6]
            chain_id = chain_id[0]
            print "Processsing file: " + file	
            pdb_names.append(pdb_id)
            chain_names.append(chain_id)

            natural_proteins = file #Load the files with the results
            designed_proteins_rosetta = "align_data_array_ordered_" + pdb_id + "_" + chain_id + "_" + str("rosetta") + ".dat"
            designed_proteins_evolved = "align_data_array_ordered_" + pdb_id + "_" + chain_id + "_" + str("evolved") + ".dat"

            split_natural_1 = "align_natural_sample1_data_array_ordered_" + pdb_id + "_" + chain_id + ".dat"
            split_natural_2 = "align_natural_sample2_data_array_ordered_" + pdb_id + "_" + chain_id + ".dat"

            raw_file_rosetta = "raw_mle_lines_ordered_rosetta.csv"
            raw_file_evolved = "raw_mle_lines_ordered_evolved.csv"
            raw_file_natural = "raw_mle_lines_ordered_natural.csv"

            #Gets the slopes/intercept data
            intercept_dict_rosetta, slope_dict_rosetta = analysis_functions.get_slope_intercept(raw_file_rosetta)
            intercept_dict_evolved, slope_dict_evolved = analysis_functions.get_slope_intercept(raw_file_evolved)
            intercept_dict_natural, slope_dict_natural = analysis_functions.get_slope_intercept_natural(raw_file_natural)

            #Calculates the quantities for comparison (ex. entropy)
            natural_distribution = analysis_functions.get_AA_distribution_KL(natural_proteins)     
            natural_entropy = analysis_functions.get_native_entropy(natural_proteins)
            natural_entropy_array = analysis_functions.make_array(natural_entropy) 
            natural_RSA = analysis_functions.get_RSA_Values(natural_proteins)
            natural_RSA_array = analysis_functions.make_array(natural_RSA)
            natural_mean_RSA_values.append(mean(natural_RSA_array)) 
            natural_mean_entropy_values.append(mean(natural_entropy_array)) 

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


            [natural_RSA_entropy_corr,p_value] = pearsonr(natural_RSA_array, natural_entropy_array)
            natural_RSA_entropy_corr = float(natural_RSA_entropy_corr)
            natural_cor_entropy_RSA_values.append(natural_RSA_entropy_corr)


            [designed_RSA_entropy_corr_rosetta,p_value] = pearsonr(designed_RSA_array_rosetta, designed_entropy_array_rosetta) 
            designed_RSA_entropy_corr_rosetta = float(designed_RSA_entropy_corr_rosetta)
            designed_cor_entropy_RSA_values_rosetta.append(designed_RSA_entropy_corr_rosetta)

            [designed_RSA_entropy_corr_evolved,p_value] = pearsonr(designed_RSA_array_evolved, designed_entropy_array_evolved) 
            designed_RSA_entropy_corr_evolved = float(designed_RSA_entropy_corr_evolved)
            designed_cor_entropy_RSA_values_evolved.append(designed_RSA_entropy_corr_evolved)

            #Calculates the KL Divergence values 
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

            #Seperates the data into three classes: buried, partially buried, and surface sites
            [buried_entropy_values_rosetta, buried_KL_values_rosetta, intermediate_entropy_values_rosetta, intermediate_KL_values_rosetta, surface_entropy_values_rosetta, surface_KL_values_rosetta] = analysis_functions.get_position_dependent_data(designed_RSA_rosetta, designed_entropy_rosetta, KL_rosetta)

            [buried_entropy_values_evolved, buried_KL_values_evolved, intermediate_entropy_values_evolved, intermediate_KL_values_evolved, surface_entropy_values_evolved, surface_KL_values_evolved] = analysis_functions.get_position_dependent_data(designed_RSA_evolved, designed_entropy_evolved, KL_evolved)

            [natural_buried_entropy_values, natural_buried_KL_values, natural_intermediate_entropy_values, natural_intermediate_KL_values, natural_surface_entropy_values, natural_surface_KL_values] = analysis_functions.get_position_dependent_data(natural_RSA, natural_entropy, split_KL)
             
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

natural_mean_buried_KL_array = array(buried_natural_mean_KL_values)
natural_mean_intermediate_KL_array = array(intermediate_natural_mean_KL_values)
natural_mean_surface_KL_array = array(surface_natural_mean_KL_values)

natural_mean_buried_entropy_array = array(buried_natural_mean_entropy_values)
natural_mean_intermediate_entropy_array = array(intermediate_natural_mean_entropy_values)
natural_mean_surface_entropy_array = array(surface_natural_mean_entropy_values)

natural_mean_graphing_filename = "graph_mean_data_ordered_natural.csv"
natural_mean_file = open(natural_mean_graphing_filename, "w")
natural_mean_file.write("PDB\tchain\tmean_RSA\tmean_entropy\tcor_entropy_RSA\tsplit_mean_KL\tintercept\tslope\n")
natural_name_length = len(pdb_names)
length_counter = 0

while(length_counter<natural_name_length):
    pdb_counter = pdb_names[length_counter]
    natural_mean_RSA = natural_mean_RSA_values[length_counter]
    natural_mean_entropy = natural_mean_entropy_values[length_counter]
    natural_cor_entropy_RSA = natural_cor_entropy_RSA_values[length_counter] 
    natural_mean_split_KL = natural_mean_split_KL_array[length_counter]
    natural_filestring = pdb_names[length_counter] +"\t" + chain_names[length_counter] + "\t" + str(natural_mean_RSA) + "\t" + str(natural_mean_entropy) + "\t" + str(natural_cor_entropy_RSA) + "\t" + str(natural_mean_split_KL) + "\t" +  str(intercept_dict_natural[pdb_counter]) + "\t" + str(slope_dict_natural[pdb_counter]) 
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
designed_mean_graphing_filename = "graph_mean_data_ordered_rosetta.csv"
designed_mean_file = open(designed_mean_graphing_filename, "w")
designed_mean_file.write("PDB\tchain\tmean_RSA\tmean_entropy\tcor_entropy_RSA\tmean_KL\tintercept\tslope\n")
designed_name_length = len(pdb_names)
while(length_counter<designed_name_length):
    pdb_counter = pdb_names[length_counter]
    designed_mean_RSA = designed_mean_RSA_values_rosetta[length_counter]
    designed_mean_entropy = designed_mean_entropy_values_rosetta[length_counter]
    designed_cor_entropy_RSA = designed_cor_entropy_RSA_values_rosetta[length_counter] 
    designed_mean_KL = designed_mean_KL_values_rosetta[length_counter]
    designed_filestring = pdb_names[length_counter] +"\t" + chain_names[length_counter] + "\t" + str(designed_mean_RSA) + "\t" + str(designed_mean_entropy) + "\t" + str(designed_cor_entropy_RSA) + "\t" + str(designed_mean_KL) +"\t" + str(intercept_dict_rosetta[pdb_counter]) + "\t" + str(slope_dict_rosetta[pdb_counter]) 
    designed_mean_file.write(designed_filestring)
    length_counter = length_counter + 1
designed_mean_file.close()

#open file name
length_counter = 0
designed_mean_graphing_filename = "graph_mean_buried_data_ordered_rosetta.csv"
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
designed_mean_graphing_filename = "graph_mean_data_ordered_evolved.csv"
designed_mean_file = open(designed_mean_graphing_filename, "w")
designed_mean_file.write("PDB\tchain\tmean_RSA\tmean_entropy\tcor_entropy_RSA\tmean_KL\tintercept\tslope\n")
designed_name_length = len(pdb_names)
while(length_counter<designed_name_length):
    pdb_counter = pdb_names[length_counter]
    designed_mean_RSA = designed_mean_RSA_values_evolved[length_counter]
    designed_mean_entropy = designed_mean_entropy_values_evolved[length_counter]
    designed_cor_entropy_RSA = designed_cor_entropy_RSA_values_evolved[length_counter] 
    designed_mean_KL = designed_mean_KL_values_evolved[length_counter]
    designed_filestring = pdb_names[length_counter] +"\t" + chain_names[length_counter] + "\t" + str(designed_mean_RSA) + "\t" + str(designed_mean_entropy) + "\t" + str(designed_cor_entropy_RSA) + "\t" + str(designed_mean_KL) + "\t" + str(intercept_dict_evolved[pdb_counter]) + "\t" + str(slope_dict_evolved[pdb_counter]) 
    designed_mean_file.write(designed_filestring)
    length_counter = length_counter + 1
designed_mean_file.close()

#open file name
length_counter = 0
designed_mean_graphing_filename = "graph_mean_buried_data_ordered_evolved.csv"
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


mean_KL_method_file = open("graph_mean_KL_all_method_data_ordered.csv", "w")
mean_KL_method_file.write("PDB\tmean_KL_method_rosetta\tmean_KL_method_evolved\tmean_KL_split_natural\n")
length_counter = 0
while(length_counter<designed_name_length):
    string = pdb_names[length_counter] + "\t" + str(designed_mean_KL_values_rosetta[length_counter]) + "\t" + str(designed_mean_KL_values_evolved[length_counter]) + "\t"+ str(natural_mean_split_KL_values[length_counter]) + "\n"
    mean_KL_method_file.write(string)
    length_counter = length_counter + 1
mean_KL_method_file.close() 

mean_KL_method_file = open("graph_mean_KL_buried_method_data_ordered.csv", "w")
mean_KL_method_file.write("PDB\tmean_KL_method_rosetta\tmean_KL_method_evolved\tmean_KL_split_natural\n")
length_counter = 0
while(length_counter<designed_name_length):
    string = pdb_names[length_counter] + "\t" + str(buried_mean_KL_values_rosetta[length_counter]) + "\t" + str(buried_mean_KL_values_evolved[length_counter]) + "\t"+ str(buried_natural_mean_KL_values[length_counter]) + "\n"
    mean_KL_method_file.write(string)
    length_counter = length_counter + 1
mean_KL_method_file.close()  

mean_KL_method_file = open("graph_mean_KL_intermediate_method_data_ordered.csv", "w")
mean_KL_method_file.write("PDB\tmean_KL_method_rosetta\tmean_KL_method_evolved\tmean_KL_split_natural\n")
length_counter = 0
while(length_counter<designed_name_length):
    string = pdb_names[length_counter] + "\t" + str(intermediate_mean_KL_values_rosetta[length_counter]) + "\t" + str(intermediate_mean_KL_values_evolved[length_counter]) + "\t"+ str(intermediate_natural_mean_KL_values[length_counter]) + "\n"
    mean_KL_method_file.write(string)
    length_counter = length_counter + 1
mean_KL_method_file.close()  

mean_KL_method_file = open("graph_mean_KL_surface_method_data_ordered.csv", "w")
mean_KL_method_file.write("PDB\tmean_KL_method_rosetta\tmean_KL_method_evolved\tmean_KL_split_natural\n")
length_counter = 0
while(length_counter<designed_name_length):
    string = pdb_names[length_counter] + "\t" + str(surface_mean_KL_values_rosetta[length_counter]) + "\t" + str(surface_mean_KL_values_evolved[length_counter]) + "\t"+ str(surface_natural_mean_KL_values[length_counter]) + "\n"
    mean_KL_method_file.write(string)
    length_counter = length_counter + 1
mean_KL_method_file.close()  






