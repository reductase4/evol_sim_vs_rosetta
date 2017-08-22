#!usr/local/bin/python
import sys, os, math, string, re, gzip, urllib, shutil, Bio, subprocess
import cStringIO 
import numpy as np
import analysis_functions as af
import random 

#Date Last Modified:  Apr 15, 2017
#Description: This file takes the designed and natural sequences for each protein at each methoderature, calculates the frequency data and prints that info to some files. 

resdict = { 'ALA': 'A', 'CYS': 'C', 'ASP': 'D', 'GLU': 'E', 'PHE': 'F',                   
            'GLY': 'G', 'HIS': 'H', 'ILE': 'I', 'LYS': 'K', 'LEU': 'L',                   
            'MET': 'M', 'ASN': 'N', 'PRO': 'P', 'GLN': 'Q', 'ARG': 'R',                   
            'SER': 'S', 'THR': 'T', 'VAL': 'V', 'TRP': 'W', 'TYR': 'Y' }          

rows = ["\"site\"", "\"RSA\"", "\"aa1\"", "\"aa2\"", "\"aa3\"", "\"aa4\"", "\"aa5\"", "\"aa6\"", "\"aa7\"", "\"aa8\"", "\"aa9\"", "\"aa10\"", "\"aa11\"", "\"aa12\"", "\"aa13\"", "\"aa14\"", "\"aa15\"", "\"aa16\"", "\"aa17\"", "\"aa18\"", "\"aa19\"", "\"aa20\""]   
PDBS = ["1ci0A", "1g58B", "1hujA", "1ibsA", "1jlwA", "1kzlA", "1m3uA", "1mozA", "1pv1A", "1qmvA", "1riiA", "1v9sB", "1w7wB", "1x1oB", "1ypiA", "1znnA", "2a84A", "2bcgY", "2br9A", "2cjmC", "2esfA", "2fliA", "2gv5D", "1b4tA", "1efvB", "1gv3A", "1hurA", "1ky2A", "1okcA", "1r6mA", "1xtdA", "1ysbA", "1zwkA", "2aiuA", "2cfeA", "2cnvA", "2eu8A", "2g0nB"]
#PDBS = ["1b4tA", "1efvB", "1gv3A", "1hurA", "1ky2A"]
methods = ["rosetta", "evolved"]

#These are the absolute PATH to our original sequence dataset files (structures, designed sequences, and natural sequences)
duncan_designed_sequence_path = "/Users/qian/Desktop/calc_score/matched_stability/sequences/designed_sequences/"
duncan_natural_sequence_path = "/Users/qian/Desktop/calc_score/matched_stability/sequences/daligned_sequences/"
duncan_structure_path = "/Users/qian/Desktop/calc_score/matched_stability/structures/"

mean_designed_protein_length = 0
designed_protein_lengths = []
for method in methods:  
    for protein in PDBS:
    	designed_length = [0]
        pdb_id = protein[0:4].upper()
        chain_id = protein[4] 
        file = pdb_id + "_" + chain_id + "_Aligned_Sequences.fasta"

        fileparts = re.split("_",file)
        print "Processsing file: " + file	
        print "PDB: " + pdb_id
        print "CHAIN: " + chain_id
        searchPDB = pdb_id + "_" + chain_id + ".pdb"  #This is the pdb file that is parsered by dssp
        pdbLocation = duncan_structure_path + searchPDB #This is the location of the pdb file
        natural_sequences =  af.get_natural_sequences_duncan(file) #Gets a list with the natural sequences  
        ancestor = natural_sequences[0] #Gets the "ancestral sequence"
        ancestor_length = len(ancestor) #Grab the length of the ancestral sequence
        #print "Ancestor Length: " + str(ancestor_length)
        counter = 0
        gaps = 0
        gap_locations = []

        #Counts the gaps within the ancestral sequence in the alignment. We must take the gaps out before counting the amino acids
        while (counter < ancestor_length):
            acid = ancestor[counter]
            if acid == '-': 
                gaps = gaps + 1
                gap_locations.append(counter) #This is an array that tracks which residues have gaps
            counter = counter + 1

        #This section takes out the ancestral gaps from all the aligned sequences and then writes the files
        #to the results_PDB_ID_CHAIN_ID.csv 
        natural_pdb_file_title = "results_natural_" + pdb_id + "_" + chain_id + ".csv"
        designed_pdb_file_title = "results_" + pdb_id + "_" + chain_id + "_" + method + ".csv" 
        natural_out_sequences = open(natural_pdb_file_title,"w") 

        #This takes all of the gaps out of the alignment
        for sequence in natural_sequences: #For each sequence in the the file with the sequences 
            new_sequence = '' 
            index = 0
            while (index < ancestor_length): #Take out all of the gaps within the aligned proteins
                if index in gap_locations:
                    index = index + 1
                    continue
                else:
                    new_sequence = new_sequence + sequence[index]
                    index = index + 1
            natural_out_sequences.write(new_sequence + "\n")
        natural_out_sequences.close()

        #list of files that contain L vs RSA data
        data = []

        #list of sequences with missing RSA values
        bad_list = []

        #Get all the RSA values using DSSP
        seq_data = af.get_values(pdb_id, chain_id)
        RSAValues = []	
        RSA = seq_data[1]
        #print "RSA Length from dssp script: " + str(len(RSA))
        index = 0

        #Open the files you just created with the gapless sequences and read in the sequences. 
        fp_natural = open(natural_pdb_file_title) 
        natural_data_arr = fp_natural.readlines()
        fp_natural.close()

        fp_designed = open(duncan_designed_sequence_path + designed_pdb_file_title) 
        designed_data_arr = fp_designed.readlines()
        fp_designed.close()

        #filepointer for the partial result
        fpW_natural = open("results_array_natural_" + pdb_id + "_" + chain_id + ".csv","w")
        fpW_designed = open("results_array_" + pdb_id + "_" + chain_id + "_" + method + ".csv","w")

        #write the RSA values:
        fpW_natural.write(af.dump_csv_line(RSA))
        fpW_designed.write(af.dump_csv_line(RSA))

        if(len(RSA) != len(natural_data_arr[0].strip())):
            print "Error !!!!!!!!!"
            #print data_arr[0].strip()
            bad_list.append(natural_pdb_file_title)

        if(len(RSA) != len(natural_data_arr[0].strip())):
            print "Error !!!!!!!!!"
            print "Length of rsa values {:0}".format(len(RSA))
            print "Length of sequence {:0}".format(len(natural_data_arr[0].strip()))

        #Re-format sequences
        for line in natural_data_arr:
            fpW_natural.write(af.dump_csv_line(line.strip()))      
        fpW_natural.close()

        for line in designed_data_arr:
            fpW_designed.write(af.dump_csv_line(line.strip()))      
        fpW_designed.close()

        #search string to use
        file = "results_array_natural_" + pdb_id + "_" + chain_id + ".csv" 
        print file 

        debug = 0

        #grab the file
        fp = open(file,"r")
        natural_data = fp.readlines()
        fp.close()

        #search string to use
        file = "results_array_" + pdb_id + "_" + chain_id +"_" + method + ".csv" 
        print file 

        debug = 0

        #find all csv files that match the search string

        #grab the file
        fp = open(file,"r")
        designed_data = fp.readlines()
        fp.close()

        #if file == "results_array_2GV5_D_0.3.csv":
        #    debug = 1

        #split the dataset by RSA
        natural_size = len(natural_data)
        for i in range(0,natural_size):
            natural_data[i] = re.split(",",natural_data[i].strip())

        if(debug):
            print len(natural_data[i])
        natural_data = np.array(natural_data, dtype = object)

        #Make arrays full of all the sequences from the files
        designed_size = len(designed_data)
        for i in range(0,designed_size):
            designed_data[i] = re.split(",",designed_data[i].strip())
            
#		if(debug):
#        	print len(designed_data[i])
        designed_data = np.array(designed_data, dtype = object)
        #print len(designed_data)

        #snag the RSA values
        RSA = natural_data[0]

        #get a list of char amino acid codes
        AA = resdict.values()
        #print "Dimensions of array {}".format(natural_data.ndim)

        n,m = natural_data.shape
        n_d, m_d = designed_data.shape 
        #print "Natural Protein Length", m
        #print "Designed Protein Length", m_d
        designed_protein_lengths.append(m_d)
        #Open .dat files that you will write the results to.
        natural_fpW = open("align_natural_data_array_" + pdb_id + "_" + chain_id + ".dat","w")
        natural_fpW.write(af.dump_csv_line2(rows))
        designed_fpW = open("align_data_array_" + pdb_id + "_" + chain_id + "_" + method + ".dat","w")
        designed_fpW.write(af.dump_csv_line2(rows))

        counter = 0
        natural_aaSum = 0
        designed_aaSum = 0
        aaCutOff = 0 #This is cut-off. If a site has less than this percent of amino acids left after taking out the gaps do not use that site in the analysis. 			  
        for i in range(0,m): #This blocks calculates the site frequency data 
            natural_aaCount = [] 
            designed_aaCount = []
            #covert back to list so we can use the "count" method
            try:
                natural_aaList = list(natural_data[1:,i])
                designed_aaList = list(designed_data[1:,i])
            except IndexError:
                print "--------------"
                print "Fudge"
                print i
                print data.ndim
                print data.shape
                print "--------------"
                quit()

            for aa in AA:
                try:
                    natural_aaCount.append(natural_aaList.count(aa))
                    designed_aaCount.append(designed_aaList.count(aa))
                except ValueError,IndexError:
                    print "-----------------------"
                    print "Blargtastic!"
                    print aaList
                    print "-----------------------"
                    quit()
            #Print the formatted frequency count data to a file
            natural_outStr = "\"" + str(counter) + "\" "+ str(counter) + " " + str(RSA[i]) + " " + af.dump_csv_line2(natural_aaCount)  
            designed_outStr = "\"" + str(counter) + "\" "+ str(counter) + " " + str(RSA[i]) + " " + af.dump_csv_line2(designed_aaCount)
            counter = counter + 1
            natural_fpW.write(natural_outStr)
            designed_fpW.write(designed_outStr)

        print "Processed: %s\n" % file
        natural_fpW.close()
        designed_fpW.close() 
print designed_protein_lengths
mean_designed_protein_length = np.mean(designed_protein_lengths)
print "Duncan Mean_Designed_Length: ", mean_designed_protein_length
