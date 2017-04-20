#!usr/local/bin/python
import sys, os, math, string, re, gzip, urllib, shutil, Bio, subprocess
import cStringIO 
import numpy as np
import analysis_functions as af
import random 

#Date Last Modified:  Apr 15, 2017
#Description: This file takes the aligned natural sequences, splits them into two samples. It then calculates the frequency data for each sample and prints the data to two different files. 

resdict = { 'ALA': 'A', 'CYS': 'C', 'ASP': 'D', 'GLU': 'E', 'PHE': 'F',                   
            'GLY': 'G', 'HIS': 'H', 'ILE': 'I', 'LYS': 'K', 'LEU': 'L',                   
            'MET': 'M', 'ASN': 'N', 'PRO': 'P', 'GLN': 'Q', 'ARG': 'R',                   
            'SER': 'S', 'THR': 'T', 'VAL': 'V', 'TRP': 'W', 'TYR': 'Y' }          

rows = ["\"site\"", "\"RSA\"", "\"aa1\"", "\"aa2\"", "\"aa3\"", "\"aa4\"", "\"aa5\"", "\"aa6\"", "\"aa7\"", "\"aa8\"", "\"aa9\"", "\"aa10\"", "\"aa11\"", "\"aa12\"", "\"aa13\"", "\"aa14\"", "\"aa15\"", "\"aa16\"", "\"aa17\"", "\"aa18\"", "\"aa19\"", "\"aa20\""]   
PDBS = ["1ci0A", "1g58B", "1hujA", "1ibsA", "1jlwA", "1kzlA", "1m3uA", "1mozA", "1pv1A", "1qmvA", "1riiA", "1v9sB", "1w7wB", "1x1oB", "1ypiA", "1znnA", "2a84A", "2bcgY", "2br9A", "2cjmC", "2esfA", "2fliA", "2gv5D", "1b4tA", "1efvB", "1gv3A", "1hurA", "1ky2A", "1okcA", "1r6mA", "1xtdA", "1ysbA", "1zwkA", "2aiuA", "2cfeA", "2cnvA", "2eu8A", "2g0nB"]

#These are the absolute PATH to our original sequence dataset files (structures, designed sequences, and natural sequences)
duncan_designed_sequence_path = "/Users/qian/Desktop/evol_sim_vs_rosetta/sequences/designed_sequences/"
duncan_natural_sequence_path = "/Users/qian/Desktop/evol_sim_vs_rosetta/sequences/aligned_sequences/"
duncan_structure_path = "/Users/qian/Desktop/evol_sim_vs_rosetta/structures/"

for protein in PDBS:
    pdb_id = protein[0:4].upper()
    chain_id = protein[4] 
    file = pdb_id + "_" + chain_id + "_Aligned_Sequences.fasta"
    fileparts = re.split("_",file)
 
    print "Processsing file: " + file	
    print "PDB: " + pdb_id
    print "CHAIN: " + chain_id
    searchPDB = pdb_id + "_" + chain_id + ".pdb"  #This is the pdb file that is parsered by dssp
    pdbLocation = duncan_structure_path + searchPDB #This is the location of the pdb file that dssp will be parsing. 
    natural_sequences = af.get_natural_sequences_duncan(file) #Gets the natural sequences  
    ancestor = natural_sequences[0]
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
    natural_out_sequences = open(natural_pdb_file_title,"w") 

    #This block of code takes out all of the natural sequence alignment
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
    print len(RSA)
    index = 0

    #Open the files you just created with the gapless sequences
    fp_natural = open(natural_pdb_file_title) 
    natural_data_arr = fp_natural.readlines()
    fp_natural.close()

    #filepointer for the partial result
    fpW_natural = open("results_array_natural_" + pdb_id + "_" + chain_id + ".csv","w")

    #write the RSA values:
    fpW_natural.write(af.dump_csv_line(RSA))

    [natural_sample1, natural_sample2] = af.split_natural_sequences(natural_data_arr) #Splits the sequences into two samples
    #filepointer for the partial result
    fpW_natural_sample1 = open("results_array_natural_sample1_" + pdb_id + "_" + chain_id + ".csv","w")
    fpW_natural_sample2 = open("results_array_natural_sample2_" + pdb_id + "_" + chain_id + ".csv","w")


    #write the RSA values:
    fpW_natural_sample1.write(af.dump_csv_line(RSA))
    fpW_natural_sample2.write(af.dump_csv_line(RSA))


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

    for line in natural_sample1:
        fpW_natural_sample1.write(af.dump_csv_line(line.strip()))      
    fpW_natural_sample1.close()

    for line in natural_sample2:
        fpW_natural_sample2.write(af.dump_csv_line(line.strip()))      
    fpW_natural_sample2.close()


    #search string to use
    file = "results_array_natural_" + pdb_id + "_" + chain_id + ".csv" 
    print file 
    debug = 0

    #find all csv files that match the search string
    #grab the file
    fp = open(file,"r")
    natural_data = fp.readlines()
    fp.close()

    debug = 0

    #find all csv files that match the search string

    #search string to use
    file = "results_array_natural_sample1_" + pdb_id + "_" + chain_id + ".csv" 
    print file 
    debug = 0

    #find all csv files that match the search string

    #grab the file
    fp = open(file,"r")
    natural_sample1_data = fp.readlines()
    fp.close()

    #search string to use
    file = "results_array_natural_sample2_" + pdb_id + "_" + chain_id + ".csv" 
    print file 

    #find all csv files that match the search string

    #grab the file
    fp = open(file,"r")
    natural_sample2_data = fp.readlines()
    fp.close()

    #split the dataset by RSA
    natural_size = len(natural_data)
    for i in range(0,natural_size):
        natural_data[i] = re.split(",",natural_data[i].strip())

    if(debug):
        print len(natural_data[i])
    natural_data = np.array(natural_data, dtype = object)

    natural_sample1_size = len(natural_sample1_data)
    for i in range(0,natural_sample1_size):
        natural_sample1_data[i] = re.split(",",natural_sample1_data[i].strip())

    if(debug):
        print len(natural_sample1_data[i])
    natural_sample1_data = np.array(natural_sample1_data, dtype = object)
    natural_sample2_size = len(natural_sample2_data)

    for i in range(0,natural_sample2_size):
        natural_sample2_data[i] = re.split(",",natural_sample2_data[i].strip())

    if(debug):
        print len(natural_sample2_data[i])
    natural_sample2_data = np.array(natural_sample2_data, dtype = object)

    #snag the RSA values
    RSA = natural_data[0]

    #get a list of char amino acid codes
    AA = resdict.values()
    #print "Dimensions of array {}".format(natural_data.ndim)

    n,m = natural_data.shape
    n_sample1, m_sample = natural_sample1_data.shape
    #Open .dat files that you will write the results to.
    natural_sample1_fpW = open("align_natural_sample1_data_array_" + pdb_id + "_" + chain_id + ".dat","w")
    natural_sample1_fpW.write(af.dump_csv_line2(rows))
    natural_sample2_fpW = open("align_natural_sample2_data_array_" + pdb_id + "_" + chain_id + ".dat","w")
    natural_sample2_fpW.write(af.dump_csv_line2(rows))

    counter = 0
    natural_sample1_aaSum = 0
    natural_sample2_aaSum = 0
		  
    for i in range(0,m): #Calculates the site frequency data
        natural_aaCount = []
        designed_aaCount = []
        natural_sample1_aaCount = []
        natural_sample2_aaCount = []
        #covert back to list so we can use the "count" method
        try:
            natural_aaList = list(natural_data[1:,i])
            natural_sample1_aaList = list(natural_sample1_data[1:,i])
            natural_sample2_aaList = list(natural_sample2_data[1:,i])
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
                natural_sample1_aaCount.append(natural_sample1_aaList.count(aa))
                natural_sample2_aaCount.append(natural_sample2_aaList.count(aa))
            except ValueError,IndexError:
                print "-----------------------"
                print "Blargtastic!"
                print aaList
                print "-----------------------"
                quit()
        #Prints formatted frequency data to a file 
        natural_sample1_outStr = "\"" + str(counter) + "\" "+ str(counter) + " " + str(RSA[i]) + " " + af.dump_csv_line2(natural_sample1_aaCount)
        natural_sample2_outStr = "\"" + str(counter) + "\" "+ str(counter) + " " + str(RSA[i]) + " " + af.dump_csv_line2(natural_sample2_aaCount)
        counter = counter + 1
        natural_sample1_fpW.write(natural_sample1_outStr)
        natural_sample2_fpW.write(natural_sample2_outStr)

    print "Processed: %s\n" % file
    natural_sample1_fpW.close()
    natural_sample2_fpW.close()

