#!/usr/bin/python
import sys, os, math, string, re, gzip, urllib, shutil, Bio, subprocess
import cStringIO 
from Bio.PDB import *
from Bio.PDB.DSSP import * 
from scipy.stats import pearsonr as pearson
import random as rnd
import numpy as np
from Bio.PDB.Polypeptide import *
from pylab import *
import matplotlib.pyplot as plt
from Bio.PDB import PDBParser

#LAST UPDATED: Apr 15, 2017
#Description: This is a series of functions that are used in the Protein Design Project Analysis

EulerGamma = 0.57721566490153286060 #Euler Gamma Constant

#This is a dictionary that relates the three letter amino acid abbreviation with its one letter abbreviation
resdict = { 'ALA': 'A', 'CYS': 'C', 'ASP': 'D', 'GLU': 'E', 'PHE': 'F', \
            'GLY': 'G', 'HIS': 'H', 'ILE': 'I', 'LYS': 'K', 'LEU': 'L', \
            'MET': 'M', 'ASN': 'N', 'PRO': 'P', 'GLN': 'Q', 'ARG': 'R', \
            'SER': 'S', 'THR': 'T', 'VAL': 'V', 'TRP': 'W', 'TYR': 'Y' }

#This is a dictionary that has the amino acid for the key and the max solvent accessibility for this amino acid
#THIS HAS BEEN UPDATED. I AM USING THE NEW THEORETICAL NUMBERS FROM THE 2013 AUSTIN, STEPHANIE, MATT, WILKE PAPER. 
residue_max_acc = {'A': 129.0, 'R': 274.0, 'N': 195.0, 'D': 193.0, \
                   'C': 158.0, 'Q': 223.0, 'E': 224.0, 'G': 104.0,  \
                   'H': 209.0, 'I': 197.0, 'L': 201.0, 'K': 237.0, \
                   'M': 218.0, 'F': 239.0, 'P': 159.0, 'S': 151.0, \
                   'T': 172.0, 'W': 282.0, 'Y': 263.0, 'V': 174.0}


#residue_max_acc = {'A': 113.0, 'R': 241.0, 'N': 158.0, 'D': 151.0, \
#		   'C': 140.0, 'Q': 189.0, 'E': 183.0, 'G': 85.0,  \
#		   'H': 194.0, 'I': 182.0, 'L': 180.0, 'K': 211.0, \
#		   'M': 204.0, 'F': 218.0, 'P': 143.0, 'S': 122.0, \
#		   'T': 146.0, 'W': 259.0, 'Y': 229.0, 'V': 160.0}


#These are the absolute PATH to our original sequence dataset files (structures, designed sequences, and natural sequences)
duncan_designed_sequence_path = "/Users/qian/Desktop/evol_sim_vs_rosetta/sequences/designed_sequences/"
duncan_natural_sequence_path = "/Users/qian/Desktop/evol_sim_vs_rosetta/sequences/aligned_sequences/"
duncan_structure_path = "/Users/qian/Desktop/evol_sim_vs_rosetta/structures/"

#These three functions take an array and then create a formatted string that is used for printing data to a file.
def dump_csv_line(line):
    new_line = ""
    size = len(line)
    for x in range(0,size):
        new_line = new_line + str(line[x])
        if(x != size-1):
            new_line = new_line + ","

    new_line += "\n"
    return new_line
 
def dump_csv_line2(line):
    new_line = ""
    size = len(line)
    for x in range(0,size):
        new_line = new_line + str(line[x])
        if(x != size-1):
            new_line = new_line + " "

    new_line += "\n"
    return new_line

def dump_csv_line3(line):
    new_line = ""
    size = len(line)
    for x in range(0,size):
        new_line = new_line + str(line[x])
        if(x != size-1):
            new_line = new_line + "\t"

    new_line += "\n"
    return new_line

#This is a function that takes a designed alignment from the UCSF dataset, and returns an array with sequences 
#with just the residues from the sequences that have a mapped residue in the corresponding natural sequences for that protein.
def get_cut_designed_sequences(file):
    all_sequences = []  
    natural_sequences = []
    designed_sequences = []
    fileparts = re.split("_",file) #These few lines stores the protein identity, pdb_id and chain information
    identity = fileparts[0]
    pdb_id = fileparts[1]
    chain_id = fileparts[2]
    index_filename = pdb_id + '_A.indices' #A string representing the filename with the index mapping information
    file_data = open(noah_designed_sequence_path + file, "r") #This is were the sequences are stored.
    seq_data = file_data.readlines()
    file_data.close()
    index_file = open(noah_structure_path + index_filename)  #Opens the file with the residue mapping information
    index_data = index_file.readlines()
    index_file.close()
    indices = []
    mod_sequences = []
    for line in index_data: #Creates an array with all of the indices of of which residues are important in the designed alignment
        parts = re.split("\t", line)
        index = int(parts[2])
        indices.append(index)
    
    #This block of code basically just removes all of the headers from the array of sequence information.
    #It creates all_sequences, which is a list of all the sequences from the natural alignment file.
    string  = ''  
    finished_sequence = ''
    for sequence in seq_data: 
        if sequence[0] == '>': #If it is a header append the last sequence that was processed
            if(string != ''):
                finished_sequence = string
                all_sequences.append(finished_sequence) #Appends the sequence into the array full of aligned sequences
                string = '' #Empty the string that represents the sequence (you don't want to append the old with the new)
        else:       
            string = string + sequence.rstrip("\n") #Strip the new line that is at the end of each sequence in alignment
    all_sequences.append(string)
    num_sequences = len(all_sequences)

    for seq in all_sequences: #For each sequence in the list representing all the sequences in the alignment
        new_seq = '' #Make an empty string representing that sequence
        for i in indices: 
            new_seq = new_seq + seq[i-1] #Only append the residues from that sequence that are in the index file to the new string
        mod_sequences.append(new_seq) #Append that new sequence with just the important residues to the new list of sequences
    return mod_sequences   #Return the list of modified sequences 


#This is a function that takes a natural alignment from the UCSF dataset, and returns an array with sequences 
#with just the residues from the sequences that have a mapped residue in the corresponding designed sequences for that protein.
def get_cut_natural_sequences(sequence_identity, all_sequences):
    natural_sequences = []
    designed_sequences = []
    file = sequence_identity + '.align.80' #Get the name of the file with the natural alignment
    pdb_id = identity_dict[sequence_identity] #Get the pdb_id corresponding to this alignemt
    chain_id = 'A' #All the pdbs are chain A
    index_filename = pdb_id + '_A.indices' #A string representing the filename with the index mapping information
    file_data = open(noah_natural_sequence_path + file, "r") #This is were the sequences are stored.
    seq_data = file_data.readlines()
    file_data.close()
    index_file = open(noah_structure_path + index_filename)  #Opens the file with the residue mapping information
    index_data = index_file.readlines()
    index_file.close()
    indices = []
    mod_sequences = []
    
    for line in index_data: #Creates an array with all of the indices of of which residues are important in the natural alignment
        parts = re.split("\t", line)
        index = int(parts[0])
        indices.append(index)

    for seq in all_sequences: #For each sequence in the list representing all the sequences in the alignment
        new_seq = '' #Make an empty string representing that sequence
        for i in indices:
            new_seq = new_seq + seq[i] #Only append the residues from that sequence that are in the index file to the new string
        mod_sequences.append(new_seq)  #Append that new sequence with just the important residues to the new list of sequences
    return mod_sequences  #Return the list of modified sequences  

#This is functions gets the RSA values that map to the residues in the UCSF natural alignments
def get_cut_RSA_values(pdb_id, RSA_values, RSA_dict):
    all_sequences = []  
    natural_sequences = []
    designed_sequences = []
    index_filename = pdb_id + '_A.indices'   #A string representing the filename with the index mapping information
    index_file = open(noah_structure_path + index_filename) #Opens the file with the residue mapping information
    index_data = index_file.readlines()
    index_file.close()
    indices = []
    mod_sequences = []
    for line in index_data:  #Creates an array with all of the indices of of which residues are important in the natural alignment
        parts = re.split("\t", line)
        #index = int(parts[3])
        index = re.split(" ",parts[3])
        index = int(index[0])
        indices.append(index)

    mod_RSA_values = []
    for i in indices:
        mod_RSA_values.append(RSA_dict[i]) #Only append the RSA values that are mapped according to the index file
    #print "Length of cut RSA values: " + str(len(mod_RSA_values))
    return mod_RSA_values  #Return RSA values
    
#This is a function that takes a list of sequencs and a new list. This list is a list of all of the sequences in the orginal that are at most 25 sequence identity when compared to every other sequence in the list
def get_dissimilar_sequences(natural_seqs):
    dissimilar_natural_seqs = []
    sequence_is_similar = False
    for seq in natural_seqs:
        seq_copy = list(natural_seqs) #This copies a new list that can be compared
        seq_copy.remove(seq) #Must take out the sequence else it will compare to itself (will be always True)
        sequence_is_similar = compare_sequences(seq, seq_copy) #Determines whether it is similar
        if sequence_is_similar == False:
            dissimilar_natural_seqs.append(seq)
    return dissimilar_natural_seqs

#This function takes a sequence and compares it to all the sequences in a list. It returns a boolean that is equal to True or False. If True then the sequnce is at least 75 percent similar to at least one sequence in the least.
def compare_sequences(seq, natural_seqs):
    cutoff = int(len(seq)*0.25) #The cutoff for sequence identity divergence
    sequence_similar = True
    for line in natural_seqs: #This compares the seq to all of the sequences in the list AA by AA
        dissimilar_count = 0
        for i in xrange(0, len(min([seq, line], key = len))): #Must take the length of the smaller seq being compared
            #print i
            natural_acid = seq[i] 
            other_acid = line[i]
            if (natural_acid != other_acid): #Checks if identity is the same at that site. 
                dissimilar_count = dissimilar_count + 1
                if(dissimilar_count>cutoff):
                    sequence_similar = False
                    return sequence_similar
                else:
                    continue
        if(dissimilar_count > cutoff):
            sequence_similar = False
    return sequence_similar 

#This functions takes a file with a bunch of natural sequences that have been aligned and then returns a list of the sequences. 
def get_natural_sequences_duncan(file):
    all_sequences = []  
    natural_sequences = []
    designed_sequences = []
    fileparts = re.split("_",file)
    pdb_id = fileparts[0]
    chain_id = fileparts[1]
    
    file_data = open(duncan_natural_sequence_path + file, "r") #This is were the sequences are stored.
    seq_data = file_data.readlines()
    file_data.close()
    
    #This block of code basically just removes all of the headers from the array of sequence information.
    #It creates all_sequences, which is a list of all the sequences from the natural alignment file.
    string  = ''
    finished_sequence = ''
    for sequence in seq_data:
        if sequence[0] == '>': #If it is a header append the last sequence that was processed
            if(string != ''):
                finished_sequence = string
                all_sequences.append(finished_sequence) #Appends the sequence into the array full of aligned sequences
                string = '' #Empty the string that represents the sequence (you don't want to append the old with the new)
        else:       
            string = string + sequence.rstrip("\n") #Strip the new line that is at the end of each sequence in alignment
    all_sequences.append(string)
    num_sequences = len(all_sequences)
    return all_sequences   

#This functions takes a file with a bunch of natural sequences that have been aligned and then returns a list of the sequences. 
def get_natural_sequences_noah(file):
    all_sequences = []  
    natural_sequences = []
    designed_sequences = []
    fileparts = re.split("_",file)
    #pdb_id = fileparts[0]
    #chain_id = fileparts[1]
    
    file_data = open(noah_natural_sequence_path + file, "r") #This is were the sequences are stored.
    seq_data = file_data.readlines()
    file_data.close()
    
    #This block of code basically just removes all of the headers from the array of sequence information.
    #It creates all_sequences, which is a list of all the sequences from the natural alignment file.
    string  = ''
    finished_sequence = ''
    for sequence in seq_data:
        if sequence[0] == '>': #If it is a header append the last sequence that was processed
            if(string != ''):
                finished_sequence = string
                all_sequences.append(finished_sequence) #Appends the sequence into the array full of aligned sequences
                string = '' #Empty the string that represents the sequence (you don't want to append the old with the new)
        else:       
            string = string + sequence.rstrip("\n") #Strip the new line that is at the end of each sequence in alignment
    all_sequences.append(string)
    num_sequences = len(all_sequences)
    return all_sequences   

#This takes a file with a bunch of aligned natural and designed sequences that are combined and returns two lists with them seperated. 
def split_merged_sequences(file):
    all_sequences = []  
    natural_sequences = []
    designed_sequences = []
    fileparts = re.split("_",file)
    pdb_id = fileparts[0]
    chain_id = fileparts[1]
    method = fileparts[2]
    file_data = open(file, "r")
    seq_data = file_data.readlines()
    file_data.close()
    string  = ''
    finished_sequence = ''
    for sequence in seq_data:
        #print sequence
        if sequence[0] == '>':
            if(string != ''):
                finished_sequence = string
                all_sequences.append(finished_sequence) #Appends the sequence into the array full of aligned sequences
                string = ''
        else:       
            string = string + sequence.rstrip("\n")
    all_sequences.append(string)
    num_sequences = len(all_sequences)
    num_designed_seq = 500
    num_natural_seq = num_sequences - num_designed_seq
    for seq_index in xrange(0,num_natural_seq):
        natural_sequences.append(all_sequences[seq_index])
    for seq_index in xrange(num_natural_seq,num_sequences):
        designed_sequences.append(all_sequences[seq_index])
    return natural_sequences, designed_sequences, method

#This takes a list of natural sequences and splits them in half. It returns two samples with the split sequences.
def split_natural_sequences(natural_sequences):
    seq_sample_1 = natural_sequences
    L = len(natural_sequences) #Gets the length of the alignment
    L_sample = int(L/2)  #Creates an int representing the value of half the length of the alignment
    seq_sample_2 = rnd.sample(natural_sequences, L_sample) #Gets a random list of sequences from the alignment. This is list L_sample long.
    for seq in seq_sample_2: #For each sequence in the randomly generated sample 
        seq_sample_1.remove(seq) #Remove it from the original alignment list
    return seq_sample_1,seq_sample_2 #Return both samples

#Gets the number of times that a count value is seen in the array and returns a list with the count, countOfCounts data 
def get_AA_counts(distribution_data):
    count_data = []
    for x in distribution_data:
        if x != 0:
            num_appearances = distribution_data.count(x)
            count_data.append((x,num_appearances))
    count_data = list(set(count_data))
    return count_data

#This function takes a list containing RSA and site AA count info and returns a list of lists that contain all of the count data for every site
def get_transformed_data(list_of_sequences):
    transformed_distribution = []
    for seq in list_of_sequences: #Removes any newlines that happen to be appended to the new of the 
        sequences = seq.strip() 
    list_of_sequences.pop(0)

    new_data = []
    for line in list_of_sequences:
        element = line.split()
        element.pop(0)  #Removes the first three elements in the line since they do not contain frequency data
        element.pop(0)
        element.pop(0)
        new_data.append(element)
    for data in new_data: #For each site  (Each site should have 20 numbers, each number representing how often each AA was seen at that site
        new_elements = []
        for count in data:  #For each number of counts representing each AA    
            new_count = int(count) #Turns the string representing the count data into an actual number
            new_elements.append(new_count) #Appends the number to the list containing the count data for that site
        transformed_distribution.append(new_elements) #Appends the count data for each site to a giant list
    return transformed_distribution

#This function takes a list containing RSA and site AA count info and returns a list of lists that contain all of the count data for every site
def get_transformed_data_KL(list_of_sequences):
    transformed_distribution = []
    for seq in list_of_sequences: #Removes any newlines that happen to be appended to the new of the 
        sequences = seq.strip() 
    list_of_sequences.pop(0)

    new_data = []
    for line in list_of_sequences:
        element = line.split()
        element.pop(0)  #Removes the first three elements in the line since they do not contain frequency data
        element.pop(0)
        element.pop(0)
        new_data.append(element)
    for data in new_data: #For each site  (Each site should have 20 numbers, each number representing how often each AA was seen at that site
        new_elements = []
        for count in data:  #For each number of counts representing each AA    
            new_count = float(count) + 0.05 #Turns the string representing the count data into an actual number
            new_elements.append(new_count) #Appends the number to the list containing the count data for that site
        transformed_distribution.append(new_elements) #Appends the count data for each site to a giant list
    return transformed_distribution

#This function takes the amino acid count data for the designed and corresponding natural sequences and returns a list of lists with all the AA frequency data at each site. 
def get_AA_distribution(proteins): #Opens the file with the frequency results (the .dat file)
    protein_distribution = []
    transformed_protein_distribution = []
    num_AA = 0
    input = open(proteins, "r")
    protein_data = input.readlines()
    input.close()
    protein_distribution = get_transformed_data(protein_data) #Gets a list of list with all of the AA counts for each site 
    for site in protein_distribution:
        new_site = site
        num_AA = sum(new_site)
        aa_probs = []
        for count in new_site:   
            if count == 0: #Turns all of the raw counts into frequencies 
                prob = 0.0
                #prob = float(1)/float(num_AA + 20)
            else:
                prob = float(count)/float(num_AA)
                #prob = float(count + 1)/float(num_AA + 20)
            aa_probs.append(prob)
        transformed_protein_distribution.append(aa_probs)
    return transformed_protein_distribution
 
 #This function takes the amino acid count data for the designed and corresponding natural sequences and returns a list of lists with all the AA frequency data at each site. 
def get_AA_distribution_KL(proteins): #Opens the file with the frequency results (the .dat file)
    protein_distribution = []
    transformed_protein_distribution = []
    num_AA = 0
    input = open(proteins, "r")
    protein_data = input.readlines()
    input.close()
    protein_distribution = get_transformed_data_KL(protein_data) #Gets a list of list with all of the AA counts for each site 
    for site in protein_distribution:
        new_site = site
        num_AA = sum(new_site)
        aa_probs = []
        for count in new_site:   
            if count == 0: #Turns all of the raw counts into frequencies 
                prob = 0.0
                #prob = float(1)/float(num_AA + 20)
            else:
                prob = float(count)/float(num_AA)
                #prob = float(count + 1)/float(num_AA + 20)
            aa_probs.append(prob)
        transformed_protein_distribution.append(aa_probs)
    return transformed_protein_distribution

#Returns a list of KL_Divergence values for two distributions. The inputs are two list of lists representing the two distributions.
def get_Kullback_Leibler(real_proteins, created_proteins):
    #print real_proteins
    #print created_proteins
    KL_Values = []
    KL_Number = 0
       
    real_array = array(real_proteins) #Array with all of the AA frequency for each site in the natural alignment
    created_array = array(created_proteins) #The array same for the designed
    created_num_residues, created_num_AA = created_array.shape
    num_residues,num_AA = real_array.shape
    
    for i in xrange(0, num_residues): #This block of code calculates the KL-Divergence for each site in the protein
        real_values = real_proteins[i]
        created_values = created_proteins[i]
        KL_Number = 0
        for j in xrange(0,20):
            value = log(float(real_values[j])/float(created_values[j]))
            value = value*float(real_values[j])
            KL_Number = KL_Number + value
             
            '''if (created_values[j] == 0.0 and real_values[j]== 0.0):
                #print "In the first cateogory", str(created_values[j]), str(real_values[j])
                continue
            elif (created_values[j] == 0.0):
                #print "In the second cateogory", str(created_values[j]), str(real_values[j])
                continue
            elif (real_values[j] == 0.0):
            	value = 0.0
                #print "In the third cateogory", str(created_values[j]), str(real_values[j])
            else:            
                value = log(float(real_values[j])/float(created_values[j]))
                value = value*float(real_values[j])
                KL_Number = KL_Number + value
                #print "In the last cateogory", str(created_values[j]), str(real_values[j]), value
            '''
            
            
        KL_Values.append(KL_Number)
    return KL_Values #Returns a list with the KL-Divergence for each site
 
#This function calculates the mean of list of data points and returns the mean
def calculate_mean(data_points):
    mean_of_data = []
    new_data_points = []
    for point in data_points:
        new_point = float(point)
        new_data_points.append(new_point)
    
    sum_of_data = sum(new_data_points)
    num_elements = len(data_points)
    mean_of_data = float(sum_of_data)/float(num_elements)
    return mean_of_data
 
#This function takes a file (the .dat file generated from the calculate_distribution python scripts) and extracts the RSA values
def get_RSA_Values(protein_file):
    input = open(protein_file, "r")
    protein_data = input.readlines() #Reads all the data from the file into an array to process 
    input.close()
    RSA = []
    for seq in protein_data:
        sequence = seq.strip() #Strips the newlines from each line in the array of file data lines
    protein_data.pop(0)
    new_data = []
    for line in protein_data:
        element = line.split()
        test = element.pop(0)
        test = element.pop(0)
        RSA_value = element.pop(0) #Grabs the RSA values from the file dat for each site
        RSA.append(RSA_value) # Appends the RSA values
        new_data.append(element)
    return RSA #Returns a list of RSA values

#The inputs to this functions are three lists with the RSA, entropy, and KL site data for a protein. This sorts the entropy, and 
#KL data by RSA into three categories: buried, partially buried, and surface. It then returns the data. 
def get_position_dependent_data(RSA_data, entropy_data, KL_data):
    num_sites = len(RSA_data)
    i = 0
    buried_KL_values = []
    buried_entropy_values = []
    buried_RSA_values = []

    intermediate_KL_values = []
    intermediate_entropy_values  = []
    intermediate_RSA_values = []

    surface_KL_values = []
    surface_entropy_values = []
    surface_RSA_values = []
    
    while i < num_sites:
        if (float(RSA_data[i])<0.05): # The RSA is less than 0.05 it is a buried site
            buried_KL_values.append(float(KL_data[i]))
            buried_entropy_values.append(float(entropy_data[i]))
            buried_RSA_values.append(float(RSA_data[i]))
            i = i + 1
        elif (0.05<=float(RSA_data[i])<=0.25):  #These are the partially buried sites
            intermediate_KL_values.append(float(KL_data[i]))
            intermediate_entropy_values.append(float(entropy_data[i]))
            intermediate_RSA_values.append(float(RSA_data[i]))
            i = i + 1
        elif (float(RSA_data[i]) > 0.25): #These are surface sites
            surface_KL_values.append(float(KL_data[i]))
            surface_entropy_values.append(float(entropy_data[i]))
            surface_RSA_values.append(float(RSA_data[i]))
            i = i + 1
        else: #This is not possible 
            print "Problem in get_core_data"
            print "RSA value is: " + str(RSA_data[i])
            i = i + 1
    return buried_entropy_values, buried_KL_values, intermediate_entropy_values, intermediate_KL_values, surface_entropy_values, surface_KL_values 

#This takes a file with site amino acid count data and then calculates the native entropy for sites. 
def get_native_entropy(protein_file): #Give it the result file (the .dat file)
    probs = get_AA_distribution(protein_file) #Gets the AA frequencies
    entropy_values = []
    effective_values = []
    entropy_number = 0
    probs_array = array(probs)
    num_residues,num_AA = probs_array.shape
    
    for i in xrange(0, num_residues): #Calculated the native entropy
        probs_values = probs_array[i,:]
        prob_sum = sum(probs_values)
        entropy_number = 0
        for j in xrange(0,20):
            if (probs_values[j] == 0.0):
                value = 0.0
            else:
                value = (float(probs_values[j])*log(float(probs_values[j])))
            entropy_number = entropy_number + value
            effective_amino_acids = exp(-entropy_number)
        entropy_values.append(-entropy_number)
        effective_values.append(effective_amino_acids)
    return effective_values

#This takes the file with the count data and calculates the entropy using a different entropy calculator. 
def get_entropy(protein_file):
    protein_distribution = []
    transformed_protein_distribution = []
    num_AA = 0
    input = open(protein_file, "r")
    protein_data = input.readlines()
    input.close()
    probs = get_transformed_data(protein_data) #Get just the raw AA appearance counts for each site
    entropy_values = []
    entropy_number = 0
    #probs_array = array(probs)
    new_probs = []
    #new_probs_array = array(probs) #Make it into an array
    for site in probs:
        new_site_array = []
        for count in site:
            new_count = count + 1
            new_site_array.append(new_count)
        new_probs.append(new_site_array)
    #num_residues,num_AA = new_probs_array.shape
    count = 0
    for site in new_probs:
        frequencies = site
        entropy_number = Entropy_H_G(frequencies)
        entropy_values.append(entropy_number)
        if (entropy_number >2.998):
            print count
            print frequencies
            print "Entropy Number: " + str(entropy_number)  
        count = count + 1
    return entropy_values

#This takes a list of values that are string values and creates an array
def make_array(list_of_values):
    new_value_list = []
    for value in list_of_values:
        new_value = float(value)
        new_value_list.append(new_value)
    value_array = array(new_value_list)
    return value_array

#Opens a file with the slope and intercept data for all the proteins for the designed proteins. 
#It returns two dictionaries that map the slope and intercept data to the PDB names. 
def get_slope_intercept(file): #Give it the raw_mle_line file as input
    slope_list = []
    intercept_list = []
    pdb_list = []
    input = open(file, "r")
    file_data = input.readlines() #Read the file into an array
    for line in file_data:
        data = re.split(",", line) #Split each line by a comma (Make an array the PDB, intercept and slope info as separate elements)
        pdb_line = data[0]  #These lines get the PDB in each line
        PDB = re.split("_", pdb_line) 
        PDB = PDB[4]
        pdb_list.append(PDB) #Append the PDB name to a list of PDB names
        intercept = data[1] #Store the intercept
        slope = data[2] #Store the slope
        intercept_list.append(intercept) #Append the intercept to a list of intercepts
        slope_list.append(slope) #Append the slopes to a list of slopes
    intercept_dict = dict(zip(pdb_list, intercept_list)) #Create a dictionary with the PDB as a key and its intercept as its value
    slope_dict = dict(zip(pdb_list, slope_list)) #Create a dictionary with the PDB as a key and its slope as its value
    return intercept_dict, slope_dict #Return the slope and intercept dictionaries

#Opens a file with the slope and intercept data for all the proteins for the natural proteins. 
#It returns two dictionaries that map the slope and intercept data to the PDB names. 
def get_slope_intercept_natural(file): #Give it the raw_mle_line file as input
    slope_list = []
    intercept_list = []
    pdb_list = []
    input = open(file, "r")
    file_data = input.readlines() #Read the file into an array
    for line in file_data:
        data = re.split(",", line) #Split each line by a comma (Make an array the PDB, intercept and slope info as separate elements)
        pdb_line = data[0] #These lines get the PDB in each line
        PDB = re.split("_", pdb_line)
        PDB = PDB[5]
        pdb_list.append(PDB) #Append the PDB name to a list of PDB names
        intercept = data[1] #Store the intercept
        slope = data[2] #Store the slope
        intercept_list.append(intercept) #Append the intercept to a list of intercepts
        slope_list.append(slope) #Append the slopes to a list of slopes
    intercept_dict = dict(zip(pdb_list, intercept_list)) #Create a dictionary with the PDB as a key and its intercept as its value
    slope_dict = dict(zip(pdb_list, slope_list))  #Create a dictionary with the PDB as a key and its slope as its value
    return intercept_dict, slope_dict #Return the slope and intercept dictionaries

#This calculates the RSA values for a PDB using DSSP. It returns a list of the Amino Acids and a list of their RSA values. 
def get_values(pdb_id, chain_id):
    searchPDB = pdb_id + "_" + chain_id + ".pdb"  #This is the pdb file that is parsered by dssp
    pdbLocation = duncan_structure_path + searchPDB #This is the location of the pdb file that dssp will be parsing.
    structure = PDBParser().get_structure(pdb_id, pdbLocation)#Creates a structure object 
    model = structure[0]
    method_file = 'new_pdb_methodfile.txt'
    Bio.PDB.Dice.extract(structure, chain_id, 0, 10000, method_file) #Isolate just the chain from the PDB, make method PDB with just the chain
    outputFile = pdb_id + "_" + chain_id + ".txt" 
    processString = 'mkdssp' + ' -i ' + '"' + method_file + '"'  + ' -o pdbOutput.txt '  #Run DSSP to and get a file with the Solvent Accessibility 
    process = subprocess.Popen(processString, shell = True, stdout = subprocess.PIPE)
    process.wait() # Wait until dssp is done processing the file and calculating the Solvent Acessiblility  values
    input = open("pdbOutput.txt" , 'r')    
    fileContents = input.readlines()	
    string = fileContents[25]
    SAValue1 = string[13]
    SAList = [] #This is is the list which will store the SA values for each site
    AAList = [] #This is the list which will store the amino acid values for each site
    index = 0
    NoRSA = 0
    for line in fileContents:
        if index<28: #This skips the first few lines that do not have the SA value data
            index = index + 1
            continue
        else:  #Goes through each line with has the SA for each amino acid in order
            string = line #This stores the current line in the string "string"
            SAValue = string[35:39] #This stores the SA value for the current 
            AA = string[13] #This stores what the amino acid type is at the current  position
            number = int(SAValue)  #This turns the string fir the SA into an int type
            if AA !=( '!' or '*'): #This takes out the missing gaps that dssp might put in
                max_acc = residue_max_acc[AA] #This uses the dictionary to find the max SA for the amino acid at the current position (site)
                SAList.append(number/max_acc) #This divides the SA value for that position by the amx SA value position for tha amino acid. This normalizes the values and gives us the Relative Solvent Accessability (RSA) value. We this appends this value to the list of RSA values
                AAList.append(AA) #This appends the amino acid to the list
            else:
           	NoRSA = NoRSA + 1 #Counts the number of residues that did not have SA values
	    index = index + 1
    input.close() #Close the file with the dssp output file
    #os.remove('pdbOutput.txt') #Deletes the dssp output file 
    return (AAList, SAList) #Return the RSA values and the SAList

def getResNumber(residue):
    id = residue.get_id()
    assert(id[0]==' ') # just make sure we're working with a properly numbered PDB file
    assert(id[2]==' ')
    return id[1]

'''
def get_wcn_values(pdb_id, chain_id):
    searchPDB = pdb_id + "_" + chain_id + ".pdb"  
    pdbLocation = duncan_structure_path + searchPDB 
    structure = PDBParser().get_structure(pdb_id, pdbLocation) 
    model = structure[0]
    for chain in model:
        wcn_values = []
        iwcn_values = []
        for r1 in chain:
            i = 0
            for r2 in chain:
                if getResNumber(r2) != getResNumber(r1):
                    i += 1/(float(r1['CA']-r2['CA'])**2)
            wcn_values.append(i)
            iwcn_values.append(1/i)
            
    return (wcn_values, iwcn_values)

def get_cn13_values(pdb_id, chain_id):
    searchPDB = pdb_id + "_" + chain_id + ".pdb"  
    pdbLocation = duncan_structure_path + searchPDB 
    structure = PDBParser().get_structure(pdb_id, pdbLocation) 
    model = structure[0]
    for chain in model:
        cn_values = []
        icn_values = []
        for r1 in chain:
            i = 0
            for r2 in chain:
                if getResNumber(r2) != getResNumber(r1):
                	distance = float(r1['CA']-r2['CA'])
                	if distance < 13:
                		i = i + 1
            cn_values.append(i)
            icn_values.append(float(1.0/i))
    return (cn_values,icn_values) 
'''

#This calculates the RSA values for a PDB using DSSP. It returns a list of the Amino Acids, a list of RSA values, and a dictionary mapping residue positions to RSA. 
def get_noah_RSA_values(pdb_id, chain_id):
    searchPDB = pdb_id + "_" + chain_id + ".pdb"  #This is the pdb file that is parsered by dssp
    pdbLocation = noah_structure_path + searchPDB #This is the location of the pdb file that dssp will be parsing. 
    structure = PDBParser().get_structure(pdb_id, pdbLocation)#Creates a structure object 
    model = structure[0]
    method_file = 'new_pdb_methodfile.txt'
    Bio.PDB.Dice.extract(structure, chain_id, 0, 10000, method_file) #Isolate just the chain from the PDB, make method PDB with just the chain
    outputFile = pdb_id + "_" + chain_id + ".txt" 
    processString = 'dssp' + ' -i ' + '"' + method_file +'"'  + ' -o pdbOutput.txt ' #Run DSSP to and get a file with the Solvent Accessibility 
    process = subprocess.Popen(processString, shell = True, stdout = subprocess.PIPE)
    process.wait() # Wait until dssp is done processing the file and calculating the Solvent Acessiblility  values
    input = open("pdbOutput.txt" , 'r')    
    fileContents = input.readlines()	
    string = fileContents[25]
    SAValue1 = string[13]
    SAList = [] #This is is the list which will store the SA values for each site
    AAList = [] #This is the list which will store the amino acid values for each site
    residue_list = [] #This is a list that stores the residue postion numbers
    index = 0
    NoRSA = 0
    for line in fileContents:
        if index<25: #This skips the first few lines that do not have the SA value data
            index = index + 1
            continue
        else:  #Goes through each line with has the SA for each amino acid in order
            string = line #This stores the current line in the string "string"
            SAValue = string[35:39] #This stores the SA value for the current 
            res_pos = string[6:10]
            AA = string[13] #This stores what the amino acid type is at the current  position
            number = int(SAValue)  #This turns the string for the SA into an int type
            if AA !=( '!' or '*'): #This takes out the missing gaps that dssp might put in
                max_acc = residue_max_acc[AA] #This uses the dictionary to find the max SA for the amino acid at the current position (site)
                SAList.append(number/max_acc) #This divides the SA value for that position by the amx SA value position for tha amino acid. This normalizes the values and gives us the Relative Solvent Accessability (RSA) value. We this appends this value to the list of RSA values
                residue = int(res_pos)
                AAList.append(AA) #This appends the amino acid to the list
                residue_list.append(residue)
            else:
           	NoRSA = NoRSA + 1 #Counts the number of residues that did not have SA values
	    index = index + 1
    RSA_dict = dict(zip(residue_list,SAList))
    input.close() #Close the file with the dssp output file
    #os.remove('pdbOutput.txt') #Deletes the dssp output file 
    return (AAList, SAList, RSA_dict) #Return the RSA values and the SAList


#This takes a file with a bunch of generated quantities (mean RSA, mean KL value, mean entropy, correlation between RSA and Entropy) 
#for a protein and returns lists with this data
def get_mean_designed_data(file_of_data):
    designed_mean_RSA_values = []
    designed_mean_entropy_values = []
    designed_cor_entropy_RSA_values = []
    designed_cor_entropy_iWCN_values = []
    designed_mean_KL_values = []
    designed_data = []
    protein_file = open(file_of_data, "r")
    designed_protein_data = protein_file.readlines()
    protein_file.close()
    header = designed_protein_data.pop(0)

    for line in designed_protein_data:  
        data = re.split("\t", line)
        designed_data.append(data)
    
    for data in designed_data:
        designed_mean_RSA_values.append(data[2])
        designed_mean_entropy_values.append(data[3])
        designed_cor_entropy_RSA_values.append(data[4])
        designed_mean_KL_values.append(data[5])
        designed_cor_entropy_iWCN_values.append(data[6])
    return designed_mean_RSA_values, designed_mean_entropy_values, designed_cor_entropy_RSA_values, designed_mean_KL_values, designed_cor_entropy_iWCN_values

def get_entropy_corr_data(file_of_data):
    natural_rosetta_corr_values = []
    natural_evolved_corr_values = []
    rosetta_evolved_corr_values = []
    designed_data = []
    protein_file = open(file_of_data, "r")
    designed_protein_data = protein_file.readlines()
    protein_file.close()
    header = designed_protein_data.pop(0)

    for line in designed_protein_data:  
        data = re.split("\t", line)
        designed_data.append(data)
    
    for data in designed_data:
        natural_rosetta_corr_values.append(data[2])
        natural_evolved_corr_values.append(data[3])
        rosetta_evolved_corr_values.append(data[4])
    return natural_rosetta_corr_values, natural_evolved_corr_values, rosetta_evolved_corr_values

#This takes a file with a bunch of generated quantities (RSA, KL value, Entropy) for a protein and returns lists with this data
def get_designed_graph_data(generated_data_file):
  designed_file = open(generated_data_file)
  data = designed_file.readlines()
  designed_data = []
  designed_file.close()
  for line in data:
      file_data = re.split("\t",line)
      designed_data.append(file_data)

  designed_entropy = designed_data[0]
  designed_entropy.pop(0)
  designed_entropy_array = make_array(designed_entropy)

  designed_RSA = designed_data[1]
  designed_RSA.pop(0)
  designed_RSA_array = make_array(designed_RSA)

  designed_KL = designed_data[2]
  designed_KL.pop(0)
  designed_KL_array = make_array(designed_KL)   
  return designed_entropy,designed_RSA,designed_KL    

#This takes a file with a bunch of generated quantities (RSA, KL value, Entropy) for a protein and returns lists with this data
def get_mean_ordered_designed_data(file_of_data):
    designed_mean_RSA_values = []
    designed_mean_entropy_values = []
    designed_cor_entropy_RSA_values = []
    designed_mean_KL_values = []
    designed_intercept_values = []
    designed_slope_values = []
    designed_data = []
    protein_file = open(file_of_data, "r")
    designed_protein_data = protein_file.readlines()
    protein_file.close()
    header = designed_protein_data.pop(0)

    for line in designed_protein_data:  
        data = re.split("\t", line)
        designed_data.append(data)
    
    for data in designed_data:
        designed_mean_RSA_values.append(data[2])
        designed_mean_entropy_values.append(data[3])
        designed_cor_entropy_RSA_values.append(data[4])
        designed_mean_KL_values.append(data[5])
        designed_intercept_values.append(data[6])
        designed_slope_values.append(data[7])

    return designed_mean_RSA_values, designed_mean_entropy_values, designed_cor_entropy_RSA_values, designed_mean_KL_values, designed_intercept_values, designed_slope_values

#This is a G() calculator that is used to calculate entropy
def Gi(n): #Function Written By Claus Wilke
    '''Helper function needed for entropy estimation.
    Defined by Grassberger 2003. http:/arxiv.org/abs/physics/0307138
    '''  
    if n == 0:
        return 0
    if n == 1:
        return -EulerGamma - math.log(2)
    if n == 2:
        return 2-EulerGamma - math.log(2)
    if (n % 2) == 1:
        return Gi( n-1 )
    return Gi(n-2) + 2./(n-1)

#This is used in the entropy calculator
def Entropy_H_G(list_of_frequency_counts): #Function Written By Claus Wilke
    '''Best entropy estimator according to Grassberget 2003,
    http:/arxiv.org/abs/physics/0307138
    '''  
    z = list_of_frequency_counts
    N = sum(z) # total number of observations
    return math.log(N) - (1./N)* sum([n*Gi(n) for n in z])
