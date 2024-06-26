# this is a modified version of "NCBI_Retrieval_Script.py" that shows the mutation frequency differences between the average frequency of the viruses vs a specific variant
# this code only has comparisons for HPV with HPV16 and HPV18

# Import necessary packages
import sys                                              # module contains methods and variables for modifying Python's runtime environment
import csv                                              # module implements classes to read and write tabular data in csv format
import numpy as np                                      # module for array creation and working with numerical data
from Bio import SeqIO                                   # module functioning as an interface to input and output fasta format files
from Bio import Entrez                                  # module to search NCBI for protein sequences with user-specified parameters
from Bio.SeqUtils.ProtParam import ProteinAnalysis      # module for analysis of protein sequences
import matplotlib.pyplot as plt                         # module to plot numerical data
import seaborn as sns                                   # module for visualization and exploratory analysis, based on matplotlib
import scipy.stats                                      # module to conduct statistical tests
import re                                               # module for regular expression matching operations
import pandas as pd                                     # module for building dataframes
from math import log10, floor                           # module used to output user-specified number of significant figures in output data
from datetime import date                               # module to pull current dates
import os                                               # module to navigate command line and avoid hard coded paths

# Create input file
input_file = open('proteinSearch.txt', 'w')             # open and write sequences to proteins text file 

# Get the current working directory
current_directory = os.getcwd()
output_file = current_directory + 'amino_acids.csv'  # CHANGE to your own directory
output_file2 = current_directory + 'prot_freq.csv' # change your directory path as well

# Protein Sequence Retrieval from NCBI based on search term
Entrez.email = input("Enter email: ")                   # user prompted to enter email (tell NCBI who you are to access sequences)

protTerm = input("Enter NCBI search term: ")            # user prompted to enter protein sequence ID
numSeqs = input("How many protein sequences would you like to extract? ") # user prompted to enter # seqs to retrieve

print("Default date range for protein sequence extraction is from 01/01/2000 - Current.") 
print("Would you like to extract protein sequences from a specified date range?") 
startDate = "2000/01/01"                                # start date for default date range
today = date.today()                                    # current date pulled using datetime module for end date for default date range
endDate = today.strftime("%d/%m/%Y")                    # end date for default date range
dateY_N = ""                                            # initialize user-specified option to opt with default date range or enter their own
while dateY_N != "N" and dateY_N != "n" and dateY_N != "Y" and dateY_N != "y": # while loop to loop through user-entered options until "Y" (yes) or "N" (no) entered
    dateY_N = input("Enter (Y/N): ")                    # user prompted to enter "Y" (yes) or "N" (no) in regard to setting their own date range
    if dateY_N != "N" and dateY_N != "n" and dateY_N != "Y" and dateY_N != "y": # if the user enters a value other than "Y" or "N" they will be asked to enter one of those options
        print("Please choose yes (Y) or no (N).")
if dateY_N == "Y" or "y" and dateY_N != "N" and dateY_N != "n": # if the user enters "Y", then they will be prompted for further inputs
    startDate == input("Using format YYYY/MM/DD, enter start date: ") # user prompted to enter start date
    endDate == input("Using format YYYY/MM/DD, enter end date: ") # user prompted to enter end date
    
    
searchResultHandle = Entrez.esearch(db = "protein", term = protTerm, retmax = numSeqs, idtype = "protein", datetype = "pdat", mindate = startDate, maxdate = endDate) # entrez search handle with user-specified options
searchResult = Entrez.read(searchResultHandle)          # read in handle with parameters set to user inputs 
ids = searchResult["IdList"]                            # list of IDs created from protein sequences retrieved

handle = Entrez.efetch(db="protein", id=ids, rettype="fasta", retmode="text") # protein sequences retrieved by IDs in fasta format
record = handle.read()                                  # record created reading in the handle containing the fasta format protein sequences           

input_file = open('proteinSearch.txt', 'w')             # input file created in write format to write protein sequences in fasta format
input_file.write(record.rstrip('\n'))                   # each fasta format protein sequence is stripped of the new line character
input_file.close()                                      # close the input file containing the fasta format protein sequences


input_file = open('proteinSearch.txt', 'r')             # reopen input file but in read format
protein_ids = []                                        # make a list of protein IDs

# Define function to extract amino acid sequences from a FASTA file
def extract_amino_acids(input_file):                    
    amino_acids = []                                    # create a list of amino acids frequencies, where each item in the list is a dictitionary of the amino acid frequencies for a protein retrieved from NCBI
    for record in SeqIO.parse(input_file, "fasta"):     # parse through each record in the input file using the SeqIO parser to detect fasta formatted sequences
        sequence = str(record.seq)                      # sequence from each record obtained
        protein = ProteinAnalysis(sequence)             # protein analysis set to be performed on each individual sequence
        aa_percentages = protein.get_amino_acids_percent() # percentages of amino acids obtained from protein sequence
        protein_ids.append(record.id)                   # each protein ID appended to the list of protein IDs
        amino_acids.append(aa_percentages)              # each amino acid percentage appended to list of amino acid percentages
    return amino_acids                                  # list of amino acid percentages returned

# Define function to write amino acid frequencies and min/max percentages to a CSV file

def write_csv(output_file, data):                     
    with open(output_file, "w", newline="") as csvfile: # open new output file as a CSV file
        writer = csv.writer(csvfile)                    # use writer function to write outputs to the CSV file
        writer.writerow(["Rank", "Amino Acid", "Frequency"]) # write the header row to the CSV file
        rank = 1                                        # initialize rank to 1
        for aa, freq in data.items():                   # iterate through each amino acid with corresponding frequency as an item in data
            writer.writerow([rank, aa, freq])           # write the rank and frequency corresponding to each amino acid to the output file
            rank += 1                                   # add 1 to rank each time an amino acid is added

# Extract amino acid sequences from input file
amino_acids = extract_amino_acids(input_file)           # amino acids functions as a list of dictionaries of the amino acids and corresponding frequencies for each protein (each protein has its own dictionary of amino acid frequencies)
input_file.close()                                      # close the input file
 
# Calculate total frequencies of each amino acid
total_aa_freqs = {}                                     # create a dictionary of amino acid frequencies
for aa in amino_acids:                                  # for loop iterates through each set of amino acid frequencies (one dictionary per protein) in the list of amino acid frequencies for the individual proteins
    for k, v in aa.items():                             # within each dictionary, k functions as the amino acid, while v functions as its corresponding frequency                             
        if k in total_aa_freqs:                         # for each amino acid, its corresponding frequency is added to total_aa_freqs
            total_aa_freqs[k] += v
        else:                                           # if the amino acid is not present in total_aa_freqs, its initial frequency is set
            total_aa_freqs[k] = v
            
sig = int(input("How many significant figures would you like to preserve for amino acid frequencies? ")) # user is prompted to enter number of desired significant figures in statistical outputs
for k, v in total_aa_freqs.items():                     # for loop iterates through each amino acid and corresponding frequency in total_aa_freqs
    total_aa_freqs[k] = round((v/int(numSeqs)),sig-int(floor(log10(abs(v/int(numSeqs)))))-1) # each frequency is rounded to the number of significant figures specified by the user
    
total_aa_freqs_sorted = dict(sorted(total_aa_freqs.items(), key = lambda item: item[1], reverse = True)) # total_aa_freqs_sorted is a sorted version of total_aa_freqs, which are the total amino acid frequencies across all proteins, sorted from most to least common
          
# Calculate Min/Max %'s of Frequencies
min_percent = min(total_aa_freqs.values())             # calculate minimum frequency
max_percent = max(total_aa_freqs.values())             # calculate maximum frequency

# Write amino acid frequencies to output file
write_csv(output_file, total_aa_freqs_sorted)

aa_outfile = open("test.out.csv", "w", newline="")     # undesired output printed in test.out.csv--ignore this file output
writer = csv.writer(aa_outfile, delimiter=",")         # writer function used to write output to both test.out.csv (ignore) and protein.csv (contains desired output)

protein_outfile = open("protein.csv","w")              # open a new outfile "protein.csv" to write amino acid frequencies for each individual protein (tabular representation of the list of dictionaries for the amino acid frequencies in each individual protein)

# Creates and writes headers
header = ["Protein ID"] + list(amino_acids[0].keys())  # write the header, which contains the protein ID in the first column and each one-letter amino acid abbreviation in subsequent columns
protein_outfile.write(",".join(header)+'\n')           # header formatting
writer.writerow(header)                                # header written to the protein.csv outfile

# Writes amino acid frequencies
for i in range(len(protein_ids)):                      # for loop iterates through each ID in the list of protein IDs
    aa_dict = amino_acids[i]                           # aa_dict is an itemized format of each item in the amino_acids list                      
    aa_freqlist = list(aa_dict.values())               # list creation from values in aa_dict
    row = protein_ids[i] + "," + ",".join([str(x) for x in aa_freqlist]) # list comprehension to write rows so that the protein ID lands in the first column, with the amino acid frequencies added to each subsequent column across each row
    protein_outfile.write(row + "\n")                  # file formatting
    writer.writerow(row)                               # write row to the protein.csv outfile

aa_outfile.close()                                     # close the test.out.csv outfile
protein_outfile.close()                                # close the protein.csv outfile

# Plot average amino acid frequencies using seaborn
sns.set_style("whitegrid")                            
plt.figure(figsize=(10, 6))
sns.barplot(x=list(total_aa_freqs.keys()), y=list(total_aa_freqs.values()), palette="Blues")
plt.title("Amino Acid Frequencies")                    # barplot title
plt.xlabel("Amino Acid")                               # barplot x-axis
plt.ylabel("Frequency")                                # barplot y-axis
plt.text(-0.6, min_percent-0.005, f"Min: {min_percent:.2%}") # displays minimum average amino acid frequency
plt.text(19.2, max_percent-0.005, f"Max: {max_percent:.2%}") # displays maximum average amino acid frequency
plt.savefig("frequency.png")                              # save outputfigure in png format
plt.show()                                             # display plot to screen

HPV_16 = "FVTVVDTTRSTNMSLCAAISTSEPTYKNTNFKEYLRHGEEYDLQFIF" # sequence of HPV16 (gencode accession: WVD73578.1)
HPV_18 = "PVPGQYDATKFKQYSRHVEEYDLQFIFQLCTITLTADVMSYIHSMNSSILEDWNFGVPPPPTTSLVDTYRFVQSVAITCQKDAAPAENKDPYDKLKFWNVDLKEKFSLDLDQ" # sequence of HPV18 (genbank: AWW87290.1)

# function that calculates amino acid frequencies
def calculate_aa_freq(sequence):
    sig = 2 # significant figure value
    all_freqs = {} # initialize dictionary with all frequencies
    total_aa_count = 0 # initialize count of all amino acids
    for aa in sequence: # loop through each amino acid in the sequence
        if aa in all_freqs:
            all_freqs[aa] += 1 # add one to each amino acid count when present
        else:
            all_freqs[aa] = 1 # initialize amino acid if this is the first occurance of it
        total_aa_count += 1
    for k, v in all_freqs.items(): # for loop iterates through each amino acid and corresponding frequency in all_freqs
        all_freqs[k] = round((v/int(total_aa_count)),sig-int(floor(log10(abs(v/int(total_aa_count)))))-1) # each frequency is rounded to the number of significant figures
    return all_freqs

HPV16_freq = calculate_aa_freq(HPV_16) # calculate aa frequency for hpv16 variant
HPV18_freq = calculate_aa_freq(HPV_18) # calculate aa frequency for hpv18 variant

# function that subtracts amino acid frequencies of a specific variant from average amino acid frequencies
def mutation_aa_freq(total_aa_freqs, variant_aa_freqs):
    mut_aa_freqs = {} # initialize variable to store difference in aa frequencies
    for aa, total_freq in total_aa_freqs.items():
        new_freq = float(total_freq) - float(variant_aa_freqs.get(aa, 0)) # subtracts the amino acid frequency of a variant from the average total frequency (if AA is not found, 0 is used)
        mut_aa_freqs[aa] = new_freq # add each new aa and frequency to the graph
    return mut_aa_freqs

HPV_16_subtract = mutation_aa_freq(total_aa_freqs, HPV16_freq) # subtract variant frequencies from average frequencies
HPV_18_subtract = mutation_aa_freq(total_aa_freqs, HPV18_freq) # subtract variant frequencies from average frequencies

# Plot HPV16 differing amino acid frequencies using seaborn
sns.set_style("whitegrid")                            
plt.figure(figsize=(10, 6))
sns.barplot(x=list(HPV_16_subtract.keys()), y=list(HPV_16_subtract.values()), palette="Blues")
plt.title("Amino Acid Frequencies in HPV16")                    # barplot title
plt.xlabel("Amino Acid")                               # barplot x-axis
plt.ylabel("Frequency")                                # barplot y-axis
plt.text(-0.6, min_percent-0.005, f"Min: {min_percent:.2%}") # displays minimum average amino acid frequency
plt.text(19.2, max_percent-0.005, f"Max: {max_percent:.2%}") # displays maximum average amino acid frequency
plt.savefig("frequency_HPV16.png")                              # save outputfigure in png format
plt.show()                                             # display plot to screen

# Plot HPV18 differing amino acid frequencies using seaborn
sns.set_style("whitegrid")                            
plt.figure(figsize=(10, 6))
sns.barplot(x=list(HPV_18_subtract.keys()), y=list(HPV_18_subtract.values()), palette="Blues")
plt.title("Amino Acid Frequencies in HPV18")                    # barplot title
plt.xlabel("Amino Acid")                               # barplot x-axis
plt.ylabel("Frequency")                                # barplot y-axis
plt.text(-0.6, min_percent-0.005, f"Min: {min_percent:.2%}") # displays minimum average amino acid frequency
plt.text(19.2, max_percent-0.005, f"Max: {max_percent:.2%}") # displays maximum average amino acid frequency
plt.savefig("frequency_HPV18.png")                              # save outputfigure in png format
plt.show()                                             # display plot to screen

# Get the amino acids as x-axis labels
x_labels = list(total_aa_freqs.keys())
# Width of each bar in bar graph
bar_width = 0.3
# Set the positions for the groups of bars
positions = range(len(x_labels))
# Plotting
sns.set_style("whitegrid")
plt.figure(figsize=(10, 6))
# Plot HPV18 amino acid frequencies
plt.bar([pos + bar_width for pos in positions], HPV_18_subtract.values(), width=bar_width, color="green", alpha=0.7, label='HPV18 Amino Acid Frequencies')
# Plot HPV16 amino acid frequencies
plt.bar([pos + 2 * bar_width for pos in positions], HPV_16_subtract.values(), width=bar_width, color="blue", alpha=0.7, label='HPV16 Amino Acid Frequencies')
plt.title("Amino Acid Frequencies")
plt.xlabel("Amino Acid")
plt.ylabel("Frequency")
plt.xticks([pos + bar_width for pos in positions], x_labels)  # Set x-axis ticks to amino acid labels
plt.legend()  # Show legend
plt.savefig("all_frequencies.png")  # Save figure
plt.show()
