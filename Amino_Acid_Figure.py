# module to create a graph
import matplotlib.pyplot as plt
# module that counts the amino acid frequency
from collections import Counter

# function that plots the average frequency of each amino acid present in a protein sequence
# sequences should be a list of amino acid sequences
def aa_freq_plot(sequences, image_file):
    # initialize counter to store total frequency of each amino acid
    counts = Counter()
    # get the frequency of each residue per sequence and add to the total
    for sequence in sequences:
        count = Counter(sequence)
        # update the total counts
        counts.update(count)
    # get the total number of AA sequences
    total_seq = len(sequences)
    #find the average frequency of each AA
    average_freq = {aa: count / total_seq for aa, count in counts.items()}
    # obtain amino acids and each frequency
    amino_acids = list(average_freq.keys())
    frequency = list(average_freq.values())
    # plot the amino acid composition graph
    plt.figure(figsize=(10,6))
    # specifiy axis
    plt.bar(amino_acids, frequency, color='blue')
    # label axis & title
    plt.xlabel('Residue')
    plt.ylabel('Average Frequency')
    plt.title('Average Frequency of Amino Acids')
    plt.grid(True)
    #save graph as image file
    plt.savefig(image_file)
    plt.close()
