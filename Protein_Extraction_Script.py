# Importing neccessary modules
# This module provides operating system functionality like reading and
# writing to the file system
import os
# This makes it easy to write user-friendly command-line interfaces 
# It parses command-line arguments and options
import argparse
# This allows for the use of the Biopython library while also allowing for NCBI access
from Bio import Entrez

# This defined function allows for command-line arguements 
def check_args(args=None):
    # Create an ArgumentParser with the following description
    parser = argparse.ArgumentParser(description='Retrieve protein sequences for a virus from NCBI')
    # Adds arguments for virus name and output directory and email 
    parser.add_argument('-v', '--virus',required=True)
    parser.add_argument('-o', '--output',required=True)
    parser.add_argument('-e', '--email',default='your@email.com')
    # Parses the arguments given 
    return parser.parse_args(args)

# This defined function retrieves the protein sequences for a given virus
def fetch_proteins(virus, output_directory, email):
     # Set the email address for NCBI Entrez requests
    Entrez.email = email 
    # Constructs the search term for the virus
    search_term = f"{virus} [ORGANISM]"
    # Performs a search in the NCBI protein database
    handle = Entrez.esearch(db="protein", term=search_term, retmax=10)
    search_result = Entrez.read(handle)
    ids = search_result["IdList"]
    # Looks over the protein IDs and retrieves the protein sequences
    for protein_id in ids:
        handle = Entrez.efetch(db="protein", id=protein_id, rettype="fasta", retmode="text")
        protein_data = handle.read()
        handle.close()
        # Defines the filename and filepath for the FASTA file 
        filename = f"{virus}_{protein_id}.fasta"
        filepath = os.path.join(output_directory, filename)
        # Writes the protein sequence data to the FASTA file 
        with open(filepath, 'w') as file:
            file.write(protein_data)
        # Prints a message indicating where the protein sequence was saved
        print(f"Protein sequence saved to: {filepath}")

# Main function used to retrieve the virus protein sequences requested
def VirusProteins():
    # Parses the command-line arguments 
    args = check_args()
    virus = args.virus
    output_directory = args.output
    email = args.email

    # Creates the output directory if it does not exist
    if not os.path.exists(output_directory):
        os.makedirs(output_directory)

    # Retrieves and saves the protein sequences
    fetch_proteins(virus, output_directory, email)

