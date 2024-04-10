
import os
import argparse
from Bio import Entrez

def check_args(args=None):
    parser = argparse.ArgumentParser(description='Retrieve protein sequences for a virus from NCBI')
    parser.add_argument('-v', '--virus',
                        help='Virus name or identifier',
                        required=True)
    parser.add_argument('-o', '--output',
                        help='Output directory for protein files',
                        required=True)
    return parser.parse_args(args)

def fetch_proteins(virus, output_directory):
    Entrez.email = "your@email.com"  # Put your email here
    search_term = f"{virus} [ORGANISM]"
    handle = Entrez.esearch(db="protein", term=search_term, retmax=10)
    search_result = Entrez.read(handle)
    ids = search_result["IdList"]
    for protein_id in ids:
        handle = Entrez.efetch(db="protein", id=protein_id, rettype="fasta", retmode="text")
        protein_data = handle.read()
        handle.close()
        filename = f"{virus}_{protein_id}.fasta"
        filepath = os.path.join(output_directory, filename)
        with open(filepath, 'w') as file:
            file.write(protein_data)
        print(f"Protein sequence saved to: {filepath}")

def VirusProteins():
    args = check_args()
    virus = args.virus
    output_directory = args.output

    if not os.path.exists(output_directory):
        os.makedirs(output_directory)

    fetch_proteins(virus, output_directory)

