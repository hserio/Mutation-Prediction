**Overview**

A virus's ability to generate genetic diversity during replication is a key factor in its effective infection of host cell and ability to grow more drug-resistant strains which ultimately results in its wide spread. The rate of mutation along with other environmental and biological factors also make it quite difficult to predict its mutation patterns as seen with the spread of SARS-CoV-2. Therefore, this pipeline here is part of a greater ML tool that predicts sequence mutations of viral proteomes based on AAC analysis. This pipeline is a command line tool that extracts NCBI records of viral proteins based on user's input, create visuals to summarize AAC distribution, and generate necessary files for further analysis by the ML tool.

**Functionality**

Fullfilling the first step of the ML prediction tool, data extraction, this pipeline operates based on minimal user input to generate neccessary files and visuals before processing the data with the ML models. Upon running this python script, the user is prompted to input the virus they are intersted in, the maximum number of records to be extracted, and an optional data range, all of which are used to refine the NCBI search query. It will also prompt the user to input the desired numer of significant figures to preserved the amino acid frequencies later. Then, it will retrieve the number of sequences specified from NCBI, generate .csv files, and create plots to summarize the distribution of the amino acid composition (AAC). 

**Features**

Data Retrieval: Easily downloads from the National Center for Biotechnology Information (NCBI) to create .csv files and visualise the difference in amino acid frequencies through boxplots 

**Software Requirements**
Linux, Python 3.12

**Dependencies**
Biopython, seaborn, matplotlib, scipy, numpy, sys, csv, re, math, datetime, os

**Installation**

1) To get started, clone the repository and install the required dependencies:

         git clone https://github.com/hserio/Mutation-Prediction.git

         cd Mutation-Prediction


**Usage**

After cloning the repository and installing any required dependencies, follow the steps below:

1) Running the Data Retrieval program (Respond to the script prompts to refine the NCBI search)

         python NCBI_Retrieval_Script.py
   
Prompts will include user email, NCBI Search Term (Ex: SARS COVID19), number of protein sequences to be extracted, a desired date range and how many signifigant figures are to be preserved for amino acid frequencies

**Example**

         Enter Email: ****@email.com
         Enter NCBI Search Term: SARS COVID19
         How many protein sequences would you like to extract? 1000
         Default fate range for protein sequence extraction is from 01/01/2000 - Current.
         Would you like to extract protein sequences from a specified date range?
         Enter (Y/N): N
         How many significant figures would you like to preserve for amino acid frequencies? 2

**Outputs**

frequency_distribution.png 
- A figure displaying the frequency distributions for each amino acid residue from the extracted sequences. 

protein.csv
- A csv file with the amino acid frequencies per extracted protein sequence.

protein_Search.txt
- A text file containing each extracted protein sequence in fasta format.

amino_acids.csv
- A csv file ranking the average amino acid frequency from all extracted sequences.

composition.png
- A figure displaying the total counts of each amino acid residue for all extracted sequences.

**Additional Information**

Another script, Strain_Mutation_Visualization.py, has been created as a template for future mutation analysis. This script outputs a bar plot (all_frequencies.png) that shows the differences in amino acid composition between a specific strain vs the average amino acid composition for the same virus. Right now it is designed to compare human papillomavirus 16 and 18 to the average amino acid frequencies and can be tested using the code below:

         python Strain_Mutation_Visualization.py

When prompted, enter these terms:
   
         Enter Email: ****@email.com
         Enter NCBI Search Term: Human papillomavirus
         How many protein sequences would you like to extract? 1000
         Default fate range for protein sequence extraction is from 01/01/2000 - Current.
         Would you like to extract protein sequences from a specified date range?
         Enter (Y/N): N
         How many significant figures would you like to preserve for amino acid frequencies? 2

We also have a 'sklearn' file in this repository containing information on retraining the machine learning models in order to aid individuals continuing the tools development.
