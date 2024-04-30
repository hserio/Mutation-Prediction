**Overview**

Viruses mutate extremely frequently, their mutation rate being higher than any other organism. Their ability to quickly evolve has led to efficient infection of their host, evading antiviral drugs. This characteristic also makes it quite difficult to estimate/quantify how accurate a prediction is to what is, since the spread of a virus, as seen with the spread of Sars-CoV-2. There is an urgent demand for accurate predictive tools for viral mutations to get ahead of any virus and stop the spread before mutation.

**Functionality**

The Viral Mutation Prediction Tool is the starting point for a comprehensive software package to extract, preprocessing, training and mutation prediction from additional viral genomic data, not just from Sars-Cov-2. This toolkit is the starting point, this will have the Python pipeline built for extraction. This will prompt the user in the command line with a query for NCBI in the form of user input for information such as the organism wanted, how many sequences, and if a specified date range is wanted. It will also prompt the user if amino acid frequencies should be preserved and if so, how many significant figures. From this, it will retrieve the number of sequences specified from NCBI, generate .csv files, and create boxplots to view amino acid frequencies. 

**Features**

Data Retrieval: Easily downloads from the National Center for Biotechnology Information (NCBI) to create .csv files and visualise the difference in amino acid frequencies through boxplots 

**Installation**

1) To get started, clone the repository and install the required dependencies:

         git clone https://github.com/hserio/Mutation-Prediction.git

         cd Mutation-Prediction

         pip install seaborn re pandas scipy

**Usage**

After cloning the repository and installing any required dependencies, follow the steps below:

1) Running the Data Retrieval program (Respond to the script prompts to refine the desired search)

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
