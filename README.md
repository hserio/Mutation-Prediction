**Overview**

Viruses mutate extremely frequently, their mutation rate being higher than any other organism. Their ability to quickly evolve has led to efficient infection of their host, evading antiviral drugs. This characteristic also makes it quite difficult to estimate mutation prediction, as seen with the spread of SARS-CoV-2. There is an urgent demand for accurate predictive tools for viral mutations to get ahead of any virus and stop the spread before mutation.

**Functionality**

The viral mutation prediction tool is being developed that will first extract protein sequences from NCBI, preprocess the data, and use machine learning models to obtain mutation predictions. This toolkit is the starting point, this will have the Python pipeline built for extraction. This will prompt the user in the command line with a query for NCBI in the form of user input for information such as the organism wanted, how many sequences should be extracted, and if a specified date range is desired. It will also prompt the user how many significant figures should be preserved for amino acid frequencies. From this, it will retrieve the number of sequences specified from NCBI, generate .csv files, and create plots to view amino acid frequencies. 

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
