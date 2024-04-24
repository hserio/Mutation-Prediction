**Overview**

Viruses mutate extremely frequently, their mutation rate being higher than any other organism. Their ability to quickly evolve has led to efficient infection of their host, evading antiviral drugs. This characteristic also makes it quite difficult to estimate/quantify how accurate a prediction is to what is, since the spread of a virus, as seen with the spread of Sars-CoV-2. There is an urgent demand for accurate predictive tools for viral mutations to get ahead of any virus and stop the spread before mutation.

**Functionality**

The Viral Mutation Prediction Tool is the starting point for a comprehensive software package to extract, preprocessing, training and mutation prediction from additional viral genomic data, not just from Sars-Cov-2. This toolkit is the starting point, this will have the Python pipeline build for extraction. This will prompt the user in the command line with a query for NCBI in the form of user-input for information such as the organism wanted, how many sequences, and if a specified date range is wanted. It will also prompt the user if amino acid frequencies should be preserved and if so, how many significant figures. From this, it will retrieve the amount of sequences specified from NCBI, generate .csv files, and create boxplots to view amino acid frequencies. 

**Features**

Data Retrieval: Easily downloads from the National Center for Biotechnology Information (NCBI) to create .csv files and visualise the difference in amino acid frequencies through boxplots 

**Installation**

To get started, clone the repository and install the required dependencies:

git clone https://github.com/hserio/Mutation-Prediction.git
cd mutation-prediction
pip install seaborn re pandas scipy
