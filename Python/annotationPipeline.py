#!/usr/bin/env python3

# import necessary packages
import re
import pandas as pd
import mysql.connector

##########################################
# Functions created for repetitive tasks #
##########################################

# Function to clean up the names
def clean_names(list):
    for item in list:
        item[0] = re.match(r'(\d+_\d+)', item[0]).group(1)
    return list

# Function to return a list of sequence ID that still need predicted products for annotation
def get_missing_ids(df1, df2):
    missing = df1[~df1.qry_id.isin(df2.qry_id)]
    missing_list = []
    for id in missing.qry_id:
        missing_list.append(id)
    return missing_list   

###########################################
# Load and parse FASTA data              #
###########################################

# Open test fasta file, extract the sequence IDs, and create a dataframe
fasta_id = []
for line in open('/Users/nickwinters/desktop/appliedBioinf/Bioinformatics/final/data/prodigal2fasta.nostars.faa'):
    line = line.rstrip()
    if line.startswith('>'):
        m = re.match(r'>(\d+_\d+)', line).group(1)
        fasta_id.append(m)

# Create a dataframe from parsed ids
fasta_col = ['qry_id']
fasta_data = fasta_id
fasta_df = pd.DataFrame(data = fasta_data, columns = fasta_col)

############################################
# Load and parse HMM data                 #
############################################

# Load hmmscan output and parse out relevant data (If the protein has a match with an E value < 1e50 keep that product)
hmm_rows = []

for line in open('/Users/nickwinters/desktop/appliedBioinf/Bioinformatics/final/data/hmmscan_slice.txt'):
    line = line.rstrip()
    values = line.split('\t')
    if float(values[2]) < 1e-50:
        hmm_rows.append(values)
    else:
        continue

# Clean up the names and create a dataframe 
hmm_rows = clean_names(hmm_rows)    
hmm_data = hmm_rows
hmm_columns = ['qry_id', 'product', 'evalue']
hmm_df = pd.DataFrame(hmm_data, columns = hmm_columns)

###########################################
# Load and parse BLAST data               #
###########################################

# Connect to MySQL server
conn = mysql.connector.connect(user = 'python', password = 'Pythonpass29!',
                                host = 'localhost', database = 'annot')

# Create a cursor object to issue SQL commands
cursor = conn.cursor()

# Select statement to look at BLAST results with an E value < 1e-50
cursor.execute("""
                SELECT qry_id, evalue, product
                FROM blast
                WHERE evalue < 1e-50
                """) 

# Fetch the results and create a dataframe
columns = [column[0] for column in cursor.description]
results = cursor.fetchall()
blast_df = pd.DataFrame(results, columns = columns)

# Commit the changes and close
conn.commit()
cursor.close()

# Filter for needed ids
filtered_blast = blast_df[blast_df['qry_id'].isin(get_missing_ids(fasta_df, hmm_df))]

# Concat filtered blast df with hmm df
hmm_blast = pd.concat([hmm_df, filtered_blast], axis = 0)

# Keep the rows with the minimum E value for each qry_id
hmm_blast_filtered = hmm_blast.loc[hmm_blast.groupby('qry_id')['evalue'].idxmin()]

# Drop duplicates and sort by sequence id
hmm_blast_df = hmm_blast_filtered.sort_values('qry_id').drop_duplicates(subset='qry_id', keep='first').sort_index()

###########################################
# Load and parse TMHMM data              #
###########################################

# Load and parse tmhmm output
tmh_rows = []

for line in open('/Users/nickwinters/desktop/appliedBioinf/Bioinformatics/final/data/tmhmm.short.slice.txt'):
    line = line.rstrip()
    values = line.split('\t')
    if re.search(r'(PredHel=)(\d+)', values[1]).group(2) > '0':
        tmh_rows.append(values)

# Clean up the names and create a dataframe
tmh_data = clean_names(tmh_rows)
tmh_columns = ['qry_id', 'product']
tmh_df = pd.DataFrame(tmh_data, columns = tmh_columns)

# Filter for needed ids, drop duplicates, and replace 'PredHel=*' with 'Putative transmembrane protein'
tmh_df = tmh_df[tmh_df['qry_id'].isin(get_missing_ids(fasta_df, hmm_blast_df))]
tmh_df['product'] = 'Putative transmembrane protein'

# Concat filtered tmh df with hmm_blast df
hmm_blast_tmh_df = pd.concat([hmm_blast_df, tmh_df], axis = 0)

########################################################
# label all remaining products as 'hypothetical protien  #
########################################################

# Identify the sequence ids that still need a product name
last_missing = fasta_df[~fasta_df.qry_id.isin(hmm_blast_tmh_df.qry_id)]

# Add these ids to the combined dataframe
output = pd.concat([hmm_blast_tmh_df, last_missing], axis = 0)

# Sort by sequence ids
output['sort_index'] = [int(re.search(r'(\d+_)(\d+)', x).group(2)) for x in output['qry_id']]
output_sorted = output.sort_values('sort_index').reset_index(drop=True)

# Drop unnecessary columns
final_output = output_sorted.drop(labels = ['evalue', 'sort_index'], axis = 1)

# Fill in the unpredicted product names with 'Hypothetical protein'
final_output.fillna('Hypothetical protein', inplace=True)

########################################
# Final output file                    # 
########################################

# Write final output to a tab delimited file called annotation_output.txt
final_output.to_csv('final_output.txt', sep='\t', index=False, header=False)
print("File created: final_output.txt")


