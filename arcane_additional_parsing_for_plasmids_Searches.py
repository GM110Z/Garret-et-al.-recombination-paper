#!/usr/bin/env python
# coding: utf-8

# In[8]:


import pandas as pd
import csv
import sys


list_of_accession = []
#input the file with a list of IDs of the plasmid genomes files
with open (sys.argv[1], 'r', encoding='utf-8-sig') as csvfile:
    efetchin=csv.reader(csvfile, delimiter = ',')
    for row in efetchin:
        list_of_accession.append(str(row[0]))
#input the efetch_output.txt file got froma arcane.py
df = pd.read_csv(sys.argv[2], sep="\t", low_memory=False)
df.columns=['ID','Source', 'Nucleotide Accession', 'Start', 'Stop', 'Strand', 'Protein', 'Protein Name', 'Organism', ' Strain', 'Assembly']
# Get names of indexes for which rows have to be dropped
indexNames = df[ df['Source'] == 'INSDC'].index
# Delete these row indexes from dataFrame
df.drop(indexNames , inplace=True)
# Get names of indexes for which rows have to be dropped
indexNames2 = df[ df['Source'] == 'PAT'].index
# Delete these row indexes from dataFrame
df.drop(indexNames2 , inplace=True)
df.to_csv('hits found anywhere.tsv','\t')
df3=df.sort_values(by = "Protein", axis=0, ascending=True, inplace=False)
df3.to_csv('neigh_analysis_input.tsv', '\t', header=False, columns = ['Assembly', 'Protein'])
#write file for stephen for plasmid hits
df2=df[df['Nucleotide Accession'].isin(list_of_accession)]
df2.to_csv('hits_on_plasmids.tsv','\t')

