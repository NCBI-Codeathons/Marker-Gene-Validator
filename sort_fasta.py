# -*- coding: utf-8 -*-
"""
Created on Fri Apr 10 10:45:44 2020

@author: mcveigh
"""

#
# script to selectively remove sequences from a fasta file based on an accession list
# default accession list is ITS_reject_seqs3.txt
# output file names are hard coded and will overwrite any existing file with the same name.
# accession match works for accession and accession.version. The vesion number is ignored for the match

import pandas as pd
import Bio
import os
import sys

inputfile = sys.argv[1]
#outputfile = sys.argv[2]

#Read in the reject list and find the accession
reject_file_name = (r'ITS_reject_seqs3.txt') 
#reject_file_name = (r'/panfs/pan1.be-md.ncbi.nlm.nih.gov/dnaorg/ITS/ITS_reject_seqs')
rejectlist_df = pd.read_csv(reject_file_name, sep='\t', index_col=None, low_memory=False, header=None, names=["accession", "type", "reason"])
rejectlist = rejectlist_df['accession']
reject_list = set(rejectlist_df['accession'].tolist())
   
#Parse the GenBank file and remove any sequences found on the reject list
from Bio import SeqIO
sequences = [] 
found = []
#missingRNA = []


for seq_record in SeqIO.parse(inputfile, "fasta"):    
    str_id = seq_record.id
    #print(seq_record.id)
    if str_id.find('.') != -1:        
        str_id = str_id[:str_id.find('.')]        
    if str_id not in reject_list:                       
        sequences.append(seq_record)
    else:
        found.append(seq_record)
        print("I found this accession on the reject list and wrote the sequence to found.fsa: ", seq_record.id)  
SeqIO.write(sequences, "stripped.fsa", "fasta")  
SeqIO.write(found, "found.fsa", "fasta")   