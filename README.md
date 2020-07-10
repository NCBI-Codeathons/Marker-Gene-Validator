# Build a marker gene reference set

## Problem description and objective 

GenBank uses Foosh pipelines to process new submissions of marker genes. An important component of these pipelines is a high quality database of reference gene sequences. Creating such a curated databases involves several manual steps.  The process starts with complicated Entrez queries.  The nucleotide/protein sequences are then manually evaluated using alignments and BLAST. Due to taxonomoic breadth this process is repeated several times to try and catch outliers and bad data. 

The goal of this project is to automate as many of these tasks as possible yet keep the process flexible enough for a curator to jump in at any step and review data. 

As a proof-of-concept, we will produce a Jupyter notebook that a curator can use to build a reference database for the gene CYTB that can be plugged into the GenBank Foosh pipeline. We chose CYTB for this because we can restrict the gene sequences to RefSeq database to keep the initial set relatively small, the gene does not have introns, and is similar to COX1 for which GenBank already has a foosh pipeline. 

##### [1] Develop query to pull nucleotide/proteins sequences of desired gene sequence from Entrez. 
This query can be as limited or as broad as the gene family of interest and its requirements 

`this works but I think it needs to be more restrictive esearch -db gene -query 'metazoa[orgn] AND cytb[gene name]' | efetch -format uid > cytb_geneids.txt`

generates 9880 geneids linking to sequences of cytb*

*We have refined the query*

##### [2] Use NCBI Datasets to download metadata and sequence data for gene IDs from step 1

	-The datasets provide the nucleotide sequence, protein sequence and associated metadata
	-From the datasets the information needed for the next steps can be extracted.
		-after extraction check the metadata table  (check query is providing what you are expecting):
			-product names
			-gene names
	-Extract protein fasta set and check fasta headers.  (may need script in case something changes in the future)
	
		>YP_003667946.1 cytochrome b (mitochondrion) [Carassius gibelio]
	
*We can extract the full data set*

*We can download expanded taxonomy information that includes tax lineages and subtrees*

*Use NCBI datasets and split into taxonomy groupings*
		
	
##### [3] Create histogram to assess length variation. Evaluate sequences at high and low ends. Create script to remove
sequences of certain lengths

*created histogram in Jupyter* 
*Going to use stats to evaluate*

##### [4] Run BLASTall to align the proteins.  From the alignment assess the following:
	
	-what can be programatically removed based on data provided
	-sort by tax groups to evaluate 'outliers'
	-Is there a way to evaulate programmatically internal frameshifts. 
	
Ran test using CLUSTALO but since it does not product tabular output this will not be useful for the next steps


##### [5] Evaluate BLASTall output

	-length
	-identity
	
Start with large dataset and then breakdown into smaller datasets based on criteria set for BLAST queries.
Repeated iterations should make it easier to detect bad sequences or outliers

*Determined that breaking down into taxonomic groups earlier is the better option*

*Reviewing smaller tax groups using BLASTall to determine cutoffs*

#### [6] Create BLAST databases (nucleotide and protein)

## Workflow

The entire workflow can be executed from a Jupyter notebook with the ability to review and tweak parameters at almost each step. At the same time, the modular nature of the components allows one to wrap the entire workflow into a single script that can be executed without user intervention. 

Individual steps of the workflow are described below:

### 1. Fetch gene data 
```
input: query
output: bdbag archive, gene list
```
In this step, user provides an Entrez query to be used with the NCBI Gene database that will be used to obtain a list of all NCBI GeneIDs and a data archive containing sequence and metadata. 

### 2. Evaluate names
```
input: bdbag archive
output: gene names list
```
The bdbag archive returned by NCBI Datasets contains a data table that will be parsed to obtain a unique list of all gene symbols from the data. A tabular output showing the gene name and the number of sequences with that name allows the curator to quickly check for outliers in gene names. 

### 3. Evaluate sequence lengths
```
input: bdbag archive
output: summary statistics table
```
Sequence length information is extracted from the data table and a set of summary statistics are presented to the curator. Curator can then set parameters that will filter out any outliers. 
	
#### Script 3a: Output tax/gene table 

#### Script 3b: Blastall

	input: bdbag, taxid list
	output: blast table, msa file
	
	msa file will be used if closer evaluation of alignments is needed. 
	
#### Script 4: Evaluate BLAST table, determine sequences to be removed

	input: BLAST table, parameters for filtering
	output: Seqids that would be removed
	
	Review removed SeqIDs file to determine if the sequences are truly bad or the 
	taxid groupings are too broad and smaller taxid groups should be removed. 
	
#### Script 5: Filter bdbag fasta file to remove SeqIds from Script 4

	input: bdbag
	output: edited bdbag
	
	Results: Final bdbag for BLAST database
	
Notes:
-Use SeqKit to add remove sequences when we do not want to rebuild from the beginning
	

	
	
