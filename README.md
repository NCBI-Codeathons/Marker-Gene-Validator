# Build a marker gene reference set

#### Problem: Efficiently create curated databases of genes for GenBank Foosh pipelines. 

Creating the curated databases involves several manual steps.  The process starts with complicated Entrez queries.  The nucleotide/protein sequences are then manually evaluated using alignments and BLAST. Due to taxonomoic breadth this process is repeated several times to try and catch outliers and bad data. 

We are going to start with the gene CYTB since we can download from the RefSeq database (keeps the initial starting sequences small), the gene does not have introns, and is similar to COX1 for which GenBank already has a foosh pipeline. 

##### [1] Develop query to pull nucleotide/proteins sequences of desired gene sequence from Entrez. 
This query can be as limited or as broad as the gene family of interest and its requirements 

`this works but I think it needs to be more restrictive esearch -db gene -query 'metazoa[orgn] AND cytb[gene name]' | efetch -format uid > cytb_geneids.txt`

generates 9880 geneids linking to sequences of cytb*

*We have refined the query*

##### [2] Use NCBI Datasets to download metadata and sequence data for gene IDs from step 1

	-The datasets provide the nucleotide sequence, protein sequence and associated metadata
	-From the datasets the information needed for the next steps can be extracted.
		-after extraction check the metadata table  (good way check query is providing what you are expecting):
			-product names
			-gene names
	-Extract protein fasta set and check fasta headers.  (may need script in case something changes in the future)
	
		>YP_003667946.1 cytochrome b (mitochondrion) [Carassius gibelio]
	
	-Split at this point into different TAXIDs due to limitations of alighments and BLAST (?) 
	
*We can extract the full data set*
*We can download the expanded taxonomy information*
*Use NCBI datasets and split into taxonomy groupings*
		
	
##### [3] Create histogram to assess length variation. Evaluate sequences at high and low ends. Create script to remove
sequences of certain lengths

*created histogram in Jupyter* 

##### [4] Run CLUSTALW to align the proteins.  From the alignment assess the following:
	
	-what can be programatically removed based on data provided
	-sort by tax groups to evaluate 'outliers'
	-Is there a way to evaulate programmatically internal frameshifts. 
	
	Ran test using CLUSTALO but since it does not product tabular output this would not work for the above steps.


##### [5] Validate the sequences with a series of BLAST queries

	-TaxID
	-length
	-identity
	
Start with large dataset and then breakdown into smaller datasets based on criteria set for BLAST queries.
Repeated iterations should make it easier to detect bad sequences or outliers

*Determined that breaking down into taxonomic groups earlier is the better option*
*Reviewing smaller tax groups using BLASTall to determine cutoffs*

#### [6] Create BLAST databases (nucleotide and protein)
	
	
