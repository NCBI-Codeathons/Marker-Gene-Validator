# Build a marker gene reference set
NCBI Datasets Codeathon Team 1

#### Problem: Efficiently create curated databases of genes for GenBank Foosh pipelines. 

Creating the curated databases involves several manual steps.  The process starts with complicated Entrez queries.  The nucleotide/protein sequences are then manually evaluated using alighnments and BLAST. Due to taxonomoic breadth this
process is repeates several times to try and catch outliers and bad data. 

We are going to start with the gene CYTB since we can download from the RefSeq database (keeps the initial starting sequences small), the gene does not have introns, and is similar to COX1 for which GenBank already has a foosh pipeline. 

##### [1] Develop query to pull nucleotide/proteins sequences of desired gene sequence from Entrez

*esearch -db gene -query 'cytb[gene name]' | efetch -format uid > cytb_geneids.txt

generates 9880 geneids linking to sequences of cytb*


	-Capture organism information
		-taxID
		-lineage 
		-organism name

	-Capture gene name/product name

	-Capture GeneIDs

##### [2] Use datasets for further analysis once GeneIDs are captured

##### [3] Validate the sequences with a series of BLAST queries

	-TaxID
	-length
	-identity
	
Start with large dataset and then breakdown into smaller datasets based on criteria set for BLAST queries.
Repeated iterations should make it easier to detect bad sequences or outliers

#### [4] Create BLAST databases (nucleotide and protein)
	-include organism name
