import sys, json

uniquesyms = {}
uniquedescs = {}

if ( len( sys.argv ) != 3 ):
	print( "Specify an input .tsv file and an output file." )
	exit()

with open( sys.argv[ 1 ], "r" ) as tsvfile:
	all_lines = tsvfile.readlines()

for x in all_lines:
	gene_info = x.rstrip().split("\t")
	id = gene_info[0]
	symbol = gene_info[1].upper()
	if (id == "gene_id"):
		continue
	if (symbol not in uniquesyms):
		uniquesyms[ symbol ] = [ id ]
	else:
		uniquesyms[ symbol ].append( id )

with open( sys.argv[ 2 ], "w" ) as outfile:
	json.dump( uniquesyms, outfile )
