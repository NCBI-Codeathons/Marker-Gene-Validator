#!/opt/python-all/bin/python3

import sys

uniquesyms = {}

if ( len( sys.argv ) != 2 ):
	print( "Specify an input .tsv file." )
	exit()

with open( sys.argv[ 1 ], "r" ) as tsvfile:
	all_lines = tsvfile.readlines()

gene_count = len( all_lines )
cutoff = gene_count * 0.1

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

for s in sorted( uniquesyms.keys() ):
	output = [ s, str( len( uniquesyms[ s ] ) ) ]
	if ( len( uniquesyms[ s ] ) <= cutoff ):
		output.append( ",".join( sorted( uniquesyms[ s ] ) ) )
	print( "\t".join( output ) )