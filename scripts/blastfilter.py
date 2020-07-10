import getopt, sys
import numpy

infile_name = ""
pident_cutoff = ""
qcov_cutoff = ""

pident_stats = {}
qcov_stats = {}
to_drop = []

try:
	opts, args = getopt.getopt(sys.argv[1:], "i:", [ "pident=", "qcov=" ])
except getopt.GetoptError as err:
	# print help information and exit:
	print(err) 
	exit()
for o, a in opts:
	if (o == "-i"):
		infile_name = a
	elif (o in ("--pident")):
		pident_cutoff = float( a )
	elif (o in ("--qcov")):
		qcov_cutoff = float( a )
	else:
		assert False, "unhandled option"
 
if (infile_name == ""):
	print("You must specify an input file using the -i parameter.")
	exit()

with open( infile_name, "r" ) as tsvfile:
	all_lines = tsvfile.readlines()

for x in all_lines:
	if (x[ 0 ] == "#" ):
		continue
	blast_info = x.rstrip().split("\t")
	qid, sid, pident, length, mismatch, gapopen, qstart, qend, sstart, send, evalue, bitscore, qcov = blast_info
	if (qid == sid):
		continue
	if (qid not in pident_stats):
		pident_stats[ qid ] = [ float( pident ) ]
	else:
		pident_stats[ qid ].append( float( pident ) ) 

	if (qid not in qcov_stats):
		qcov_stats[ qid ] = [ float( qcov ) ]
	else:
		qcov_stats[ qid ].append( float( qcov ) ) 

# If a percent identity cutoff was specified, add all the query seq ids to the "to_drop" list whose average is under the cutoff and that haven't already been added to "to_drop"
if (pident_cutoff != ""):
	for qid in pident_stats:
		if ( ( qid not in to_drop ) and ( numpy.mean( pident_stats[ qid ] ) <= pident_cutoff ) ):
			to_drop.append( qid )

# If a query coverage cutoff was specified, add all the query seq ids to the "to_drop" list whose qcov average is under the cutoff and that haven't already been added to "to_drop"
if (qcov_cutoff != ""):
	for qid in qcov_stats:
		if ( ( qid not in to_drop ) and ( numpy.mean( qcov_stats[ qid ] ) <= qcov_cutoff ) ):
			to_drop.append( qid )
 
 
print( "\n".join( to_drop ) )












