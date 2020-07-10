#!/bin/bash

set -euo pipefail

usage_text="Usage: blast_all.sh -b <bdbag> -a <acclist> -t <threads>"

while getopts :b:a:t:h option ;
do
  case "${option}" in
    b ) bdbag=${OPTARG} ;;
    a ) acclist=${OPTARG} ;;
    t ) threads=${OPTARG} ;;
    h ) echo "${usage_text}" 
    exit 1 ;;
    \?) echo "${usage_text}"
    exit 1 ;;
  esac
done

prot_fasta="${bdbag}/data/protein.faa"

if [ ! -f "${prot_fasta}" ]; then
    echo "${prot_fasta} does not exist."
    exit 1
fi

txid_list=$(cut -f1 ${acclist} | sort -u)
echo -e "Number of taxids in ${acclist}: $(echo $txid_list | sed 's/ /\n/g' | wc -l)"

for txid in ${txid_list} ; do 
    echo -e "Processing ${txid}"

    echo -e "$(date) Filtering protein fasta..." ;
    seqkit grep -f <(grep -w ${txid} ${acclist} | cut -f2) ${prot_fasta} > ${txid}_input.fa ;

    echo -e "$(date) Create a BLAST database..." ;
    makeblastdb -in ${txid}_input.fa -parse_seqids -dbtype prot -out ${txid}_blastdb ;

    echo -e "$(date) Running all-vs-all blast..." ;
    blastp -db ${txid}_blastdb -parse_deflines -num_threads ${threads} -query ${txid}_input.fa -outfmt 11 -out ${txid}_archive.asn ;

    echo -e "$(date) Generating blast tabular output..."
    blast_formatter -archive ${txid}_archive.asn -outfmt '7 std qcovs' -out ${txid}_output.tsv 

    echo -e "$(date) Generating blast seq-align asn..."
    blast_formatter -archive ${txid}_archive.asn -outfmt 8 -out ${txid}_output.asn
done