#!/bin/bash

set -euo pipefail

usage_text="Usage: make_finaldb.sh -b <bdbag> -a <acclist> -p <filename_prefix> -t <threads>"

while getopts :b:a:p:t:h option ;
do
  case "${option}" in
    b ) bdbag=${OPTARG} ;;
    a ) acclist=${OPTARG} ;;
    p ) filepref=${OPTARG} ;;
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

echo -e "$(date) Filtering protein fasta..." ;
seqkit grep -v -f ${acclist} ${prot_fasta} > ${filepref}.fasta 

echo -e "$(date) Create a BLAST database..." 
makeblastdb -in ${filepref}.fasta -parse_seqids -dbtype prot -out ${filepref}_blastdb 
