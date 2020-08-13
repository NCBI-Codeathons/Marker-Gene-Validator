import argparse
import logging
import json
import subprocess
from Bio import Entrez

gene_ids_file = "gene_ids.txt"

logging.basicConfig(filename='process_markers.log', level=logging.DEBUG)
logger = logging.getLogger(f"{__file__}:{__name__}")


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--output-file', required=True, type=str, dest='output_file',
                        help='output zipfile name for dataset')
    parser.add_argument('--email', required=True, type=str, dest='email',
                        help='email to provide biopython')
    parser.add_argument('--query', required=False, type=str, dest='query',
                        default='Hominoidea [ORGN] AND cytb [GENE] AND source mitochondrion [PROP] NOT rnatype mrna [PROP] NOT srcdb pdb [PROP] NOT uncultured NOT unverified',
                        help='query to collect genes for')
    args = parser.parse_args()
    populate_gene_ids_file(args.query, args.email, gene_ids_file)
    json_data = format_file_data_into_json(gene_ids_file)
    obtain_gene_datasets(json_data, args.output_file)

    
def log_and_raise_exception(error_message):
    """input: error_message (error message string)
       logs error and raises exception
    """ 
    logger.error(error_message)
    raise Exception(error_message)

    
def populate_gene_ids_file(query, email, gene_ids_file):
    """input: query (query to give esearch), email (email to provide biopython),
              gene_ids_file (where to output gene ids)
       Calls esearch to obtain ids for gene search and writes the ids to gene_ids_file
    """
    Entrez.email = email
    esearch_record = Entrez.read(Entrez.esearch(db="gene", term=query))
    total_record_number = int(esearch_record["Count"])
    ids = esearch_record["IdList"]
    if total_record_number > len(ids):
        ids = Entrez.read(Entrez.esearch(db="gene", term=query, retmax=total_record_number))["IdList"]
    result_count_message = f"Gene search for query '{query}' returned {total_record_number} results"
    print(result_count_message)
    logger.info(result_count_message)
    if len(ids) == 0:
        log_and_raise_exception("The provided search query gave yielded no results in esearch")
    with open(gene_ids_file, "w+") as f:
        for result in ids:
            f.write(f"{result}\n")
            

def obtain_gene_datasets(gene_json, output_file):
    """input: gene_json (json string request data), output_file (file to write datasets to)
       Calls datasets api to get gene datasets
    """
    # Fix later to not use curl and handle non-200 responses
    return_code = subprocess.check_call(f'curl -X POST "https://api.ncbi.nlm.nih.gov/datasets/v1alpha/download/gene?filename=ncbi_dataset.zip" -H "accept: application/zip" -H "Content-Type: application/json" --data \'{gene_json}\' > {output_file}', shell=True)
    if return_code != 0:
        log_and_raise_exception("Command to obtain datasets failed")


def format_file_data_into_json(data_file):
    """input: data_file (sequence ids separated by \n)
       output: json request structure of gene ids to pass to datasets api
    """
    with open(data_file, "r") as f:
        content = f.read()
    genes = content.strip().split("\n")
    return json.dumps({'gene_ids': [int(gene) for gene in genes], 'include_annotation_type': ['FASTA_ALL']})


if __name__ == "__main__":
    main()
