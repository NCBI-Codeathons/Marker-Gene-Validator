import argparse
import json
import subprocess
from Bio import Entrez

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--output-file', required=True, type=str, dest='output_file',
                        help='output zipfile name for dataset')
    parser.add_argument('--email', required=True, type=str, dest='email',
                        help='email to provide biopython')
    parser.add_argument('--query', required=True, type=str, dest='query',
                        help='query to collect genes for')
    gene_ids_file = "gene_ids.txt"
    args = parser.parse_args()
    populate_gene_ids_file(args.query, args.email, gene_ids_file)
    json_data = format_file_data_into_json(gene_ids_file)
    obtain_gene_datasets(json_data, args.output_file)
    
def populate_gene_ids_file(query, email, gene_ids_file):
    Entrez.email = email
    esearch_results = Entrez.read(Entrez.esearch(db="gene", term=query, retmax=10000))["IdList"]
    if len(esearch_results ) == 0:
        raise Exception("The provided search query gave yielded no results in esearch")
    with open(gene_ids_file, "w+") as f:
        for result in esearch_results:
            f.write(f"{result}\n")

def obtain_gene_datasets(gene_json, output_file):
    # Fix later to not use curl and handle non-200 responses
    return_code = subprocess.check_call(f'curl -X POST "https://api.ncbi.nlm.nih.gov/datasets/v1alpha/download/gene?filename=ncbi_dataset.zip" -H "accept: application/zip" -H "Content-Type: application/json" --data \'{gene_json}\' > {output_file}', shell=True)
    if return_code != 0:
        raise Exception("Command to obtain datasets failed")


def format_file_data_into_json(data_file):
    with open(data_file, "r") as f:
        content = f.read()
    genes = content.strip().split("\n")
    return json.dumps({'gene_ids': [int(gene) for gene in genes]})


if __name__ == "__main__":
    main()
