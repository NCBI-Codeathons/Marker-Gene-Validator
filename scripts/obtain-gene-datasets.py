import argparse
import json
import subprocess

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--gene-id-file', required=True, type=str, dest='gene_id_file',
                        help='file of gene ids to download datasets for')
    parser.add_argument('--output-file', required=True, type=str, dest='output_file',
                        help='output zipfile name for dataset')
    args = parser.parse_args()
    json_data = format_file_data_into_json(args.gene_id_file)
    obtain_gene_datasets(json_data, args.output_file)

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
