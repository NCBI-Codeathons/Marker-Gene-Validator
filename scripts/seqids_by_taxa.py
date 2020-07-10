import argparse
import json
import os
import subprocess
from tax_tree import TaxTree
from report_reader import DatasetsReportReader

UNASSIGNED_SEQIDS = "-1"


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--bdbag', required=True, type=str,
                        help='toplevel bdbag directory, e.g. path to ncbi_datasets')
    parser.add_argument('--taxids', required=True, type=str,
                        help='comma-separated list of taxids (at any level of tax tree)')
    parser.add_argument('--output', required=True, type=str,
                        help='Output file for the 2 column table that mapping seqids to taxids')
    parser.add_argument('--email', required=True, type=str,
                        help='email to provide bipython')
    args = parser.parse_args()

    taxtree = TaxTree()
    taxids = args.taxids.split(',')
    taxmap = {}
    tax_to_seqid = {}
    unassigned_seqids = []

    # Retrieve data report from the bd bag
    gene_data_report = get_data_report(args.bdbag)
    # Find any taxids in the report that are not alredy in the tree and add them using Entrez
    add_missing_taxids(taxtree, gene_data_report, args.email)

    # Assign all taxids to the tax group to which they belong. If the tax groups overlap
    # taxids in the overlapping area will be assigned to the last taxid group to which they belong
    for group_taxid in taxids:
        tax_to_seqid[group_taxid] = []
        child_ids = taxtree.get_subtree_taxids(group_taxid)
        if not child_ids:
            print(f"Warning - No taxids found for {group_taxid} subtree")
            continue
        for child_id in child_ids:
            taxmap[child_id] = group_taxid

    # add a taxid group to capture seqids which were not found to belong to any group
    tax_to_seqid[UNASSIGNED_SEQIDS] = []

    # Go over all transcripts and assign the seqids (based on their taxid) to the group
    # to which they belong, or the UNASSSIGNED group
    for gene in gene_data_report.genes:
        taxgroup = taxmap.get(str(gene.tax_id), "")
        for transcript in gene.transcripts:
            if not taxgroup:
                tax_to_seqid[UNASSIGNED_SEQIDS].append(transcript.accession_version)
            else:
                tax_to_seqid[taxgroup].append(transcript.accession_version)

    # write out the seqids with their corresponding tax group
    with open(args.output, "w") as f:
        for taxgroup,taxids in tax_to_seqid.items():
            for taxid in taxids:
                f.write(f'{taxgroup}\t{taxid}\n')


def get_data_report(bdbag):
    report_reader = DatasetsReportReader()
    report_file = os.path.join(bdbag, "data/data_report.yaml")
    if not os.path.isfile(report_file):
        print(f'Error opening report file. File not found: {report_file}')
        return 1

    return report_reader.gene_report_from_file(report_file)


def add_missing_taxids(taxtree, gene_data_report, email):
    missing_taxa = set()
    for gene in gene_data_report.genes:
        if not taxtree.get_org_if_exists(str(gene.tax_id)):
            missing_taxa.add(str(gene.tax_id))
    uids = ",".join(missing_taxa)
    taxtree.add_entrez_taxa(uids, email)

if __name__ == "__main__":
    main()
