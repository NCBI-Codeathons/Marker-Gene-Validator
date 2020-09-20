import logging
from Bio import Entrez
import ncbi.datasets
from ncbi.datasets.rest import ApiException

logging.basicConfig(filename='process_markers.log', level=logging.DEBUG)
logger = logging.getLogger(__name__)


'''
Construct a taxtree using the ncbi datasets GeneApi service.  
Tree is constructed using a root tax id and tax info can be returned
based on tax-id.
'''
class TaxTree:
    # defaults to metazoan
    def __init__(self, tax_id="33208"):
        self.node_map = {}
        self.root = {}
        with ncbi.datasets.ApiClient() as api_client:
            # api_instance = ncbi.datasets.TaxTreeApi(api_client)
            api_instance = ncbi.datasets.GeneApi(api_client)
            try:
                # Retrieve tax tree by taxonomy ID 
                api_response = api_instance.gene_tax_tree(tax_id)
                tree = api_response.to_dict()
                self.node_map = {}
                self._build_node_map(tree)
                self.root = self.node_map[str(tax_id)]
            except ApiException as e:
                logger.error(f"Exception when calling GeneApi->gene_tax_tree: e")

    def _build_node_map(self, node):
        self.node_map[node['tax_id']] =  { 
            'tax_id': node.get('tax_id', None),
            'parent_tax_id': node.get('parent_tax_id', None),
            'sci_name': node.get('sci_name', ""),
            'common_name': node.get('common_name', ""),
            'rank': node.get('rank', ""),
            'children': [org['tax_id'] for org in (node['children'] or [])],
        }

        if node['children']:
            for org in node.get('children', []):
                self._build_node_map(org)

    def add_tax_node(self, tax_id, tax_lineage, sci_name, common_name, rank):
        current = self.get_org_if_exists(tax_id)
        if current:
            # this can happen if an added node was an ancestor of previously added node, so we just update in that case
            current['sci_name'] = sci_name
            current['common_namme'] = common_name
            current['rank'] = rank
            return True

        new_nodes = []
        prev_node = tax_lineage[0]
        merge_node_found = False
        for parent in tax_lineage[1:]:
            new_nodes.insert(0, (prev_node, parent['TaxId']))
            if self.get_org_if_exists(parent['TaxId']):
                merge_node_found = True
                break
            prev_node = parent

        if merge_node_found:
            # Add 1 or more nodes (this includes intermediate nodes between the node being added and
            # the lowest common ancestor in the tax tree).  If there are intermediate nodes, we add
            # in order from highest-ranked down to the new node so that the parent nodes exist
            # when a new child node is created
            for node in new_nodes:
                self.node_map[node[0]['TaxId']] =  { 
                    'tax_id': node[0]['TaxId'],
                    'parent_tax_id': node[1],
                    'sci_name': node[0]['ScientificName'],
                    'children': [],
                }
                self.node_map[node[1]]['children'].append(node[0]['TaxId']) 

                if node[0] == tax_id:
                    self.node_map[tax_id]['sci_name'] = common_name
                    self.node_map[tax_id]['sci_name'] = rank
            return True
        else:
            logger.warning(f"No connection found for {[org['TaxId'] for org in tax_lineage]} to existing tree")
            return False

    def add_entrez_taxa(self, taxids, email):
        Entrez.email = email

        handle = Entrez.efetch(db="taxonomy", id=taxids, mode="text", rettype="xml")
        records = Entrez.read(handle)
        for taxon in records:
            taxid = str(taxon["TaxId"])
            name = taxon["ScientificName"]
            rank = taxon.get("Rank", "")
            common_name = taxon.get("OtherNames", {}).get("ComonName", "")
            lineage = []
            for t in taxon["LineageEx"]:
                lineage.insert(0, t)
            lineage.insert(0, {'TaxId': taxid, 'ScientificName': name})
            self.add_tax_node(str(taxid), lineage, name, common_name, rank)

    def get_lineage(self, tax_id):
        lineage = []
        if self.node_map.get(tax_id):
            lineage.append((tax_id, self.node_map[tax_id]['sci_name']))
            parent = self.node_map[tax_id].get('parent_tax_id')
            while parent and self.node_map.get(parent):
                lineage.append((parent, self.node_map[parent]['sci_name']))
                parent = self.node_map[parent].get('parent_tax_id')
        return lineage

    def get_org_if_exists(self, tax_id):
        return self.node_map.get(tax_id)

    def _all_children(self, node, result):
        if not node:
            return []

        result.append(node['tax_id'])
        if not node['children'] or (node['children'] is None):
            return

        for child_taxid in node['children']:
            self._all_children(self.get_org_if_exists(child_taxid), result)

    def get_subtree_taxids(self, parent_tax_id):
        parent_node = self.get_org_if_exists(parent_tax_id)
        results = []
        if not parent_node:
            logger.warning(f"parent node not found: {parent_tax_id}")
        else: 
            self._all_children(parent_node, results)
        return results
