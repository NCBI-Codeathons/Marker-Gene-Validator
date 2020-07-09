from __future__ import print_function
import time
import ncbi.datasets
from ncbi.datasets.rest import ApiException


'''
Construct a taxtree using the ncbi datasets TaxTreeApi service.  
Tree is constructed using a root tax id and tax info can be returned
based on tax-id.
Tree also supports separate lineage function and subtree function. The subtree
returns all nodes which do not themselves have children, or are at the species level.

Examples:
    taxtree = TaxTree(33208)
    org = taxtree.get_org_if_exists("9606")
    print("human: ", org)
    lineage = taxtree.get_lineage("9606")
    print("human lineage: ", lineage)
    all_children = taxtree.get_all_children_for("207598")
    print("children of 207598:", all_children)
'''
class TaxTree:
    # defaults to metazoan
    def __init__(self, tax_id=33208):
        self.tree = {}
        self.node_map = {}
        with ncbi.datasets.ApiClient() as api_client:
            api_instance = ncbi.datasets.TaxTreeApi(api_client)
            try:
                # Retrieve tax tree by taxonomy ID 
                api_response = api_instance.tax_tree_by_tax_id(tax_id)
                self.tree = api_response.to_dict()
                self.node_map = {}
                self._build_node_map(self.tree)
            except ApiException as e:
                print("Exception when calling TaxTreeApi->tax_tree_by_tax_id: %s\n" % e)

    def _build_node_map(self, node):
        self.node_map[node['tax_id']] =  { 
            'parent_tax_id': node.get('parent_tax_id', None),
            'sci_name': node.get('sci_name', ""),
            'common_name': node.get('common_name', ""),
            'rank': node.get('rank', "")
        }

        if node['children']:
            for org in node.get('children', []):
                self._build_node_map(org)

    def get_tree(self, tax_id=None):
        if not tax_id:
            return self.tree
        return self.node_map.get(tax_id)

    def get_node(self, tax_id):
        return self.get_org_if_exists(tax_id)

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

    def get_subtree(self, node, tax_id):
        if tax_id == node['tax_id']:
            return node
        if node['children']:
            for org in node.get('children', []):
                result = self.get_subtree(org, tax_id)
                if result:
                    return result
        return None

    def all_children(self, node, result):
        if not node:
            return []

        if not node['children'] or (node['children'] is None):
            result.append(self.node_map[node['tax_id']])
            return

        for child in node['children']:
            self.all_children(child, result)

    def get_all_children_for(self, parent_tax_id):
        parent_node = self.get_subtree(self.tree, parent_tax_id)
        result = []
        if not parent_node:
            print("parent node not found: ", parent_tax_id)
        else: 
            self.all_children(parent_node, result)
        return result
