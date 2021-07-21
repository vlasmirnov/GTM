'''
Created on Jul 8, 2021

@author: Vlad
'''

import dendropy

def populateEdgeDesc(edge):
    edge.desc = {}
    edge.descUnassigned = 0
    if len(edge.head_node.child_edges()) == 0:
        if edge.head_node.taxon.taxonSubtree is not None:
            edge.desc = {edge.head_node.taxon.taxonSubtree : edge.head_node.taxon.taxonBitmask}
        else:
            edge.descUnassigned = 1
    else:
        edge.desc = combineDescMaps([child.desc for child in edge.head_node.child_edges()])
        edge.descUnassigned = sum(child.descUnassigned for child in edge.head_node.child_edges())
        
def populateEdgeMaps(tree):
    tree.edgeMap = {}
    tree.subEdgeMap = {}    
    for edge in tree.preorder_edge_iter():
        bitmask = edge.desc[tree]
        tree.edgeMap[buildEdgeKey(tree, bitmask)] = (edge, edge.length)
        tree.subEdgeMap[bitmask & tree.startTreeBitmask] = edge 

def buildEdgeKey(tree, bitmask):
    return min(bitmask, ~bitmask & tree.taxon_namespace.all_taxa_bitmask())

def combineDescMaps(descMaps):
    newDesc = {}
    for desc in descMaps:
        for subtree, bitmask in desc.items():
            newDesc[subtree] = newDesc.get(subtree, 0) | bitmask
    return newDesc    

def resolvePolytomy(node, children):
    newNode = dendropy.Node()
    node.add_child(newNode)
    for child in children:
        newNode.add_child(node.remove_child(child))
    populateEdgeDesc(newNode.edge)
    return newNode    

def checkBitmaskViolatesSubtree(bitmask, subtree):
    return bitmask & subtree.startTreeBitmask not in subtree.subEdgeMap and \
        ~bitmask & subtree.startTreeBitmask not in subtree.subEdgeMap