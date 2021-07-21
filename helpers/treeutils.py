'''
Created on Nov 29, 2020

@author: Vlad
'''

import random
import dendropy

def loadTree(treePath):
    tree = dendropy.Tree.get(path=treePath, schema="newick", preserve_underscores=True)
    tree.is_rooted = False
    #tree.resolve_polytomies(limit=2)
    tree.collapse_basal_bifurcation()
    return tree

def writeTree(tree, treePath):
    treeString = tree.as_string(schema="newick", suppress_rooting = True)
    with open(treePath, "w") as f:
        f.write(treeString) 
        
def createRandomTree(labels):
    #leafs = list(nameSpace.bitmask_taxa_list(nameSpace.all_taxa_bitmask()))
    random.shuffle(labels)
    taxa = [dendropy.Taxon(label) for label in labels]   
    
    center = dendropy.Node()
    nodes = [dendropy.Node(taxon = taxa[0]), dendropy.Node(taxon = taxa[1]), dendropy.Node(taxon = taxa[2])]
    center.set_child_nodes((nodes[0], nodes[1], nodes[2]))
    tree = dendropy.Tree(seed_node=center)
    
    for leaf in taxa[3:]:
        newNode = dendropy.Node(taxon = leaf)
        edge = random.choice(nodes).edge
        n = edge.tail_node.new_child()
        n.add_child(edge.tail_node.remove_child(edge.head_node))
        n.add_child(newNode)
        nodes.append(n)
        nodes.append(newNode)
    
    tree.is_rooted = False
    return tree  

def collapseEdges(headNodes):
    for n in headNodes:
        childs = n.child_nodes()
        for child in childs:
            n.remove_child(child)
            n.parent_node.add_child(child)
        n.parent_node.remove_child(n)    

def attachNodetoNode(sNode, cNode):
    if cNode.parent_node is not None:
        cNode.parent_node.remove_child(cNode)
    sNode.add_child(cNode)

def extractSpanningTree(tree, numLeaves, edgeSet):
    spanningTree = tree.extract_tree()
    edgeFunc = lambda e: e.head_node.extraction_source.edge in edgeSet
    toRemove = []
    for edge in spanningTree.preorder_internal_edge_iter(exclude_seed_edge=True):
        if not edgeFunc(edge):
            toRemove.append(edge.head_node)
    collapseEdges(toRemove)
    
    leafClusters = [[n] for n in spanningTree.leaf_nodes() if edgeFunc(n.edge)]
    for node in spanningTree.preorder_internal_node_iter():
        sCluster = [nbr for nbr in node.incident_edges() if edgeFunc(nbr)]
        if len(sCluster) < 3:
            cluster = [c for c in node.child_nodes() if c.is_leaf() and not edgeFunc(c.edge)]
            if len(cluster) > 0:
                leafClusters.append(cluster)
    
    leafClusters.sort(key = lambda c : len(c))
    leavesRemaining = numLeaves        
    leaves = []
    print("Extracting spanning tree with maximum {} leaves over {} clusters..".format(numLeaves, len(leafClusters)))
    for i, cluster in enumerate(leafClusters):
        if numLeaves is not None and numLeaves > 0:
            n = min(len(cluster), int(leavesRemaining / (len(leafClusters) - i)) )
            leavesRemaining = leavesRemaining - n
        else:
            n = 1
        leaves.extend(cluster[:n])
        
    
    spanningTree = tree.extract_tree_with_taxa([l.taxon for l in leaves])
    return spanningTree
