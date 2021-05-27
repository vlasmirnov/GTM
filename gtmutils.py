'''
Created on May 24, 2021

@author: Vlad
'''

import dendropy
import treeutils

def extractSpanningTree(startTree, constraintTrees, numLeaves):
    spanningTree = startTree.extract_tree()
    toRemove = []
    for edge in spanningTree.preorder_internal_edge_iter(exclude_seed_edge=True):
        if edge.head_node.label not in ("scaffold", "scaffold_adjacent"):
            toRemove.append(edge.head_node)
    treeutils.collapseEdges(toRemove)
    
    leafClusters = [[n] for n in spanningTree.leaf_nodes() if n.label in ("scaffold", "scaffold_adjacent")]
    for node in spanningTree.preorder_internal_node_iter():
        sCluster = [nbr for nbr in node.incident_edges() if nbr.label in ("scaffold", "scaffold_adjacent")]
        if len(sCluster) < 3:
            cluster = [c for c in node.child_nodes() if c.is_leaf() and c.label not in ("scaffold", "scaffold_adjacent")]
            if len(cluster) > 0:
                leafClusters.append(cluster)
    
    leafClusters.sort(key = lambda c : len(c))
    leavesRemaining = numLeaves        
    leaves = []
    for i, cluster in enumerate(leafClusters):
        n = min(len(cluster), int(leavesRemaining / (len(leafClusters) - i)) )
        leaves.extend(cluster[:n])
        leavesRemaining = leavesRemaining - n
    
    spanningTree = startTree.extract_tree_with_taxa([l.taxon for l in leaves])
    return spanningTree

def getSpanningTreeEdgeSources(startTree, spanningTree):
    edgeSources = set()       
    for edge in startTree.preorder_edge_iter():
        if edge.tail_node is None:
            continue  
        splits = [s for s, b in edge.desc.items() if s is not None and b != s.taxon_namespace.all_taxa_bitmask()]
        if edge.head_node.label in ("scaffold", "scaffold_adjacent") and len(splits) == 1:
            edgeSources.add((edge, splits[0]))
    return edgeSources

def applySpanningTreeBranchLengths(startTree, spanningTree):
    treeutils.annotateTrees(startTree, [spanningTree])
    treeutils.rerootConstraintTrees(startTree, [spanningTree])    
    edgeMap = treeutils.buildEdgeLengthMap(startTree)
    edgeSources = getSpanningTreeEdgeSources(startTree, spanningTree)
    treeutils.applyEdgeLengths(edgeSources, edgeMap)
    
