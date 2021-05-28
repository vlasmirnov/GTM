'''
Created on May 24, 2021

@author: Vlad
'''

import dendropy
import treeutils
import os
import shutil
import sequenceutils
import random
from tools import external_tools

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

def computeSpanningTree(tree, alignmentPath, treePath):
    baseName = os.path.basename(treePath).split(".")[0]
    unoptimizedPath = os.path.join(os.path.dirname(treePath), "unoptimized_{}.tre".format(baseName))
    workingDir = os.path.join(os.path.dirname(treePath), "raxml_{}".format(baseName))
    outModelPath = os.path.join(os.path.dirname(treePath), "model_{}.txt".format(baseName))
    reducedAlignPath = os.path.join(os.path.dirname(treePath), "alignment_{}.txt".format(baseName))
    
    taxa = [n.taxon.label for n in tree.leaf_nodes()]
    align = sequenceutils.readFromFasta(alignmentPath)
    sequenceutils.writeFasta(align, reducedAlignPath, taxa)
    
    tree.resolve_polytomies()
    for edge in tree.edges():
        edge.length = None
    treeutils.writeTree(tree, unoptimizedPath)
    
    if os.path.exists(workingDir):
        shutil.rmtree(workingDir)
    os.makedirs(workingDir)    
    task = external_tools.runRaxmlNgOptimize(reducedAlignPath, None, unoptimizedPath, workingDir, outModelPath, treePath)
    external_tools.runCommand(**task)

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

def createRandomTree(taxa):
    #leafs = list(nameSpace.bitmask_taxa_list(nameSpace.all_taxa_bitmask()))
    random.shuffle(taxa)   
    
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
    
