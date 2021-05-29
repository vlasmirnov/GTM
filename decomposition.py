'''
Created on May 27, 2021

@author: Vlad
'''

import os
import shutil
import random
import time
import dendropy
import sequenceutils
from configs import Configs
    
'''
def decomposeSequences(context):
    time1 = time.time()    
    time2 = time.time()  
    Configs.log("Decomposed {} into {} subsets in {} sec..".format(context.sequencesPath, len(context.subsetPaths), time2-time1))    
'''


def randomDecomposition(subsetsDir, sequences, numSubsets):
    allTaxa = list(sequences.keys())
    random.shuffle(allTaxa)
    
    taxonSubsets = [allTaxa[i :: numSubsets] for i in range(numSubsets)]
    subsetPaths = []
    for n, subset in enumerate(taxonSubsets):
        subsetPath = os.path.join(subsetsDir, "subset_{}.txt".format(n+1))
        subsetPaths.append(subsetPath)                    
        sequenceutils.writeFasta(sequences, subsetPath, subset) 
    return subsetPaths



def decomposeGuideTree(subsetsDir, sequencesPath, guideTreePath, maxSubsetSize, maxNumSubsets, strategy):
    sequences = sequenceutils.readFromFasta(sequencesPath, removeDashes = False)
    guideTree = dendropy.Tree.get(path=guideTreePath, schema="newick", preserve_underscores=True)
    guideTree.collapse_basal_bifurcation()
    
    for edge in guideTree.postorder_edge_iter():
        if len(edge.head_node.child_edges()) > 0:
            edge.childs = sum([e.childs for e in edge.head_node.child_edges()])
        else:
            edge.childs = 1
    guideTree.childs = guideTree.seed_node.edge.childs       
    trees = decomposeTree(guideTree, maxSubsetSize, maxNumSubsets, strategy)
    
    taxonSubsets = []
    for tree in trees:
        keep = [n.taxon.label for n in tree.leaf_nodes()]
        taxonSubsets.append(keep)
    
    return sequenceutils.writeSubsetsToDir(subsetsDir, sequencesPath, taxonSubsets)

def decomposeTree(tree, maxSubsetSize, numSubsets, strategy):
    trees = [tree]
    while len(trees) < numSubsets:
        largestTree = max(trees, key=lambda t : t.childs)
        
        if maxSubsetSize is not None and largestTree.childs <= maxSubsetSize:
            return trees
        else:
            numChilds = largestTree.childs
            if strategy == "centroid":
                e = getCentroidEdge(largestTree)
            elif strategy == "randomcentroid":
                e = getCentroidEdgeRandom(largestTree, int(maxSubsetSize/3))
            t1, t2 = bipartitionByEdge(largestTree, e)
            Configs.log("Decomposing a tree with {} leaves into {} and {}..".format(numChilds, t1.childs, t2.childs))
            trees.remove(largestTree)
            trees = trees + [t1, t2]
    return trees

def bipartitionByEdge(tree, edge):
    newRoot = edge.head_node
    tail = edge.tail_node
    tail.remove_child(newRoot)
    newTree = dendropy.Tree(seed_node=newRoot, taxon_namespace = tree.taxon_namespace)
    newTree.childs = edge.childs
    tree.childs = tree.childs - newTree.childs
    
    curNode = tail
    while curNode is not None:
        curNode.edge.childs = curNode.edge.childs - edge.childs
        curNode = curNode.edge.tail_node
        
    newTree.collapse_basal_bifurcation()
    if tail == tree.seed_node:
        tree.collapse_basal_bifurcation()
    elif len(tail.child_nodes()) == 1:
        childs = tail.child_nodes()
        for child in childs:
            tail.remove_child(child)
            tail.parent_node.add_child(child)
        tail.parent_node.remove_child(tail)
    
    return tree, newTree

def getCentroidEdge(tree):
    bestBalance = float('inf')
    for edge in tree.postorder_edge_iter():
        if edge.tail_node is None or edge.head_node.label == "scaffold":
            continue
        balance = abs(tree.childs/2 - edge.childs)
        if balance < bestBalance:
            bestBalance = balance
            bestEdge = edge  
    return bestEdge

def getCentroidEdgeRandom(tree, minBound = 5):
    candidates = []
    for edge in tree.postorder_internal_edge_iter():
        if edge.tail_node is None: # or edge.head_node.label == "scaffold":
            continue
        if min(edge.childs, tree.childs - edge.childs) >= minBound:
            candidates.append(edge)    
    return random.choice(candidates)


