'''
Created on May 27, 2021

@author: Vlad
'''

import os
import shutil
import random
import time
import dendropy
import heapq
import sequenceutils
from configs import Configs
from platform import node
    
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
    bestBalance = float('inf')
    for edge in tree.postorder_internal_edge_iter():
        if edge.tail_node is None: # or edge.head_node.label == "scaffold":
            continue
        if min(edge.childs, tree.childs - edge.childs) >= minBound:
            candidates.append(edge) 
        balance = abs(tree.childs/2 - edge.childs)
        if balance < bestBalance:
            bestBalance = balance
            bestEdge = edge     
    return random.choice(candidates) if len(candidates) > 0 else bestEdge

def treemer(tree, skeletonSize):
    pairHeap = []
    usedPairs = set()
    labelNodeMap = {}
    assigns = {}
    
    for node in tree.postorder_node_iter():
        node.shortestDistance = None
        node.childDescMap = {}
        node.descMap = {}  
        node.edge.length = max(node.edge.length,0) if node.edge.length is not None else 1
              
        if node.taxon is not None:
            node.descMap = {node.taxon.label : (0, None, 0)}
            labelNodeMap[node.taxon.label] = node
        else:
            updateTreemerNode(node)
            updateTreemerNodeShortestPath(node)
            #if (node.shortestDistance[1], node.shortestDistance[2]) not in usedPairs:
            if node.shortestDistance is not None:
                heapq.heappush(pairHeap, node.shortestDistance)
    
    while len(labelNodeMap) > skeletonSize:
        dist, tax1, tax2 = heapq.heappop(pairHeap)        
        if tax1 in assigns or tax2 in assigns:
            continue
        print(dist, tax1, tax2)
        assigns[tax1] = tax2  
        leaf = labelNodeMap.pop(tax1)
        parnt = leaf.parent_node
        parnt.remove_child(leaf)
        
        node = parnt
        parnt = parnt.parent_node
        if len(node.child_nodes()) == 1:
            child = node.child_nodes()[0]
            node.remove_child(child)
            if parnt is None:
                #print("rooting..")
                tree.reseed_at(child, collapse_unrooted_basal_bifurcation = False, suppress_unifurcations = False)
            else:
                child.edge.length = child.edge.length + node.edge.length
                parnt.remove_child(node)
                parnt.add_child(child)
                #print("fixing..")
                updateTreemerNodeRemoveTaxon(parnt, tax1, pairHeap)
        else:
            updateTreemerNodeRemoveTaxon(node, tax1, pairHeap)
    
    #print(assigns)
    for taxon in assigns:
        recurseAssign(assigns, taxon)
        
    return tree, labelNodeMap, assigns

def recurseAssign(assigns, taxon):
    stack = []
    while assigns[taxon] in assigns:
        stack.append(taxon)
        taxon = assigns[taxon]
    for tax in stack:
        assigns[tax] = assigns[taxon]
    
def updateTreemerNode(node, childs = None):
    node.childDescMap = {}
    node.descMap = {}  
    childs = node.child_nodes() if childs is None else childs
    for child in childs:
        node.childDescMap[child] = set()
        for taxon, value in child.descMap.items():
            dist, desc, path = value
            if path < 2:
                node.descMap[taxon] = (dist + child.edge.length, child, path + 1)
                node.childDescMap[child].add(taxon) 
        if len(node.childDescMap[child]) == 0:
            node.childDescMap.pop(child)        
    
def updateTreemerNodeShortestPath(node):
    node.shortestDistance = None
    sortedNearestChilds = [min([(node.descMap[taxon][0], taxon) for taxon in s]) for c,s in node.childDescMap.items()]
    sortedNearestChilds.sort()
    if len(sortedNearestChilds) > 1:
        first, second = sortedNearestChilds[0], sortedNearestChilds[1]
        node.shortestDistance = (first[0] + second[0], min(first[1], second[1]), max(first[1], second[1]))

def updateTreemerNodeRemoveTaxon(node, taxon, pairHeap):
    updateTreemerNode(node)
        
    oldShortest = node.shortestDistance
    updateTreemerNodeShortestPath(node)
    if node.shortestDistance is not None and node.shortestDistance != oldShortest:
        heapq.heappush(pairHeap, node.shortestDistance)
                   
    if node.parent_node is not None and taxon in node.parent_node.descMap:
        updateTreemerNodeRemoveTaxon(node.parent_node, taxon, node, pairHeap)
           
       
