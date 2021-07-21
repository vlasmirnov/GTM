'''
Created on May 27, 2021

@author: Vlad
'''

import os
import random
import dendropy
import heapq
from helpers import sequenceutils
from helpers import treeutils
import math
from configs import Configs
from tools import methods
    
'''
def decomposeSequences(context):
    time1 = time.time()    
    time2 = time.time()  
    Configs.log("Decomposed {} into {} subsets in {} sec..".format(context.sequencesPath, len(context.subsetPaths), time2-time1))    
'''

def decomposeSequences(context, workingDir):
    if context.decompositionStrategy == "fulltreemer":
        subsetPaths = treemerDecomposition(workingDir, context.alignmentPath, context.startTreePath)
    elif context.decompositionStrategy == "random":
        subsetPaths = randomDecomposition(workingDir, context.alignmentPath, context.startTreePath)
    elif context.decompositionStrategy == "random_postorder":
        subsetPaths = randomPostorderDecomposition(workingDir, context.alignmentPath, context.startTreePath)
    elif context.decompositionStrategy == "random_postorder_stagger":
        subsetPaths = randomPostorderStaggeredDecomposition(workingDir, context.alignmentPath, context.startTreePath)
    elif context.decompositionStrategy == "random_postorder_stagger_skip":
        subsetPaths = randomPostorderStaggeredSkipDecomposition(workingDir, context.alignmentPath, context.startTreePath)
    elif context.decompositionStrategy == "random_p_pss":
        if context.iterations % 2 == 0:
            subsetPaths = randomPostorderDecomposition(workingDir, context.alignmentPath, context.startTreePath)
        elif context.iterations % 2 == 1:
            subsetPaths = randomPostorderStaggeredSkipDecomposition(workingDir, context.alignmentPath, context.startTreePath)
    elif context.decompositionStrategy == "random_p_rc":
        if context.iterations % 2 == 0:
            subsetPaths = randomPostorderDecomposition(workingDir, context.alignmentPath, context.startTreePath)
        elif context.iterations % 2 == 1:
            subsetPaths = decomposeGuideTree(context.subsetsDir, context.alignmentPath, context.startTreePath, 
                                             Configs.decompositionMaxSubsetSize, Configs.decompositionMaxNumSubsets, "randomcentroid")
    elif context.decompositionStrategy == "random_rc_p":
        if context.iterations % 2 == 0:  
            subsetPaths = decomposeGuideTree(context.subsetsDir, context.alignmentPath, context.startTreePath, 
                                             Configs.decompositionMaxSubsetSize, Configs.decompositionMaxNumSubsets, "randomcentroid")
        elif context.iterations % 2 == 1:
            subsetPaths = randomPostorderDecomposition(workingDir, context.alignmentPath, context.startTreePath)
    else:                
        subsetPaths = decomposeGuideTree(context.subsetsDir, context.alignmentPath, context.startTreePath, Configs.decompositionMaxSubsetSize, 
                                                       Configs.decompositionMaxNumSubsets, context.decompositionStrategy)
    return subsetPaths

def randomDecompositionOld(subsetsDir, sequences, numSubsets):
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
           
def treemerDecomposition(workingDir, alignmentPath, startTreePath, numSubsets = None):
    subsetsDir = os.path.join(workingDir, "subsets")
    tree = treeutils.loadTree(startTreePath)
    numSubsets = numSubsets if numSubsets is not None else math.ceil(len(tree.leaf_nodes())/Configs.decompositionMaxSubsetSize)
    tree, labelNodeMap, assigns = treemer(tree, numSubsets)
    subsets = {label : [label] for label in labelNodeMap}
    for tax in assigns:
        subsets[assigns[tax]].append(tax)
    subsets = list(subsets.values())
    return sequenceutils.writeSubsetsToDir(subsetsDir, alignmentPath, subsets)     

def randomDecomposition(workingDir, alignmentPath, startTreePath, numSubsets = None):
    subsetsDir = os.path.join(workingDir, "subsets")
    taxa = sequenceutils.readTaxaFromFasta(alignmentPath)
    numSubsets = numSubsets if numSubsets is not None else math.ceil(len(taxa)/Configs.decompositionMaxSubsetSize)
    random.shuffle(taxa)
    subsets = [taxa[i :: numSubsets] for i in range(numSubsets)]
    return sequenceutils.writeSubsetsToDir(subsetsDir, alignmentPath, subsets)

def randomPostorderDecomposition(workingDir, alignmentPath, startTreePath, numSubsets = None):
    subsetsDir = os.path.join(workingDir, "subsets")
    tree = treeutils.loadTree(startTreePath)
    taxa = [n.taxon.label for n in tree.postorder_node_iter() if n.taxon is not None]
    numSubsets = numSubsets if numSubsets is not None else math.ceil(len(taxa)/Configs.decompositionMaxSubsetSize)
    rotator = random.randint(0, len(taxa)-1)
    taxa = taxa[rotator:] + taxa[:rotator]
    subsets = []
    pos, leavesRemaining = 0, len(taxa)        
    for i in range(numSubsets):
        n = int(leavesRemaining / (numSubsets - i))
        subsets.append(taxa[pos : pos + n])
        leavesRemaining = leavesRemaining - n
        pos = pos + n
    return sequenceutils.writeSubsetsToDir(subsetsDir, alignmentPath, subsets)

def randomPostorderStaggeredDecomposition(workingDir, alignmentPath, startTreePath, numSubsets = None):
    subsetsDir = os.path.join(workingDir, "subsets")
    tree = treeutils.loadTree(startTreePath)
    taxa = [n.taxon.label for n in tree.postorder_node_iter() if n.taxon is not None]
    numSubsets = numSubsets if numSubsets is not None else math.ceil(len(taxa)/Configs.decompositionMaxSubsetSize)
    subsets = [taxa[i :: numSubsets] for i in range(numSubsets)]
    return sequenceutils.writeSubsetsToDir(subsetsDir, alignmentPath, subsets)

def randomPostorderStaggeredSkipDecomposition(workingDir, alignmentPath, startTreePath, numSubsets = None):
    subsetsDir = os.path.join(workingDir, "subsets")
    tree = treeutils.loadTree(startTreePath)
    taxa = [n.taxon.label for n in tree.postorder_node_iter() if n.taxon is not None]
    numSubsets = numSubsets if numSubsets is not None else math.ceil(len(taxa)/Configs.decompositionMaxSubsetSize)
    rotator = random.randint(0, len(taxa)-1)
    taxa = taxa[rotator:] + taxa[:rotator]
    oddTaxa, evenTaxa = taxa[0 :: 2], taxa[1 :: 2]
    numSubsets = int(numSubsets / 2)
    subsets = []
    for taxaList in [oddTaxa, evenTaxa]:
        pos, leavesRemaining = 0, len(taxaList)        
        for i in range(numSubsets):
            n = int(leavesRemaining / (numSubsets - i))
            subsets.append(taxaList[pos : pos + n])
            leavesRemaining = leavesRemaining - n
            pos = pos + n
    return sequenceutils.writeSubsetsToDir(subsetsDir, alignmentPath, subsets)

def guideTreeSpan(workingDir, guideTreePath, **kwargs):
    subsets = [s for p, s in kwargs["subsetPaths"].items()]
    subsets.sort(key = lambda c : len(c))
    leavesRemaining = Configs.spanningTreeSize        
    leaves = []
    for i, cluster in enumerate(subsets):
        n = min(len(cluster), int(leavesRemaining / (len(subsets) - i)) )
        leaves.extend(cluster[:n])
        leavesRemaining = leavesRemaining - n
    alignPath = os.path.join(workingDir, "guide_tree_align.txt")
    sequenceutils.writeSubsetsToFiles(kwargs["alignmentPath"], {alignPath : leaves})
    startTreePath = kwargs.get("startTreePath")
    model = kwargs.get("model")
    if startTreePath is not None:
        startTreePath = methods.extractInducedTree(workingDir, startTreePath, alignPath)        
    methods.buildTree(Configs.treeMethod, guideTreePath, alignmentPath = alignPath, startTreePath = startTreePath, model = model)
    
def guideTreeSpanDiverse(workingDir, guideTreePath, **kwargs):
    subtrees = list(kwargs["subtreePaths"].keys())
    leavesRemaining = Configs.decompositionMaxSubsetSize       
    leaves = []
    for i, path in enumerate(subtrees):
        subTree = treeutils.loadTree(path)
        taxa = [n.taxon.label for n in subTree.postorder_node_iter() if n.taxon is not None]
        nmax = min(len(taxa), int(leavesRemaining / (len(subtrees) - i)))
        n = math.ceil(len(taxa) / nmax)        
        sample = taxa[random.randint(0,n-1) :: n]
        leaves.extend(sample)
        leavesRemaining = leavesRemaining - len(sample)
        #print(len(sample))
    alignPath = os.path.join(workingDir, "guide_tree_align.txt")
    sequenceutils.writeSubsetsToFiles(kwargs["alignmentPath"], {alignPath : leaves})
    startTreePath = kwargs.get("startTreePath")
    model = kwargs.get("model")
    if startTreePath is not None:
        startTreePath = methods.extractInducedTree(workingDir, startTreePath, alignPath)    
    methods.buildTree(Configs.treeMethod, guideTreePath, alignmentPath = alignPath, startTreePath = startTreePath, model = model)