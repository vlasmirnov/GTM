'''
Created on Jul 9, 2021

@author: Vlad
'''

import time
from helpers import treeutils, gtmutils
from gtm import gtm_old
from configs import Configs

def runGtm(gtmContext):
    startTime = time.time()
    Configs.log("Starting GTM..")
    
    Configs.log("Guide tree path: {}".format(gtmContext.startTreePath))    
    Configs.log("Constraint tree paths: {}".format([p for p in gtmContext.constraintTreePaths]))
    
    if gtmContext.mode == "old":
        resultTree = gtm_old.gtmMerge(gtmContext.constraintTrees, gtmContext.startTree)
    
    elif gtmContext.mode in ("convex", "fp"):
        annotateTrees(gtmContext)  
        collapseViolatingEdges(gtmContext)   
        rerootConstraintTrees(gtmContext)
        distributeNodes(gtmContext)
        applyEdgeLengths(gtmContext, None)
        applyNodeLabels(gtmContext, None)
        resultTree = gtmContext.startTree     
    
    endTime = time.time()
    Configs.log("Finished GTM in {} seconds..".format(endTime-startTime))  
    return resultTree   

def annotateTrees(context):
    Configs.log("Annotating trees..")
    context.startTree.taxonMap = {n.taxon.label : n for n in context.startTree.leaf_nodes()}
    for label, node in context.startTree.taxonMap.items():
        node.taxon.taxonSubtree, node.taxon.taxonBitmask = None, 0
    
    for tree in context.constraintTrees:
        tree.rootBipartition = None
        taxa = [n.taxon for n in tree.leaf_nodes()]
        tree.numTaxa = len(taxa)
        startTreeTaxa = [t.label for t in taxa if t.label in context.startTree.taxonMap]
        tree.startTreeBitmask = tree.taxon_namespace.taxa_bitmask(labels = startTreeTaxa)
        
        for taxon in taxa:
            taxon.taxonSubtree = tree
            taxon.taxonBitmask = tree.taxon_namespace.taxon_bitmask(taxon)
            if taxon.label in context.startTree.taxonMap:
                context.startTree.taxonMap[taxon.label].taxon = taxon        
        
        for edge in tree.postorder_edge_iter():
            gtmutils.populateEdgeDesc(edge)
        gtmutils.populateEdgeMaps(tree)   
    
    for edge in context.startTree.postorder_edge_iter():
        gtmutils.populateEdgeDesc(edge)

def collapseViolatingEdges(context):
    Configs.log("Collapsing violations..")
    toRemove = []
    for edge in context.startTree.postorder_edge_iter():
        violatingEdge = False
        subtreeList = []
        
        for subtree, bitmask in edge.desc.items():
            if gtmutils.checkBitmaskViolatesSubtree(bitmask, subtree):
                violatingEdge = True
                break
            elif bitmask != subtree.startTreeBitmask:
                subtreeList.append(subtree)
                if context.mode == "convex" and len(subtreeList) > 1:
                    violatingEdge = True
                    break
        if violatingEdge:
            toRemove.append(edge.head_node)   
             
    treeutils.collapseEdges(toRemove)

def rerootConstraintTrees(context):
    Configs.log("Re-rooting constraint trees..")
    for edge in context.startTree.postorder_edge_iter():
        subtreeList = [s for s, b in edge.desc.items() if b != s.startTreeBitmask]
        for subtree in subtreeList:
            subtree.rootBipartition = edge.desc[subtree]    
    
    for tree in context.constraintTrees:        
        if tree.numTaxa < 3:
            gtmutils.populateEdgeMaps(tree)
            continue
        
        rootEdge = None
        if tree.rootBipartition is None and len(tree.seed_node.child_edges()) > 0:
            rootEdge = tree.seed_node.child_edges()[0]
        elif tree.rootBipartition is not None:
            for edge in tree.edges():
                bitmask = edge.desc[tree]
                if bitmask & tree.startTreeBitmask == tree.rootBipartition or ~bitmask & tree.startTreeBitmask == tree.rootBipartition:
                    rootEdge = edge
                    break
        
        if rootEdge is not None:    
            newRoot = rootEdge.tail_node
            tree.reseed_at(newRoot, update_bipartitions=False,
            collapse_unrooted_basal_bifurcation=False,
            suppress_unifurcations=False)
            for edge in tree.postorder_edge_iter():
                gtmutils.populateEdgeDesc(edge)
            gtmutils.populateEdgeMaps(tree) 

def distributeNodes(context):
    for i, tree in enumerate(context.constraintTrees):
        Configs.log("Applying constraint topology {}..".format(i+1))
        recurseDistributeNodes(context, context.startTree.seed_node, tree.seed_node, tree) 
    context.startTree.suppress_unifurcations()
    context.startTree.collapse_basal_bifurcation()     
        
def recurseDistributeNodes(context, sNode, cNode, cTree):
    stack = [cNode]
    while stack:
        cNode = stack.pop()       
        cMask = cNode.edge.desc[cTree] & cTree.startTreeBitmask
        maskEdges = {e.desc[cTree] & cMask : e for e in sNode.child_edges() if cTree in e.desc and e.desc[cTree] & cMask != 0}
        
        if cMask == 0:
            treeutils.attachNodetoNode(sNode, cNode)
        
        elif cNode.is_leaf():
            continue
        
        elif cNode == cTree.seed_node:
            stack.extend(cNode.child_nodes())
        
        elif cMask in maskEdges and (len(maskEdges[cMask].desc) > 1 or maskEdges[cMask].descUnassigned > 0 or maskEdges[cMask].desc[cTree] != cMask):
            recurseDistributeNodes(context, maskEdges[cMask].head_node, cNode, cTree)
            
        elif not any(len(maskEdges[mask].desc) > 1 or maskEdges[mask].descUnassigned > 0 for mask in maskEdges):
            for mask in maskEdges:
                sNode.remove_child(maskEdges[mask].head_node)
            treeutils.attachNodetoNode(sNode, cNode)    
            
        elif len(maskEdges) < len(sNode.child_edges()):
            newDesc = gtmutils.combineDescMaps([maskEdges[mask].desc for mask in maskEdges])
            if not any(gtmutils.checkBitmaskViolatesSubtree(b, s) for s, b in newDesc.items()):
                newNode = gtmutils.resolvePolytomy(sNode, [maskEdges[m].head_node for m in maskEdges])
                recurseDistributeNodes(context, newNode, cNode, cTree)
            else:
                stack.extend(cNode.child_nodes())     
        else:
            stack.extend(cNode.child_nodes())

def buildTreeKeyEdges(context):
    if context.treeKeyEdges is None:
        for edge in context.startTree.postorder_edge_iter():
            gtmutils.populateEdgeDesc(edge)
        
        context.treeKeyEdges = {}    
        for edge in context.startTree.preorder_edge_iter():
            if edge.tail_node is None:
                continue
            splits = [(s, b) for s, b in edge.desc.items() if b != s.taxon_namespace.all_taxa_bitmask()]
            for tree, bitmask in splits:
                key = gtmutils.buildEdgeKey(tree, bitmask)
                if key in tree.edgeMap:
                    context.treeKeyEdges[tree, key] = context.treeKeyEdges.get((tree, key), []) + [edge]

def applyEdgeLengths(context, edges = None):
    Configs.log("Applying constraint tree edge lengths..")
    buildTreeKeyEdges(context)
    
    for edge in context.startTree.preorder_edge_iter():
        if edges is None or edge in edges:
            edge.length = None
    
    for treeKey, edgeList in context.treeKeyEdges.items():
        tree, key = treeKey
        if len(edgeList) == 1:
            edge = edgeList[0]
            if edges is None or edge in edges:
                sourceEdge, length, label = tree.edgeMap[key]
                edge.length = length

def applyNodeLabels(context, edges = None):
    Configs.log("Applying constraint tree node labels..")
    buildTreeKeyEdges(context)
    
    for edge in context.startTree.preorder_edge_iter():
        if edges is None or edge in edges:
            edge.head_node.label = None
    
    for treeKey, edgeList in context.treeKeyEdges.items():
        tree, key = treeKey
        if len(edgeList) == 1:
            edge = edgeList[0]
            if edges is None or edge in edges:
                sourceEdge, length, label = tree.edgeMap[key]
                edge.head_node.label = label

