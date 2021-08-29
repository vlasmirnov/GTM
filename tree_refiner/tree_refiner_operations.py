'''
Created on Jul 15, 2021

@author: Vlad
'''

from helpers import treeutils, sequenceutils
import os
import time
from tools import methods
from configs import Configs
from gtm import gtm_operations
from gtm.gtm_context import GtmContext


def refineTree(treePath, subtreePaths, workingDir, polytomyStrategy, branchLengthStrategy, estimateBS, alignmentPath, model):
    startTime = time.time()
    Configs.log("Starting Tree Refiner..")
    
    if model is not None and os.path.exists(model):
        model, score = methods.estimateRaxmlModel(model, alignmentPath)
    
    if polytomyStrategy is not None and polytomyStrategy.lower() != "none":
        polytomyWorkingDir = os.path.join(workingDir, "resolve_polytomies")
        polytomyOutputPath = os.path.join(workingDir, "resolve_polytomies.tre")
        resolvePolytomies(treePath, polytomyStrategy, polytomyWorkingDir, 
                          model, alignmentPath, polytomyOutputPath)
        treePath = polytomyOutputPath
    
    if branchLengthStrategy is not None and branchLengthStrategy.lower() != "none":
        branchLengthWorkingDir = os.path.join(workingDir, "branch_lengths")
        branchLengthOutputPath = os.path.join(workingDir, "resolve_branch_lengths.tre")                
        resolveBranchLengths(treePath, branchLengthStrategy, branchLengthWorkingDir, 
                            model, alignmentPath, branchLengthOutputPath)
        treePath = branchLengthOutputPath
    
    if estimateBS:
        bsWorkingDir = os.path.join(workingDir, "bootstrap")
        bsOutputPath = os.path.join(workingDir, "resolve_bootstrap.tre")       
        resolveBootstrapSupport(treePath, subtreePaths, bsWorkingDir, 
                                model, alignmentPath, bsOutputPath)         
        treePath = bsOutputPath
    
    endTime = time.time()
    Configs.log("Finished Tree Refiner in {} seconds..".format(endTime-startTime))
    return treePath

def resolvePolytomies(inputTreePath, strategy, workingDir, model, alignmentPath, outputTreePath):
    Configs.log("Resolving polytomies using method {}..".format(strategy))
    tree = treeutils.loadTree(inputTreePath)
    incidentEdges = lambda node : [e for e in node.incident_edges() if e.tail_node is not None]
    
    polytomies = [node for node in tree.preorder_internal_node_iter() if len(incidentEdges(node)) > 3]
    Configs.log("Tree has {} leaves and {} polytomies..".format(len(tree.leaf_nodes()), len(polytomies)))
    while len(polytomies) > 0:
        if strategy.lower() == "random":
            tree.resolve_polytomies()
        else:
            if not os.path.exists(workingDir):
                os.makedirs(workingDir)
            
            polyEdges = set()
            for node in polytomies:
                for edge in node.incident_edges():
                    polyEdges.add(edge)    
            Configs.log("Found {} polytomies and {} polytomy-incident edges..".format(len(polytomies), len(polyEdges)))  
            spanningTree = treeutils.extractSpanningTree(tree, Configs.polytomyTreeSize, polyEdges, spanInternalNodes = False)
            unresolvedPath = os.path.join(workingDir, "unresolved_spanning_tree.tre")
            resolvedPath = os.path.join(workingDir, "resolved_spanning_tree.tre")
            spanningAlignPath = os.path.join(workingDir, "alignment_spanning_tree.txt") 
            Configs.log("Using polytomy resolution starting tree with {} leaves and {} internal edges..".format(
                len(spanningTree.leaf_nodes()), len(spanningTree.internal_edges(exclude_seed_edge = True))))
            
            treeutils.writeTree(spanningTree, unresolvedPath)                        
            taxa = [n.taxon.label for n in spanningTree.leaf_nodes()]
            sequenceutils.writeSubsetsToFiles(alignmentPath, {spanningAlignPath : taxa})   
            methods.buildTree(strategy, resolvedPath, alignmentPath = spanningAlignPath, constraintTreePath = unresolvedPath, model = model)
            spanningTree = treeutils.loadTree(resolvedPath)
                        
            gtmContext = GtmContext(startTree = tree, constraintTrees = [spanningTree])
            gtm_operations.annotateTrees(gtmContext)
            gtm_operations.rerootConstraintTrees(gtmContext)  
            gtm_operations.distributeNodes(gtmContext)
            #Configs.log("New tree has {} leaves and {} internal edges..".format(len(tree.leaf_nodes()), len(tree.internal_edges(exclude_seed_edge = True))))
            edges = set(e for e in tree.preorder_edge_iter() if e.length is None or e.length <= 0)
            gtm_operations.applyEdgeLengths(gtmContext, edges)
            
            polytomies = [node for node in tree.preorder_internal_node_iter() if len(incidentEdges(node)) > 3]
            Configs.log("New tree has {} leaves and {} polytomies..".format(len(tree.leaf_nodes()), len(polytomies)))
            
    treeutils.writeTree(tree, outputTreePath)
        
def resolveBranchLengths(inputTreePath, strategy, workingDir, model, alignmentPath, outputTreePath):
    Configs.log("Resolving missing branch lengths using method {}..".format(strategy))
    tree = treeutils.loadTree(inputTreePath)
    
    missingLengths = set(e for e in tree.preorder_edge_iter() if e.tail_node is not None and (e.length is None or e.length <= 0))
    Configs.log("Tree has {} edges with missing or non-positive lengths..".format(len(missingLengths)))
    while len(missingLengths) > 0 and strategy is not None:
        if os.path.exists(strategy):
            spanningTree = treeutils.loadTree(strategy)
        else:
            if not os.path.exists(workingDir):
                os.makedirs(workingDir)
            
            for edge in list(missingLengths):
                for nbr in edge.get_adjacent_edges():
                    missingLengths.add(nbr)
            spanningTree = treeutils.extractSpanningTree(tree, Configs.branchLengthTreeSize, missingLengths)
            unoptimizedPath = os.path.join(workingDir, "unoptimized_spanning_tree.tre")
            optimizedPath = os.path.join(workingDir, "optimized_spanning_tree.tre")
            Configs.log("Using branch length resolution starting tree with {} leaves and {} internal edges..".format(
                len(spanningTree.leaf_nodes()), len(spanningTree.internal_edges(exclude_seed_edge = True))))
            
            treeutils.writeTree(spanningTree, unoptimizedPath)
            methods.raxmlEvaluateModelParameters(unoptimizedPath, alignmentPath, model, optimizedPath)
            spanningTree = treeutils.loadTree(optimizedPath)
            
        applyTreeBranchLengths(tree, [spanningTree])
        missingLengths = set(e for e in tree.preorder_edge_iter() if e.tail_node is not None and (e.length is None or e.length <= 0))
        Configs.log("New tree has {} edges with missing or non-positive lengths..".format(len(missingLengths)))
    treeutils.writeTree(tree, outputTreePath)  

def resolveBootstrapSupport(inputTreePath, subtreePaths, workingDir, model, alignmentPath, outputTreePath):
    Configs.log("Estimating bootstrap support values..")
    tree = treeutils.loadTree(inputTreePath)    
    if subtreePaths is not None and len(subtreePaths) > 0:
        Configs.log("Applying bootstrap support values from {} provided subtrees..".format(len(subtreePaths)))
        pathMap = getBootstrapSupport(subtreePaths, alignmentPath, model, Configs.bootstrapTrees)
        constraintBSTrees = [treeutils.loadTree(bsPath) for bsPath in pathMap.values()]
        applyTreeBootstrapSupportValues(tree, constraintBSTrees)
    
    missingLabels = set(e for e in tree.preorder_edge_iter() if e.is_internal() and e.tail_node is not None and e.head_node.label is None)
    Configs.log("Tree has {} edges with missing labels (BS support values)..".format(len(missingLabels)))
    while len(missingLabels) > 0:
        if not os.path.exists(workingDir):
            os.makedirs(workingDir)
        
        for edge in list(missingLabels):
            for nbr in edge.get_adjacent_edges():
                missingLabels.add(nbr)
        spanningTree = treeutils.extractSpanningTree(tree, Configs.bootstrapTreeSize, missingLabels)
        unoptimizedPath = os.path.join(workingDir, "unoptimized_spanning_tree.tre")
        Configs.log("Using bootstrap support starting tree with {} leaves and {} internal edges..".format(
            len(spanningTree.leaf_nodes()), len(spanningTree.internal_edges(exclude_seed_edge = True))))
        
        treeutils.writeTree(spanningTree, unoptimizedPath)
        bsPath = methods.estimateRaxmlBootstrapSupport(unoptimizedPath, alignmentPath, model, Configs.bootstrapTrees)
        spanningTree = treeutils.loadTree(bsPath)
            
        applyTreeBootstrapSupportValues(tree, [spanningTree])
        missingLabels = set(e for e in tree.preorder_edge_iter() if e.is_internal() and e.tail_node is not None and e.head_node.label is None)
        Configs.log("New tree has {} edges with missing labels (BS support values)..".format(len(missingLabels)))
    treeutils.writeTree(tree, outputTreePath)  
   
def getBootstrapSupport(treePaths, alignmentPath, model, bsTrees, edgeCollapseThreshold = 0):
    pathMap = {}
    for treePath in treePaths:
        bsPath = methods.estimateRaxmlBootstrapSupport(treePath, alignmentPath, model, bsTrees)
        pathMap[treePath] = bsPath
    
    return pathMap

def applyTreeBranchLengths(startTree, trees):
    gtmContext = GtmContext(startTree = startTree, constraintTrees = trees)
    gtm_operations.annotateTrees(gtmContext)
    gtm_operations.rerootConstraintTrees(gtmContext)    
    edges = set(e for e in startTree.preorder_edge_iter() if e.length is None or e.length <= 0)
    gtm_operations.applyEdgeLengths(gtmContext, edges)

def applyTreeBootstrapSupportValues(startTree, trees):
    gtmContext = GtmContext(startTree = startTree, constraintTrees = trees)
    gtm_operations.annotateTrees(gtmContext)
    gtm_operations.rerootConstraintTrees(gtmContext)    
    edges = set(e for e in startTree.preorder_edge_iter() if e.head_node.label is None)
    gtm_operations.applyNodeLabels(gtmContext, edges)
    
'''
def labelScaffold(tree):
    for edge in tree.preorder_edge_iter(): 
        edge.head_node.label = None
        
    for edge in tree.preorder_edge_iter():
        if edge.tail_node is None:
            continue  
        splits = [s for s, b in edge.desc.items() if b != s.taxon_namespace.all_taxa_bitmask()]
        if len(splits) != 1:
            edge.head_node.label = "scaffold"
            for e in edge.get_adjacent_edges():
                if e.head_node.label != "scaffold":
                    e.head_node.label = "scaffold_adjacent"
'''