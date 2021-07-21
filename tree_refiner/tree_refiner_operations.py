'''
Created on Jul 15, 2021

@author: Vlad
'''

from helpers import treeutils, sequenceutils
import os
import shutil
import time
from tools import methods
from configs import Configs
from gtm import gtm_operations
from gtm.gtm_context import GtmContext


def refineTree(treePath, workingDir, polytomyStrategy, branchLengthStrategy, alignmentPath, model):
    startTime = time.time()
    Configs.log("Starting Tree Refiner..")
    
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
        
    endTime = time.time()
    Configs.log("Finished Tree Refiner in {} seconds..".format(endTime-startTime))
    return treePath

def resolvePolytomies(inputTreePath, strategy, workingDir, model, alignmentPath, outputTreePath):
    Configs.log("Resolving polytomies using method {}..".format(strategy))
    tree = treeutils.loadTree(inputTreePath)
    incidentEdges = lambda node : [e for e in node.incident_edges() if e.tail_node is not None]
    polytomies = [node for node in tree.preorder_internal_node_iter() if len(incidentEdges(node)) > 3]
    Configs.log("Tree has {} leaves and {} internal edges..".format(len(tree.leaf_nodes()), len(tree.internal_edges(exclude_seed_edge = True))))
    if len(polytomies) > 0:
        if strategy.lower() == "random":
            tree.resolve_polytomies()
        else:
            if os.path.exists(workingDir):
                shutil.rmtree(workingDir)
            os.makedirs(workingDir)
            
            polyEdges = set()
            for node in polytomies:
                for edge in node.incident_edges():
                    polyEdges.add(edge)    
            Configs.log("Found {} polytomies and {} polytomy-incident edges..".format(len(polytomies), len(polyEdges)))  
            spanningTree = treeutils.extractSpanningTree(tree, Configs.spanningTreeSize, polyEdges)
            unresolvedPath = os.path.join(workingDir, "unresolved_spanning_tree.tre")
            resolvedPath = os.path.join(workingDir, "resolved_spanning_tree.tre")
            spanningAlignPath = os.path.join(workingDir, "alignment_spanning_tree.txt") 
            Configs.log("Using polytomy resolution starting tree with {} leaves and {} internal edges..".format(
                len(spanningTree.leaf_nodes()), len(spanningTree.internal_edges(exclude_seed_edge = True))))
            
            treeutils.writeTree(spanningTree, unresolvedPath)                        
            taxa = [n.taxon.label for n in spanningTree.leaf_nodes()]
            sequenceutils.writeSubsetsToFiles(alignmentPath, {spanningAlignPath : taxa})   
            if model is not None and os.path.exists(model):
                model = methods.estimateRaxmlModel(model, alignmentPath)
            methods.buildTree(strategy, resolvedPath, alignmentPath = spanningAlignPath, constraintTreePath = unresolvedPath, model = model)
            spanningTree = treeutils.loadTree(resolvedPath)
            Configs.log("Resolved polytomy resolution starting tree has {} leaves and {} internal edges..".format(
                len(spanningTree.leaf_nodes()), len(spanningTree.internal_edges(exclude_seed_edge = True))))
            
            gtmContext = GtmContext(startTree = tree, constraintTrees = [spanningTree])
            gtm_operations.annotateTrees(gtmContext)
            gtm_operations.rerootConstraintTrees(gtmContext)  
            gtm_operations.distributeNodes(gtmContext)
            Configs.log("Resolved tree has {} leaves and {} internal edges..".format(len(tree.leaf_nodes()), len(tree.internal_edges(exclude_seed_edge = True))))
            
            edges = set(e for e in tree.preorder_edge_iter() if e.length is None or e.length <= 0)
            gtm_operations.applyEdgeLengths(gtmContext, edges)
            
    treeutils.writeTree(tree, outputTreePath)
        
def resolveBranchLengths(inputTreePath, strategy, workingDir, model, alignmentPath, outputTreePath):
    Configs.log("Resolving missing branch lengths using method {}..".format(strategy))
    tree = treeutils.loadTree(inputTreePath)
    missingLengths = set(e for e in tree.preorder_edge_iter() if e.tail_node is not None and (e.length is None or e.length <= 0))
    Configs.log("Tree has {} edges with missing or non-positive lengths..".format(len(missingLengths)))
    if len(missingLengths) > 0 and strategy is not None:
        if os.path.exists(strategy):
            spanningTree = treeutils.loadTree(strategy)
        else:
            if os.path.exists(workingDir):
                shutil.rmtree(workingDir)
            os.makedirs(workingDir)
            
            for edge in list(missingLengths):
                for nbr in edge.get_adjacent_edges():
                    missingLengths.add(nbr)
            spanningTree = treeutils.extractSpanningTree(tree, Configs.spanningTreeSize, missingLengths)
            unoptimizedPath = os.path.join(workingDir, "unoptimized_spanning_tree.tre")
            optimizedPath = os.path.join(workingDir, "optimized_spanning_tree.tre")
            Configs.log("Using branch length resolution starting tree with {} leaves and {} internal edges..".format(
                len(spanningTree.leaf_nodes()), len(spanningTree.internal_edges(exclude_seed_edge = True))))
            
            treeutils.writeTree(spanningTree, unoptimizedPath)
            if model is not None and os.path.exists(model):
                model = methods.estimateRaxmlModel(model, alignmentPath)
            methods.raxmlEvaluateModelParameters(unoptimizedPath, alignmentPath, model, optimizedPath)
            spanningTree = treeutils.loadTree(optimizedPath)
            
        applySpanningTreeBranchLengths(tree, spanningTree)
        missingLengths = set(e for e in tree.preorder_edge_iter() if e.tail_node is not None and (e.length is None or e.length <= 0))
        Configs.log("Refined tree has {} edges with missing or non-positive lengths..".format(len(missingLengths)))
    treeutils.writeTree(tree, outputTreePath)  

def applySpanningTreeBranchLengths(startTree, spanningTree):
    gtmContext = GtmContext(startTree = startTree, constraintTrees = [spanningTree])
    gtm_operations.annotateTrees(gtmContext)
    gtm_operations.rerootConstraintTrees(gtmContext)    
    edges = set(e for e in startTree.preorder_edge_iter() if e.length is None or e.length <= 0)
    gtm_operations.applyEdgeLengths(gtmContext, edges)

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