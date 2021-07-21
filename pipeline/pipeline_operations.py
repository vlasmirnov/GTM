'''
Created on Jul 9, 2021

@author: Vlad
'''

import shutil
from helpers import treeutils, sequenceutils
import os
from tools import methods, decomposition
from configs import Configs
from gtm.gtm_context import GtmContext
from gtm import gtm_operations

def iterateGtm(context, workingDir, outputPath):
    if os.path.exists(workingDir):
        shutil.rmtree(workingDir)
    os.makedirs(workingDir)
    
    if len(context.constraintTreePaths) == 0:
        buildConstraintTrees(context, workingDir)
    
    if context.guideTreeStrategy == "hierarchal":
        buildHierarchalGuideTrees(context, workingDir)  
        resultTree = treeutils.loadTree(list(context.subtreesList[0].keys())[0])
        i = 1
        for i in range(1, len(context.subtreesList)):
            gtmContext = GtmContext(startTree = resultTree, constraintTreePaths = context.subtreesList[i], mode = context.mode)
            resultTree = gtm_operations.runGtm(gtmContext)
    else:
        if context.guideTreePath is None:
            buildGuideTree(context, workingDir)             
        gtmContext = GtmContext(startTreePath = context.guideTreePath, constraintTreePaths = context.constraintTreePaths, mode = context.mode)
        resultTree = gtm_operations.runGtm(gtmContext)
    
    context.outputPath = outputPath    
    treeutils.writeTree(resultTree, context.outputPath)
    context.startTreePath = context.outputPath
    context.currentTreePath = context.outputPath
    context.constraintTreePaths = []
    context.guideTreePath = None
    context.iterations = context.iterations + 1  

def buildConstraintTrees(context, workingDir):
    if context.startTreePath is None:
        Configs.log("Building starting tree..")
        context.startTreePath = os.path.join(workingDir, "start_tree.tre")
        updateModel(context)
        methods.buildTree(context.startTreeMethod, context.startTreePath, alignmentPath = context.alignmentPath, model = context.model)
        context.currentTreePath = context.startTreePath
    
    Configs.log("Building constraint trees..")    
    context.subsetsDir = os.path.join(workingDir, "subsets")
    context.subtreesDir = os.path.join(workingDir, "subtrees")
    context.subsetPaths = decomposition.decomposeSequences(context, workingDir)        
    updateModel(context)
    context.constraintTreePaths = buildSubtrees(context.subtreesDir, context.getStartTreeForML(), 
                                                context.subsetPaths, Configs.treeMethod, context.model)
        
def buildSubtrees(subtreesDir, startTreePath, subsetPaths, method, model):
    if os.path.exists(subtreesDir):
        shutil.rmtree(subtreesDir)
    os.makedirs(subtreesDir)
    
    subtreePaths = {}
    for subsetPath, subset in subsetPaths.items():
        baseName = os.path.basename(subsetPath).split(".")[0]
        subTreePath = os.path.join(subtreesDir, "subtree_{}.tre".format(baseName))
        subtreePaths[subTreePath] = subset   
        inducedStartTree = None
        if startTreePath is not None:
            inducedStartTree = methods.extractInducedTree(subtreesDir, startTreePath, subsetPath)
        methods.buildTree(method, subTreePath, alignmentPath = subsetPath, startTreePath = inducedStartTree, model = model)
    
    return subtreePaths  

def buildHierarchalGuideTrees(context, workingDir):
    context.subtreesList = [context.constraintTreePaths]
    subsetsDir, subtreesDir, subsetPaths = context.subsetsDir, context.subtreesDir, context.subsetPaths
    constraintTreePaths = context.constraintTreePaths
    s = 1
    while len(subsetPaths) > 1:
        orderedSubtrees = [os.path.join(subtreesDir, "subtree_subset_{}.tre".format(n+1)) for n in range(len(subsetPaths))]
        subsetsDir = os.path.join(workingDir, "subsets_{}".format(s))
        subtreesDir = os.path.join(workingDir, "subtrees_{}".format(s))
        subsets = []
        if len(subsetPaths) % 2 == 1:
            subsets.append(constraintTreePaths[orderedSubtrees[-1]])
        for i in range(int(len(subsetPaths)/2)):
            tree1, tree2 = treeutils.loadTree(orderedSubtrees[2*i]), treeutils.loadTree(orderedSubtrees[2*i+1])
            taxa1 = [n.taxon.label for n in tree1.postorder_node_iter() if n.taxon is not None]
            taxa2 = [n.taxon.label for n in tree2.postorder_node_iter() if n.taxon is not None]
            subset = taxa1[0 :: 2] + taxa2[0 :: 2]
            subsets.append(subset)
        subsetPaths = sequenceutils.writeSubsetsToDir(subsetsDir, context.alignmentPath, subsets)
        constraintTreePaths = buildSubtrees(subtreesDir, context.startTreePath, subsetPaths, Configs.treeMethod, context.model)
        context.subtreesList.append(constraintTreePaths)
        s = s + 1
    context.subtreesList = list(reversed(context.subtreesList))  

def buildGuideTree(context, workingDir):
    Configs.log("Building guide tree..")
    if context.guideTreeStrategy is None:
        Configs.log("No guide tree or strategy passed in, using the starting tree..")
        context.guideTreePath = context.startTreePath
    elif context.guideTreeStrategy == "span":
        context.guideTreePath = os.path.join(workingDir, "guide_tree.tre")
        decomposition.guideTreeSpan(workingDir, context.guideTreePath, startTreePath = context.getStartTreeForML(), model = context.model,
                                            alignmentPath = context.alignmentPath, subsetPaths = context.subsetPaths)
    elif context.guideTreeStrategy == "span_diverse":
        context.guideTreePath = os.path.join(workingDir, "guide_tree.tre")
        decomposition.guideTreeSpanDiverse(workingDir, context.guideTreePath, startTreePath = context.getStartTreeForML(), model = context.model,
                                                   alignmentPath = context.alignmentPath, subtreePaths = context.constraintTreePaths)     

def updateModel(context):
    if context.modelSourcePath == "estimate":
        if context.currentTreePath is not None:
            Configs.log("Updating global RAxML model from {}..".format(context.currentTreePath))
            context.model = methods.estimateRaxmlModel(context.currentTreePath, context.alignmentPath, Configs.raxmlModelLimit)
    elif context.modelSourcePath is not None and context.model is None:
        Configs.log("Updating global RAxML model from {}..".format(context.modelSourcePath))
        context.model = methods.estimateRaxmlModel(context.modelSourcePath, context.alignmentPath)
