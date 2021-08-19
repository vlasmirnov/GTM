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
from pipeline.pipeline_context import PipelineContext

def iterateGtm(context, workingDir, outputPath):
    if not os.path.exists(workingDir):
        os.makedirs(workingDir)
    #    shutil.rmtree(workingDir)
    
    if len(context.constraintTreePaths) == 0:
        buildConstraintTrees(context, workingDir)
    
    if context.guideTreePath is None:
        buildGuideTree(context, workingDir)             
    gtmContext = GtmContext(startTreePath = context.guideTreePath, constraintTreePaths = context.constraintTreePaths, mode = Configs.mode)
    resultTree = gtm_operations.runGtm(gtmContext)
    
    context.currentTreePath = outputPath
    treeutils.writeTree(resultTree, outputPath)

def buildConstraintTrees(context, workingDir):
    if context.startTreePath is None:
        updateModel(context)
        Configs.log("Building starting tree..")
        context.startTreePath = os.path.join(workingDir, "start_tree.tre")
        methods.buildTree(Configs.startTreeMethod, context.startTreePath, alignmentPath = context.alignmentPath, model = context.model)
    
    if context.topLevelStartTreePath is None:
        context.topLevelStartTreePath = context.startTreePath
    
    if context.currentTreePath != context.startTreePath:
        context.currentTreePath = context.startTreePath  
        updateModel(context)        
        logMLScore(context, clean = True)
    else:
        updateModel(context)
    
    Configs.log("Building constraint trees..")        
    context.subsetPaths = decomposition.decomposeSequences(context, workingDir)  
    context.constraintTreePaths = buildSubtrees(context, workingDir)

def buildSubtrees(context, workingDir):
    subtreesDir = os.path.join(workingDir, "subtrees")  
    if os.path.exists(subtreesDir):
        shutil.rmtree(subtreesDir)
    os.makedirs(subtreesDir)
    
    subtreePaths = {}
    for subsetPath, subset in context.subsetPaths.items():
        baseName = os.path.basename(subsetPath).split(".")[0]
        subTreePath = os.path.join(subtreesDir, "subtree_{}.tre".format(baseName))
        subtreePaths[subTreePath] = subset   
        inducedStartTree = None
        if context.getStartTreeForML() is not None:
            inducedStartTree = methods.extractInducedTree(subtreesDir, context.getStartTreeForML(), subsetPath)
        methods.buildTree(Configs.treeMethod, subTreePath, alignmentPath = subsetPath, 
                          startTreePath = inducedStartTree, model = context.model)
    
    return subtreePaths  

def buildGuideTree(context, workingDir):
    Configs.log("Building guide tree..")
    if Configs.guideTreeStrategy is None:
        buildFullGuideTree(context, workingDir)
    elif Configs.guideTreeStrategy == "span":
        buildPartialGuideTree(context, workingDir)
    elif Configs.guideTreeStrategy == "recursive":
        buildRecursiveGuideTree(context, workingDir) 

def buildFullGuideTree(context, workingDir):
    Configs.log("Building guide tree over full set of sequences..")
    if Configs.guideTreeMethod is None:
        Configs.log("No guide tree or method passed in, using the starting tree..")
        context.guideTreePath = context.startTreePath
    else:
        context.guideTreePath = os.path.join(workingDir, "guide_tree.tre")
        methods.buildTree(Configs.guideTreeMethod, context.guideTreePath, alignmentPath = context.alignmentPath, 
                      startTreePath = context.getStartTreeForML(), model = context.model)
        
def buildRecursiveGuideTree(context, workingDir):
    Configs.log("Building recursive guide tree..")
    context.guideTreePath = os.path.join(workingDir, "guide_tree.tre")
    
    taxa = []       
    taxaSets = []
    for path in context.constraintTreePaths:
        subTree = treeutils.loadTree(path)
        taxonSet = [n.taxon.label for n in subTree.postorder_node_iter() if n.taxon is not None]
        taxa.extend(taxonSet)
        taxaSets.append(taxonSet)
    numGuideLeaves = max(Configs.guideTreeRecursionBaseSize, int(len(taxa) / Configs.guideTreeRecursionFactor))
    guideLeaves = decomposition.sampleEvenly(taxa, numGuideLeaves, values = context.getTaxaLengths())
    #guideLeaves = decomposition.sampleLongestRandom6(taxaSets, numGuideLeaves, values = context.getTaxaLengths())
    Configs.log("Building recursive guide tree over {} leaves..".format(len(guideLeaves)))
    print(sum([context.getTaxaLengths()[t] for t in taxa])/len(taxa))
    print(sum([context.getTaxaLengths()[t] for t in guideLeaves])/len(guideLeaves))
    
    if len(guideLeaves) <= Configs.guideTreeRecursionBaseSize:
        if Configs.guideTreeMethod is None:
            Configs.log("No guide tree or method passed in, using the starting tree..")
            context.guideTreePath = context.topLevelStartTreePath
        else:
            alignPath = os.path.join(workingDir, "guide_tree_align.txt")
            sequenceutils.writeSubsetsToFiles(context.alignmentPath, {alignPath : guideLeaves})
            startTreePath = context.getStartTreeForML()
            if startTreePath is not None:
                startTreePath = methods.extractInducedTree(workingDir, startTreePath, alignPath)    
            methods.buildTree(Configs.guideTreeMethod, context.guideTreePath, alignmentPath = alignPath, 
                              startTreePath = startTreePath, model = context.model)
    else:
            newWorkingDir = os.path.join(workingDir, "guide_tree")
            newStartTreePath = os.path.join(newWorkingDir, "start_tree.tre")
            newContext = PipelineContext(workingDir = newWorkingDir,
                                         alignmentPath = context.alignmentPath,
                                         startTreePath = newStartTreePath,
                                         model = context.model,
                                         taxaLengths = context.taxaLengths,
                                         topLevelStartTreePath = context.topLevelStartTreePath)
            startTree = treeutils.loadTree(context.startTreePath)
            newStartTree = startTree.extract_tree_with_taxa_labels(guideLeaves)
            treeutils.writeTree(newStartTree, newStartTreePath)
            iterateGtm(newContext, newWorkingDir, context.guideTreePath)

def buildPartialGuideTree(context, workingDir):
    Configs.log("Building partial guide tree..")
    context.guideTreePath = os.path.join(workingDir, "guide_tree.tre")
    
    taxa = []       
    for path in context.constraintTreePaths:
        subTree = treeutils.loadTree(path)
        taxa.extend([n.taxon.label for n in subTree.postorder_node_iter() if n.taxon is not None])
    leaves = decomposition.sampleEvenly(taxa, Configs.decompositionMaxSubsetSize, values = context.getTaxaLengths())
    alignPath = os.path.join(workingDir, "guide_tree_align.txt")
    sequenceutils.writeSubsetsToFiles(context.alignmentPath, {alignPath : leaves})
    startTreePath = context.getStartTreeForML()
    if startTreePath is not None:
        startTreePath = methods.extractInducedTree(workingDir, startTreePath, alignPath)    
    methods.buildTree(Configs.guideTreeMethod, context.guideTreePath, alignmentPath = alignPath, 
                      startTreePath = startTreePath, model = context.model)

def updateModel(context):
    if context.modelSourcePath == "estimate":
        if context.currentTreePath is not None:
            Configs.log("Updating global RAxML model from {}..".format(context.currentTreePath))
            context.model, score = methods.estimateRaxmlModel(context.currentTreePath, context.alignmentPath, Configs.raxmlModelLimit)
    elif context.modelSourcePath is not None and context.model is None:
        Configs.log("Updating global RAxML model from {}..".format(context.modelSourcePath))
        context.model, score = methods.estimateRaxmlModel(context.modelSourcePath, context.alignmentPath)
        
def logMLScore(context, clean = False):
    if context.trackMLScores:
        Configs.log("Scoring current tree: {}".format(context.currentTreePath))
        #score = methods.raxmlGetScore(context.currentTreePath, context.alignmentPath, context.model, clean)
        model, score = methods.estimateRaxmlModel(context.currentTreePath, context.alignmentPath, Configs.raxmlModelLimit)
        Configs.log("Current ML score: {}".format(score))

        
        
'''
def buildSubtrees(context, subtreesDir, subsetPaths):
    if os.path.exists(subtreesDir):
        shutil.rmtree(subtreesDir)
    os.makedirs(subtreesDir)
    
    subtreePaths = {}
    for subsetPath, subset in subsetPaths.items():
        baseName = os.path.basename(subsetPath).split(".")[0]
        subTreePath = os.path.join(subtreesDir, "subtree_{}.tre".format(baseName))
        subtreePaths[subTreePath] = subset   
        
        taxa = sequenceutils.readTaxaFromFasta(subsetPath)
        if len(taxa) <= Configs.decompositionMaxSubsetSize:
            inducedStartTree = None
            if context.getStartTreeForML() is not None:
                inducedStartTree = methods.extractInducedTree(subtreesDir, context.getStartTreeForML(), subsetPath)
            methods.buildTree(Configs.treeMethod, subTreePath, alignmentPath = subsetPath, startTreePath = inducedStartTree, model = context.model)
        else:
            newWorkingDir = os.path.join(subtreesDir, "subtree_{}".format(baseName))
            newStartTreePath = os.path.join(newWorkingDir, "start_tree.tre")
            newContext = PipelineContext(workingDir = newWorkingDir,
                                         alignmentPath = subsetPath,
                                         startTreePath = newStartTreePath,
                                         guideTreeStrategy = context.guideTreeStrategy,
                                         decompositionStrategy = context.decompositionStrategy,
                                         mode = context.mode,
                                         model = context.model,
                                         useInducedStartTreeForML = context.useInducedStartTreeForML,
                                         guideTreeRecursionFactor = context.guideTreeRecursionFactor,
                                         taxaLengths = context.taxaLengths)
            startTree = treeutils.loadTree(context.startTreePath)
            newStartTree = startTree.extract_tree_with_taxa_labels(taxa)
            treeutils.writeTree(newStartTree, newStartTreePath)
            iterateGtm(newContext, newWorkingDir, subTreePath)
    
    return subtreePaths
'''
