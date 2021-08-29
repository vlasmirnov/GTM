'''
Created on May 27, 2021

@author: Vlad
'''

import os
import shutil
from tools import external_tools
from configs import Configs
from helpers import treeutils
from helpers import sequenceutils
from tools import decomposition
import random
import dendropy
import math

def buildTree(method, outputPath, **kwargs):
    workingDir = os.path.join(os.path.dirname(outputPath), "temp_{}_{}".format(method, os.path.basename(outputPath).split(".")[0]))
    if os.path.exists(workingDir):
        shutil.rmtree(workingDir)
    
    Configs.log("Building {} tree {}".format(method.upper(), outputPath))
    
    if method.lower() == "clustal":
        task = external_tools.runClustalOmegaGuideTree(kwargs["alignmentPath"], workingDir, outputPath, Configs.numCores)
        external_tools.runCommand(**task)
        
    elif method.lower() == "fasttree":
        task = external_tools.runFastTree(kwargs["alignmentPath"], workingDir, outputPath, kwargs.get("mode","normal"), kwargs.get("startTreePath"))
        external_tools.runCommand(**task)
    
    elif method.lower() == "iqtree":
        task = external_tools.runIqTree(kwargs["alignmentPath"], kwargs.get("model"), workingDir, kwargs.get("startTreePath"), 
                                         kwargs.get("constraintTreePath"), outputPath, Configs.numCores)
        external_tools.runCommand(**task)
        
    elif method.lower() == "raxml":
        task = external_tools.runRaxmlNg(kwargs["alignmentPath"], kwargs.get("model"), workingDir, kwargs.get("startTreePath"), 
                                         kwargs.get("constraintTreePath"), outputPath, Configs.numCores)
        external_tools.runCommand(**task)
        
    elif method.lower() == "random":
        taxa = sequenceutils.readTaxaFromFasta(kwargs["alignmentPath"])
        tree = treeutils.createRandomTree(taxa)
        treeutils.writeTree(tree, outputPath)
        
    elif method.lower() == "treemer":
        buildTreemerTree(workingDir, kwargs["alignmentPath"], 10, outputPath)
    else:
        raise Exception("Tree estimation method {} not recognized..".format(method))
    
    
def raxmlEvaluateModelParameters(treePath, alignmentPath, model, outputPath):
    baseName = os.path.basename(treePath).split(".")[0]
    workingDir = os.path.join(os.path.dirname(treePath), "model_estimation_{}".format(baseName))
    if os.path.exists(workingDir):
        shutil.rmtree(workingDir)
    os.makedirs(workingDir)    
    
    reducedAlignPath = os.path.join(workingDir, "alignment_{}.txt".format(baseName))    
    unoptimizedPath = os.path.join(workingDir, "unoptimized_{}".format(os.path.basename(treePath)))
    modelPath = os.path.join(workingDir, "model_{}.txt".format(baseName))
    logPath = os.path.join(workingDir, "raxmllog_{}.txt".format(baseName))
    
    tree = treeutils.loadTree(treePath)
    tree.resolve_polytomies()
    tree.suppress_unifurcations()
    for edge in tree.edges():
        edge.length = None
    treeutils.writeTree(tree, unoptimizedPath)
    
    taxa = [n.taxon.label for n in tree.leaf_nodes()]
    sequenceutils.writeSubsetsToFiles(alignmentPath, {reducedAlignPath : taxa})   
    
    task = external_tools.runRaxmlNgOptimize(reducedAlignPath, model, unoptimizedPath, workingDir, 
                                             logPath, modelPath, outputPath, Configs.numCores)
    external_tools.runCommand(**task)
    return modelPath, logPath

def raxmlGetBootstrap(treePath, alignmentPath, model, bsTrees, outputPath):
    baseName = os.path.basename(treePath).split(".")[0]
    workingDir = os.path.join(os.path.dirname(treePath), "bootstrapping_{}".format(baseName))
    if os.path.exists(workingDir):
        shutil.rmtree(workingDir)
    os.makedirs(workingDir)
    
    task = external_tools.runRaxmlNgBootstrap(alignmentPath, model, treePath, bsTrees, workingDir, outputPath, Configs.numCores)
    external_tools.runCommand(**task)
    return outputPath

def raxmlGetScore(treePath, alignmentPath, model, clean = False):
    if clean:        
        tree = treeutils.loadTree(treePath)
        tree.resolve_polytomies()
        tree.suppress_unifurcations()
        for edge in tree.edges():
            if edge.length is None or edge.length < 0:
                edge.length = 0
        treePath = os.path.join(os.path.dirname(treePath), "clean_{}".format(os.path.basename(treePath)))
        treeutils.writeTree(tree, treePath)
    
    task = external_tools.runRaxmlNgScore(alignmentPath, model, treePath, os.path.dirname(treePath), Configs.numCores)
    runner = external_tools.runCommand(**task)
    lines = runner.stdout.splitlines()
    for line in lines:
        if line.startswith("Final LogLikelihood: "):
            score = float(line.split("Final LogLikelihood: ")[1])
            return score

def raxmlReadModel(modelPath):
    with open(modelPath) as file:
        firstLine = file.readline().strip()
        modelToken = firstLine.split(",")[0]
        return modelToken

def raxmlReadScore(logPath):
    with open(logPath) as file:
        for line in file:
            if line.startswith("Final LogLikelihood: "):
                score = float(line.split("Final LogLikelihood: ")[1])
                return score

def estimateRaxmlBootstrapSupport(treePath, alignmentPath, model, bsTrees = Configs.bootstrapTrees):
    Configs.log("Estimating bootstrap supports with {} replicates for {}".format(bsTrees, treePath))
    baseName = os.path.basename(treePath).split(".")[0]
    workingDir = os.path.join(os.path.dirname(treePath), "bootstrapping_{}".format(baseName))
    if os.path.exists(workingDir):
        shutil.rmtree(workingDir)
    os.makedirs(workingDir)
    
    bsTreesPath = os.path.join(workingDir, "bstrees_{}".format(os.path.basename(treePath)))
    outputPath = os.path.join(workingDir, "support_{}".format(os.path.basename(treePath)))
    reducedAlignPath = os.path.join(workingDir, "alignment_{}.txt".format(baseName))    
    unoptimizedPath = os.path.join(workingDir, "unoptimized_{}".format(os.path.basename(treePath)))
    
    tree = treeutils.loadTree(treePath)
    tree.resolve_polytomies()
    tree.suppress_unifurcations()
    for edge in tree.edges():
        edge.length = None
    treeutils.writeTree(tree, unoptimizedPath)
    taxa = [n.taxon.label for n in tree.leaf_nodes()]
    sequenceutils.writeSubsetsToFiles(alignmentPath, {reducedAlignPath : taxa}) 
    
    task = external_tools.runRaxmlNgBootstrap(reducedAlignPath, model, unoptimizedPath, bsTrees, 
                                              workingDir, bsTreesPath, Configs.numCores)
    external_tools.runCommand(**task)
    
    task = external_tools.runRaxmlNgBootstrapSupport(reducedAlignPath, model, unoptimizedPath, bsTreesPath, 
                                              workingDir, outputPath, Configs.numCores)
    external_tools.runCommand(**task)
    
    return outputPath

def estimateRaxmlModel(treePath, alignmentPath, maxLeaves = None):
    try:
        tree = treeutils.loadTree(treePath)
        Configs.log("Found valid newick in {}".format(treePath))
    except:
        with open(treePath) as f:
            firstLine = f.readline()
            model = firstLine.split(',')[0]
        Configs.log("Found model {} in file {}".format(model, treePath))
        return model, None
    
    if maxLeaves is not None:
        taxa = [n.taxon.label for n in tree.postorder_node_iter() if n.taxon is not None]
        if len(taxa) > maxLeaves:
            n = math.ceil(len(taxa) / maxLeaves)        
            sample = taxa[random.randint(0,n-1) :: n]
            tree = tree.extract_tree_with_taxa_labels(sample)
        treePath = os.path.join(os.path.dirname(treePath), "reduced_{}".format(os.path.basename(treePath)))
        treeutils.writeTree(tree, treePath)
    
    Configs.log("Estimating RAxML model from {}".format(treePath))
    modelTreePath = os.path.join(os.path.dirname(treePath), "optimized_{}".format(os.path.basename(treePath)))
    modelPath, logPath = raxmlEvaluateModelParameters(treePath, alignmentPath, None, modelTreePath)
    model = raxmlReadModel(modelPath)
    score = raxmlReadScore(logPath)
    Configs.log("Estimated model {} with score {}".format(model, score))
    return model, score

def extractInducedTree(workingDir, startTreePath, alignmentPath, inducedTreePath = None):
    if inducedTreePath is None:
        baseName = "induced_{}_{}".format(os.path.basename(alignmentPath).replace(".", "_"), os.path.basename(startTreePath))
        inducedTreePath = os.path.join(workingDir, baseName)
    taxa = sequenceutils.readTaxaFromFasta(alignmentPath)
    tree = treeutils.loadTree(startTreePath)
    inducedTree = tree.extract_tree_with_taxa_labels(taxa)
    inducedTree.resolve_polytomies()
    for edge in inducedTree.edges():
        edge.length = None
    treeutils.writeTree(inducedTree, inducedTreePath)
    return inducedTreePath

def buildTreemerTree(workingDir, alignmentPath, numSubsets, outputPath):
    if os.path.exists(workingDir):
        shutil.rmtree(workingDir)
    os.makedirs(workingDir)
    subsetsDir = os.path.join(workingDir, "subsets")
    subtreesDir = os.path.join(workingDir, "subtrees")
    os.makedirs(subtreesDir)
    
    taxa = sequenceutils.readTaxaFromFasta(alignmentPath)
    random.shuffle(taxa)
    subsets = [taxa[i :: numSubsets] for i in range(numSubsets)]
    subsetPaths = sequenceutils.writeSubsetsToDir(subsetsDir, alignmentPath, subsets)
    subtreePaths = {}
    for subsetPath, subset in subsetPaths.items():
        baseName = os.path.basename(subsetPath).split(".")[0]
        subTreePath = os.path.join(subtreesDir, "subtree_{}.tre".format(baseName))
        subtreePaths[subTreePath] = subset        
        buildTree("raxml", subTreePath, alignmentPath = subsetPath)
    
    assigns = {}
    skeletonTaxa = {}
    for treePath in subtreePaths:
        tree = treeutils.loadTree(treePath)
        tree, labelNodeMap, subAssigns = decomposition.treemer(tree, int(Configs.decompositionMaxSubsetSize / numSubsets))
        skeletonTaxa.update(labelNodeMap)
        assigns.update(subAssigns)
    
    skeletonPath = os.path.join(workingDir, "skeleton_taxa.txt")
    skeletonTreePath = os.path.join(workingDir, "skeleton.tre")
    sequenceutils.writeSubsetsToFiles(alignmentPath, {skeletonPath : skeletonTaxa})
    buildTree("raxml", skeletonTreePath, alignmentPath = skeletonPath)
    skeletonTree = treeutils.loadTree(skeletonTreePath)
    
    skeletonLabelNodes = {node.taxon.label : node for node in skeletonTree.leaf_nodes()}       
    for tax1, tax2 in assigns.items():
        newNode = dendropy.Node(taxon = dendropy.Taxon(tax1))
        skeletonLabelNodes[tax2].parent_node.add_child(newNode)
    treeutils.writeTree(skeletonTree, outputPath)
        
    