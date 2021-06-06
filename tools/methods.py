'''
Created on May 27, 2021

@author: Vlad
'''

import os
import shutil
from tools import external_tools
from configs import Configs
import treeutils
import sequenceutils
import decomposition
import random
import dendropy

def buildTree(method, outputPath, **kwargs):
    workingDir = os.path.join(os.path.dirname(outputPath), "temp_{}_{}".format(method, os.path.basename(outputPath).split(".")[0]))
    
    if method.lower() == "clustal":
        task = external_tools.runClustalOmegaGuideTree(kwargs["alignmentPath"], workingDir, outputPath, Configs.numCores)
        external_tools.runCommand(**task)
        
    elif method.lower() == "fasttree":
        task = external_tools.runFastTree(kwargs["alignmentPath"], workingDir, outputPath, kwargs.get("mode","fast"), kwargs.get("startTreePath"))
        external_tools.runCommand(**task)
        
    elif method.lower() == "raxml":
        model = kwargs.get("model", getRaxmlModel())        
        task = external_tools.runRaxmlNg(kwargs["alignmentPath"], model, workingDir, kwargs.get("startTreePath"), 
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
    
    tree = treeutils.loadTree(treePath)
    tree.resolve_polytomies()
    for edge in tree.edges():
        edge.length = None
    treeutils.writeTree(tree, unoptimizedPath)
    
    taxa = [n.taxon.label for n in tree.leaf_nodes()]
    sequenceutils.writeSubsetsToFiles(alignmentPath, {reducedAlignPath : taxa})   
    
    task = external_tools.runRaxmlNgOptimize(reducedAlignPath, model, unoptimizedPath, workingDir, modelPath, outputPath, Configs.numCores)
    external_tools.runCommand(**task)
    return modelPath

def raxmlGetScore(treePath, alignmentPath, model):
    task = external_tools.runRaxmlNgScore(alignmentPath, model, treePath, os.path.dirname(treePath), Configs.numCores)
    runner = external_tools.runCommand(**task)
    lines = runner.stdout.decode("utf-8").splitlines()
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

def getRaxmlModel():
    if Configs.raxmlModel is None and Configs.raxmlModelSourcePath is not None and Configs.raxmlModelSourcePath != "estimate":
        Configs.log("Estimating RAxML model from {}".format(Configs.raxmlModelSourcePath))
        modelTreePath = os.path.join(os.path.dirname(Configs.raxmlModelSourcePath), "optimized_{}".format(os.path.basename(Configs.raxmlModelSourcePath)))
        modelPath = raxmlEvaluateModelParameters(Configs.raxmlModelSourcePath, Configs.alignmentPath, None, modelTreePath)
        Configs.raxmlModel = raxmlReadModel(modelPath)
        Configs.log("Setting RAxML model to {}".format(Configs.raxmlModel))
    return Configs.raxmlModel

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
        
    