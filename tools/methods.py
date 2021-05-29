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

def buildTree(method, outputPath, **kwargs):
    workingDir = os.path.join(os.path.dirname(outputPath), "method_" + os.path.basename(outputPath).split(".")[0])
    
    if method.lower() == "clustal":
        task = external_tools.runClustalOmegaGuideTree(kwargs["alignmentPath"], workingDir, outputPath, Configs.numCores)
        
    elif method.lower() == "fasttree":
        task = external_tools.runFastTree(kwargs["alignmentPath"], workingDir, outputPath, kwargs.get("mode","fast"), kwargs.get("startTreePath"))
    
    elif method.lower() == "raxml":
        model = kwargs.get("model", getRaxmlModel())        
        task = external_tools.runRaxmlNg(kwargs["alignmentPath"], model, workingDir, kwargs.get("startTreePath"), 
                                         kwargs.get("constraintTreePath"), outputPath, Configs.numCores)
    else:
        raise Exception("Tree estimation method {} not recognized..".format(method))
    
    external_tools.runCommand(**task)

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
    
    task = external_tools.runRaxmlNgOptimize(reducedAlignPath, model, unoptimizedPath, workingDir, modelPath, outputPath)
    external_tools.runCommand(**task)
    return modelPath

def raxmlReadModel(modelPath):
    with open(modelPath) as file:
        firstLine = file.readline().strip()
        modelToken = firstLine.split(",")[0]
        return modelToken

def getRaxmlModel():
    if Configs.raxmlModel is None and Configs.raxmlModelSourcePath is not None and Configs.raxmlModelSourcePath != "estimate":
        Configs.log("Estimating RAxML model from {}".format(Configs.raxmlModelSourcePath))
        modelTreePath = os.path.join(os.path.dirname(Configs.raxmlModelSourcePath), "optimized_{}".format(os.path.basename(Configs.raxmlModelSourcePath)))
        modelPath = raxmlEvaluateModelParameters(Configs.raxmlModelSourcePath, Configs.alignmentPath, None, modelTreePath)
        Configs.raxmlModel = raxmlReadModel(modelPath)
        Configs.log("Setting RAxML model to {}".format(Configs.raxmlModel))
    return Configs.raxmlModel
    