'''
Created on May 25, 2021

@author: Vlad
'''

import shutil
import treeutils
import sequenceutils
import argparse
import dendropy
import time
import os
import gtm
import gtmutils
import decomposition
from tools import external_tools, methods
from configs import buildConfigs, Configs

        
def main(args):   
    startTime = time.time()
    
    buildConfigs(args)
    
    alignmentPath = Configs.alignmentPath
    startTreeArg = Configs.startTree
    constraintTreePaths = Configs.constraintTrees
    guideTreeArg = Configs.guideTree
    for i in range(Configs.iterations):
        workingDir = os.path.join(Configs.workingDir, "iteration_{}".format(i+1))
        outputPath = os.path.join(workingDir, "gtm_tree.tre")
        iterateGtm(alignmentPath, startTreeArg, constraintTreePaths, guideTreeArg, workingDir, outputPath)
        startTreeArg = outputPath
        guideTreeArg = outputPath
        constraintTreePaths = []
    
    #gtmResultPath = os.path.join(Configs.workingDir, "gtm_result.tre")
    #shutil.copyfile(outputPath, gtmResultPath)
    
    branchLengthsTree = Configs.branchLengthsTree
    if branchLengthsTree is not None:  
        gtmResultPath = outputPath  
        workingDir = os.path.join(Configs.workingDir, "branch_lengths")
        outputPath = os.path.join(workingDir, "optimized_result.tre")
        updateScaffoldBranchLengths(gtmResultPath, alignmentPath, branchLengthsTree, workingDir, outputPath)  
        
    shutil.copyfile(outputPath, Configs.outputPath)
    
    endTime = time.time()
    print("Finished GTM pipeline in {} seconds..".format(endTime-startTime))

def iterateGtm(alignmentPath, startTreeArg, constraintTreePaths, guideTreeArg, workingDir, outputTreePath):
    if os.path.exists(workingDir):
        shutil.rmtree(workingDir)
    os.makedirs(workingDir)
    
    startTreePath = startTreeArg
    if len(constraintTreePaths) == 0:
        if not os.path.exists(startTreeArg):
            startTreePath = os.path.join(workingDir, "start_tree.tre")
            methods.buildTree(startTreeArg, startTreePath, alignmentPath = alignmentPath)
        subsetsDir = os.path.join(workingDir, "subsets")        
        subsetPaths = decomposition.decomposeGuideTree(subsetsDir, alignmentPath, startTreePath, Configs.decompositionMaxSubsetSize, 
                                                       Configs.decompositionMaxNumSubsets, Configs.decompositionStrategy)
        subtreesDir = os.path.join(workingDir, "subtrees")
        constraintTreePaths = buildSubtrees(subtreesDir, subsetPaths, Configs.treeMethod)
    constraintTrees = [treeutils.loadTree(path) for path in constraintTreePaths]
    
    guideTreePath = guideTreeArg
    if guideTreeArg is None:
        guideTreePath = startTreePath
    guideTree = treeutils.loadTree(guideTreePath)
    
    print(guideTreePath)
    print(constraintTreePaths.keys())
    resultTree = gtm.runGtm(constraintTrees, guideTree, Configs.mode)
    treeutils.writeTree(resultTree, outputTreePath)  

def buildSubtrees(subtreesDir, subsetPaths, method):
    if os.path.exists(subtreesDir):
        shutil.rmtree(subtreesDir)
    os.makedirs(subtreesDir)
    
    subtreePaths = {}
    for subsetPath, subset in subsetPaths.items():
        baseName = os.path.basename(subsetPath).split(".")[0]
        subTreePath = os.path.join(subtreesDir, "subtree_{}.tre".format(baseName))
        subtreePaths[subTreePath] = subset        
        methods.buildTree(method, subTreePath, alignmentPath = subsetPath)
    
    return subtreePaths

def updateScaffoldBranchLengths(inputTreePath, alignmentPath, branchLengthsTreePath, workingDir, outputTreePath):
    if os.path.exists(workingDir):
        shutil.rmtree(workingDir)
    os.makedirs(workingDir)
    
    tree = treeutils.loadTree(inputTreePath)
    if os.path.exists(branchLengthsTreePath):
        spanningTree = treeutils.loadTree(branchLengthsTreePath)
    else:
        spanningTree = gtmutils.extractSpanningTree(tree, None, Configs.spanningTreeSize)
        branchLengthsTreePath = os.path.join(workingDir, "spanning_tree.tre")
        gtmutils.computeSpanningTree(spanningTree, alignmentPath, branchLengthsTreePath)
        spanningTree = treeutils.loadTree(branchLengthsTreePath)
    
    gtmutils.applySpanningTreeBranchLengths(tree, spanningTree)
    treeutils.writeTree(tree, outputTreePath)     


if __name__ == "__main__":
    parser = argparse.ArgumentParser()

    parser.add_argument("-a", "--alignment", type=str,
                        help="Path to sequence alignment fasta", required=False, default=None)

    parser.add_argument("-s", "--start", type=str,
                        help="Starting tree file or method", required=False, default = "clustal")

    parser.add_argument("-t", "--trees", type=str, nargs="+",
                        help="Input subset tree files or method", required=False, default = ["raxml"])
    
    parser.add_argument("-g", "--guide", type=str,
                        help="Base method for building the decomposition tree", required=False, default=None)
    
    parser.add_argument("-d", "--directory", type=str,
                        help="Path to working directory", required=False, default=None)

    parser.add_argument("-o", "--output", type=str,
                        help="Output tree",
                        required=True)
    
    parser.add_argument("-m", "--mode", type=str,
                        help="Desired GTM mode (convex|fp|old)",
                        required=False, default="new")
    
    parser.add_argument("--branchlengthstree", type=str,
                        help="Path to tree for updating branch lengths", required=False, default=None)  
    
    parser.add_argument("--iterations", type=int,
                        help="Number of GTM iterations", required=False, default=1)
    
    parser.add_argument("--maxsubsetsize", type=int,
                        help="Maximum subset size for GTM",
                        required=False, default=500)
    
    parser.add_argument("--maxnumsubsets", type=int,
                        help="Maximum number of subsets for GTM",
                        required=False, default=1000)
    
    parser.add_argument("--decompstrategy", type=str,
                        help="Decomposition strategy (centroid or randomcentroid)",
                        required=False, default="centroid")

    main(parser.parse_args())