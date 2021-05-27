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
from tools import external_tools

        
def main(args):   
    startTime = time.time()
    
    startTreePath = args.start
    outputPath = args.output
    mode = args.mode
    
    constraintTreePaths = []
    for p in args.trees:
        path = os.path.abspath(p)
        if os.path.isdir(path):
            for filename in os.listdir(path):
                constraintTreePaths.append(os.path.join(path, filename))
        else:
            constraintTreePaths.append(path)
    
    workingDir = args.directory if args.directory is not None else os.path.dirname(outputPath)
    branchLengthTreePath = args.branchlengthstree
    alignmentPath = args.alignment
    
    constraintTrees = [treeutils.loadTree(path) for path in constraintTreePaths]
    startTree = treeutils.loadTree(startTreePath)    
    resultTree = gtm.runGtm(constraintTrees, startTree, mode)
    treeutils.writeTree(resultTree, outputPath)       
    
    constraintTrees = [treeutils.loadTree(path) for path in constraintTreePaths]
    startTree = treeutils.loadTree(outputPath)
    if branchLengthTreePath is not None:
        spanningTree = treeutils.loadTree(branchLengthTreePath)
    else:
        spanningTree = gtmutils.extractSpanningTree(resultTree, constraintTrees, 1000)
        branchLengthTreePath = os.path.join(workingDir, "spanning_tree.tre")
        computeSpanningTree(spanningTree, alignmentPath, branchLengthTreePath)
        spanningTree = treeutils.loadTree(branchLengthTreePath)
    
    optimizedPath = os.path.join(os.path.dirname(outputPath), "lengths_" + os.path.basename(outputPath))    
    gtmutils.applySpanningTreeBranchLengths(startTree, spanningTree)
    treeutils.writeTree(startTree, optimizedPath)
    
    endTime = time.time()
    print("Finished GTM pipeline in {} seconds..".format(endTime-startTime))

def computeSpanningTree(tree, alignmentPath, treePath):
    baseName = os.path.basename(treePath).split(".")[0]
    unoptimizedPath = os.path.join(os.path.dirname(treePath), "unoptimized_{}.tre".format(baseName))
    workingDir = os.path.join(os.path.dirname(treePath), "raxml_{}".format(baseName))
    outModelPath = os.path.join(os.path.dirname(treePath), "model_{}.txt".format(baseName))
    reducedAlignPath = os.path.join(os.path.dirname(treePath), "alignment_{}.txt".format(baseName))
    
    taxa = [n.taxon.label for n in tree.leaf_nodes()]
    align = sequenceutils.readFromFasta(alignmentPath)
    sequenceutils.writeFasta(align, reducedAlignPath, taxa)
    
    treeutils.writeTree(tree, unoptimizedPath)
    if os.path.exists(workingDir):
        shutil.rmtree(workingDir)
    os.makedirs(workingDir)    
    task = external_tools.runRaxmlNgOptimize(reducedAlignPath, None, unoptimizedPath, workingDir, outModelPath, treePath)
    external_tools.runCommand(**task)

if __name__ == "__main__":
    parser = argparse.ArgumentParser()

    parser.add_argument("-a", "--alignment", type=str,
                        help="Path to sequence alignment fasta", required=False, default=None)

    parser.add_argument("-s", "--start", type=str,
                        help="Input starting tree file", required=True)

    parser.add_argument("-t", "--trees", type=str, nargs="+",
                        help="Input subset tree files", required=True)
    
    parser.add_argument("-d", "--directory", type=str,
                        help="Path to working directory", required=False, default=None)

    parser.add_argument("-o", "--output", type=str,
                        help="Output tree",
                        required=True)
    
    parser.add_argument("-m", "--mode", type=str,
                        help="Desired GTM mode (convex|fp|old)",
                        required=False, default="convex")
    
    parser.add_argument("--branchlengthstree", type=str,
                        help="Path to tree for updating branch lengths", required=False, default=None)

    main(parser.parse_args())