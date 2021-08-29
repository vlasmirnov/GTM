'''
Created on Jul 8, 2021

@author: Vlad
'''

import argparse
import time
import shutil
import os
from tree_refiner import tree_refiner_operations
from configs import Configs


def main(args):       
    treePath = args.start
    outputPath = args.output
    workingDir = os.path.join(os.path.dirname(outputPath), "tree_refiner_working_dir")
    workingDir = os.path.abspath(args.directory) if args.directory is not None else workingDir
    
    subtrees = []
    for p in args.trees:
        path = os.path.abspath(p)
        if os.path.isdir(path):
            for filename in os.listdir(path):
                subtrees.append(os.path.join(path, filename))
        else:
            subtrees.append(path)
    
    polytomyStrategy = args.polytomies
    branchLengthStrategy = args.branchlengths
    bootstrap = args.bootstrap.lower() == "true"
    alignmentPath = args.alignment
    model = args.model   
    
    Configs.polytomyTreeSize = args.polytomytreesize
    Configs.branchLengthTreeSize = args.branchlengthtreesize
    Configs.bootstrapTrees = args.bootstraptrees
    Configs.bootstrapTreeSize = args.bootstraptreesize
    
    treePath = tree_refiner_operations.refineTree(treePath, subtrees, workingDir, polytomyStrategy, 
                                                  branchLengthStrategy, bootstrap, alignmentPath, model)
    shutil.copyfile(treePath, outputPath)

if __name__ == "__main__":
    parser = argparse.ArgumentParser()

    parser.add_argument("-a", "--alignment", type=str,
                        help="Path to sequence alignment fasta", required=False, default=None)

    parser.add_argument("-s", "--start", type=str,
                        help="Input starting tree file", required=True)
    
    parser.add_argument("-t", "--trees", type=str, nargs="+",
                        help="Input subset tree files", required=False, default = [])

    parser.add_argument("-o", "--output", type=str,
                        help="Output tree",
                        required=True)
    
    parser.add_argument("-d", "--directory", type=str,
                        help="Path to working directory", required=False, default=None)
    
    parser.add_argument("-p", "--polytomies", type=str,
                        help="Resolve polytomies (tree method or none)",
                        required=False, default="raxml")
    
    parser.add_argument("-l", "--branchlengths", type=str,
                        help="Resolve branch lengths (tree method or none)",
                        required=False, default="raxml")
    
    parser.add_argument("-b", "--bootstrap", type=str,
                        help="Estimate bootstrap support",
                        required=False, default="false")
    
    parser.add_argument("--model", type=str,
                        help="Model argument for RAxML/IQTree calls. <Model>, 'estimate', or omit (optimize separately for each run)",
                        required=False, default=None)
    
    parser.add_argument("--polytomytreesize", type=int,
                        help="Maximum subtree size for Polytomy Tree Refiner",
                        required=False, default=200)
    
    parser.add_argument("--branchlengthtreesize", type=int,
                        help="Maximum subtree size for Branch Length Tree Refiner",
                        required=False, default=50000)
    
    parser.add_argument("--bootstraptreesize", type=int,
                        help="Maximum subtree size for Bootstrap Tree Refiner",
                        required=False, default=200)
    
    parser.add_argument("--bootstraptrees", type=int,
                        help="Number of bootstrap replicates to use",
                        required=False, default=10)

    main(parser.parse_args())