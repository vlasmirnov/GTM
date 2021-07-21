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
    
    polytomyStrategy = args.polytomies
    branchLengthStrategy = args.branchlengths
    alignmentPath = args.alignment
    model = args.model   
    Configs.spanningTreeSize = args.refinertreesize 
    
    treePath = tree_refiner_operations.refineTree(treePath, workingDir, polytomyStrategy, branchLengthStrategy, alignmentPath, model)
    shutil.copyfile(treePath, outputPath)

if __name__ == "__main__":
    parser = argparse.ArgumentParser()

    parser.add_argument("-a", "--alignment", type=str,
                        help="Path to sequence alignment fasta", required=False, default=None)

    parser.add_argument("-s", "--start", type=str,
                        help="Input starting tree file", required=True)

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
    
    parser.add_argument("--model", type=str,
                        help="Model argument for RAxML/IQTree calls. <Model>, 'estimate', or omit (optimize separately for each run)",
                        required=False, default=None)
    
    parser.add_argument("--refinertreesize", type=int,
                        help="Maximum subtree size for Tree Refiner",
                        required=False, default=500)

    main(parser.parse_args())