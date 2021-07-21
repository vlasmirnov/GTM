'''
Created on Jun 24, 2019

@author: Vlad
'''

from helpers import treeutils
from gtm.gtm_context import GtmContext
from gtm import gtm_operations
import argparse
import os
      

def main(args):   
    startTreePath = args.start
    outputPath = args.output
    mode = args.mode
    
    #if mode not in ["convex", "old", "fp"]:
    #    raise Exception("Mode {} not recognized.".format(mode))
    
    constraintTreePaths = []
    for p in args.trees:
        path = os.path.abspath(p)
        if os.path.isdir(path):
            for filename in os.listdir(path):
                constraintTreePaths.append(os.path.join(path, filename))
        else:
            constraintTreePaths.append(path)
    
    gtmContext = GtmContext(startTreePath = startTreePath, constraintTreePaths = constraintTreePaths, mode = mode)
    resultTree = gtm_operations.runGtm(gtmContext)
    treeutils.writeTree(resultTree, outputPath)       

#Following the arguments convention from NJMerge and TreeMerge (Erin Molloy)
if __name__ == "__main__":
    parser = argparse.ArgumentParser()

    parser.add_argument("-s", "--start", type=str,
                        help="Input starting tree file", required=True)

    parser.add_argument("-t", "--trees", type=str, nargs="+",
                        help="Input subset tree files", required=True)

    parser.add_argument("-o", "--output", type=str,
                        help="Output tree",
                        required=True)
    
    parser.add_argument("-m", "--mode", type=str,
                        help="Desired GTM mode (convex|fp|old)",
                        required=False, default="convex")

    main(parser.parse_args())