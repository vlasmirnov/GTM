'''
Created on May 25, 2021

@author: Vlad
'''

import shutil
import argparse
import time
import os
from configs import buildConfigs, Configs
from pipeline.pipeline_context import PipelineContext
from pipeline import pipeline_operations
from tree_refiner import tree_refiner_operations

        
def main(args):   
    startTime = time.time()
    
    buildConfigs(args)    
    Configs.log("Starting GTM pipeline..")
    regularPipeline()
    
    endTime = time.time()
    Configs.log("Finished GTM pipeline in {} seconds..".format(endTime-startTime))

def regularPipeline():
    context = PipelineContext(workingDir = Configs.workingDir,
                              alignmentPath = Configs.alignmentPath, 
                              startTreePath = Configs.startTreePath, startTreeMethod = Configs.startTreeMethod,
                              constraintTreePaths = Configs.constraintTreePaths, 
                              guideTreePath = Configs.guideTreePath, guideTreeStrategy = Configs.guideTreeStrategy,
                              decompositionStrategy = Configs.decompositionStrategy,
                              mode = Configs.mode,
                              model = Configs.model, modelSourcePath = Configs.modelSourcePath,
                              useInducedStartTreeForML = Configs.useInducedStartTreeForML)
    
    for i in range(Configs.iterations):
        startTime = time.time()
        Configs.log("Starting iteration {} of {}..".format(i+1, Configs.iterations))
        workingDir = os.path.join(Configs.workingDir, "iteration_{}".format(i + 1))
        outputPath = os.path.join(workingDir, "gtm_tree.tre")
        pipeline_operations.iterateGtm(context, workingDir, outputPath)
        endTime = time.time()
        Configs.log("Finished iteration {} in {} seconds..".format(i+1, endTime-startTime))
    
    outputPath = tree_refiner_operations.refineTree(context.outputPath, Configs.workingDir, Configs.polytomyStrategy, 
                                       Configs.branchLengthStrategy, context.alignmentPath, context.model)

    shutil.copyfile(outputPath, Configs.outputPath)



if __name__ == "__main__":
    parser = argparse.ArgumentParser()

    parser.add_argument("-a", "--alignment", type=str,
                        help="Path to sequence alignment fasta", required=False, default=None)

    parser.add_argument("-s", "--start", type=str,
                        help="Starting tree file or method", required=False, default = "clustal")

    parser.add_argument("-t", "--trees", type=str, nargs="+",
                        help="Input subset tree files or method", required=False, default = ["raxml"])
    
    parser.add_argument("-g", "--guide", type=str,
                        help="Base method for building the guide tree", required=False, default=None)
    
    parser.add_argument("-d", "--directory", type=str,
                        help="Path to working directory", required=False, default=None)

    parser.add_argument("-o", "--output", type=str,
                        help="Output tree",
                        required=True)
    
    parser.add_argument("-m", "--mode", type=str,
                        help="Desired GTM mode (convex|fp|old)",
                        required=False, default="convex")
    
    parser.add_argument("-p", "--polytomies", type=str,
                        help="Resolve polytomies (tree method or none)",
                        required=False, default="raxml")
    
    parser.add_argument("-l", "--branchlengths", type=str,
                        help="Resolve branch lengths (tree method or none)",
                        required=False, default="raxml")
    
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
    
    parser.add_argument("--refinertreesize", type=int,
                        help="Maximum subtree size for Tree Refiner",
                        required=False, default=500)
    
    parser.add_argument("--model", type=str,
                        help="Model argument for RAxML/IQTree calls. <Model>, 'estimate', or None (optimize separately for each run)",
                        required=False, default=None)
    
    parser.add_argument("--useinducedstarttreeforml", type=str,
                        help="Pass current working tree as a starting tree for ML runs",
                        required=False, default="false")

    main(parser.parse_args())