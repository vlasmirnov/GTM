'''
Created on May 27, 2021

@author: Vlad
'''
import os
import sequenceutils

class Configs:
    treeMethods = ("clustal", "raxml", "fasttree")
    
    startTree = None
    outputPath = None
    mode = None
    constraintTrees = []   
    treeMethod = "fasttree"
    guideTree = None
    workingDir = None
    branchLengthsTree = None
    alignmentPath = None    
    
    dataType = None
    
    iterations = 1
    spanningTreeSize = 200
    decompositionMaxNumSubsets = 500
    decompositionMaxSubsetSize = 1000
    decompositionStrategy = "centroid"
    
    numCores = os.cpu_count()

    @staticmethod
    def log(msg, path = None):
        print(msg)
        #path = Configs.logPath if path is None else path
        #Configs.writeMsg(msg, path)

    @staticmethod
    def inferDataType(sequencesFile):
        if Configs.dataType is None:
            Configs.dataType = sequenceutils.inferDataType(sequencesFile)
            Configs.log("Data type wasn't specified. Inferred data type {} from {}".format(Configs.dataType.upper(), sequencesFile))
        return Configs.dataType 

def buildConfigs(args):
    if args.start is not None:
        Configs.startTree = os.path.abspath(args.start) if os.path.exists(os.path.abspath(args.start)) else args.start   
    if args.guide is not None:
        Configs.guideTree = os.path.abspath(args.guide) if os.path.exists(os.path.abspath(args.guide)) else args.guide    
    
    Configs.outputPath = os.path.abspath(args.output)
    Configs.mode = args.mode
    
    Configs.workingDir = os.path.join(os.path.dirname(Configs.outputPath), "gtm_working_dir")
    Configs.workingDir = os.path.abspath(args.directory) if args.directory is not None else Configs.workingDir
    if not os.path.exists(Configs.workingDir):
        os.makedirs(Configs.workingDir)
    
    if len(args.trees) == 1 and args.trees[0].lower() in Configs.treeMethods:
        Configs.constraintTrees = []
        Configs.treeMethod = args.trees[0]
    else:
        for p in args.trees:
            path = os.path.abspath(p)
            if os.path.isdir(path):
                for filename in os.listdir(path):
                    Configs.constraintTrees.append(os.path.join(path, filename))
            else:
                Configs.constraintTrees.append(path)
    
    if args.branchlengthstree is not None:
        Configs.branchLengthsTree = os.path.abspath(args.branchlengthstree) if os.path.exists(os.path.abspath(args.branchlengthstree)) else args.branchlengthstree
    if args.alignment is not None:
        Configs.alignmentPath = os.path.abspath(args.alignment) if os.path.exists(os.path.abspath(args.alignment)) else args.alignment
        
    Configs.iterations = args.iterations
    Configs.decompositionMaxSubsetSize = args.maxsubsetsize
    Configs.decompositionMaxNumSubsets = args.maxnumsubsets
    Configs.decompositionStrategy = args.decompstrategy