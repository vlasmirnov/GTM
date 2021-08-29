'''
Created on May 27, 2021

@author: Vlad
'''
import os
import time
from helpers import sequenceutils

class Configs:
    #treeMethods = ("clustal", "raxml", "fasttree", "iqtree")
    
    workingDir = None
    alignmentPath = None   
    outputPath = None 
    
    startTreePath = None
    startTreeMethod = None
        
    constraintTreePaths = []   
    treeMethod = "fasttree"
    
    guideTreePath = None
    guideTreeStrategy = None  
    guideTreeMethod = None
    
    polytomyStrategy = None
    branchLengthStrategy = None
    
    mode = None
    dataType = None
    model = None
    modelSourcePath = None
    useInducedStartTreeForML = False
    trackMLScores = False
    
    iterations = 1
    decompositionMaxNumSubsets = 1000
    decompositionMaxSubsetSize = 500
    decompositionStrategy = "centroid"
    raxmlModelLimit = 10000
    polytomyTreeSize = 200
    branchLengthTreeSize = 50000
    bootstrapTreeSize = 200
    guideTreeRecursionFactor = 2
    guideTreeRecursionBaseSize = 1000
    
    useBootstrap = True
    bootstrapTrees = 10
    
    numCores = os.cpu_count()
    logPath = None
    errorPath = None
    debugPath = None

    @staticmethod
    def log(msg, path = None):
        print(msg)
        path = Configs.logPath if path is None else path
        Configs.writeMsg(msg, path)
    
    @staticmethod
    def error(msg, path = None):
        Configs.log(msg)
        path = Configs.errorPath if path is None else path
        Configs.writeMsg(msg, path)
    
    @staticmethod
    def debug(msg, path = None):
        print(msg)
        path = Configs.debugPath if path is None else path
        Configs.writeMsg(msg, path)
    
    @staticmethod
    def writeMsg(msg, path):
        if path is not None:
            with open(path, 'a') as logFile:
                logFile.write("{}    {}\n".format(time.strftime("%Y-%m-%d %H:%M:%S"), msg))

    @staticmethod
    def inferDataType(sequencesFile):
        if Configs.dataType is None:
            Configs.dataType = sequenceutils.inferDataType(sequencesFile)
            Configs.log("Data type wasn't specified. Inferred data type {} from {}".format(Configs.dataType.upper(), sequencesFile))
        return Configs.dataType 
    
def buildConfigs(args):
    if args.start is not None and os.path.exists(os.path.abspath(args.start)):
        Configs.startTreePath = args.start
    elif args.start is not None:
        Configs.startTreeMethod = args.start

    if args.guide is not None and os.path.exists(os.path.abspath(args.guide)):
        Configs.guideTreePath = os.path.abspath(args.guide)
    elif args.guide is not None:
        Configs.guideTreeMethod = args.guide   
    Configs.guideTreeStrategy = args.guidetreestrategy 
    
    Configs.outputPath = os.path.abspath(args.output)
    if args.alignment is not None:
        Configs.alignmentPath = os.path.abspath(args.alignment) if os.path.exists(os.path.abspath(args.alignment)) else args.alignment
    Configs.workingDir = os.path.join(os.path.dirname(Configs.outputPath), "gtm_working_dir")
    Configs.workingDir = os.path.abspath(args.directory) if args.directory is not None else Configs.workingDir
    if not os.path.exists(Configs.workingDir):
        os.makedirs(Configs.workingDir)
    
    if len(args.trees) == 1 and not os.path.exists(os.path.abspath(args.trees[0])):
        Configs.treeMethod = args.trees[0]
    else:
        for p in args.trees:
            path = os.path.abspath(p)
            if os.path.isdir(path):
                for filename in os.listdir(path):
                    Configs.constraintTreePaths.append(os.path.join(path, filename))
            else:
                Configs.constraintTreePaths.append(path)
    
    Configs.polytomyStrategy = args.polytomies
    if args.branchlengths is not None:
        Configs.branchLengthStrategy = os.path.abspath(args.branchlengths) if os.path.exists(os.path.abspath(args.branchlengths)) else args.branchlengths
        
    Configs.iterations = args.iterations
    Configs.decompositionMaxSubsetSize = args.maxsubsetsize
    Configs.polytomyTreeSize = args.polytomytreesize
    Configs.branchLengthTreeSize = args.branchlengthtreesize
    Configs.bootstrapTreeSize = args.bootstraptreesize
    Configs.decompositionMaxNumSubsets = args.maxnumsubsets
    Configs.decompositionStrategy = args.decompstrategy
    Configs.guideTreeRecursionFactor = args.guidetreerecursionfactor
    Configs.guideTreeRecursionBaseSize = args.guidetreerecursionbasesize
    
    if args.model == "estimate":
        Configs.modelSourcePath = "estimate"
        if Configs.startTreePath is not None and os.path.exists(Configs.startTreePath):
            Configs.modelSourcePath = Configs.startTreePath       
    elif args.model is not None and os.path.exists(os.path.abspath(args.model)):
        Configs.modelSourcePath = os.path.abspath(args.model)
    elif args.model is not None:
        Configs.model = args.model
    
    Configs.mode = args.mode
    Configs.useInducedStartTreeForML = args.useinducedstarttreeforml.lower() == "true"
    Configs.trackMLScores = args.trackmlscores.lower() == "true"
    
    Configs.useBootstrap = args.bootstrap.lower() == "true"
    Configs.bootstrapTrees = args.bootstraptrees
        
    Configs.logPath = os.path.join(Configs.workingDir, "log.txt")    
    Configs.errorPath = os.path.join(Configs.workingDir, "log_errors.txt")
    Configs.debugPath = os.path.join(Configs.workingDir, "log_debug.txt") 
    