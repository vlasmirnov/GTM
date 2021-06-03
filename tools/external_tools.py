'''
Created on May 25, 2021

@author: Vlad
'''

import os
import subprocess
import shutil
from configs import Configs

CLUSTAL_PATH = os.path.join(os.path.dirname(os.path.abspath(__file__)), "clustal/clustalo")
FASTTREE_PATH = os.path.join(os.path.dirname(os.path.abspath(__file__)), "fasttree/FastTree")
RAXML_PATH = os.path.join(os.path.dirname(os.path.abspath(__file__)), "raxmlng/raxml-ng")

def runCommand(**kwargs):
    if not os.path.exists(kwargs["workingDir"]):
        os.makedirs(kwargs["workingDir"])
    command = kwargs["command"]
    print("Running an external tool, command: {}".format(command))
    runner = subprocess.run(command, shell = True, cwd = kwargs["workingDir"], universal_newlines = True, stdout=subprocess.PIPE, stderr=subprocess.STDOUT)
    try:    
        runner.check_returncode()
    except:
        print("Command encountered error: {}".format(command))
        print("Exit code: {}".format(runner.returncode))
        print("Output: {}".format(runner.stdout))
        raise
    for srcPath, destPath in kwargs.get("fileCopyMap", {}).items():
        shutil.move(srcPath, destPath)
    return runner

def runClustalOmegaGuideTree(fastaPath, workingDir, outputPath, threads = 1):
    tempPath = os.path.join(os.path.dirname(outputPath), "temp_{}".format(os.path.basename(outputPath)))
    args = [CLUSTAL_PATH]
    args.extend(["-i", fastaPath, "--max-hmm-iterations=-1", "--guidetree-out={}".format(tempPath)])
    args.extend(["--threads={}".format(threads)])
    taskArgs = {"command" : subprocess.list2cmdline(args), "fileCopyMap" : {tempPath : outputPath}, "workingDir" : workingDir, "outputFile" : outputPath}
    return taskArgs

def runFastTree(fastaFilePath, workingDir, outputPath, mode = "normal", intree = None):
    tempPath = os.path.join(os.path.dirname(outputPath), "temp_{}".format(os.path.basename(outputPath)))
    
    args = [FASTTREE_PATH]
    if Configs.inferDataType(fastaFilePath) == "protein":
        args.extend(["-lg"])
    else:
        args.extend(["-nt", "-gtr"])
    
    if intree is not None:
        args.extend(["-intree", intree])
    
    if mode == "fast":
        args.extend(["-fastest", "-nosupport"]) 
    elif mode == "faster":
        args.extend(["-fastest", "-nosupport", "-mlnni", "4" ]) 
    elif mode == "noml":
        args.extend(["-fastest", "-nosupport", "-noml"])
    
    args.extend([fastaFilePath, ">", tempPath])
    taskArgs = {"command" : subprocess.list2cmdline(args), "fileCopyMap" : {tempPath : outputPath}, "workingDir" : workingDir, "outputFile" : outputPath}
    return taskArgs

def runRaxmlNg(fastaFilePath, model, workingDir, startTreePath, constraintTreePath, outputPath, threads = 8):
    # raxml-ng --msa prim.phy --model GTR+G --prefix T4 --threads 2 --seed 2 --tree pars{25},rand{25}
    baseName = os.path.basename(outputPath).replace(".","")
    raxmlFile = os.path.join(workingDir, "{}.raxml.bestTree".format(baseName))
    #seed = random.randint(1, 1000000)
    seed = 74
    args = [RAXML_PATH,
            "--msa", fastaFilePath,
            "--prefix", baseName,
            "--threads", str(threads),
            "--seed", str(seed)]
    
    if model is not None:
        args.extend(["--model", model])
    elif Configs.inferDataType(fastaFilePath) == "protein":
        args.extend(["--model", "LG+G"])
    else:
        args.extend(["--model", "GTR+G"])
    
    if constraintTreePath is not None:
        args.extend(["--tree-constraint", constraintTreePath])
        args.extend(["--tree", "rand{{{}}}".format(1)])
    elif startTreePath is not None:
        args.extend(["--tree", startTreePath])
    else:
        args.extend(["--tree", "pars{{{}}}".format(1)])
        
    
    taskArgs = {"command" : subprocess.list2cmdline(args), "fileCopyMap" : {raxmlFile : outputPath}, "workingDir" : workingDir, "outputFile" : outputPath}
    return taskArgs

def runRaxmlNgOptimize(fastaFilePath, model, treePath, workingDir, outputModelPath, outputPath, threads = 8):
    #raxml-ng --evaluate --msa prim.phy --threads 2 --model GTR+G+FC --tree T3.raxml.bestTree --prefix E4
    baseName = os.path.basename(outputPath).replace(".","")
    logFile = os.path.join(workingDir, "{}.raxml.log".format(baseName))
    modelFile = os.path.join(workingDir, "{}.raxml.bestModel".format(baseName))
    raxmlFile = os.path.join(workingDir, "{}.raxml.bestTree".format(baseName))
    args = [RAXML_PATH,
            "--evaluate",
            "--msa", fastaFilePath,
            "--prefix", baseName,
            "--threads", str(threads),
            "--tree", treePath]  
    
    if model is not None:
        args.extend(["--model", model])
    elif Configs.inferDataType(fastaFilePath) == "protein":
        args.extend(["--model", "LG+G"])
    else:
        args.extend(["--model", "GTR+G"])
        
    taskArgs = {"command" : subprocess.list2cmdline(args), "fileCopyMap" : {modelFile : outputModelPath, raxmlFile : outputPath}, "workingDir" : workingDir, "outputFile" : outputPath}
    return taskArgs

def runRaxmlNgScore(fastaFilePath, model, treePath, workingDir, threads = 8):
    args = [RAXML_PATH,
            "--loglh",
            "--msa", fastaFilePath,
            "--threads", str(threads),
            "--tree", treePath]  
    
    if model is not None:
        args.extend(["--model", model])
    elif Configs.inferDataType(fastaFilePath) == "protein":
        args.extend(["--model", "LG+G"])
    else:
        args.extend(["--model", "GTR+G"])
        
    taskArgs = {"command" : subprocess.list2cmdline(args), "workingDir" : workingDir}
    return taskArgs