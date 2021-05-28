'''
Created on May 27, 2021

@author: Vlad
'''

import os
from tools import external_tools
from configs import Configs

def buildTree(method, outputPath, **kwargs):
    workingDir = os.path.join(os.path.dirname(outputPath), os.path.basename(outputPath).split(".")[0])
    if method.lower() == "clustal":
        task = external_tools.runClustalOmegaGuideTree(kwargs["alignmentPath"], workingDir, outputPath, Configs.numCores)
    elif method.lower() == "fasttree":
        task = external_tools.runFastTree(kwargs["alignmentPath"], workingDir, outputPath, kwargs.get("mode","fast"), kwargs.get("startTreePath"))
    elif method.lower() == "raxml":
        task = external_tools.runRaxmlNg(kwargs["alignmentPath"], workingDir, kwargs.get("startTreePath"), 
                                         kwargs.get("constraintTreePath"), outputPath, Configs.numCores)
    else:
        raise Exception("Tree estimation method {} not recognized..".format(method))
    
    external_tools.runCommand(**task)