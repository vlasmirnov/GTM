'''
Created on Jul 8, 2021

@author: Vlad
'''

import dendropy
import os
from helpers import treeutils, gtmutils


class PipelineContext:
        
    def __init__(self, **kwargs):
        self.workingDir = None
        self.alignmentPath = None
        self.currentTreePath = None
        self.outputFile = None        
        
        self.startTreePath = None
        self.startTreeMethod = None
        self.startTree = None        
        
        self.constraintTreePaths = None
        self.constraintTrees = None
        
        self.guideTreePath = None
        self.guideTreeStrategy = None        
        self.guideTree = None
        
        self.subsetsDir = None
        self.subtreesDir = None
        self.subsetPaths = None
        
        self.iterations = 0
        self.decompositionStrategy = None
        self.mode = None
        self.model = None
        self.modelSourcePath = None
        self.useInducedStartTreeForML = False
        
        
        for attr in kwargs:
            vars(self)[attr] = kwargs.get(attr)
        
        if not os.path.exists(self.workingDir):
            os.makedirs(self.workingDir)
        
        self.initialize()    
    
    def initialize(self):
        if self.startTree is None and self.startTreePath is not None:
            self.startTree = treeutils.loadTree(self.startTreePath)
        if self.constraintTrees is None:
            self.constraintTrees = [treeutils.loadTree(path) for path in self.constraintTreePaths]
    
    def getStartTreeForML(self):
        return self.startTreePath if self.useInducedStartTreeForML else None
            