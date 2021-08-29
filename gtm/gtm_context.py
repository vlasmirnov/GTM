'''
Created on Jul 8, 2021

@author: Vlad
'''

from helpers import treeutils

class GtmContext:
        
    def __init__(self, **kwargs):
        self.startTreePath = None
        self.constraintTreePaths = None
        
        self.startTree = None
        self.constraintTrees = None
        self.guideTree = None     
        
        self.mode = None
        
        self.treeKeyEdges = None
        
        #self.outputFile = None
        #self.workingDir = None
        
        for attr in kwargs:
            vars(self)[attr] = kwargs.get(attr)
        
        #if not os.path.exists(self.workingDir):
        #    os.makedirs(self.workingDir)
        
        self.initialize()    
    
    def initialize(self):
        if self.startTree is None:
            self.startTree = treeutils.loadTree(self.startTreePath)
        if self.constraintTrees is None:
            self.constraintTrees = [treeutils.loadTree(path) for path in self.constraintTreePaths]
        
    
            
 