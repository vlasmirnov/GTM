'''
Created on Jul 8, 2021

@author: Vlad
'''

import os
from helpers import sequenceutils
from configs import Configs


class PipelineContext:
        
    def __init__(self, **kwargs):
        self.workingDir = None
        self.alignmentPath = None
        self.currentTreePath = None
        
        self.startTreePath = None        
        self.subsetPaths = None
        self.constraintTreePaths = []
        self.guideTreePath = None
        
        self.subsetsDir = None
        self.subtreesDir = None
        
        self.iterations = 0
        self.model = None
        self.modelSourcePath = None
        self.trackMLScores = False
        
        self.taxaLengths = None
        self.topLevelStartTreePath = None
        
        for attr in kwargs:
            vars(self)[attr] = kwargs.get(attr)
        
        if not os.path.exists(self.workingDir):
            os.makedirs(self.workingDir)
    
    def getStartTreeForML(self):
        return self.startTreePath if Configs.useInducedStartTreeForML else None
    
    def getTaxaLengths(self):
        if self.taxaLengths is None and self.alignmentPath is not None:
            sequences = sequenceutils.readFromFasta(self.alignmentPath, removeDashes = True)
            self.taxaLengths = {t : len(s.seq) for t, s in sequences.items()}
        return self.taxaLengths