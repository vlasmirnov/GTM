'''
Created on Nov 29, 2020

@author: Vlad
'''

import dendropy

def loadTree(treePath):
    tree = dendropy.Tree.get(path=treePath, schema="newick", preserve_underscores=True)
    tree.is_rooted = False
    #tree.resolve_polytomies(limit=2)
    tree.collapse_basal_bifurcation()
    return tree

def writeTree(tree, treePath):
    treeString = tree.as_string(schema="newick", suppress_rooting = True)
    with open(treePath, "w") as f:
        f.write(treeString) 

def annotateTrees(startTree, constraintTrees):
    startTree.taxonMap = {n.taxon.label : n for n in startTree.leaf_nodes()}
    for label, node in startTree.taxonMap.items():
        node.taxon.taxonSubtree, node.taxon.taxonBitmask = None, 0
    
    for tree in constraintTrees:
        tree.rootBipartition = None
        taxa = [n.taxon for n in tree.leaf_nodes()]
        tree.numTaxa = len(taxa)
        startTreeTaxa = [t.label for t in taxa if t.label in startTree.taxonMap]
        tree.startTreeBitmask = tree.taxon_namespace.taxa_bitmask(labels = startTreeTaxa)
        
        for taxon in taxa:
            taxon.taxonSubtree = tree
            taxon.taxonBitmask = tree.taxon_namespace.taxon_bitmask(taxon)
            if taxon.label in startTree.taxonMap:
                startTree.taxonMap[taxon.label].taxon = taxon        
        
        #tree.update_bipartitions()
        for edge in tree.postorder_edge_iter():
            populateEdgeDesc(edge)
        populateEdgeMaps(tree)   
        populateEdgeLengths(tree, False)  
    
    for edge in startTree.postorder_edge_iter():
        populateEdgeDesc(edge)
    populateEdgeLengths(startTree, True)

def populateEdgeDesc(edge):
    if len(edge.head_node.child_edges()) == 0:
        edge.desc = {edge.head_node.taxon.taxonSubtree : edge.head_node.taxon.taxonBitmask}
    else:
        edge.desc = combineDescMaps([child.desc for child in edge.head_node.child_edges()])
    
    edge.subtree = None
    for subtree, bitmask in edge.desc.items():
        if subtree is None:
            continue
        if bitmask != subtree.startTreeBitmask and edge.subtree is None:
            edge.subtree = subtree
        elif bitmask != subtree.startTreeBitmask and edge.subtree is not None:
            edge.subtree = None
            break   
        
def populateEdgeMaps(tree):
    tree.edgeMap = {}
    tree.subEdgeMap = {}    
    for edge in tree.preorder_edge_iter():
        bitmask = edge.desc[tree]
        tree.edgeMap[bitmask] = edge
        tree.subEdgeMap[bitmask & tree.startTreeBitmask] = edge 

def populateEdgeLengths(tree, isStartTree):
    tree.edgeLengthMap = {}
    for edge in tree.preorder_edge_iter():
        tree.edgeLengthMap[buildDescKey(edge.desc, isStartTree)] = edge.length

def buildDescKey(desc, isStartTree):
    if isStartTree:
        keyList = [(s, b & s.startTreeBitmask) for s,b in desc.items() if s is not None]
        keyList = [(s,b) for s,b in keyList if b != 0]
    else:
        keyList = [(s, min(b, ~b & s.taxon_namespace.all_taxa_bitmask())) for s,b in desc.items()]
    return frozenset(keyList)

def combineDescMaps(descMaps):
    newDesc = {}
    for desc in descMaps:
        for subtree, bitmask in desc.items():
            newDesc[subtree] = newDesc.get(subtree, 0) | bitmask
    return newDesc

def distributeNodes(startTree, constraintTrees):
    for tree in constraintTrees:
        recurseDistributeNodes(startTree.seed_node, tree.seed_node, tree) 
    startTree.suppress_unifurcations()
    #startTree.collapse_basal_bifurcation()     
            
def recurseDistributeNodes(sNode, cNode, cTree):
    stack = [cNode]
    while stack:
        cNode = stack.pop()        
        cMask = cNode.edge.desc[cTree] & cTree.startTreeBitmask
        maskEdges = {e.desc[cTree] : e for e in sNode.child_edges() if cTree in e.desc and e.desc[cTree] & cMask == e.desc[cTree]}
        
        if cMask == 0:
            if len(sNode.adjacent_nodes())>2:
                resolvePolytomy(sNode, sNode.child_nodes())
            attachNodetoNode(sNode, cNode)
        elif cMask in maskEdges and len(maskEdges[cMask].desc) > 1:
            recurseDistributeNodes(maskEdges[cMask].head_node, cNode, cTree)
        else:            
            if not any(len(maskEdges[mask].desc) > 1 for mask in maskEdges):
                for mask in maskEdges:
                    sNode.remove_child(maskEdges[mask].head_node)
                attachNodetoNode(sNode, cNode) 
            elif len(maskEdges) < len(sNode.child_edges()):
                newDesc = combineDescMaps([maskEdges[mask].desc for mask in maskEdges])
                if not any(s is not None and b not in s.subEdgeMap and b ^ s.startTreeBitmask not in s.subEdgeMap for s, b in newDesc.items()):
                    newNode = resolvePolytomy(sNode, [maskEdges[m].head_node for m in maskEdges])
                    recurseDistributeNodes(newNode, cNode, cTree)
                else:
                    stack.extend(cNode.child_nodes())
            else:
                stack.extend(cNode.child_nodes())

def resolvePolytomy(node, children):
    newNode = dendropy.Node()
    node.add_child(newNode)
    for child in children:
        newNode.add_child(node.remove_child(child))
    populateEdgeDesc(newNode.edge)
    return newNode    

def attachNodetoNode(sNode, cNode):
    if cNode.parent_node is not None:
        cNode.parent_node.remove_child(cNode)
    sNode.add_child(cNode)

def collapseViolatingEdges(startTree, removeConvexityViolations):
    toRemove = []
    for edge in startTree.postorder_edge_iter():
        violatingEdge = False
        subtreeList = []
        
        for subtree, bitmask in edge.desc.items():
            if subtree is None:
                continue
            if bitmask not in subtree.subEdgeMap and bitmask ^ subtree.startTreeBitmask not in subtree.subEdgeMap: 
                violatingEdge = True
                break
            elif bitmask != subtree.startTreeBitmask:
                subtreeList.append(subtree)
                if removeConvexityViolations and len(subtreeList) > 1:
                    violatingEdge = True
                    break
        if violatingEdge:
            toRemove.append(edge.head_node)   
             
    collapseEdges(toRemove)
        
def collapseEdges(headNodes):
    for n in headNodes:
        childs = n.child_nodes()
        for child in childs:
            n.remove_child(child)
            n.parent_node.add_child(child)
        n.parent_node.remove_child(n)    
            
def rerootConstraintTrees(startTree, constraintTrees):
    for edge in startTree.postorder_edge_iter():
        subtreeList = [s for s, b in edge.desc.items() if s is not None and b != s.startTreeBitmask]
        for subtree in subtreeList:
            subtree.rootBipartition = edge.desc[subtree]    
    
    for tree in constraintTrees:        
        if tree.numTaxa < 3:
            populateEdgeMaps(tree)
            continue
        
        rootEdge = None
        if tree.rootBipartition is None and len(tree.seed_node.child_edges()) > 0:
            rootEdge = tree.seed_node.child_edges()[0]
        elif tree.rootBipartition is not None:
            for edge in tree.edges():
                bitmask = edge.desc[tree]
                if bitmask & tree.startTreeBitmask == tree.rootBipartition or ~bitmask & tree.startTreeBitmask == tree.rootBipartition:
                    rootEdge = edge
                    break
        
        if rootEdge is not None:    
            joinPoint = rootEdge.tail_node.new_child()
            joinPoint.add_child(rootEdge.tail_node.remove_child(rootEdge.head_node))
            tree.reseed_at(joinPoint, update_bipartitions=False,
            collapse_unrooted_basal_bifurcation=False,
            suppress_unifurcations=False)
            for edge in tree.postorder_edge_iter():
                populateEdgeDesc(edge)
            populateEdgeMaps(tree)

def dealWithEdgeLengths(startTree): 
    for edge in startTree.postorder_edge_iter():
        populateEdgeDesc(edge)
    if len(startTree.seed_node.child_edges()) == 2:
        childs = startTree.seed_node.child_edges()
        if buildDescKey(childs[0].desc, True) not in startTree.edgeLengthMap:  
            collapseEdges([childs[0].head_node])
        else:
            collapseEdges([childs[1].head_node])
    
    labelScaffold(startTree)    
    scaffoldEdgeSources, constraintEdgeSources = getEdgeSources(startTree)
    edgeMap = buildEdgeLengthMap(startTree)
    #applyEdgeLengths(scaffoldEdgeSources, edgeMap)
    applyEdgeLengths(constraintEdgeSources, edgeMap)

def labelScaffold(startTree):
    for edge in startTree.preorder_edge_iter(): 
        edge.head_node.label = None
        
    for edge in startTree.preorder_edge_iter():
        if edge.tail_node is None:
            continue  
        splits = [s for s, b in edge.desc.items() if s is not None and b != s.taxon_namespace.all_taxa_bitmask()]
        if len(splits) != 1:
            edge.head_node.label = "scaffold"
            for e in edge.get_adjacent_edges():
                if e.head_node.label != "scaffold":
                    e.head_node.label = "scaffold_adjacent"
        else:
            edge.subtree = splits[0]

def getEdgeSources(startTree): 
    scaffoldEdgeSources = set()
    constraintEdgeSources = set()
    for edge in startTree.preorder_edge_iter():
        if edge.tail_node is None:
            continue  
        if edge.head_node.label == "scaffold":
            scaffoldEdgeSources.add((edge, startTree))
        else:
            constraintEdgeSources.add((edge, edge.subtree))
    return scaffoldEdgeSources, constraintEdgeSources

def buildEdgeLengthMap(startTree):
    edgeMap = {}    
    for edge in startTree.preorder_edge_iter():
        edgeMap[edge] = {}
        if edge.tail_node is None:
            continue
        key = buildDescKey(edge.desc, True)
        edgeMap[startTree, key] = edgeMap.get((startTree, key), []) + [edge]
        edgeMap[edge][startTree] = key
        
        splits = [(s, b) for s, b in edge.desc.items() if s is not None and b != s.taxon_namespace.all_taxa_bitmask()]
        for subtree, bitmask in splits:
            key = buildDescKey({subtree : bitmask}, False)
            edgeMap[subtree, key] = edgeMap.get((subtree, key), []) + [edge]
            edgeMap[edge][subtree] = key
    return edgeMap

def applyEdgeLengths(edgeSources, edgeMap):
    keySet = set((sourceTree, edgeMap[edge][sourceTree]) for edge, sourceTree in edgeSources)    
    for sourceTree, key in keySet:        
        length = sourceTree.edgeLengthMap[key]           
        edgeGroup = edgeMap[sourceTree, key]        
        updateEdges = [e for e in edgeGroup if (e, sourceTree) in edgeSources]
        skipEdges = [e for e in edgeGroup if (e, sourceTree) not in edgeSources]
        try:
            sSum = sum([e.length for e in skipEdges])
            for edge in updateEdges:                
                edge.length = (length - sSum) / len(updateEdges)
        except:
            pass
       
def mapConstraintTreeNodes(startTree, constraintTrees):
    startTree.nodeMap = {}
    bipartSets = {tree : set() for tree in constraintTrees}
    for edge in startTree.preorder_edge_iter():
        for subtree, bitmask in edge.desc.items():
            bipartSets[subtree].add(bitmask)
            childMasks = [child.desc[subtree] for child in edge.head_node.child_edges() if child.desc.get(subtree, bitmask) != bitmask]
            for childMask in childMasks:
                startTree.nodeMap[subtree, bitmask, childMask] = edge.head_node
    
    toRemove = []
    for tree in constraintTrees:
        for edge in tree.postorder_edge_iter():
            edge.desc = {}
            bitmask = edge.bipartition.leafset_bitmask & tree.startTreeBitmask
            if bitmask != 0:
                if bitmask not in bipartSets[tree] and bitmask ^ tree.startTreeBitmask not in bipartSets[tree]:
                    toRemove.append(edge.head_node)
    collapseEdges(toRemove)
    
    for tree in constraintTrees:    
        tree.fullSubEdgeMap = {}
        tree.nodeMap = {}
            
        for edge in tree.postorder_edge_iter():
            bitmask = edge.bipartition.leafset_bitmask & tree.startTreeBitmask
            if bitmask != 0:
                for child in edge.head_node.child_edges():
                    childMask = child.bipartition.leafset_bitmask & tree.startTreeBitmask
                    if childMask == 0:
                        if bitmask not in tree.fullSubEdgeMap:
                            tree.fullSubEdgeMap[bitmask] = [child]
                        else:
                            tree.fullSubEdgeMap[bitmask].append(child)  
                    elif childMask != bitmask and bitmask != tree.startTreeBitmask:
                        if (tree, bitmask, childMask) in startTree.nodeMap:
                            tree.nodeMap[edge.head_node] = startTree.nodeMap[tree, bitmask, childMask]
                        else:
                            tree.nodeMap[edge.head_node] = startTree.nodeMap[tree, tree.startTreeBitmask, childMask]


