'''
Created on Nov 29, 2020

@author: Vlad
'''

import dendropy

def loadTree(treePath):
    tree = dendropy.Tree.get(path=treePath, schema="newick", preserve_underscores=True)
    tree.is_rooted = False
    tree.resolve_polytomies(limit=2)
    tree.collapse_basal_bifurcation()
    return tree

def writeTree(tree, treePath):
    treeString = tree.as_string(schema="newick", suppress_rooting = True)
    with open(treePath, "w") as f:
        f.write(treeString) 

def annotateTrees(startTree, constraintTrees):
    startTree.taxaSubtrees = {}
    startTree.taxaBitmasks = {}
    
    for tree in constraintTrees:
        tree.startTreeTaxa = set()
        taxa = [n.taxon for n in tree.leaf_nodes()]
        tree.numTaxa = len(taxa)
        for taxon in taxa:
            startTree.taxaSubtrees[taxon.label] = tree
            startTree.taxaBitmasks[taxon.label] = tree.taxon_namespace.taxon_bitmask(taxon)
    
    for n in startTree.leaf_nodes():
        taxon = n.taxon.label
        subtree  = startTree.taxaSubtrees[taxon]
        subtree.startTreeTaxa.add(taxon)
        
    for tree in constraintTrees:
        tree.rootBipartition = None
        tree.startTreeBitmask = tree.taxon_namespace.taxa_bitmask(labels = tree.startTreeTaxa)
        tree.update_bipartitions()
        populateEdgeMap(tree)
        populateEdgeLengthMap(tree)
    
    populateEdgeDesc(startTree)
    populateStartTreeEdgeLengthMap(startTree)

def populateEdgeDesc(startTree):
    for edge in startTree.postorder_edge_iter():
        edge.desc = {}
        
        if len(edge.head_node.child_edges()) == 0:
            taxon = edge.head_node.taxon.label
            subtree = startTree.taxaSubtrees[taxon]
            edge.desc = {subtree : startTree.taxaBitmasks[taxon]}
        else:
            for child in edge.head_node.child_edges():
                for subtree, bitmask in child.desc.items():
                    edge.desc[subtree] = edge.desc.get(subtree, 0) | bitmask

def distributeNodes(startTree, constraintTrees):
    for tree in constraintTrees:
        recurseDistributeNodes(startTree.seed_node, tree.seed_node, tree)        
            
def recurseDistributeNodes(sNode, cNode, cTree):
    internalChilds = [e for e in sNode.child_edges() if len(e.desc) > 1 and cTree in e.desc]
    treeChilds = [e for e in sNode.child_edges() if len(e.desc) == 1 and cTree in e.desc]
    internalMasks = {e.desc[cTree] : e for e in internalChilds}
    treeMasks = {e.desc[cTree] : e for e in treeChilds}
    
    stack = [cNode]
    while stack:
        cNode = stack.pop()        
        cMask = cNode.edge.bipartition.leafset_bitmask & cTree.startTreeBitmask
        cNode.edge.desc = {cTree : cMask} if cMask > 0 else {}
        if cMask == 0:
            attachNode(sNode, cNode)
        elif cMask in internalMasks:
            recurseDistributeNodes(internalMasks[cMask].head_node, cNode, cTree)
        else:
            subInternalMasks = [mask for mask in internalMasks if mask & cMask == mask]
            subTreeMasks = [mask for mask in treeMasks if mask & cMask == mask]
            if len(subTreeMasks) > 0 and len(subInternalMasks) == 0:
                for mask in subTreeMasks:
                    sNode.remove_child(treeMasks[mask].head_node)
                attachNode(sNode, cNode)
            else:
                stack.extend(cNode.child_nodes())

def attachNode(sNode, cNode):
    if cNode.parent_node is not None:
        cNode.parent_node.remove_child(cNode)
    sNode.add_child(cNode)

def dealWithEdgeLengths(startTree):   
    populateEdgeDesc(startTree)
    edgeMap = {}
    sEdgeMap = {}
    
    for edge in startTree.preorder_edge_iter():
        if edge.tail_node is None:
            continue
        edge.length = None
        
        splits = [(subtree, bitmask) for subtree, bitmask in edge.desc.items() if bitmask != subtree.taxon_namespace.all_taxa_bitmask()]
        if len(splits) != 1:
            edge.head_node.label = "scaffold"
            for e in edge.get_adjacent_edges():
                if e.head_node.label != "scaffold":
                    e.head_node.label = "scaffold_adjacent"
            
            desc = {subtree : subtree.startTreeBitmask & bitmask for subtree, bitmask in edge.desc.items()}
            key = frozenset(desc)
            sEdgeMap[startTree, key] = sEdgeMap.get((startTree, key), []) + [edge]
            
        for subtree, bitmask in splits:
            bitmask = bitmask if bitmask in subtree.edgeLengthMap else bitmask ^ subtree.taxon_namespace.all_taxa_bitmask()
            edgeMap[subtree, bitmask] = edgeMap.get((subtree, bitmask), []) + [edge]
    
    for key, edges in sEdgeMap.items():
        subtree, bitmask = key
        for edge in edges:
            try:
                edge.length = subtree.edgeLengthMap[bitmask] / len(edges)
            except:
                pass
        
    for key, edges in edgeMap.items():
        subtree, bitmask = key
        cEdges = [e for e in edges if e.length is None]
        sEdges = [e for e in edges if e.length is not None]
        sSum = sum([e.length for e in sEdges])
        for edge in cEdges:
            try:
                edge.length = (subtree.edgeLengthMap[bitmask] - sSum) / len(cEdges)
            except:
                pass
            
def populateEdgeMap(tree):
    tree.edgeMap = {}
    tree.subEdgeMap = {}
    for edge in tree.postorder_edge_iter():
        bitmask = edge.bipartition.leafset_bitmask
        tree.edgeMap[bitmask] = edge
        tree.subEdgeMap[bitmask & tree.startTreeBitmask] = edge 

def populateEdgeLengthMap(tree):
    tree.edgeLengthMap = {}
    for edge in tree.preorder_edge_iter():
        bitmask = edge.bipartition.leafset_bitmask
        tree.edgeLengthMap[bitmask] = edge.length

def populateStartTreeEdgeLengthMap(tree):
    tree.edgeLengthMap = {}
    for edge in tree.preorder_edge_iter():
        key = frozenset(edge.desc)
        tree.edgeLengthMap[key] = edge.length

def collapseViolatingEdges(startTree, removeConvexityViolations):
    toRemove = []
    for edge in startTree.postorder_edge_iter():
        violatingEdge = False
        edge.subtree = None

        subtreeList = []
        for subtree, bitmask in edge.desc.items():
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
        else:
            for subtree in subtreeList:
                subtree.rootBipartition = edge.desc[subtree]
            if len(subtreeList) == 1:
                edge.subtree = subtreeList[0]
    
    collapseEdges(toRemove)
        
def collapseEdges(headNodes):
    for n in headNodes:
        childs = n.child_nodes()
        for child in childs:
            n.remove_child(child)
            n.parent_node.add_child(child)
        n.parent_node.remove_child(n)    
            
def rerootConstraintTrees(trees):
    for tree in trees:        
        if tree.numTaxa < 3:
            populateEdgeMap(tree)
            continue
        
        rootEdge = None
        if tree.rootBipartition is None and len(tree.seed_node.child_edges()) > 0:
            rootEdge = tree.seed_node.child_edges()[0]
        elif tree.rootBipartition is not None:
            for edge in tree.edges():
                bitmask = edge.bipartition.leafset_bitmask
                if bitmask & tree.startTreeBitmask == tree.rootBipartition or ~bitmask & tree.startTreeBitmask == tree.rootBipartition:
                    rootEdge = edge
                    break
        
        if rootEdge is not None:    
            joinPoint = rootEdge.tail_node.new_child()
            joinPoint.add_child(rootEdge.tail_node.remove_child(rootEdge.head_node))
            tree.reseed_at(joinPoint, update_bipartitions=True,
            collapse_unrooted_basal_bifurcation=False,
            suppress_unifurcations=False)
            populateEdgeMap(tree)

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


