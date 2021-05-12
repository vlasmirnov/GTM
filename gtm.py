'''
Created on Jun 24, 2019

@author: Vlad
'''

import treeutils
import gtm_old
import argparse
import dendropy
import time
import os
from platform import node


def runGtm(constraintTrees, startTree, mode):
    if mode == "old":
        return gtm_old.gtmMerge(constraintTrees, startTree)
    
    if mode == "new":
        treeutils.annotateTrees(startTree, constraintTrees)  
        treeutils.collapseViolatingEdges(startTree, True)   
        treeutils.rerootConstraintTrees(constraintTrees)
        mapTreeNodes(startTree)
        return assembleTree(startTree, constraintTrees)
    
    treeutils.annotateTrees(startTree, constraintTrees)  
    treeutils.collapseViolatingEdges(startTree, mode == "convex")   
    treeutils.rerootConstraintTrees(constraintTrees)
        
    if mode == "convex":
        return joinConvexSubtrees(startTree)
    
    treeutils.mapConstraintTreeNodes(startTree, constraintTrees)    
    return resolveTree(startTree)                 

def mapTreeNodes(startTree):
    edgeMapping = {}
    for node in startTree.preorder_node_iter():
        if node.edge.subtree is None:
            node.mappedNode = node
            if node != startTree.seed_node:
                node.label = "scaffold"
        else:
            subtree, bitmask = node.edge.subtree, node.edge.desc[node.edge.subtree]
            cEdge = subtree.subEdgeMap[bitmask]
            node.mappedNode = cEdge.head_node
            for child in node.child_edges():
                if child.desc.get(subtree) == bitmask:
                    edgeMapping[cEdge] = edgeMapping.get(cEdge, []) + [node]
    
    for edge, nodes in edgeMapping.items():
        for node in nodes:
            joinPoint = edge.tail_node.new_child()
            joinPoint.add_child(edge.tail_node.remove_child(edge.head_node))
            try:
                joinPoint.edge.length = edge.length / (len(nodes) + 1)
            except:
                pass
            node.mappedNode = joinPoint
        try:
            edge.length = edge.length / (len(nodes) + 1)
        except:
            pass

def assembleTree(startTree, constraintTrees):
    joinPoints = set()
    unjoins = []
    for node in startTree.preorder_node_iter():
        if node.edge.subtree is None and len(node.child_edges()) == 0:
            subtree = list(node.edge.desc.keys())[0]
            joinPoints.update(getJoinPoints(node, subtree))
            continue
               
        for child in node.child_edges():
            if child.subtree != node.edge.subtree:
                unjoins.append(child.head_node)
                if child.subtree != None:
                    joinPoints.update(getJoinPoints(node, child.subtree))
                    #joinPoints.add((node.mappedNode, child.subtree.seed_node))
                else:
                    joinPoints.add((node.mappedNode, child.head_node))
                
    for child in unjoins:
        child.parent_node.remove_child(child)
    
    for tree in constraintTrees:
        if tree not in startTree.seed_node.edge.desc:
            attachPoint = startTree.seed_node
            joinPoints.update(getJoinPoints(attachPoint, tree))
                
    result = glueJoinPoints(startTree, joinPoints)
    result.collapse_basal_bifurcation()
    return result     

def getJoinPoints(node, subtree):
    if len(subtree.seed_node.child_edges()) == 0:
        return [(node.mappedNode, subtree.seed_node)]
    else:
        return [(node.mappedNode, child.head_node) for child in subtree.seed_node.child_edges()]

def glueJoinPoints(startTree, joinPoints):  
    for parent, child in joinPoints:
        if child.parent_node is not None:
            child.parent_node.remove_child(child)
        parent.add_child(child)
    return dendropy.Tree(seed_node = startTree.seed_node)  

def joinConvexSubtrees(startTree):
    joinPoints = {}
    
    for node in startTree.preorder_node_iter():                
        if node.edge.subtree is None and len(node.child_edges()) == 0:
            subtree = list(node.edge.desc.keys())[0]
            joinPoints[node] = [subtree.seed_node]
            continue
        
        for adj in node.child_edges():
            if adj.subtree != node.edge.subtree:
                joinPoints[node] = []  
                subtrees = set()
                if node.edge.subtree is not None:
                    subtrees.add(node.edge.subtree)
                    bitmask = node.edge.desc[node.edge.subtree]                            
                    joinEdge = node.edge.subtree.subEdgeMap[bitmask]
                    joinPoint = joinEdge.tail_node.new_child()
                    joinPointChild = joinEdge.tail_node.remove_child(joinEdge.head_node)
                    joinPoint.add_child(joinPointChild)
                    joinPoints[node].append(joinPoint)
                    node.edge.subtree.subEdgeMap[bitmask] = joinPointChild.edge
                for child in node.child_edges():
                    if child.subtree is None or child.subtree in subtrees:
                        continue
                    subtrees.add(child.subtree)
                    joinPoints[node].append(child.subtree.seed_node)
                break
    
    root = startTree.seed_node
    for node, points in joinPoints.items():
        point = points[0]
        if node.edge.subtree is None:
            if node.parent_node is None:
                root = point
            elif node.parent_node is not None:
                node.parent_node.add_child(point)
                node.parent_node.remove_child(node)
        for child in node.child_nodes():
            if child.edge.subtree is None:
                node.remove_child(child)
                point.add_child(child)
        for child in points[1:]:
            for gchild in child.child_nodes():
                child.remove_child(gchild)
                point.add_child(gchild)


    return dendropy.Tree(seed_node = root)    

def resolveTree(startTree):
    used = set()
    for edge in startTree.postorder_edge_iter():
        head, tail = edge.head_node, edge.tail_node
        if edge.subtree is None and head != startTree.seed_node:
            head.label = "scaffold"
        joins = [head]
        
        for subtree, bitmask in edge.desc.items():
            if (subtree, bitmask) in used:
                continue
            used.add((subtree, bitmask))
            used.add((subtree, bitmask ^ subtree.startTreeBitmask))
            edgeSet = subtree.fullSubEdgeMap.get(bitmask, [])
            edgeSet.extend(reversed(subtree.fullSubEdgeMap.get(bitmask ^ subtree.startTreeBitmask, [])))
            
              
            for branch in edgeSet:
                if branch.tail_node in subtree.nodeMap:
                    joinPoint = subtree.nodeMap[branch.tail_node]
                else:
                    joinPoint = dendropy.Node()
                    if tail is not None:
                        tail.remove_child(head)
                        tail.add_child(joinPoint)
                        joins.append(joinPoint)
                    
                    joinPoint.add_child(head)
                    head = joinPoint
                    joinPoint.edge.desc = dict(edge.desc)
                    
                branch.tail_node.remove_child(branch.head_node)         
                joinPoint.add_child(branch.head_node)
            
        try:
            if len(joins) > 1:
                el = edge.length / (len(joins) + 1)
                for nd in joins:
                    nd.edge.length = el
        except:
            pass
                        
    return startTree
        
def main(args):   
    startTime = time.time()
    
    startTreePath = args.start
    outputPath = args.output
    mode = args.mode
    
    #if mode not in ["convex", "old", "fp"]:
    #    raise Exception("Mode {} not recognized.".format(mode))
    
    constraintTreePaths = []
    for p in args.trees:
        path = os.path.abspath(p)
        if os.path.isdir(path):
            for filename in os.listdir(path):
                constraintTreePaths.append(os.path.join(path, filename))
        else:
            constraintTreePaths.append(path)
    
    constraintTrees = [treeutils.loadTree(path) for path in constraintTreePaths]
    startTree = treeutils.loadTree(startTreePath)    
    resultTree = runGtm(constraintTrees, startTree, mode)
    treeutils.writeTree(resultTree, outputPath)       
    
    endTime = time.time()
    print("Finished GTM in {} seconds..".format(endTime-startTime))

#Following the arguments convention from NJMerge and TreeMerge (Erin Molloy)
if __name__ == "__main__":
    parser = argparse.ArgumentParser()

    parser.add_argument("-s", "--start", type=str,
                        help="Input starting tree file", required=True)

    parser.add_argument("-t", "--trees", type=str, nargs="+",
                        help="Input subset tree files", required=True)

    parser.add_argument("-o", "--output", type=str,
                        help="Output tree",
                        required=True)
    
    parser.add_argument("-m", "--mode", type=str,
                        help="Desired GTM mode (convex|fp|old)",
                        required=False, default="convex")

    main(parser.parse_args())