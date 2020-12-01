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


def runGtm(constraintTrees, startTree, mode):
    if mode == "old":
        return gtm_old.gtmMerge(constraintTrees, startTree)
    
    treeutils.annotateTrees(startTree, constraintTrees)
    treeutils.collapseViolatingEdges(startTree, mode == "convex")   
    treeutils.rerootConstraintTrees(constraintTrees)
        
    if mode == "convex":
        return joinConvexSubtrees(startTree)
    
    treeutils.mapConstraintTreeNodes(startTree, constraintTrees)    
    return resolveTree(startTree)                 

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
                    
                    joinPoint.add_child(head)
                    head = joinPoint
                    joinPoint.edge.desc = dict(edge.desc)
                    
                branch.tail_node.remove_child(branch.head_node)         
                joinPoint.add_child(branch.head_node)
                        
    return startTree
        
def main(args):   
    startTime = time.time()
    
    startTreePath = args.start
    outputPath = args.output
    mode = args.mode
    
    if mode not in ["convex", "old", "fp"]:
        raise Exception("Mode {} not recognized.".format(mode))
    
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