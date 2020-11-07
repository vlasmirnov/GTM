# GTM
Guide Tree Merger - a tool for merging disjoint trees.

- - - -

## Dependencies
GTM requires
* Python 3
* dendropy Python package

- - - -

## Getting Started
Please navigate to the "example" directory to get started.\
GTM's commandline arguments are as follows:  

**-s** specifies the path to the guide tree (AKA starting tree)\
**-t** specifies the list of subtrees being merged\
**-o** specifies the path to the output tree

We can run our example merge with the following command, which will save the result as gtm_tree.tre:

**python3 ../gtm.py -s guide_tree.tre -t subtrees/subset_tree_0.tre subtrees/subset_tree_1.tre subtrees/subset_tree_10.tre subtrees/subset_tree_11.tre subtrees/subset_tree_12.tre subtrees/subset_tree_13.tre subtrees/subset_tree_2.tre subtrees/subset_tree_3.tre subtrees/subset_tree_4.tre subtrees/subset_tree_5.tre subtrees/subset_tree_6.tre subtrees/subset_tree_7.tre subtrees/subset_tree_8.tre subtrees/subset_tree_9.tre -o gtm_tree.tre**

- - - -

## Things to Keep in Mind

* dendropy's memory usage starts to become excessive when approaching one million taxa. Future work will optimize GTM to deal with this.
