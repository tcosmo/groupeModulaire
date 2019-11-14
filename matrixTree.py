import numpy as np
import matplotlib.pyplot as plt
plt.style.use('dark_background')

import cypari

class treeNode(object):
    def __init__(self, matrix=None, address="", 
                       leftChild=None, rightChild=None):
        """ A node of the Farey binary tree.
            Args:
            
            matrix: np 2x2 matrix (dtype=np.int64)
            address: string of 'L' and 'R' giving the address in the tree.
            leftChild, rightChild: children of that node.
        """
        
        self.matrix = matrix
        self.address = address
        self.leftChild = leftChild
        self.rightChild = rightChild
        
    def get_coord(self):
        """ Transforms the matrix in coordinates.
        """
        return sum(self.matrix[0]), sum(self.matrix[1])
    
    def get_trace(self):
        return self.matrix.trace()

def tree_size(node):
    """ Returns the number of nodes in the tree of a given root.
    """
    if node is None:
        return 0
    return 1 + tree_size(node.leftChild) + tree_size(node.rightChild)

def list_of_descendants(node):
    """ Returns the list of descendants of a node including this node.
    """
    return get_nodes_with_property(node, lambda x: True)

def get_nodes_with_property(root, prop):
    """ Returns all nodes of the tree having the given property.
        Args:
        
        root: root of the tree (type treeNode)
        prop: predicate treeNode -> bool
    """
    
    if root is None:
        return []
    
    if prop(root):
        l = [root]
    else:
        l = []
    
    return l + get_nodes_with_property(root.leftChild, prop) + \
               get_nodes_with_property(root.rightChild, prop)

def get_trace_property(trace):
    """ Returns a predicate treeNode -> bool, 
        choosing nodes which have a given trace.
    """
    def trace_property(node):
        return node.get_trace() == trace
    return trace_property

def construct_tree(include_node_criterion, 
                   left_coord=(1,0), right_coord=(0,1), 
                   last_address = "", last_move = "", height=0):
    """ Construct the Farey binary tree.
        Args:
        
        include_node_criterion: function that takes a node and the height 
                                and decide whether or not to include this node in the tree (thus stops the construction).
        left_coord: starting at (1,0)
        right_coord: starting at (0,1)
        last_addresss: address of the father
        last_move: last move 'L' or 'R' to get to the current node
        height: current height in the tree
    """
    node = treeNode()
    node.matrix = np.array([[left_coord[0],right_coord[0]], [left_coord[1],right_coord[1]]], dtype=np.int64)
    node.address = last_address + last_move
    
    if not include_node_criterion(node, height):
        return None
    
    node.leftChild = construct_tree(include_node_criterion,
                                    left_coord, node.get_coord(), 
                                    node.address, "L", height+1)
    
    node.rightChild = construct_tree(include_node_criterion,
                                    node.get_coord(), right_coord,
                                    node.address, "R", height+1)
    
    return node

def get_nodes_with_trace(t):
    """ Construct the tree from root to find those nodes.
        We skim the branches that are too large in trace.
    """
    criterion = lambda node, height: node.get_trace() <= t and height <= t
    root = construct_tree(criterion)
    all_nodes = list_of_descendants(root)

    return list(filter(get_trace_property(t),all_nodes))

def get_max_height_criterion(max_height):
    """ Return a stopping criterion for a given height.
    """
    def max_height_criterion(node, height):
        return height <= max_height
    return max_height_criterion

def get_matrix_of_address(address):
    """ Returns the matrix corresponding to a given 'L'/'R' address.
    """
    T = np.array([[1,1],[0,1]])
    L = np.array([[1,0],[1,1]])
    prod = np.array([[1,0],[0,1]])
    for c in address:
        if c == 'L':
            prod = prod.dot(T)
        else:
            prod = prod.dot(L)
    return prod

def get_circular_shifts(word):
    """ Returns all the circular shifts of the given input word.
    """
    n = len(word)
    return ["".join([word[i - j] for i in range(n)]) for j in range(n)]

def get_circular_rep(word):
    """ Returns the canonical representant of the set of circular shifts
        of a word. We choose: min (lex order) of the circular permutations 
        of 'word'.
    """
    if word == '':
        return ''
    return min(get_circular_shifts(word))

def compress(word):
    """ Returns compact representation of word.
        Example: LRRRR -> L R4 
    """

    letters = ['L', 'R']
    first_letter = word[0]
    current_letter = word[0]
    counts = [0]
    for c in word:
        if c == current_letter:
            counts[-1] += 1
        else:
            current_letter = c
            counts.append(1)

    compress = ""
    for i,c in enumerate(counts):
        choose = i%2 if first_letter == 'L' else 1-i%2
        compress += letters[choose]+ ("" if c == 1 else str(c)) + " "

    return compress[:-1]

def get_conjugaison_classes(t):
    """ Return the conjugaison classes of trace t.

        Input:
            t: a trace (different than 2)

        Output:
            1. The number of classes
            2. np.array [coord[0], cood[1], class_number]
            3. Array of representant of each classes
    """
    nodes = get_nodes_with_trace(t)
    classes = {}
    for n in nodes:
        rep = get_circular_rep(n.address)
        if not rep in classes:
            classes[rep] = []
        classes[rep].append(n)
    
    sorted_class = sorted(list(classes.items()), key=lambda x: len(x[0]))
    
    to_return = []
    reps = []
    for i,(rep,class_) in enumerate(sorted_class):
        for n in class_:
            to_return.append([n.get_coord()[0], n.get_coord()[1],i])
        reps.append(compress(rep))
            
    return len(classes), np.array(to_return), reps 