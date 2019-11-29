## BIBLIS
 
import numpy as np
import matplotlib.pyplot as plt
plt.style.use('dark_background')
 
import cypari
 
## CLASS treeNode
 
class treeNode(object):
    def __init__(self, matrix=None, address="",
                       leftChild=None, rightChild=None):
        """ A node of the Farey binary tree.
           Args:
           
           matrix: np 2x2 matrix (dtype=np.int64)
           address: string of 'L' and 'T' giving the address in the tree.
           leftChild, rightChild: children of that node.
       """
       
        self.matrix = matrix
        self.address = address
        self.leftChild = leftChild
        self.rightChild = rightChild
       
    def get_coord(self, right_column=False):
        """ Transforms the matrix in coordinates.
        """
        if right_column:
            return self.matrix[:,1]

        return sum(self.matrix[0]), sum(self.matrix[1])
   
    def get_trace(self):
        return self.matrix.trace()

### J'AI ECHANGE L ORDRE DES FONCTIONS
 
## FAREY TREE AND ITS DESCENDANCE
 
# construct_tree(node_criterion) : d'abord on construit l'arbre
# tree_size(node) : on peut calculer sa taille
 
# list_of_descendants(node) : on peut lister la descendance
# get_noes_with_property(node,prop) : on peut lister la descendance filtrée selon un critere
 
# get_node_with_trace(t) : liste des nodes de trace t
# get_trace_property(t) : renvoie une fonction node --> bool (testant si trace == t)
 
# get_max_height_criterion(max_height) : renvoie fonction
 
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
       last_move: last move 'L' or 'T' to get to the current node
       height: current height in the tree
   """
    node = treeNode()
    node.matrix = np.array([[left_coord[0],right_coord[0]], [left_coord[1],right_coord[1]]], dtype=np.int64)
    node.address = last_address + last_move
   
    if not include_node_criterion(node, height):
        return None
   
    node.leftChild = construct_tree(include_node_criterion,
                                    left_coord, node.get_coord(),
                                    node.address, "T", height+1)
   
    node.rightChild = construct_tree(include_node_criterion,
                                    node.get_coord(), right_coord,
                                    node.address, "L", height+1)
   
    return node
 
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
 
    l=prop(root)*[root] ### ICI JE ME SUIS PERMIS UNE MODIF DE CODE
   
    return l + get_nodes_with_property(root.leftChild, prop) + \
               get_nodes_with_property(root.rightChild, prop)
 
def get_trace_property(trace):
    """ Returns a predicate treeNode -> bool,
       choosing nodes which have a given trace.
   """
    def trace_property(node):
        return node.get_trace() == trace
    return trace_property
 
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
 
## REPRESENTATION DES MATRICES ET DES CLASSES DE CONJUGAISON
 
def get_matrix_of_address(address):
    """ Returns the matrix corresponding to a given 'L'/'T' address.
   """
    T = np.array([[1,1],[0,1]])
    L = np.array([[1,0],[1,1]])
    prod = np.array([[1,0],[0,1]])
    for c in address:
        if c == 'T':
            prod = prod.dot(T)
        elif c=='L':
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
    return min(get_circular_shifts(word)) ### J'AI LAISSE LE MIN CAR L<T ET C CE QUE JE VEUX

def compress(word):
    """ Returns compact representation of word.
       Example: LTTTT -> L T4
   """
 
    letters = ['L', 'T']
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
        choose = i%2 if first_letter == 'L' else 1-i%2 ### ICI JE N'AI PAS CHANGE LE L : JE NE SAIS PAS SI JE DOIS
        compress += letters[choose]+ ("" if c == 1 else str(c)) + " "
 
    return compress[:-1]
 
def get_conjugaison_classes(t, right_column=False):
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
            to_return.append([n.get_coord(right_column=right_column)[0], n.get_coord(right_column=right_column)[1],i])
        reps.append(rep)
           
    return len(classes), np.array(to_return), reps
 
## CALCUL DE L'ENLACEMENT
 
def get_word_base(max_word_length, current_pair=('LT', 'TL')):
    """ Returns the word base for the computation of enl.
    # The set of pairs form a binary tree which we fill in à la Pascal
    # The root pair is ('LT','TL')
    # The pairs to the extreme left are ('Ln T', 'T Ln')
    # The pairs on the extreme right are ('L Tn', 'Tn L')
    # The children of (G, D) are
    # to the left (P, Q) with P=G[:-1]+'LT' (on enlève 'T' on rajoute 'LT') and Q=D+'L'
    # to the right (R, S) with R=G+'T' and S=D[:-1]+'TL' (on enlève 'L' on rajoute 'TL')
    # COUCOU!!!
    # Note the properties : 'G' ends with 'T' and 'D' ends with 'L' are preserved
    # which is why G[:-1] = G - 'T' and D[:-1] = D - 'L'
    """
 
    if len(current_pair[0]) > max_word_length:
        return []
 
    pair_left = current_pair[0][:-1]+'LT', current_pair[1]+'L'
    pair_right = current_pair[0]+'T', current_pair[1][:-1]+'TL'
 
    return [current_pair] + get_word_base(max_word_length, pair_left) +\
                            get_word_base(max_word_length, pair_right)
 
def occ(P,A):
    """ Returns the number of times P appears at the begining of circular
       shifts of A.
   """
    shifts = get_circular_shifts(A)
    counter = 0
    n = len(P)
 
    if n > len(A):
        return counter
 
    for shift in shifts:
        if shift[:n] == P:
            counter += 1
 
    return counter
 
def scal_PQ(P,Q,A,B):
    return 0.5*(occ(P,A)*occ(Q,B)+ occ(P,B)*occ(Q,A))
 
def enl(A,B):
    """ Returns the enl metric on the words A and B in the L/T alphabet.
   """
    word_base = get_word_base(max(len(A),len(B)))
    return sum([scal_PQ(P,Q,A,B) for P,Q in word_base])
 
## CALCUL DE LA TRACE
def get_quantum(word):
    Lq = cypari.pari("[q,0;1,1/q]")
    Tq = cypari.pari("[q,1;0,1/q]")

    if len(word) == 0:
        return cypari.pari("[1,0;0,1]")
    result = Lq if word[0] == "L" else Tq

    for c in word[1:]:
        result *= (Lq if c == "L" else Tq)

    return result 

def get_Fricke(word):
    M = get_quantum(word)
    return M[0][0]+M[1][1]

def are_Fricke_equiv(word1,word2):
    return get_Fricke(word1) == get_Fricke(word2) 


## FORME QUADRATIQUE ASSOCIEE
 
## INVERSE DE GAUSS
def Gauss_inverse(word):
    return word[::-1]

def is_ambiguous(word):
    return get_circular_rep(word) == get_circular_rep(Gauss_inverse(word))

def one_complent(word):
    to_return = ""
    for c in word:
        if c == 'L':
            to_return += 'T'
        else:
            to_return += 'L'

    return to_return

def is_inert(word):
    return get_circular_rep(word) == get_circular_rep(one_complent(word))

def is_reciprocal(word):
    return get_circular_rep(word) == get_circular_rep(Gauss_inverse(one_complent(word)))

## Rademacher
def Rademacher(word):
    return word.count("L")-word.count("T")

## GENUS

def is_trival_genus(word):
    rep = get_circular_rep(word)
    return (compress(rep).count('L') == 1) and (compress(rep).count('T') == 1)