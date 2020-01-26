## BIBLIS
 
import numpy as np
import matplotlib.pyplot as plt
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
 
## FAREY TREE AND ITS DESCENDANCE
 
# construct_tree(node_criterion) : d'abord on construit l'arbre
# tree_size(node) : on peut calculer sa taille
 
# list_of_descendants(node) : on peut lister la descendance
# list_of_nodes_with_prop(node,prop) : on peut lister la descendance filtrée selon un critere
# list_of_nodes_with_trace_equal(t) : liste des nodes de trace t
# list_of_descendants(node)
# list_of_adresses(node)

# get_prop_trace_equal(t) : renvoie une fonction node --> bool (testant si trace == t) 
# get_prop_max_height(max_height) : renvoie fonction
 
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
    # Returns the number of nodes in the tree of a given root.
    if node is None:
        return 0
    return 1 + tree_size(node.leftChild) + tree_size(node.rightChild)
 
def list_of_descendants(node): 
    #Returns the list of descendants of a node including this node.
    return list_of_nodes_with_prop(node, lambda x: True)

def list_of_addresses(root, prop):
    all_nodes = list_of_nodes_with_prop(root, prop)
    return [n.address for n in all_nodes]

def list_of_nodes_with_prop(root, prop):
    """ Returns all nodes of the tree having the given property.
       Args:
       
       root: root of the tree (type treeNode)
       prop: predicate treeNode -> bool
    """
   
    if root is None:
        return []
 
    l=prop(root)*[root]
   
    return l + list_of_nodes_with_prop(root.leftChild, prop) + \
               list_of_nodes_with_prop(root.rightChild, prop)
 
def get_prop_trace_equal(trace):
    """ Returns a predicate treeNode -> bool,
       choosing nodes which have a given trace.
   """
    def trace_equal(node):
        return node.get_trace() == trace
    return trace_equal
 
def list_of_nodes_with_trace_equal(t):
    """ Construct the tree from root to find those nodes.
       We skim the branches that are too large in trace.
   """
    criterion = lambda node, height: node.get_trace() <= t and height <= t
    root = construct_tree(criterion)
    all_nodes = list_of_descendants(root)
 
    return list(filter(get_prop_trace_equal(t),all_nodes))
 
def get_prop_max_height(max_height):
    # Return a stopping criterion for a given height.
    def prop_max_height(node, height):
        return height <= max_height
    return prop_max_height
 
## REPRESENTATION DES MATRICES ET DES CLASSES DE CONJUGAISON

# matrix_of_address
# list_of_circular_shifts
# circular_min_rpz
# compress
# conj_class_with_trace


def matrix_of_address(address):
    # Returns the matrix corresponding to a given 'L'/'T' address.
    T = np.array([[1,1],[0,1]])
    L = np.array([[1,0],[1,1]])
    prod = np.array([[1,0],[0,1]])
    for c in address:
        if c == 'T':
            prod = prod.dot(T)
        elif c=='L':
            prod = prod.dot(L)
    return prod
 
def list_of_circular_shifts(word):
    # Returns all the circular shifts of the given input word.
    n = len(word)
    return ["".join([word[i - j] for i in range(n)]) for j in range(n)]
 
def circular_min_rpz(word):
    """ Returns the canonical representant of the set of circular shifts
       of a word. We choose: min (lex order) of the circular permutations
       of 'word'.
   """
    if word == '':
        return ''
    return min(list_of_circular_shifts(word)) ### J'AI LAISSE LE MIN CAR L<T ET C CE QUE JE VEUX

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
 
def conj_class_with_trace(t, right_column=False):
    """ Return the conjugaison classes of trace t.
 
       Input:
           t: a trace (different than 2)
 
       Output:
           1. The number of classes
           2. np.array [coord[0], coord[1], class_number]
           3. Array of representant of each classes
   """
    nodes = list_of_nodes_with_trace_equal(t)
    classes = {}
    for n in nodes:
        rep = circular_min_rpz(n.address)
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
 
 
## INVERSE DE GAUSS

def is_invertible(word): # peut être mieux is_invertible(node) ?
    A = matrix_of_address(word)
    a,b,c,d = A[0,0], A[0,1], A[1,0], A[1,1]
    return np.gcd.reduce(c,d-a,b)==1

def Gauss_inverse(word):
    return word[::-1]

def is_ambiguous(word):
    return circular_min_rpz(word) == circular_min_rpz(Gauss_inverse(word))

def one_complent(word):
    to_return = ""
    for c in word:
        if c == 'L':
            to_return += 'T'
        else:
            to_return += 'L'

    return to_return

def is_inert(word):
    return circular_min_rpz(word) == circular_min_rpz(one_complent(word))

def is_reciprocal(word):
    return circular_min_rpz(word) == circular_min_rpz(Gauss_inverse(one_complent(word)))

## Rademacher
def Rademacher(word):
    return word.count("L")-word.count("T")

## GENUS

def is_principal_genus(word):
    rep = circular_min_rpz(word)
    return (compress(rep).count('L') == 1) and (compress(rep).count('T') == 1)

## CALCUL DE LA TRACE QUANTIQUE

def quantize(word):
    Lq = cypari.pari("[q,0;1,1/q]")
    Tq = cypari.pari("[q,1;0,1/q]")

    if len(word) == 0:
        return cypari.pari("[1,0;0,1]")
    result = Lq if word[0] == "L" else Tq

    for c in word[1:]:
        result *= (Lq if c == "L" else Tq)

    return result 

def Fricke(word):
    M = quantize(word)
    return M[0][0]+M[1][1]

def are_Fricke_equiv(word1,word2):
    return Fricke(word1) == Fricke(word2)    

## INVERSE MATRICIEL

def matrix_sl_inverse(mat):
    a,b = mat[0][0],mat[1][0]
    c,d = mat[0][1],mat[1][1]
    return np.array([[d,-c],[-b,a]])


## FORME QUADRATIQUE ET MATRICE

def quadratic_form_of_matrix(matrix):
    if matrix.trace()<0:
        matrix=-matrix
    l =  matrix[0][1] 
    m =  matrix[1][1]-matrix[0][0]
    u = -matrix[1][0]
    return (l,m,u)

def matrix_of_quadratic_form(quad):
    l, m, u = quad(0), quad(1), quad(2)
    t = np.sqrt(m**2-4*l*u+4)
    a, d = (t-m)/2, (t+m)/2
    b, c = -u, l
    mat = np.array([[a,c],[b,d]])
    return mat

## DISCRIMINANT ET PRODUIT SCALAIRE DE KILLING DE DEUX FORMES QUADRATIQUES ET RESULTANT

def discriminant(mat):
    return (mat.trace()**2-4)

def killing_form(matrix_A,matrix_B):
    #m_A = matrix_of_address(word_a) 
    #m_B = matrix_of_address(word_b)
    scal = 2*(matrix_A * matrix_B).trace()-matrix_A.trace()*matrix_B.trace()
    return scal

def resultant(mat_A,mat_B):
    la,ma,ua = quadratic_form_of_matrix(mat_A)
    lb,mb,ub = quadratic_form_of_matrix(mat_B)
    res = (la*ua-lb*ub)**2-ma*mb*(la*ua+lb*ub)+(la*mb**2*ua)+(lb*ma**2*ub)
    return res


## FRACTIONS CONTIUES


def continued_fraction(x,long, expansion=[], A = np.matrix([[1,0],[0,1]]), colonnes=[] ):
    """ 
    entrer : un réel x et un entier long 
    retourne :
    expansion : 
    colonnes : liste des 'long' premières approximations par fractions continuées partielles de 'x'
    """
    
    if long <= 0:
        return(expansion, colonnes, A)
    
    n=int(np.floor(x))
    expansion.append(n)
    
    A=A*np.matrix([[n,1],[1,0]])
    colonnes.append("{} /{}".format(int(A[0,1]), int(A[1,1]) ))
    
    return continued_fraction(1/(x-n),long-1, expansion, A, colonnes)


def continued_fraction_hj(x,long, expansion=[], A = np.matrix([[1,0],[0,1]]), colonnes=[] ):
    """ Idem mais pour les fractions continues d'Hirzebruch Youg
    """
    if long <= 0:
        return expansion#, colonnes, A
    
    n=int(np.ceil(x))
    expansion.append(n)
    
    A=A*np.matrix([[n,-1],[1,0]])
    colonnes.append("{} /{}".format(int(A[0,1]), int(A[1,1]) ))
    
    return continued_fraction_hj(1/(n-x),long-1, expansion, A, colonnes)


