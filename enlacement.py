## BIBLIS
 
import numpy as np
import matplotlib.pyplot as plt
import cypari

import matrixTree as mT

## CALCUL DE L'ENLACEMENT

# linking_patterns
# occ
# scal_PQ
# enl
 
def linking_patterns(max_word_length, current_pair=('LT', 'TL')):
    """ Returns the word base for the computation of enl.
    
    # The set of pairs form a binary tree which we fill in à la Pascal
    
    # The root pair is ('LT','TL')
    # The pairs to the extreme left are ('Ln T', 'T Ln')
    # The pairs on the extreme right are ('L Tn', 'Tn L')

    # The children of (G, D) are
    # to the left (P, Q) with P=G[:-1]+'LT' (on enlève 'T' on rajoute 'LT') and Q=D+'L'
    # to the right (R, S) with R=G+'T' and S=D[:-1]+'TL' (on enlève 'L' on rajoute 'TL')
    
    # Note the properties : 'G' ends with 'T' and 'D' ends with 'L' are preserved
    # which is why G[:-1] = G - 'T' and D[:-1] = D - 'L'
    """
 
    if len(current_pair[0]) > max_word_length:
        return []
 
    pair_left = current_pair[0][:-1]+'LT', current_pair[1]+'L'
    pair_right = current_pair[0]+'T', current_pair[1][:-1]+'TL'
 
    return [current_pair] + linking_patterns(max_word_length, pair_left) +\
                            linking_patterns(max_word_length, pair_right)
 
def occ(P,A):
    # Returns the number of times P appears at the begining of circular shifts of A.
    shifts = mT.list_of_circular_shifts(A)
    counter = 0
    n,l = len(P), len(A)
 
    for shift in shifts:
        power = shift + shift * (n//l)
        if power[:n] == P:
            counter += 1
 
    return counter
 
def scal_PQ(P,Q,A,B):
    return (occ(P,A)*occ(Q,B)+ occ(P,B)*occ(Q,A))
 
def enl(A,B):
    # Returns the enl metric on the words A and B in the L/T alphabet.
    patterns = linking_patterns(len(A)+len(B)+1)
    return sum([scal_PQ(P,Q,A,B) for P,Q in patterns])

 
def linking_matrix(root, size):
    # renvoie matrice d'enlacement entre classes de long combinatoire size
    addresses = mT.list_of_addresses(root, lambda n: len(n.address) == size)
    conj_class = {}
    
    for a in addresses:
        rep = mT.circular_min_rpz(a)
        conj_class[rep] = True
        
    print(len(addresses), len(conj_class))
    
    all_reps = sorted(conj_class.keys())
    #all_reps = conj_class.keys()
    
    
    enl_mat = []
    for rep1 in conj_class:
        enl_mat.append([])
        for rep2 in conj_class:
            enl_mat[-1].append(cross_word(rep1,rep2))
            
    enl_mat = np.array(enl_mat)
    return all_reps, enl_mat

## L'algo de Pierre pour le crossing number

def orderlex(l1, l2, n1, n2):
    """ lexicographic order on infinite words, returns : 
    0 if the n1-th shift of l1^infty is smaller than the n2-th shift of l2^infty, 1 if it is larger, and 2 if the two words coincide 
    (this can be determined by looking only length(l1)+length(l2) letters since u^infty = v^infty iff uv = vu).
    """
    i = 0
    while( (l1[(i+n1)%len(l1)] == l2[(i+n2)%len(l2)]) & (i <= (len(l1)+len(l2))) ):
        i = i+1
    if l1[(i+n1)%len(l1)] < l2[(i+n2)%len(l2)]:
        return 0
    else:
        if l1[(i+n1)%len(l1)] > l2[(i+n2)%len(l2)]:
            return 1
        else:
            return 2

def cross(l1, l2):
    """ The function cross computes the crossing number of l1 and l2 on the template.
    It relies on the observation that a crossing occurs 
    when an arc coming from the left ear (the 0) goes to the right 
    of an arc coming from the right ear (the 1).
    """
    c = 0
    for i in range(len(l1)):
        for j in range(len(l2)):
            if ( (l1[i] == 0) & (l2[j] == 1) & (orderlex(l1,l2,i+1,j+1) == 1) ):
                c = c+1
            if ( (l1[i] == 1) & (l2[j] == 0) & (orderlex(l1,l2,i+1,j+1) == 0) ):
                c = c+1
    return c

def cross_word(w1, w2):
    """ Same as before with L and T words """
    c = 0
    for i in range(len(w1)):
        for j in range(len(w2)):
            if ( (w1[i] == 'L') & (w2[j] == 'T') & (orderlex(w1,w2,i+1,j+1) == 1) ):
                c = c+1
            if ( (w1[i] == 'T') & (w2[j] == 'L') & (orderlex(w1,w2,i+1,j+1) == 0) ):
                c = c+1
    return c

## On teste l'égalité entre ma formule pour enl et l'algo de Pierre

def bin_list_from_word(word):
    # We translate our L and T words into lists of 0 and 1
    liste = []
    for letter in word:
        if letter == 'L':
            liste.append(0)
        elif letter == 'T':
            liste.append(1)
        else :
            raise ValueError("Not an L & T word")
    return liste

def cross_word_p(w1,w2):
	return cross(bin_list_from_word(w1), bin_list_from_word(w2))

bin_list_from_word('LLTLT')

# We check that cross = enl on all_reps, it is always true : youpii !!

"""
for i,rep1 in enumerate(all_reps):
    for j,rep2 in enumerate(all_reps):
        cross, enlac = cross(bin_list_from_word(rep1),bin_list_from_word(rep2)), 2*mT.enl(rep1,rep2)
        print(rep1,rep2,cross, int(enlac), enlac == cross)

"""

# Another (better) test

"""
k1='LLLLLLLTTLT'
k2='LLLLLTLLLLT'
l1 = lk.bin_list_from_word(k1)
l2 = lk.bin_list_from_word(k2)

height = 4
root = mT.construct_tree(mT.get_prop_max_height(height))
words_bound_length = [node.address for node in mT.list_of_nodes_with_prop(root, lambda node : node.get_trace()>2)]

for k0 in words_bound_length:
    l0 = lk.bin_list_from_word(k0)
    print(k0, (lk.cross(l0,l1), lk.enl(k0,k1)))#, mT.Rademacher(k0),mT.Rademacher(k1)))
    #print(k0, (lk.cross(l0,l2), lk.enl(k0,k2)))#, mT.Rademacher(k0),mT.Rademacher(k2)))

"""