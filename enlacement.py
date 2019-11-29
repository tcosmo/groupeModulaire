

## BIBLIS
 
import numpy as np
import matplotlib.pyplot as plt
import cypari

import matrixTree as mT

## Du code de Tristan

def get_addresses(root, predicate):
    all_nodes = mT.get_nodes_with_property(root, predicate)
    return [n.address for n in all_nodes]

# matrice d'enlacement 
def get_enl_matrix(root, size):
    addresses = get_addresses(root, lambda n: len(n.address) == size)
    conj_class = {}
    
    for a in addresses:
        rep = mT.get_circular_rep(a)
        conj_class[rep] = True
        
    print(len(addresses), len(conj_class))
    
    all_reps = sorted(conj_class.keys())
    
    enl_mat = []
    for rep1 in conj_class:
        enl_mat.append([])
        for rep2 in conj_class:
            enl_mat[-1].append(mT.enl(rep1,rep2))
            
    enl_mat = np.array(enl_mat)
    return all_reps, enl_mat

all_reps, enl_mat = get_enl_matrix(root,7)

# un test
for i,rep1 in enumerate(all_reps):
    for j,rep2 in enumerate(all_reps):
        if rep1.count('L') >= 1 and rep2.count('L') >= 1 and rep1.count('T') >= 1 and rep2.count('T') >= 1:
            print("enl({},{})={}".format(rep1,rep2,enl_mat[i,j]))

## L'algo de Pierre

"""
Orbits are coded by sequences of 0's and 1's (instead of L and T).
Here are three orbits or period 5, 5 and 7.
"""
gamma1 = [0,0,1,0,1]
gamma2 = [0,1,0,1,1]
gamma3 = [0,0,1,0,1,0,1]
reste_euclid_div = 10%3

"""
The function ordrelex implements the lexicographic order on infinite words: it says 
0 if the n1-th shift of l1^infty is smaller than the n2-th shift of l2^infty, 
1 if it is larger, and
2 if the two words coincide 
(this can be determined only by looking length(l1)+length(l2) letters since u^infty = v^infty iff uv = vu).
"""

def ordrelex(l1, l2, n1, n2):
    i = 0
    while( (l1[(i+n1)%len(l1)] == l2[(i+n2)%len(l2)]) & (i <= (len(l1)+len(l2))) ):
        i = i+1;
    if l1[(i+n1)%len(l1)] < l2[(i+n2)%len(l2)]:
        return 0
    else:
        if l1[(i+n1)%len(l1)] > l2[(i+n2)%len(l2)]:
            return 1
        else:
            return 2


# We check that gamma1 and gamma1 coincice, while gamma1 is smaller than gamma2.
print(ordrelex(gamma1,gamma1,0,0))
print(ordrelex(gamma1,gamma2,0,0))
print(ordrelex(gamma2,gamma1,0,0))

"""
The function cro computes the crossing number of l1 and l2 on the template.
It relies on the observation that a crossing occurs 
when an arc coming from the left ear (the 0) goes to the right 
of an arc coming from the right ear (the 1).
"""

def cro(l1, l2):
    c = 0
    for i in range(len(l1)):
        for j in range(len(l2)):
            if ( (l1[i] == 0) & (l2[j] == 1) & (ordrelex(l1,l2,i+1,j+1) == 1) ):
                c = c+1
            if ( (l1[i] == 1) & (l2[j] == 0) & (ordrelex(l1,l2,i+1,j+1) == 0) ):
                c = c+1
    return c


## On teste l'égalité entre ma formule pour enl et l'algo de Pierre

""" 
We translate our L and T words into lists of 0 and 1
"""

def bin_list_from_word(word):
    liste = []
    for letter in word:
        if letter == 'L':
            liste.append(0)
        elif letter == 'T':
            liste.append(1)
        else :
            raise ValueError("Not an L & T word")
    return liste

bin_list_from_word('LLTLT')

# the folowing should give 4 and 8
mT.enl('LLTLT','LTLTT'), cro(bin_list_from_word('LLTLT'), bin_list_from_word('LTLTT'))

# we check that cro = 2*enl on all_reps, it is always true : youpii !!

for i,rep1 in enumerate(all_reps):
    for j,rep2 in enumerate(all_reps):
        cross, enlac = cro(bin_list_from_word(rep1),bin_list_from_word(rep2)), 2*mT.enl(rep1,rep2)
        print(rep1,rep2,cross, int(enlac), enlac == cross)


