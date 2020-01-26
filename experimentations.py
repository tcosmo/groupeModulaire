## BIBLIS
 
import numpy as np
import matplotlib.pyplot as plt
import cypari

import matrixTree as mT
import enlacement as lk
##


def group_trace_words_by_Fricke_Rad(trace):
	nodes_by_trace=mT.list_of_nodes_with_trace_equal(trace)
	dic_Fricke = dict()
	for node in nodes_by_trace :
	    node_Fricke = mT.Fricke(node.address)
	    if  node_Fricke not in dic_Fricke:
	        dic_Fricke[node_Fricke]={mT.circular_min_rpz(node.address)}
	    else:
	        dic_Fricke[node_Fricke].add(mT.circular_min_rpz(node.address))
	for key in dic_Fricke:
	    print("\n", key)
	    for word in dic_Fricke[key]:
	        print(word, "\t", mT.compress(word), "\t", mT.Rademacher(word))


from importlib import reload
reload(mT)

def compare_linking_Rademacher_resultant_Killing(height):
	root = mT.construct_tree(mT.get_prop_max_height(height))
	hyp_nodes = mT.list_of_nodes_with_prop(root, lambda node : node.get_trace()>2)
	for node_0 in hyp_nodes:
	    word_0 = node_0.address
	    mat_0  = node_0.matrix
	    lst_0 = lk.bin_list_from_word(word_0)
	    for node_1 in hyp_nodes:
	        word_1 = node_1.address
	        mat_1  = node_1.matrix
	        lst_1 = lk.bin_list_from_word(word_1)
	        #if not mT.circular_min_rpz(word_0) == mT.circular_min_rpz(word_1) :
	        print(word_0, word_1)
	        print('link=', lk.cross(lst_0, lst_1), '\t', \
	              'kill=', mT.killing_form(mat_0, mat_1), '\t', \
	              'res=', mT.resultant(mat_0,mat_1))
	        #print(mat_0.trace(), mat_1.trace(), (mat_0*mat_1).trace())
	        #comut = mat_0 * mat_1 * mT.matrix_sl_inverse(mat_0) * mT.matrix_sl_inverse(mat_1)
	        #print(comut.trace())#, (mat_0*mat_1).trace())
	        print('rad0=', mT.Rademacher(word_0), '\t', 'rad1=', mT.Rademacher(word_1))

"""
###### Negative results ######

Same Fricke does NOT implie :
    same element (not even up to V4-symetries)
    same Rademacher (not even up to sign)
    same knot
    same set of exponents

Beign Fricke equivalent to an ambiguous does not implie ambiguous.

Those two different words of trace 40 have same Fricke, but 
different knots T(2,9) & T(2,3) ; different |Rademacher| ; different set of exponents one is ambiguous and not the other :
LLLLLTLLLLT 	 L5 T L4 T 	 Rad = 7
LLLLLLLTTLT 	 L7 T2 L T 	 Rad = 5

Same genus (and even in the principal genus) does not imply same linking functions.
For instance L38 T and L19 T 2 are in the principal genus for trace 40 
but here their cross functions against a sample of small words :
[('TTL', 2), ('TL', 2), ('TLT', 2), ('TLL', 4), ('LT', 2), ('LTT', 2), ('LTL', 4), ('LLT', 4)]
[('TTL', 4), ('TL', 4), ('TLT', 4), ('TLL', 6), ('LT', 4), ('LTT', 4), ('LTL', 6), ('LLT', 6)]

###### Positive results ######

Fricke is bi-unitaire (q^{2d}+...+1)/q^d = q^d+...+q^-d
Degree of Fricke is d = lenght of word L&T

Test if : roots of Fricke are the roots of reduced representatives of quadratic forms
Impossible since Fricke does not imply same class
Test if : linking is resultant of Fricke
Why not, but now its less compelling...

###### Questions ######

Does same Fricke imply :
"same genus" ?
"same number of syllables" ? I think so.

Have not been invalidated :
"same Fricke and |Rademacher| implies equal up to V4-symmetry"
"same Fricke and |Rademacher| implies same genus"


###### Experimental data ######

Those two have same Fricke, trace = 40, but not same exponents :
LLLLLLLTTLT 	 L7 T2 L T 	 Rad = 5
LLLLLTLLLLT 	 L5 T L4 T 	 Rad = 7

These two have same Fricke, trace = 128, but not same exponents :
LLLLTLLLLTLLT 	 L4 T L4 T L2 T 	 Rad = 7
LLLLLLTTLTLLT 	 L6 T2 L T L2 T 	 Rad = 5

These two have same Fricke, trace = 130, but not same exponents :
LLLLLLLTLTTTLT 	 L7 T L T3 L T 	 Rad = 4
LLLLLTLLLLLTLT 	 L5 T L5 T L T 	 Rad = 8

These two have same Fricke, trace = 130, but not same exponents :
LLLLTLLLTLLLT 	 L4 T L3 T L3 T 	 Rad = 7
LLLLLLTLTLLTT 	 L6 T L T L2 T2 	 Rad = 5

Those two have same Fricke, trace = 205, but not same exponents :
LLLLLTTTTLLTLT 	 L5 T4 L2 T L T 	 Rad = 2
LLLLTTTTLTTTLT 	 L4 T4 L T3 L T 	 Rad = -2
"""