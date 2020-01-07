## BIBLIS
 
import numpy as np
import matplotlib.pyplot as plt
import cypari

import matrixTree as mT
import enlacement as lk
##


def group_trace_words_by_fricke_Rad(trace):
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

###### Questions ######

Does same Fricke imply :
"same genus" ? This is algebraically plausible.
"same number of syllables" ? I think so.
"form a subgroup or a submonoid ?"

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