import numpy as np
import matplotlib.pyplot as plt
plt.style.use('default')

import matrixTree as mT
import plotting_routines as pR

inversibles = ["L38 ", "L3 R3 L2 R", "L R L R12", "L3 R L2 R3", "L19 R2", "L2 R7 L R", "L5 R L4 R", "L2 R L R7"]
#pR.plot_matrix_per_trace(40)
#pR.plot_matrix_per_trace_with_conj_class(40)
pR.plot_matrix_per_trace_with_conj_class(40, inversibles=inversibles, 
                                          only_inverse=False, figsize=(20,20))