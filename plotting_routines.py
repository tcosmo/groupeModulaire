import numpy as np
import matplotlib.pyplot as plt
import matrixTree as mT

def plot_matrix_per_trace(t, figsize=(10,10), right_column=False, save=True):
    """ Plot all matrices with trace t.
    """
    all_nodes = mT.get_nodes_with_trace(t)
    all_coords = np.array([n.get_coord(right_column=right_column) for n in all_nodes])
    plt.figure(figsize=figsize)
    plt.scatter(all_coords[:,0], all_coords[:,1])
    plt.title("All matrices with trace t={}".format(t))
    if save:
        plt.savefig("images/matrix_scatter_per_trace/matrices_for_t={}.png".format(t))
    plt.show()

def plot_matrix_per_trace_with_conj_class(t, right_column=False, inversibles=[], only_inverse=False, 
                                            figsize=(10,10), save=True, special_cases=[]):
    """ Plots all matrices with traces t, differenting them by their 
        conjugaison classes. If inversible are given, their label is changed.
        If only_inverse is true, plots only inversible.
    """

    nb_classes, conj_class, reps = mT.get_conjugaison_classes(t,right_column)
    
    nb_i = 0
    plt.figure(figsize=figsize)
    for class_nb in range(nb_classes):
        filter_ = conj_class[:,2] == class_nb
        style = "o"

        is_inversible = (mT.compress(reps[class_nb]) in inversibles)
        if is_inversible:
            nb_i += 1
            style = "*"

        for case,marker in special_cases:
            if case(reps[class_nb]):
                style=marker

        if (not only_inverse) or is_inversible:
            plt.scatter(conj_class[filter_,0],conj_class[filter_,1], marker=style, label="Class of {}".format(mT.compress(reps[class_nb])))
    
    s = ""
    if nb_i != 0:
        s = ", {} inverses, ".format(nb_i)

    plt.title("The {} conj classes".format(nb_classes)+ s 
              +" of all matrices with trace t={}".format(t))
    plt.legend()
    if save:
        s = ""
        if only_inverse:
            s = "only_inverse_"
        plt.savefig("images/matrix_scatter_per_trace_and_conj/matrices_conj_"+s+"for_t={}.png".format(t))
    plt.show()