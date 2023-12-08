import itertools
import dendropy
import msprime
import numpy as np


def find_confidence_interval(L_arr, bp=None, CI=0.95):
    # print("start finding")
    if bp is None:
        max_s = np.argmax(L_arr)
    else:
        max_s = bp
    l_pointer = int(max_s)
    r_pointer = int(max_s + 1)
    while True:
        if np.sum(L_arr[l_pointer:r_pointer]) >= CI or (l_pointer <= 0 and 999 <= r_pointer):
            break
        elif np.sum(L_arr[l_pointer:r_pointer]) < CI:
            if 0 < l_pointer and r_pointer < 999:
                if L_arr[l_pointer-1] > L_arr[r_pointer+1]:
                    l_pointer -= 1
                elif L_arr[l_pointer-1] < L_arr[r_pointer+1]:
                    r_pointer += 1
                else:
                    l_pointer -= 1
                    r_pointer += 1
            elif l_pointer <= 0:
                r_pointer += 1
            elif 999 <= r_pointer:
                l_pointer -= 1
    # print("end finding")

    # bp_CI = r_pointer - l_pointer

    return l_pointer, r_pointer


def merge_intervals(intervals):

    intervals.sort(key = lambda x: x[0])
    intervals = list(intervals for intervals,_ in itertools.groupby(intervals))
    p = intervals[0]
    i = 1
    while i < len(intervals):
        c = intervals[i]
        if c[0] <= p[1]:
            intervals[i-1][1] = max(p[1], c[1])
            del intervals[i]
        else:
            p = c
            i += 1
    return intervals


def recomb_segments(intervals, genome_length):
    intervals.sort(key = lambda x:x[0])
    recomb_intervals = []
    if intervals[0][0] > 0:
        recomb_intervals.append([0, intervals[0][0]])
    i = 0
    while i < len(intervals) - 1:
        if intervals[i + 1][0] > intervals[i][1]:
            recomb_intervals.append([intervals[i][1], intervals[i+1][0]])
        i += 1
    if intervals[-1][1] < genome_length:
        recomb_intervals.append([intervals[-1][1] , genome_length])
    return recomb_intervals


def search_mapped_interval(interval, interval_list, left=True, right=True):
    idx = 0
    while idx < len(interval_list):
        temp_interval = interval_list[idx]
        if temp_interval[0] <= interval[0] <= interval[1] <= temp_interval[1]:
            return idx
        idx += 1

    return -1


def generate_a_random_tree(sample_size, length):

    # Simulate an ancestral history for 3 diploid samples under the coalescent
    # with recombination on a 5kb region with human-like parameters.
    ts = msprime.simulate(sample_size=sample_size, Ne=0.5, length=length, record_full_arg=True)
    # Visualise the simulated ancestral history.
    for tree in ts.trees():
        return tree.as_newick()


def tree_dist(tree1, tree2):

    tree_list = [tree1, tree2]
    tns = dendropy.TaxonNamespace()
    dendropy_tree_list = list(map(lambda x: dendropy.Tree.get(data=x, schema="newick",
                                                              taxon_namespace=tns), tree_list))
    rf_dist = dendropy.calculate.treecompare.symmetric_difference(dendropy_tree_list[0], dendropy_tree_list[1])
    return rf_dist

