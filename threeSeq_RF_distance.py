# -*- coding: utf-8 -*-
"""
-------------------------------------------------
   File Name：     threeSeq_RF_distance
   Description :
   Author :       shicen
   date：          6/23/23
-------------------------------------------------
   Change Activity:
                   6/23/23:
-------------------------------------------------
"""
import math

import numpy as np
from dendropy.calculate import treecompare
import dendropy
from analysis_treeSeq import analysis_treeSeq
from pathlib import Path
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
import os
import sys
import subprocess
import copy
import msprime
from IPython.display import SVG, display
import itertools
import utility


def try_run_cmd(cmd):
    try:
        output = subprocess.check_output(cmd, shell=True, stderr=subprocess.STDOUT)
        # sys.stdout.write(output.decode("UTF-8"))
    except subprocess.CalledProcessError:
        print('Execution of "%s" failed!\n' % cmd)
        sys.exit(1)


def all_intervals(exp_id, iter_i):

    threeSeq_file = f"./{exp_id}/iter_{iter_i}/{exp_id}.3s.rec"
    intervals = []
    with open(threeSeq_file, 'r') as f:
        l = f.readline()
        l = f.readline()
        while l:
            columns = l.strip().split('\t')
            bp_lst = columns[12:]
            for bps in bp_lst:  # bps:  "1-9 & 631-631"
                bp_ranges = bps.split(" & ")  # bp_ranges = ['1-9', '631-631']
                for i in range(2):
                    bp_range = bp_ranges[i]
                    lower_b = min(max(int(bp_range.split('-')[0]), 0), 1000)
                    upper_b = min(int(bp_range.split('-')[1]), 1000)
                    intervals.append([lower_b, upper_b])
            l = f.readline()

    return intervals


def generate_trees(exp_id, iter_i, interval, symbol, local):

    Path(exp_id, f"iter_{iter_i}", "consv").mkdir(parents=True, exist_ok=True)
    for p in Path(exp_id, f"iter_{iter_i}").glob(f"threeseq_{symbol}_alignment_*"):
        p.unlink(missing_ok=True)
    for p in Path(exp_id, f"iter_{iter_i}").glob(f"threeseq_T_alignment_*"):
        p.unlink(missing_ok=True)
    for p in Path(exp_id, f"iter_{iter_i}").glob(f"True_T_alignment_*"):
        p.unlink(missing_ok=True)

    for i in range(len(interval)):
        seq_file = Path(exp_id, f"iter_{iter_i}", "concatenated_seq.fasta")
        seq_rec_list = list(SeqIO.parse(seq_file, "fasta"))
        seq_list = []
        for seq_rec in seq_rec_list:
            seq = seq_rec.seq[interval[i][0]:interval[i][1]]
            seq_list.append(SeqRecord(seq, id=seq_rec.id, description=""))

        if symbol == "R":
            SeqIO.write(seq_list, Path(exp_id, f"iter_{iter_i}", "consv", f"threeseq_{symbol}_alignment_{i}.fasta"), 'fasta')

        if local:
            if symbol == "R":
                cmd_raxml = '~/Downloads/raxml-ng_v1.2.0_macos_M1/raxml-ng --msa ' + str(
                    Path(exp_id, f"iter_{iter_i}", "consv", f"threeseq_{symbol}_alignment_{i}.fasta")) + ' --model GTR+G --redo'
                try_run_cmd(cmd_raxml)
        else:
            if symbol == "R":
                cmd_raxml = '~/raxml-ng/bin/raxml-ng --msa ' + str(
                    Path(exp_id, f"iter_{iter_i}", "consv", f"threeseq_{symbol}_alignment_{i}.fasta")) + ' --model GTR+G --redo'
                try_run_cmd(cmd_raxml)


def generate_a_random_tree(sample_size, length):

    # Simulate an ancestral history for 3 diploid samples under the coalescent
    # with recombination on a 5kb region with human-like parameters.
    ts = msprime.simulate(sample_size=sample_size, Ne=0.5, length=length, record_full_arg=True)
    # Visualise the simulated ancestral history.
    for tree in ts.trees():
        return (tree.as_newick())


def cal_RF_distance_random(exp_id, iter_i, breakpoint_list, recomb_intervals, T_non_recomb, uninformative_intervals,
                    genome_length, sample_size=10,
                    symbol_R="R", symbol_T="T"):
    weighted_RF_distance = 0
    weighted_RF_distance_by_random = 0
    weighted_RF_distance_np = 0
    # subalignment_length = []
    subalignment_RF_dist_weighted = []
    subalignment_RF_dist = []
    tmp_RF_dist_weighted = 0
    tmp_RF_dist = 0
    flag_in_uninfo = 1
    for i in range(len(breakpoint_list) - 1):
        c_interval = [breakpoint_list[i], breakpoint_list[i + 1]]

        inf_recomb_interval_idx = utility.search_mapped_interval(c_interval, recomb_intervals)
        true_recomb_interval_idx = utility.search_mapped_interval(c_interval, T_non_recomb)

        print(c_interval, inf_recomb_interval_idx, true_recomb_interval_idx)

        with open(
                Path(exp_id, f"iter_{iter_i}", f"tree_{true_recomb_interval_idx}.tre")) as f:
            t_tree = f.readline().strip()

        if inf_recomb_interval_idx != -1:
            with open(Path(exp_id, f"iter_{iter_i}", "consv",
                           f"threeseq_{symbol_R}_alignment_{inf_recomb_interval_idx}.fasta.raxml.bestTree")) as f:
                r_tree = f.readline().strip()
                tree_list = [t_tree, r_tree]
            flag = 1
        elif utility.search_mapped_interval(c_interval, uninformative_intervals) != -1:
            random_newick = generate_a_random_tree(sample_size, length=breakpoint_list[i + 1] - breakpoint_list[i])
            tree_list = [t_tree, random_newick]
            flag = 0
        else:
            raise Exception("didn't find interval")

        # print(tree_list)
        tns = dendropy.TaxonNamespace()
        dendropy_tree_list = list(map(lambda x: dendropy.Tree.get(data=x, schema="newick",
                                                                  taxon_namespace=tns), tree_list))
        rf_dist = dendropy.calculate.treecompare.symmetric_difference(dendropy_tree_list[0], dendropy_tree_list[1])

        # curr_RF_distance = rf_dist * (breakpoint_list[i + 1] - breakpoint_list[i]) / genome_length
        curr_RF_distance = (rf_dist / (2 * (sample_size - 2))) * (breakpoint_list[i + 1] - breakpoint_list[i]) / genome_length
        weighted_RF_distance += curr_RF_distance
        print(curr_RF_distance)
        if flag == 0:
            weighted_RF_distance_by_random += curr_RF_distance
            if i != 0:
                if flag_in_uninfo == 0:
                    subalignment_RF_dist_weighted.append(tmp_RF_dist_weighted)
                    subalignment_RF_dist.append(tmp_RF_dist)
                    tmp_RF_dist_weighted = 0
                    tmp_RF_dist = 0
                    flag_in_uninfo = 1
        else:
            flag_in_uninfo = 0
            curr_RF_distance2 = (rf_dist / (2 * (sample_size - 2)))
            weighted_RF_distance_np += curr_RF_distance2
            tmp_RF_dist_weighted += curr_RF_distance
            tmp_RF_dist += curr_RF_distance2

    if flag_in_uninfo == 0:
        subalignment_RF_dist_weighted.append(tmp_RF_dist_weighted)
        subalignment_RF_dist.append(tmp_RF_dist)

    print(f"weighted_RF_distance: {weighted_RF_distance}, weighted_RF_distance_np: {weighted_RF_distance_np},\
            subalignment_RF_dist_weighted: {subalignment_RF_dist_weighted}, subalignment_RF_dist: {subalignment_RF_dist}")
    return [weighted_RF_distance, weighted_RF_distance_np,
            subalignment_RF_dist_weighted, subalignment_RF_dist]


def cal_RF_distance_neighbor(exp_id, iter_i, breakpoint_list, recomb_intervals, T_non_recomb, uninformative_intervals,
                    genome_length, sample_size=10,
                    symbol_R="R"):
    weighted_RF_distance = 0
    weighted_RF_distance_by_neighbor = 0
    weighted_RF_distance_np = 0
    subalignment_RF_dist_weighted = []
    subalignment_RF_dist = []
    tmp_RF_dist_weighted = 0
    tmp_RF_dist = 0
    flag_in_uninfo = 1
    most_recent_inf_idx = -1
    most_recent_true_idx = -1

    for i in range(len(breakpoint_list) - 1):
        c_interval = [breakpoint_list[i], breakpoint_list[i + 1]]

        inf_recomb_interval_idx = utility.search_mapped_interval(c_interval, recomb_intervals)
        true_recomb_interval_idx = utility.search_mapped_interval(c_interval, T_non_recomb)
        most_recent_inf_idx = inf_recomb_interval_idx if inf_recomb_interval_idx != -1 else most_recent_inf_idx

        print(c_interval, inf_recomb_interval_idx, true_recomb_interval_idx)

        with open(
                Path(exp_id, f"iter_{iter_i}", f"tree_{true_recomb_interval_idx}.tre")) as f:
            t_tree = f.readline().strip()

        if inf_recomb_interval_idx != -1:
            with open(Path(exp_id, f"iter_{iter_i}", "consv",
                           f"threeseq_{symbol_R}_alignment_{inf_recomb_interval_idx}.fasta.raxml.bestTree")) as f:
                r_tree = f.readline().strip()

            rf_dist = utility.tree_dist(t_tree, r_tree)

            curr_dist = (rf_dist / (2 * (sample_size - 2))) * (c_interval[1] - c_interval[0] + 1) / genome_length
            curr_dist2 = (rf_dist / (2 * (sample_size - 2)))
            print(curr_dist)
            weighted_RF_distance += curr_dist
            weighted_RF_distance_np += curr_dist2

            flag_in_uninfo = 0
            tmp_RF_dist_weighted += curr_dist
            tmp_RF_dist += curr_dist2

        else:
            uninform_interval_idx = utility.search_mapped_interval(c_interval, uninformative_intervals)
            if uninform_interval_idx == -1:
                raise Exception("didn't find interval")

            if most_recent_inf_idx == -1:
                with open(Path(exp_id, f"iter_{iter_i}", "consv", f"threeseq_{symbol_R}_alignment_0.fasta.raxml.bestTree")) as f:
                    r_tree = f.readline().strip()

                rf_dist = utility.tree_dist(t_tree, r_tree)
                # curr_dist = rf_dist * (c_interval[1] - c_interval[0] + 1) / genome_length
                curr_dist = (rf_dist / (2 * (sample_size - 2))) * (c_interval[1] - c_interval[0] + 1) / genome_length
                print(curr_dist)
                weighted_RF_distance += curr_dist
                weighted_RF_distance_by_neighbor += curr_dist
                continue

            if flag_in_uninfo == 0:
                subalignment_RF_dist_weighted.append(tmp_RF_dist_weighted)
                subalignment_RF_dist.append(tmp_RF_dist)
                tmp_RF_dist_weighted = 0
                tmp_RF_dist = 0
                flag_in_uninfo = 1

            if most_recent_inf_idx == len(recomb_intervals) - 1:
                with open(Path(exp_id, f"iter_{iter_i}", "consv",
                               f"threeseq_{symbol_R}_alignment_{len(recomb_intervals) - 1}.fasta.raxml.bestTree")) as f:
                    r_tree = f.readline().strip()

                rf_dist = utility.tree_dist(t_tree, r_tree)
                # curr_dist = rf_dist * (c_interval[1] - c_interval[0] + 1) / genome_length
                curr_dist = (rf_dist / (2 * (sample_size - 2))) * (c_interval[1] - c_interval[0] + 1) / genome_length
                print(curr_dist)
                weighted_RF_distance += curr_dist
                weighted_RF_distance_by_neighbor += curr_dist
                continue

            curr_uninformative_interval = uninformative_intervals[uninform_interval_idx]
            curr_uninformative_interval_start = curr_uninformative_interval[0]
            curr_uninformative_interval_end = curr_uninformative_interval[1]
            curr_uninformative_interval_mid = (curr_uninformative_interval_start + curr_uninformative_interval_end) / 2
            no_break_in_window_flag = 0

            if c_interval[0] == recomb_intervals[most_recent_inf_idx][1] and \
                    c_interval[1] == recomb_intervals[most_recent_inf_idx + 1][0]:
                assert c_interval[0] == curr_uninformative_interval_start and \
                       c_interval[1] == curr_uninformative_interval_end
                with open(Path(exp_id, f"iter_{iter_i}", "consv",
                               f"threeseq_{symbol_R}_alignment_{most_recent_inf_idx}.fasta.raxml.bestTree")) as f:
                    r_tree = f.readline().strip()
                rf_dist = utility.tree_dist(t_tree, r_tree)
                curr_dist = (rf_dist / (2 * (sample_size - 2))) * (
                        c_interval[1] - c_interval[0] + 1) / genome_length / 2
                weighted_RF_distance += curr_dist
                weighted_RF_distance_by_neighbor += curr_dist
                with open(Path(exp_id, f"iter_{iter_i}","consv",
                               f"threeseq_{symbol_R}_alignment_{most_recent_inf_idx + 1}.fasta.raxml.bestTree")) as f:
                    r_tree = f.readline().strip()
                rf_dist = utility.tree_dist(t_tree, r_tree)
                # curr_dist = rf_dist * (c_interval[1] - c_interval[0] + 1) / genome_length / 2
                curr_dist = (rf_dist / (2 * (sample_size - 2))) * (
                        c_interval[1] - c_interval[0] + 1) / genome_length / 2
                print(curr_dist)
                weighted_RF_distance += curr_dist
                weighted_RF_distance_by_neighbor += curr_dist
            else:
                if curr_uninformative_interval_start <= c_interval[0] <= c_interval[1] <= \
                        curr_uninformative_interval_mid:
                    with open(Path(exp_id, f"iter_{iter_i}", "consv",
                                   f"threeseq_{symbol_R}_alignment_{max(most_recent_inf_idx,0)}.fasta.raxml.bestTree")) as f:
                        r_tree = f.readline().strip()
                    rf_dist = utility.tree_dist(t_tree, r_tree)
                    # curr_dist = rf_dist * (c_interval[1] - c_interval[0] + 1) / genome_length / 2
                    curr_dist = (rf_dist / (2 * (sample_size - 2))) * (
                            c_interval[1] - c_interval[0] + 1) / genome_length
                    print(curr_dist)
                    weighted_RF_distance += curr_dist
                    weighted_RF_distance_by_neighbor += curr_dist
                elif curr_uninformative_interval_start <= c_interval[0] <= \
                        curr_uninformative_interval_mid <= c_interval[1]:
                    with open(Path(exp_id, f"iter_{iter_i}","consv",
                                   f"threeseq_{symbol_R}_alignment_{max(most_recent_inf_idx,0)}.fasta.raxml.bestTree")) as f:
                        r_tree = f.readline().strip()
                    rf_dist = utility.tree_dist(t_tree, r_tree)
                    # curr_dist = rf_dist * (c_interval[1] - c_interval[0] + 1) / genome_length / 2
                    curr_dist = (rf_dist / (2 * (sample_size - 2))) * (
                            curr_uninformative_interval_mid - c_interval[0] + 1) / genome_length
                    print(curr_dist)
                    weighted_RF_distance += curr_dist
                    weighted_RF_distance_by_neighbor += curr_dist
                    with open(Path(exp_id, f"iter_{iter_i}", "consv",
                                   f"threeseq_{symbol_R}_alignment_{min(most_recent_inf_idx,len(recomb_intervals)-1)}.fasta.raxml.bestTree")) as f:
                        r_tree = f.readline().strip()
                    rf_dist = utility.tree_dist(t_tree, r_tree)
                    # curr_dist = rf_dist * (c_interval[1] - c_interval[0] + 1) / genome_length / 2
                    curr_dist = (rf_dist / (2 * (sample_size - 2))) * (
                            c_interval[1] - curr_uninformative_interval_mid + 1) / genome_length
                    print(curr_dist)
                    weighted_RF_distance += curr_dist
                    weighted_RF_distance_by_neighbor += curr_dist
                elif curr_uninformative_interval_start <= curr_uninformative_interval_mid <= c_interval[0] <= \
                         c_interval[1] <= curr_uninformative_interval_end:
                    with open(Path(exp_id, f"iter_{iter_i}", "consv",
                                   f"threeseq_{symbol_R}_alignment_{min(most_recent_inf_idx,len(recomb_intervals)-1)}.fasta.raxml.bestTree")) as f:
                        r_tree = f.readline().strip()
                    rf_dist = utility.tree_dist(t_tree, r_tree)
                    # curr_dist = rf_dist * (c_interval[1] - c_interval[0] + 1) / genome_length / 2
                    curr_dist = (rf_dist / (2 * (sample_size - 2))) * (
                            c_interval[1] - c_interval[0] + 1) / genome_length
                    print(curr_dist)
                    weighted_RF_distance += curr_dist
                    weighted_RF_distance_by_neighbor += curr_dist
                else:
                    raise Exception("Interesting case!")

    if flag_in_uninfo == 0:
        subalignment_RF_dist_weighted.append(tmp_RF_dist_weighted)
        subalignment_RF_dist.append(tmp_RF_dist)

    print(f"weighted_RF_distance: {weighted_RF_distance}, weighted_RF_distance_np: {weighted_RF_distance_np},\
            subalignment_RF_dist_weighted: {subalignment_RF_dist_weighted}, subalignment_RF_dist: {subalignment_RF_dist}")
    return [weighted_RF_distance, weighted_RF_distance_np,
            subalignment_RF_dist_weighted, subalignment_RF_dist]


def threeSeq_RF_distance(exp_id, iter_i, genome_length, bp_list, sample_size, local, rege_tree):

    intervals = all_intervals(exp_id, iter_i)
    merged_intervals = utility.merge_intervals(intervals)

    left_incl_right_excl_merged_intervals = copy.deepcopy(merged_intervals)
    for mi in left_incl_right_excl_merged_intervals:
        if mi[1] + 1 <= 1000:
            mi[1] += 1

    recomb_intervals = utility.recomb_segments(left_incl_right_excl_merged_intervals, genome_length)
    subalignment_length_list = []
    no_missed_bp_list = []
    for ri in recomb_intervals:
        subalignment_length_list.append(ri[1] - ri[0])
        tmp_no_missed_bp = 0
        for bp in bp_list:
            if ri[0] <= bp < ri[1]:
                tmp_no_missed_bp += 1
        no_missed_bp_list.append(tmp_no_missed_bp)

    T_non_recomb = [[0, int(bp_list[0] + 1)]]
    for i in range(len(bp_list) - 1):
        T_non_recomb.append([int(bp_list[i]) + 1, int(bp_list[i + 1]) + 1])
    T_non_recomb.append([int(bp_list[-1]+1), genome_length])

    breakpoint_set = set()

    for i in left_incl_right_excl_merged_intervals:
        breakpoint_set.add(i[0])
        breakpoint_set.add(i[1])

    for i in T_non_recomb:
        breakpoint_set.add(i[0])
        breakpoint_set.add(i[1])

    breakpoint_list = sorted(list(breakpoint_set))

    # generate_trees(exp_id, iter_i, T_non_recomb, "T", local)
    if rege_tree:
        generate_trees(exp_id, iter_i, recomb_intervals, "R", local)

    RF_distance_random = cal_RF_distance_random(exp_id, iter_i, breakpoint_list, recomb_intervals, T_non_recomb,
                                      left_incl_right_excl_merged_intervals, genome_length,
                                      sample_size, symbol_R="R")

    RF_distance_neighbor = cal_RF_distance_neighbor(exp_id, iter_i, breakpoint_list, recomb_intervals, T_non_recomb,
                                      left_incl_right_excl_merged_intervals, genome_length,
                                      sample_size, symbol_R="R")

    # print([RF_distance_random, RF_distance_neighbor], subalignment_length_list, no_missed_bp_list)
    assert all(list(map(lambda x, y: math.isclose(x, y), RF_distance_neighbor[3], RF_distance_random[3])))
    assert math.isclose(RF_distance_neighbor[1], RF_distance_random[1])
    assert len(subalignment_length_list) == len(no_missed_bp_list) == len(RF_distance_neighbor[3]) \
            == len(RF_distance_random[3])

    return [RF_distance_random, RF_distance_neighbor], subalignment_length_list, no_missed_bp_list
