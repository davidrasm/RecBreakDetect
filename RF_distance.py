# -*- coding: utf-8 -*-
"""
-------------------------------------------------
   File Name：     RF_distance
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
    except subprocess.CalledProcessError as e:
        print('Execution of "%s" failed!\n' % cmd)
        print(e.output)
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
                    # print(lower_b, upper_b)
                    # breakpoint_list.add(lower_b)
                    # breakpoint_list.add(upper_b + 1)
                    # threeSeq_mid_list.add(int((lower_b + upper_b) / 2 + 0.5))
                    intervals.append([lower_b, upper_b])
                    # seq_array[lower_b] = seq_array[upper_b] = 1
                    # if upper_b >= (lower_b + 2):
                    #     seq_array[lower_b + 1:upper_b] = -1

            l = f.readline()

    return intervals


def generate_trees(exp_id, iter_i, interval, method, symbol, window_size, local):
    Path(exp_id, f"iter_{iter_i}", method, window_size).mkdir(parents=True, exist_ok=True)
    for p in Path(exp_id, f"iter_{iter_i}").glob(f"{method}_{symbol}_alignment_*"):
        p.unlink(missing_ok=True)
    for p in Path(exp_id, f"iter_{iter_i}").glob(f"{method}_T_alignment_*"):
        p.unlink(missing_ok=True)

    for i in range(len(interval)):
        seq_file = Path(exp_id, f"iter_{iter_i}", "concatenated_seq.fasta")
        seq_rec_list = list(SeqIO.parse(seq_file, "fasta"))
        seq_list = []
        for seq_rec in seq_rec_list:
            seq = seq_rec.seq[interval[i][0]:interval[i][1] + 1]
            seq_list.append(SeqRecord(seq, id=seq_rec.id, description=""))

        # print(R_seq_list)
        SeqIO.write(seq_list,
                    Path(exp_id, f"iter_{iter_i}", method, window_size, f"{method}_{symbol}_alignment_{i}.fasta"),
                    'fasta')

        if local:
            if method == "True":
                raise Exception
                # cmd_raxml = '~/Downloads/raxml-ng_v1.2.0_macos_M1/raxml-ng --msa ' + str(
                #     Path(exp_id, f"iter_{iter_i}", f"{method}_{symbol}_alignment_{i}.fasta")) + ' --model GTR+G'
            else:
                cmd_raxml = '~/Downloads/raxml-ng_v1.2.0_macos_M1/raxml-ng --msa ' + str(
                    Path(exp_id, f"iter_{iter_i}", method, window_size,
                         f"{method}_{symbol}_alignment_{i}.fasta")) + ' --model GTR+G --redo'
        else:
            if method == "True":
                raise Exception
                # cmd_raxml = '~/raxml-ng/bin/raxml-ng --msa ' + str(
                #     Path(exp_id, f"iter_{iter_i}", f"{method}_{symbol}_alignment_{i}.fasta")) + ' --model GTR+G'
            else:
                cmd_raxml = '~/raxml-ng/bin/raxml-ng --msa ' + str(
                    Path(exp_id, f"iter_{iter_i}", method, window_size,
                         f"{method}_{symbol}_alignment_{i}.fasta")) + ' --model GTR+G --redo'
        try_run_cmd(cmd_raxml)

        # with open(Path(exp_id, f"iter_{iter_i}", f"threeseq_{symbol}_alignment_{i}.fasta.raxml.bestTree")) as f:
        #     r_tree = f.readline().strip()

        # print(r_tree)


def cal_RF_distance(exp_id, iter_i, method, merged_intervals, breakpoint_list, genome_length, sample_size,
                    R_non_recomb, T_non_recomb, window_size, R_symbol="R", T_symbol="T",
                    T_name="True"):
    weighted_RF_distance_adjacent = 0
    weighted_RF_distance_by_random = 0
    weighted_RF_distance_np = 0

    most_recent_inf_idx = -1
    most_recent_true_idx = -1

    subalignment_RF_dist_weighted = []
    subalignment_RF_dist = []
    tmp_RF_dist_weighted = 0
    tmp_RF_dist = 0
    flag_in_uninfo = 1
    for i in range(len(breakpoint_list) - 1):
        c_interval = [breakpoint_list[i], breakpoint_list[i + 1] - 1]
        inf_recomb_interval_idx = utility.search_mapped_interval(c_interval, R_non_recomb)
        true_recomb_interval_idx = utility.search_mapped_interval(c_interval, T_non_recomb)
        most_recent_inf_idx = inf_recomb_interval_idx if inf_recomb_interval_idx != -1 else most_recent_inf_idx
        # most_recent_true_idx = true_recomb_interval_idx if true_recomb_interval_idx != -1 else most_recent_true_idx

        print(c_interval, inf_recomb_interval_idx, true_recomb_interval_idx)

        with open(
                Path(exp_id, f"iter_{iter_i}", f"tree_{true_recomb_interval_idx}.tre")) as f:
            t_tree = f.readline().strip()

        if inf_recomb_interval_idx != -1:
            with open(Path(exp_id, f"iter_{iter_i}", method, window_size,
                           f"{method}_{R_symbol}_alignment_{inf_recomb_interval_idx}.fasta.raxml.bestTree")) as f:
                r_tree = f.readline().strip()

            rf_dist = utility.tree_dist(t_tree, r_tree)

            # curr_dist = rf_dist * (c_interval[1] - c_interval[0] + 1) / genome_length
            curr_dist = (rf_dist / (2 * (sample_size - 2))) * (c_interval[1] - c_interval[0] + 1) / genome_length
            curr_dist2 = (rf_dist / (2 * (sample_size - 2)))

            print(curr_dist)
            weighted_RF_distance_adjacent += curr_dist
            weighted_RF_distance_by_random += curr_dist
            weighted_RF_distance_np += curr_dist2

            flag_in_uninfo = 0
            tmp_RF_dist_weighted += curr_dist
            tmp_RF_dist += curr_dist2

        else:
            uninform_interval_idx = utility.search_mapped_interval(c_interval, merged_intervals)

            if uninform_interval_idx == -1:
                raise Exception("didn't find interval")

            print("add a random distance!")
            random_newick = utility.generate_a_random_tree(sample_size, length=c_interval[1] - c_interval[0] + 1)
            # tree_list = [t_tree, random_newick]
            rf_dist = utility.tree_dist(t_tree, random_newick)

            # curr_dist = rf_dist * (c_interval[1] - c_interval[0] + 1) / genome_length
            curr_dist = (rf_dist / (2 * (sample_size - 2))) * (c_interval[1] - c_interval[0] + 1) / genome_length
            print(curr_dist)
            weighted_RF_distance_by_random += curr_dist
            flag = 0

            print("use adjacent tree!")

            if most_recent_inf_idx == -1:
                with open(Path(exp_id, f"iter_{iter_i}", method, window_size,
                               f"{method}_{R_symbol}_alignment_0.fasta.raxml.bestTree")) as f:
                    r_tree = f.readline().strip()
                # tree_list = [t_tree, r_tree]

                rf_dist = utility.tree_dist(t_tree, r_tree)

                curr_dist = (rf_dist / (2 * (sample_size - 2))) * (c_interval[1] - c_interval[0] + 1) / genome_length
                print(curr_dist)
                weighted_RF_distance_adjacent += curr_dist
                continue

            if flag_in_uninfo == 0:
                subalignment_RF_dist_weighted.append(tmp_RF_dist_weighted)
                subalignment_RF_dist.append(tmp_RF_dist)
                tmp_RF_dist_weighted = 0
                tmp_RF_dist = 0
                flag_in_uninfo = 1

            if most_recent_inf_idx == len(R_non_recomb) - 1:
                with open(Path(exp_id, f"iter_{iter_i}", method, window_size,
                               f"{method}_{R_symbol}_alignment_{len(R_non_recomb) - 1}.fasta.raxml.bestTree")) as f:
                    r_tree = f.readline().strip()
                # tree_list = [t_tree, r_tree]

                rf_dist = utility.tree_dist(t_tree, r_tree)

                # curr_dist = rf_dist * (c_interval[1] - c_interval[0] + 1) / genome_length
                curr_dist = (rf_dist / (2 * (sample_size - 2))) * (c_interval[1] - c_interval[0] + 1) / genome_length
                print(curr_dist)
                weighted_RF_distance_adjacent += curr_dist
                continue

            curr_uninformative_interval = merged_intervals[uninform_interval_idx]
            curr_uninformative_interval_start = curr_uninformative_interval[0]
            curr_uninformative_interval_end = curr_uninformative_interval[1]
            curr_uninformative_interval_mid = (curr_uninformative_interval_start + curr_uninformative_interval_end) / 2
            no_break_in_window_flag = 0

            if c_interval[0] - 1 == R_non_recomb[most_recent_inf_idx][1] and \
                    c_interval[1] + 1 == R_non_recomb[most_recent_inf_idx + 1][0]:
                assert c_interval[0] == curr_uninformative_interval_start and c_interval[
                    1] == curr_uninformative_interval_end - 1
                with open(Path(exp_id, f"iter_{iter_i}", method, window_size,
                               f"{method}_{R_symbol}_alignment_{most_recent_inf_idx}.fasta.raxml.bestTree")) as f:
                    r_tree = f.readline().strip()
                rf_dist = utility.tree_dist(t_tree, r_tree)
                curr_dist = (rf_dist / (2 * (sample_size - 2))) * (
                        c_interval[1] - c_interval[0] + 1) / genome_length / 2
                weighted_RF_distance_adjacent += curr_dist
                with open(Path(exp_id, f"iter_{iter_i}", method, window_size,
                               f"{method}_{R_symbol}_alignment_{most_recent_inf_idx + 1}.fasta.raxml.bestTree")) as f:
                    r_tree = f.readline().strip()
                rf_dist = utility.tree_dist(t_tree, r_tree)
                # curr_dist = rf_dist * (c_interval[1] - c_interval[0] + 1) / genome_length / 2
                curr_dist = (rf_dist / (2 * (sample_size - 2))) * (
                        c_interval[1] - c_interval[0] + 1) / genome_length / 2
                print(curr_dist)
                weighted_RF_distance_adjacent += curr_dist
            else:
                if curr_uninformative_interval_start <= c_interval[0] <= c_interval[1] <= \
                        curr_uninformative_interval_mid:
                    with open(Path(exp_id, f"iter_{iter_i}", method, window_size,
                                   f"{method}_{R_symbol}_alignment_{max(most_recent_inf_idx,0)}.fasta.raxml.bestTree")) as f:
                        r_tree = f.readline().strip()
                    rf_dist = utility.tree_dist(t_tree, r_tree)
                    # curr_dist = rf_dist * (c_interval[1] - c_interval[0] + 1) / genome_length / 2
                    curr_dist = (rf_dist / (2 * (sample_size - 2))) * (
                            c_interval[1] - c_interval[0] + 1) / genome_length
                    print(curr_dist)
                    weighted_RF_distance_adjacent += curr_dist
                elif curr_uninformative_interval_start <= c_interval[0] <= \
                        curr_uninformative_interval_mid <= c_interval[1]:
                    with open(Path(exp_id, f"iter_{iter_i}", method, window_size,
                                   f"{method}_{R_symbol}_alignment_{max(most_recent_inf_idx,0)}.fasta.raxml.bestTree")) as f:
                        r_tree = f.readline().strip()
                    rf_dist = utility.tree_dist(t_tree, r_tree)
                    # curr_dist = rf_dist * (c_interval[1] - c_interval[0] + 1) / genome_length / 2
                    curr_dist = (rf_dist / (2 * (sample_size - 2))) * (
                            curr_uninformative_interval_mid - c_interval[0] + 1) / genome_length
                    print(curr_dist)
                    weighted_RF_distance_adjacent += curr_dist
                    with open(Path(exp_id, f"iter_{iter_i}", method, window_size,
                                   f"{method}_{R_symbol}_alignment_{min(most_recent_inf_idx,len(R_non_recomb)-1)}.fasta.raxml.bestTree")) as f:
                        r_tree = f.readline().strip()
                    rf_dist = utility.tree_dist(t_tree, r_tree)
                    # curr_dist = rf_dist * (c_interval[1] - c_interval[0] + 1) / genome_length / 2
                    curr_dist = (rf_dist / (2 * (sample_size - 2))) * (
                            c_interval[1] - curr_uninformative_interval_mid + 1) / genome_length
                    print(curr_dist)
                    weighted_RF_distance_adjacent += curr_dist
                elif curr_uninformative_interval_start <= curr_uninformative_interval_mid <= c_interval[0] <= \
                         c_interval[1] <= curr_uninformative_interval_end:
                    with open(Path(exp_id, f"iter_{iter_i}", method, window_size,
                                   f"{method}_{R_symbol}_alignment_{min(most_recent_inf_idx,len(R_non_recomb)-1)}.fasta.raxml.bestTree")) as f:
                        r_tree = f.readline().strip()
                    rf_dist = utility.tree_dist(t_tree, r_tree)
                    # curr_dist = rf_dist * (c_interval[1] - c_interval[0] + 1) / genome_length / 2
                    curr_dist = (rf_dist / (2 * (sample_size - 2))) * (
                            c_interval[1] - c_interval[0] + 1) / genome_length
                    print(curr_dist)
                    weighted_RF_distance_adjacent += curr_dist
                else:
                    raise Exception("Interesting case!")


    if flag_in_uninfo == 0:
        subalignment_RF_dist_weighted.append(tmp_RF_dist_weighted)
        subalignment_RF_dist.append(tmp_RF_dist)

    print(f"weighted_RF_distance_by_random: {weighted_RF_distance_by_random}, "
          f"weighted_RF_distance_adjacent: {weighted_RF_distance_adjacent}, "
          f"weighted_RF_distance_np: {weighted_RF_distance_np},\
            subalignment_RF_dist_weighted: {subalignment_RF_dist_weighted}, subalignment_RF_dist: {subalignment_RF_dist}")
    return [weighted_RF_distance_by_random, weighted_RF_distance_adjacent, weighted_RF_distance_np,
            subalignment_RF_dist_weighted, subalignment_RF_dist]


def RF_distance(exp_id, iter_i, method, genome_length, bp_list, sample_size, local, rege_tree,
                can_bp_list=[]):
    if method == "threeSeq":
        intervals = all_intervals(exp_id, iter_i)
        merged_intervals = utility.merge_intervals(intervals)

        can_bp_list = list(map(lambda x: int((x[0] + x[1]) / 2 + 0.5), merged_intervals))
    else:
        can_bp_list = sorted(can_bp_list)

    # no need if already having trees
    # T_non_recomb = [[0, int(bp_list[0] + 1)]]
    # for i in range(len(bp_list) - 1):
    #     T_non_recomb.append([int(bp_list[i]) + 1, int(bp_list[i + 1]) + 1])
    # T_non_recomb.append([int(bp_list[-1] + 1), genome_length])

    # generate_trees(exp_id, iter_i, T_non_recomb, "True", "T", local)

    T_non_recomb = [[0, int(bp_list[0])]]
    for i in range(len(bp_list) - 1):
        T_non_recomb.append([int(bp_list[i]) + 1, int(bp_list[i + 1])])
    T_non_recomb.append([int(bp_list[-1] + 1), genome_length - 1])

    window_size = [1, 5, 10, 25, 50]
    subalignment_length_list = []
    RF_dist_list_1 = [[-1, -1, -1, -1, -1] for x in window_size]
    no_missed_bp_list = []

    for cnt, ws in enumerate(window_size):
        print(ws)
        R_non_recomb = []
        intervals = []
        subalignment_length_window_list = []
        no_missed_bp_window_list = []

        if int(can_bp_list[0]) - 1 - ws >= 0:
            R_non_recomb.append([0, int(can_bp_list[0]) - 1 - ws])
        for i in range(len(can_bp_list) - 1):
            if int(can_bp_list[i + 1] - 1 - ws) >= int(can_bp_list[i]) + ws:
                R_non_recomb.append(
                    [int(can_bp_list[i]) + ws, int(can_bp_list[i + 1] - 1 - ws)])
            intervals.append(
                [max(can_bp_list[i] - ws, 0), min(can_bp_list[i] + ws, genome_length - 1)])
        if int(can_bp_list[-1] + ws) < 999:
            R_non_recomb.append([int(can_bp_list[-1] + ws), 999])
        intervals.append([max(can_bp_list[-1] - ws, 0), min(can_bp_list[-1] + ws, genome_length - 1)])

        if len(R_non_recomb) == 0:
            RF_dist_list_1[cnt] = [-1, -1, -1, -1, -1]
            continue

        R_non_recomb = utility.merge_intervals(R_non_recomb)
        merged_intervals = utility.merge_intervals(intervals)

        for rn in R_non_recomb:
            subalignment_length_window_list.append(rn[1] - rn[0] + 1)
            tmp_no_missed_bp = 0
            for bp in bp_list:
                if rn[0] < int(bp) < rn[1] - 1:
                    tmp_no_missed_bp += 1
            no_missed_bp_window_list.append(tmp_no_missed_bp)

        subalignment_length_list.append(subalignment_length_window_list)
        no_missed_bp_list.append(no_missed_bp_window_list)

        if rege_tree:
            generate_trees(exp_id, iter_i, R_non_recomb, method, "R", str(ws), local)

        breakpoint_set = set()
        for i in T_non_recomb:
            breakpoint_set.add(i[0])

        for i in R_non_recomb:
            breakpoint_set.add(i[0])
            breakpoint_set.add(i[1] + 1)

        breakpoint_list = sorted(list(breakpoint_set))

        method_RF_distance = cal_RF_distance(exp_id, iter_i, method, merged_intervals, breakpoint_list,
                                             genome_length, sample_size, R_non_recomb, T_non_recomb, str(ws))

        RF_dist_list_1[cnt] = method_RF_distance

        print(method_RF_distance, subalignment_length_window_list, no_missed_bp_window_list)
        assert len(subalignment_length_window_list) == len(no_missed_bp_window_list) == len(method_RF_distance[3]) \
               == len(method_RF_distance[4])

    RF_dist_list_2 = []

    '''
    seq_file = str(Path(exp_id, f"iter_{iter_i}", 'concatenated_seq.fasta'))
    record_dict = SeqIO.to_dict(SeqIO.parse(seq_file, 'fasta'))

    seqs = [x.seq for x in record_dict.values()]
    n_SNP = sum([len(set([x[i] for x in seqs])) > 1 for i in range(len(seqs[0]))])
    window_size_factor = [0.01, 0.02, 0.05, 0.1, 0.2, 0.5, 1]
    window_size = [int(x * n_SNP) for x in window_size_factor]

    for ws in window_size:
        print(ws)
        R_non_recomb = []
        intervals = []
        if int(can_bp_list[0]) - 1 - ws >= 0:
            R_non_recomb.append([0, int(can_bp_list[0]) - 1 - ws])
        for i in range(len(can_bp_list) - 1):
            if int(can_bp_list[i + 1] - 1 - ws) >= int(can_bp_list[i]) + ws:
                R_non_recomb.append(
                    [int(can_bp_list[i]) + ws, int(can_bp_list[i + 1] - 1 - ws)])
            intervals.append(
                [max(can_bp_list[i] - ws, 0), min(can_bp_list[i] + ws, genome_length - 1)])
        if 999 > min(int(can_bp_list[-1] + ws), 999):
            R_non_recomb.append([int(can_bp_list[-1] + ws), 999])
        intervals.append([max(can_bp_list[-1] - ws, 0), min(can_bp_list[-1] + ws, genome_length - 1)])

        if len(R_non_recomb) == 0:
            seq_tree_file = Path(exp_id, f"iter_{iter_i}",
                                 "concatenated_seq.fasta.raxml.bestTree")
            with open(seq_tree_file) as f:
                seq_tree = f.readline().strip()
            random_newick = utility.generate_a_random_tree(sample_size, length=genome_length)
            tree_list = [seq_tree, random_newick]
            tns = dendropy.TaxonNamespace()
            dendropy_tree_list = list(map(lambda x: dendropy.Tree.get(data=x, schema="newick",
                                                                      taxon_namespace=tns), tree_list))
            rf_dist = dendropy.calculate.treecompare.symmetric_difference(dendropy_tree_list[0], dendropy_tree_list[1])
            curr_dist = (rf_dist / (2 * (sample_size - 2)))
            RF_dist_list_2.append([curr_dist, -1])
            continue

        R_non_recomb = utility.merge_intervals(R_non_recomb)
        merged_intervals = utility.merge_intervals(intervals)

        generate_trees(exp_id, iter_i, R_non_recomb, method, "R", local)

        breakpoint_set = set()
        for i in T_non_recomb:
            breakpoint_set.add(i[0])

        for i in R_non_recomb:
            breakpoint_set.add(i[0])
            breakpoint_set.add(i[1] + 1)

        breakpoint_list = sorted(list(breakpoint_set))

        method_RF_distance = cal_RF_distance(exp_id, iter_i, method, merged_intervals, breakpoint_list,
                                             genome_length, sample_size, R_non_recomb, T_non_recomb)

        RF_dist_list_2.append(method_RF_distance)
        '''

    return RF_dist_list_1, subalignment_length_list, no_missed_bp_list
