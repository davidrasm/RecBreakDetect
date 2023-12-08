# -*- coding: utf-8 -*-
"""
-------------------------------------------------
   File Name：     threeSeq
   Description :
   Author :       shicen
   date：          11/22/22
-------------------------------------------------
   Change Activity:
                   11/22/22:
-------------------------------------------------
"""
import os
import subprocess
import sys
from itertools import accumulate
from pathlib import Path
import re

import numpy as np
from Bio import SeqIO

from matplotlib import pyplot as plt
import tskit
from analysis_treeSeq import analysis_treeSeq
import time
from scipy.special import comb
import check_topo
import utility
from threeSeq_RF_distance import threeSeq_RF_distance
from RF_distance import RF_distance


def run_cmd(cmd, input=None):
    # try:
    #     output = subprocess.check_output(cmd, shell=True, stderr=subprocess.STDOUT)
    #     sys.stdout.write(output.decode("UTF-8"))
    # except subprocess.CalledProcessError:
    #     print('Execution of "%s" failed!\n' % cmd)
    #     sys.exit(1)
    output = subprocess.run(cmd, shell=True, text=True, input=input, stdout=subprocess.DEVNULL,
                            stderr=subprocess.STDOUT)


def threeSeq_multiproc_task(i, exp_id, rerun, local, p_dir=""):
    recomb_rate = float(re.findall(r'recombRate_(.*?)_mutRate_(.*)', exp_id)[0][0])
    mut_rate = float(re.findall(r'recombRate_(.*?)_mutRate_(.*)', exp_id)[0][1])

    tree_list = []

    # giving directory name to Path() function
    paths = sorted(list(Path(exp_id, f"iter_{i}").glob('*.tre')))

    # iterating over all files
    for path in paths:
        # print(path)  # printing file name

        with open(str(path), 'r') as f:
            tree = f.readline()
            tree_list.append(tree)

    treeSeq_info, _ = analysis_treeSeq(exp_id, i, rewrite=False, tree_list=tree_list)

    # res_file_path = Path(exp_id, "iter_" + str(i), "result_table_threeSeq.txt")

    print(f"threeSeq: {i}th iteration for {exp_id}")

    try:
        if type(treeSeq_info['bp']) != list and treeSeq_info['sample_size'] == 3:
            res_arr = threeSeq_detection(treeSeq_info, exp_id, i, rerun, local)
        else:
            res_arr = threeSeq_detection_2(treeSeq_info, exp_id, i, rerun, local)
        if res_arr[0] >= 0:
            print(res_arr)
            res_file_path = Path(p_dir, exp_id, f"iter_{i}", "threeSeq_result.txt")
            with open(res_file_path, 'w') as f:
                f.write("\t".join(map(str, res_arr)) + f'\t{recomb_rate}\t{mut_rate}\t{i}\tthreeSeq\n')
    except Exception as e:
        print(e)
        exc_type, exc_obj, exc_tb = sys.exc_info()
        fname = os.path.split(exc_tb.tb_frame.f_code.co_filename)[1]
        print(exc_type, fname, exc_tb.tb_lineno)
        res_arr = np.array([-1, 0, 0, 0, 0, 0, 0, 0, 0, 0])
        res_file_path = Path(p_dir, exp_id, f"iter_{i}", "threeSeq_result.txt")
        with open(res_file_path, 'w') as f:
            f.write("\t".join(map(str, res_arr)) + f'\t{recomb_rate}\t{mut_rate}\t{i}\tthreeSeq\n')


def threeSeq_detection(treeSeq_info, exp_id, iter_n, rerun, local, p_dir=""):
    genome_length = treeSeq_info['genome_length']
    sample_size = treeSeq_info['sample_size']

    n_combination = int(sample_size * (sample_size - 1) / 2)

    power_arr = np.zeros(n_combination)
    t_child_flag_arr = np.zeros(n_combination)
    bp_dist = np.repeat(999.0, n_combination)
    bp_precision = np.empty(n_combination)
    bp_CI = np.zeros(n_combination)
    n_informative_sites_arr = np.zeros(n_combination)
    d1_arr = np.zeros(n_combination)
    d2_arr = np.zeros(n_combination)
    p_arr = np.ones(n_combination)
    m_arr = np.empty(n_combination)
    n_arr = np.empty(n_combination)
    bp_result = ["" for i in range(n_combination)]
    seq_arange = np.arange(genome_length)
    seqCode_arr = np.array(
        [np.zeros(genome_length) for i in range(n_combination)]
    )
    no_weighted_info_sites_arr = np.zeros(n_combination)
    max_accuracy_n = 0

    # treeSeq_info = self.ts[exp_id][iter_n]
    bp = treeSeq_info['bp']
    recomb_node = treeSeq_info['recomb_node']
    print("The recombination node is {0}".format(recomb_node))

    seq_dir = str(Path(exp_id, "iter_" + str(iter_n)))

    if local:
        cmd = 'cd %s \n~/Downloads/3seq/3seq -f ./concatenated_seq.fasta -id %s -ptable ~/Downloads/3seq/myPvalueTable1000 ' % (
            seq_dir, exp_id)
    else:
        cmd = 'cd %s \n~/recomb/3seq/3seq -f ./concatenated_seq.fasta -id %s -ptable ~/recomb/3seq/myPvalueTable1000 ' % (
            seq_dir, exp_id)

    if rerun:
        file_list = [Path(p_dir, exp_id, "iter_" + str(iter_n), exp_id + ".3s." + ext)
                     for ext in ["rec", "log", "pvalHist", "longRec"]]
        for f in file_list:
            f.unlink(missing_ok=True)
        run_cmd(cmd, 'Y')
        time.sleep(5)

    rec_file = str(Path(p_dir, exp_id, "iter_" + str(iter_n), exp_id + '.3s.rec'))
    seq_file = str(Path(p_dir, exp_id, "iter_" + str(iter_n), 'concatenated_seq.fasta'))
    record_dict = SeqIO.to_dict(SeqIO.parse(seq_file, 'fasta'))
    # plt.clf()

    with open(rec_file, 'r') as f:

        l = f.readline()
        l = f.readline()
        n_res = 0
        # new logic
        while l:
            power_arr[n_res] = 1
            columns = l.strip().split('\t')
            bp_lst = columns[12:]
            m = int(columns[3])
            n = int(columns[4])
            k = int(columns[5])
            m_arr[n_res] = m
            n_arr[n_res] = n
            p_arr[n_res] = float(columns[8])
            p_label = str(columns[0])
            q_label = str(columns[1])
            c_label = str(columns[2])

            n_informative_sites_arr[n_res] = m + n
            if int(c_label.replace("n", "")) in recomb_node:
                t_child_flag_arr[n_res] = 1

            for bps in bp_lst:  # bps:  "1-9 & 631-631"
                bp_ranges = bps.split(" & ")  # bp_ranges = ['1-9', '631-631']

                for i in range(2):
                    bp_range = bp_ranges[i]
                    # for bp_range in bp_ranges:
                    lower_b = int(bp_range.split('-')[0])
                    upper_b = int(bp_range.split('-')[1])

                    mid = (lower_b + upper_b) / 2.0
                    if abs(mid - bp) < bp_dist[n_res]:
                        bp_dist[n_res] = abs(mid - bp)

                        bp_result[n_res] = [bps, i]
                        # bp_precision[n_res] = 1 / (upper_b - lower_b)
                        if lower_b <= bp <= upper_b:
                            bp_CI[n_res] = 1

            p = record_dict[p_label]
            q = record_dict[q_label]
            c = record_dict[c_label]

            seqpair_list = [list(p), list(q), list(c)]

            seqCode = list(map(lambda p, q, c: (c == p and c != q) or (c != p and c == q),
                               seqpair_list[0], seqpair_list[1], seqpair_list[2]))
            d1_arr[n_res] = np.inner(abs(seq_arange - bp), seqCode)
            d2_arr[n_res] = np.inner((seq_arange - bp) ** 2, seqCode)

            weighted_info_sites = list(
                map(lambda c, p: 1 - abs(p - treeSeq_info['bp']) / genome_length if c else 0,
                    seqCode, seq_arange))
            no_weighted_info_sites = sum(weighted_info_sites)
            no_weighted_info_sites_arr[n_res] = no_weighted_info_sites

            seqCode_pq = np.array(
                list(map(lambda p, q, c: 1 if (c == p and c != q) else -1 if (c != p and c == q) else 0,
                         seqpair_list[0], seqpair_list[1], seqpair_list[2])))
            seqCode_arr[n_res] = seqCode_pq

            n_res += 1
            l = f.readline()

    d1_arr[0:n_res] = d1_arr[0:n_res] / n_informative_sites_arr[0:n_res]
    d2_arr[0:n_res] = np.sqrt(d2_arr[0:n_res]) / n_informative_sites_arr[0:n_res]

    # original
    max_accuracy_n = np.argmin(bp_dist)

    site_probability = np.zeros(genome_length)

    if power_arr[max_accuracy_n] and (1 not in treeSeq_info['type']):

        # find site-specific probability

        m = m_arr[max_accuracy_n]
        n = n_arr[max_accuracy_n]

        final_bp_regions = bp_result[max_accuracy_n][0]
        final_bp_intervals = final_bp_regions.split(" & ")

        bp_position_idx = bp_result[max_accuracy_n][1]

        if bp_position_idx == 0:
            seqCode_res = seqCode_arr[max_accuracy_n]
        else:
            m, n = n, m
            seqCode_res = - seqCode_arr[max_accuracy_n]

        curr_h = 0
        for i in range(genome_length):
            curr_h += seqCode_res[i]
            if curr_h <= 0:
                site_probability[i] = 0
            else:
                P_k = (comb(m + n, n + curr_h + 1)) / (comb(m + n, n))

                site_probability[i] = np.max([1 - P_k, 0])

        site_probability = site_probability / np.sum(site_probability)

        plt.plot(np.arange(genome_length), site_probability)
        plt.savefig(str(Path(p_dir, exp_id, "iter_" + str(iter_n), f'threeSeq_likelihood_{exp_id}.png')))

        Likelihood_file = \
            Path(p_dir, exp_id, "iter_" + str(iter_n), "threeSeq_Likelihood.txt")

        with open(str(Likelihood_file), 'w') as f:
            f.write("\t".join(map(str, site_probability)) + '\tthreeSeq\n')

    # find variance
    pos_seq = np.arange(genome_length)
    pos_expectation = np.inner(pos_seq, site_probability)
    pos_var = np.inner(np.power(pos_seq - bp, 2), site_probability)
    bp_precision2 = 1 / pos_var

    # find CI
    if np.sum(site_probability) > 0:
        bp_CI_l, bp_CI_r = utility.find_confidence_interval(site_probability)
        bp_width = bp_CI_r - bp_CI_l
        bp_precision = 1 / bp_width
    else:
        bp_precision = 0

    if power_arr[max_accuracy_n] > 0:
        with open(str(Path(p_dir, exp_id, "iter_" + str(iter_n), 'informative_sites.txt')), 'a') as f2:
            L_str = '\t'.join(map(str, seqCode_arr[max_accuracy_n].astype(int)))
            L_str += '\tthreeSeq\n'
            f2.write(L_str)

        threeSeq_step = np.array(list(accumulate(seqCode_arr[max_accuracy_n], lambda x, y: x + y)))
        np.savetxt(str(Path(p_dir, exp_id, "iter_" + str(iter_n), 'threeSeq_step.txt')),
                   threeSeq_step.astype(int), fmt="%d")

    res = np.array([power_arr[max_accuracy_n], t_child_flag_arr[max_accuracy_n],
                    bp_dist[max_accuracy_n], bp_CI[max_accuracy_n], bp_precision, bp_precision2,
                    n_informative_sites_arr[max_accuracy_n], no_weighted_info_sites_arr[max_accuracy_n],
                    d1_arr[max_accuracy_n], d2_arr[max_accuracy_n]])

    return res


def threeSeq_detection_2(treeSeq_info, exp_id, iter_n, rerun, local, p_dir=""):
    genome_length = treeSeq_info['genome_length']
    sample_size = treeSeq_info['sample_size']

    bp_list = treeSeq_info['bp']
    recomb_node = treeSeq_info['recomb_node']
    print("The recombination node is {0}".format(recomb_node))

    seq_dir = str(Path(p_dir, exp_id, "iter_" + str(iter_n)))

    if local:
        cmd = 'cd %s \n~/Downloads/3seq/3seq -f ./concatenated_seq.fasta -id %s -ptable ~/Downloads/3seq/myPvalueTable1000 ' % (
            seq_dir, exp_id)
    else:
        cmd = 'cd %s \n~/recomb/3seq/3seq -f ./concatenated_seq.fasta -id %s -ptable ~/recomb/3seq/myPvalueTable1000 ' % (
            seq_dir, exp_id)

    if rerun:
        file_list = [Path(p_dir, exp_id, "iter_" + str(iter_n), exp_id + ".3s." + ext)
                     for ext in ["rec", "log", "pvalHist", "longRec"]]
        for f in file_list:
            f.unlink(missing_ok=True)
        run_cmd(cmd, 'Y')
        time.sleep(60)

    rec_file = str(Path(p_dir, exp_id, "iter_" + str(iter_n), exp_id + '.3s.rec'))
    seq_file = str(Path(p_dir, exp_id, "iter_" + str(iter_n), 'concatenated_seq.fasta'))

    if not Path(rec_file).is_file():
        time.sleep(10)

    power = 0
    n_bp = len(bp_list)

    with open(rec_file, 'r') as f:
        res_list = f.readlines()

    n_res = len(res_list) - 1
    print(n_res)
    bp_dist = np.repeat(genome_length * 1.0, n_res)
    can_bp_list = np.repeat(-genome_length - 1 * 1.0, n_res)
    can_bp_recombtype_list = []
    bp_width = np.empty(n_res)
    cover_list = np.zeros(n_res)
    threeSeq_conserv_RF_distance = []
    threeSeq_window_RF_distance1 = []
    subalignment_length_list = []
    no_missed_bp_list = []
    no_consistent_info_sites_list = []
    no_weighted_consistent_info_sites_list = []
    no_informative_sites_list = []

    ts_file = str(Path(p_dir, exp_id, "iter_" + str(iter_n), 'ts.txt'))
    ts = tskit.load(ts_file)

    plt.clf()

    n_cnt = 0

    with open(rec_file, 'r') as f:

        l = f.readline()
        l = f.readline()

        while l:
            power = 1
            columns = l.strip().split('\t')

            p_node = int(str(columns[0]).replace("n", ""))
            q_node = int(str(columns[1]).replace("n", ""))
            c_node = int(str(columns[2]).replace("n", ""))

            p_label = str(columns[0])
            q_label = str(columns[1])
            c_label = str(columns[2])

            record_dict = SeqIO.to_dict(SeqIO.parse(seq_file, 'fasta'))

            p = record_dict[p_label]
            q = record_dict[q_label]
            c = record_dict[c_label]

            seqpair_list = [list(p), list(q), list(c)]

            seqCode = list(map(lambda p, q, c: (c == p and c != q) or (c != p and c == q),
                               seqpair_list[0], seqpair_list[1], seqpair_list[2]))
            no_informative_sites_list.append(sum(seqCode))

            print(c_node, q_node, p_node)
            simplified_ts = ts.simplify([p_node, q_node, c_node])

            n_trees = len(simplified_ts.trees())
            simplified_ts_bps = list(simplified_ts.breakpoints())
            relevant_breakpoints = []
            relevant_breakpoints_recombtype = -1

            for i in range(n_trees - 1):
                flag = False
                tree1 = simplified_ts.at_index(i)
                tree2 = simplified_ts.at_index(i + 1)
                tree1_t = check_topo.print_topology(tree1.as_newick())[0]
                tree2_t = check_topo.print_topology(tree2.as_newick())[0]
                if tree1_t != tree2_t:
                    # flag = True
                    relevant_breakpoints.append(simplified_ts_bps[i + 1])
                    relevant_breakpoints_recombtype = 3
                    continue
                else:
                    tree1_root = tree1.root
                    tree2_root = tree2.root
                    if tree1_root == tree2_root and (tree1.as_newick() != tree2.as_newick()):
                        # flag = True
                        relevant_breakpoints.append(simplified_ts_bps[i + 1])
                        relevant_breakpoints_recombtype = 2
                        continue

            can_bp_recombtype_list.append(relevant_breakpoints_recombtype)

            bp_pred_lst = columns[12:]
            if len(relevant_breakpoints) > 0:
                for bps in bp_pred_lst:  # bps:  1-9 & 631-631
                    bp_ranges = bps.split(" & ")  # bp_ranges = ['1-9', '631-631']
                    for bp_range in bp_ranges:
                        lower_b = int(bp_range.split('-')[0])
                        upper_b = int(bp_range.split('-')[1])

                        mid = (lower_b + upper_b) / 2.0

                        for idx, rbp in enumerate(relevant_breakpoints):
                            # print(mid, rbp, abs(mid - rbp), abs(can_bp_list[n_cnt] - rbp))
                            if abs(mid - rbp) < bp_dist[n_cnt]:
                                can_bp_list[n_cnt] = mid
                                bp_dist[n_cnt] = abs(mid - rbp)
                                # print(mid, rbp, abs(mid - rbp), abs(can_bp_list[n_cnt] - rbp))
                                bp_width[n_cnt] = upper_b - lower_b
                                if lower_b <= rbp <= upper_b:
                                    cover_list[n_cnt] = 1

                n_cnt += 1
            l = f.readline()
            # print(bp_dist)

    threeSeq_conserv_RF_distance_random = []
    threeSeq_conserv_RF_distance_neighbor = []
    subalignment_length = 0
    no_missed_bp = -1
    if power:
        threeSeq_conserv_RF_distance, subalignment_length, no_missed_bp = threeSeq_RF_distance(exp_id, iter_n,
                                                                                               genome_length,
                                                                                               bp_list, sample_size,
                                                                                               local=local,
                                                                                               rege_tree=True)
        threeSeq_conserv_RF_distance_neighbor = threeSeq_conserv_RF_distance[1]
        threeSeq_conserv_RF_distance_random = threeSeq_conserv_RF_distance[0]
        threeSeq_window_RF_distance1, subalignment_length_list, no_missed_bp_list \
            = RF_distance(exp_id, iter_n, "threeSeq", genome_length, bp_list, sample_size, local, rege_tree=True)


    res = [power, can_bp_list.tolist(), bp_dist.tolist(),
           cover_list.tolist(), no_informative_sites_list, can_bp_recombtype_list, [threeSeq_conserv_RF_distance_random,
                                                                                    threeSeq_conserv_RF_distance_neighbor,
                                                                                    subalignment_length, no_missed_bp],
           threeSeq_window_RF_distance1, subalignment_length_list, no_missed_bp_list]

    return res


if __name__ == "__main__":
    exp_id = sys.argv[1]
    iter_i = int(float(sys.argv[2]))

    threeSeq_multiproc_task(iter_i, exp_id, False, False)
