# -*- coding: utf-8 -*-
"""
-------------------------------------------------
   File Name：     maxChi
   Description :
   Author :       shicen
   date：          11/2/22
-------------------------------------------------
   Change Activity:
                   11/2/22:
-------------------------------------------------
"""
import functools
import sys
from itertools import combinations
import tskit
from scipy.stats import chi2_contingency
from pathlib import Path
import pandas as pd
from matplotlib import pyplot as plt
import numpy as np
from Bio import SeqIO
import check_topo
import utility
import re
from analysis_treeSeq import analysis_treeSeq
import subprocess
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


def maxChi_multiproc_task(i, exp_id, rerun, local, max_bp, p_dir=""):
    recomb_rate = float(re.findall(r'recombRate_(.*?)_mutRate_(.*)', exp_id)[0][0])
    mut_rate = float(re.findall(r'recombRate_(.*?)_mutRate_(.*)', exp_id)[0][1])

    paths = sorted(list(Path(exp_id, f"iter_{i}").glob('*.tre')))

    tree_list = []

    # iterating over all files
    for path in paths:
        # print(path)  # printing file name

        with open(str(path), 'r') as f:
            tree = f.readline()
            tree_list.append(tree)

    treeSeq_info, _ = analysis_treeSeq(exp_id, i, rewrite=False, tree_list=tree_list)

    print(f"maxChi: {i}th iteration for {exp_id}")
    # n_bp = len(treeSeq_info['bp'])
    try:
        if type(treeSeq_info['bp']) != list and treeSeq_info['sample_size'] == 3:
            res_arr = maxchi_detection(treeSeq_info, exp_id, i, rerun, local)
        else:
            res_arr = maxchi_detection_2(treeSeq_info, exp_id, i, rerun, local,
                                         max_bp=max_bp)
        if res_arr[0] >= 0:
            print(res_arr)
            res_file_path = Path(p_dir, exp_id, f"iter_{i}", "maxChi_result.txt")
            with open(res_file_path, 'w') as f:
                f.write("\t".join(map(str, res_arr)) + f'\t{recomb_rate}\t{mut_rate}\t{i}\tmaxChi\n')
    except Exception as e:
        print(e)
        res_arr = np.array([-1, 0, 0, 0, 0, 0, 0, 0, 0, 0])
        res_file_path = Path(p_dir, exp_id, f"iter_{i}", "maxChi_result.txt")
        with open(res_file_path, 'w') as f:
            f.write("\t".join(map(str, res_arr)) + f'\t{recomb_rate}\t{mut_rate}\t{i}\tmaxChi\n')


def maxchi_detection(treeSeq_info, exp_id, i, rerun, local, left=None, right=None, p_dir=""):
    genome_length = treeSeq_info['genome_length']
    sample_size = treeSeq_info['sample_size']

    if left is None:
        left = 0
    if right is None:
        right = genome_length

    window = False
    if window:
        print("window open")
    else:
        print("window close")

    seqfile = str(Path(p_dir, exp_id) / ("iter_" + str(i)) / "concatenated_seq.fasta")

    n_comb = int(sample_size * (sample_size - 1) / 2)

    power_arr = np.zeros(n_comb)
    t_child_flag_arr = np.zeros(n_comb)
    bp_dist_arr = np.repeat(999.0, n_comb)
    bp_width_arr = np.zeros(n_comb)
    bp_precision_arr = np.zeros(n_comb)
    bp_precision_arr2 = np.zeros(n_comb)
    bp_CI_arr = np.empty(n_comb)
    chi_arr = np.zeros(n_comb)
    n_informative_sites_arr = np.zeros(n_comb)
    d1_arr = np.zeros(n_comb)
    d2_arr = np.zeros(n_comb)
    no_weighted_info_sites_arr = np.zeros(n_comb)
    can_pos_arr = np.empty(n_comb)
    Likelihood_arr = np.array(
        [np.empty(genome_length) for i in range(n_comb)]
    )
    info_site_code_list = np.array(
        [np.empty(genome_length) for i in range(n_comb)]
    )

    bp = treeSeq_info['bp']
    recomb_node = treeSeq_info['recomb_node']
    print("The recombination node is {0}".format(recomb_node))

    seq_dict = SeqIO.to_dict(SeqIO.parse(seqfile, 'fasta'))

    # method 1

    # seq_list = list(map(list, seq_list))
    pairIter = combinations(range(sample_size), 2)
    # P_threshold = 1

    n_res = 0
    seq_arange = np.arange(genome_length)

    maxChi_output_file = Path(p_dir, exp_id, "iter_" + str(i), exp_id + "_maxChi.txt")
    Likelihood_file = \
        Path(p_dir, exp_id, "iter_" + str(i), "maxChi_Likelihood.txt")

    plt.clf()

    if rerun:
        Likelihood_file.unlink(missing_ok=True)
        maxChi_output_file.unlink(missing_ok=True)

        for item in pairIter:
            p_list = np.ones(int(genome_length))
            chi_list = np.zeros(int(genome_length))
            seqpair_list = [list(seq_dict["n" + str(node)]) for node in item]

            seqCode = list(map(lambda x, y: x != y, seqpair_list[0], seqpair_list[1]))
            pair_list = sorted(list(item))
            pair_name = "_".join("n" + str(node) for node in pair_list)
            n_SNP = sum(seqCode)

            # for k in range(1 + left, right):
            for k in range(n_SNP, genome_length - n_SNP):
                if not window:
                    left_seq = seqCode[0:k]
                    right_seq = seqCode[k:genome_length]
                    contingency_table = np.array([[sum(left_seq), sum(right_seq)],
                                                  [k - sum(left_seq), (genome_length - k) - sum(right_seq)]])
                else:
                    left_seq = seqCode[k - n_SNP:k]
                    right_seq = seqCode[k:k + n_SNP]
                    contingency_table = np.array([[sum(left_seq), sum(right_seq)],
                                                      [n_SNP - sum(left_seq), n_SNP - sum(right_seq)]])
                if any(contingency_table.sum(axis=1) == 0) or any(contingency_table.sum(axis=0) == 0):
                    p_list[k] = 1
                    chi_list[k] = 0
                else:
                    tmp_chi, tmp_chi_pval, _, _ = chi2_contingency(contingency_table)
                    p_list[k] = tmp_chi_pval
                    chi_list[k] = tmp_chi

            likelihood_list = np.zeros(genome_length)
            likelihood_list[n_SNP:genome_length - n_SNP + 1] = (1 / p_list[n_SNP:genome_length - n_SNP + 1]) / \
                                                              np.sum(1 / p_list[n_SNP:genome_length - n_SNP + 1])

            with open(str(maxChi_output_file), 'a') as f:
                chi_str = '\t'.join(list(map(str, chi_list)))
                chi_str += f'\t{pair_name}\tchi\n'
                f.write(chi_str)
                p_str = '\t'.join(list(map(str, p_list)))
                p_str += f'\t{pair_name}\tp_value\n'
                f.write(p_str)
                likelihood_str = '\t'.join(list(map(str, likelihood_list)))
                likelihood_str += f'\t{pair_name}\tlikelihood\n'
                f.write(likelihood_str)

    if maxChi_output_file.is_file():
        pairIter = combinations(range(sample_size), 2)
        for item in pairIter:
            pair_list = sorted(list(item))
            pair_name = "_".join("n" + str(node) for node in pair_list)
            seqpair_list = [list(seq_dict["n" + str(node)]) for node in item]
            seqCode = list(map(lambda x, y: x != y, seqpair_list[0], seqpair_list[1]))
            n_SNP = sum(seqCode)

            maxChi_table = pd.read_csv(str(maxChi_output_file), header=None, sep='\t')

            chi_list = np.array(
                maxChi_table[(maxChi_table.iloc[:, genome_length] == pair_name) &
                             (maxChi_table.iloc[:, genome_length + 1] == "chi")].iloc[:, 0:genome_length]
            ).reshape(-1)

            p_list = np.array(
                maxChi_table[(maxChi_table.iloc[:, genome_length] == pair_name) &
                             (maxChi_table.iloc[:, genome_length + 1] == "p_value")].iloc[:,
                0:genome_length]
            ).reshape(-1)

            likelihood_list = np.array(
                maxChi_table[(maxChi_table.iloc[:, genome_length] == pair_name) &
                             (maxChi_table.iloc[:, genome_length + 1] == "likelihood")].iloc[:,
                0:genome_length]
            ).reshape(-1)

            if len(chi_list) == 0:
                continue

            can_pos = np.argmax(chi_list)
            can_p = p_list[can_pos]

            if not window:
                P_threshold = 1 / n_comb / genome_length  # Bonferroni correction
            else:
                P_threshold = 1 / n_comb / (2 * n_SNP)  # Bonferroni correction
            if can_p <= P_threshold:

                pos_seq = np.arange(genome_length)
                pos_expectation = np.inner(pos_seq, likelihood_list)
                pos_var = np.inner(np.power(pos_seq - bp, 2), likelihood_list)
                bp_precision = 1 / pos_var
                bp_precision_arr2[n_res] = bp_precision

                likelihood_entropy = - np.sum(np.nan_to_num(np.log2(likelihood_list)) * likelihood_list)

                can_pos_arr[n_res] = can_pos
                Likelihood_arr[n_res] = likelihood_list

                power_arr[n_res] = 1
                bp_dist_arr[n_res] = abs(can_pos - bp)
                if recomb_node & set(pair_list):
                    t_child_flag_arr[n_res] = 1
                bp_CI_l, bp_CI_r = utility.find_confidence_interval(likelihood_list, 0.95)
                bp_CI_arr[n_res] = 1 if bp_CI_l <= bp <= bp_CI_r else 0

                # add new precision
                bp_precision_arr[n_res] = 1 / (bp_CI_r - bp_CI_l)

                chi_arr[n_res] = chi_list[can_pos]

                d1_arr[n_res] = np.inner(abs(seq_arange - bp), seqCode)
                d2_arr[n_res] = np.inner((seq_arange - bp) ** 2, seqCode)
                n_informative_sites_arr[n_res] = sum(seqCode)
                informative_sites_code = np.array(list(map(
                    lambda x, y: 1 if (x and y < bp) else
                    -1 if (x and y >= bp) else 0, seqCode, seq_arange
                )))
                weighted_info_sites = list(map(
                    lambda c, p: 1 - abs(p - treeSeq_info['bp']) / genome_length if c else 0,
                    seqCode, seq_arange
                ))
                no_weighted_info_sites = sum(weighted_info_sites)
                no_weighted_info_sites_arr[n_res] = no_weighted_info_sites

                info_site_code_list[n_res] = informative_sites_code

                print("_".join("n" + str(node) for node in pair_list))
                print([can_pos, can_p, t_child_flag_arr[n_res]])
                n_res += 1

    d1_arr[0:n_res] = d1_arr[0:n_res] / n_informative_sites_arr[0:n_res]
    d2_arr[0:n_res] = np.sqrt(d2_arr[0:n_res]) / n_informative_sites_arr[0:n_res]

    max_accuracy_n = np.argmin(bp_dist_arr)

    Path(p_dir, exp_id, "iter_" + str(i), f'maxchi_likelihood_{exp_id}.png').unlink(missing_ok=True)

    if power_arr[max_accuracy_n] and (1 not in treeSeq_info['type']):
        plt.plot(np.arange(genome_length), Likelihood_arr[max_accuracy_n])
        plt.savefig(str(Path(p_dir, exp_id, "iter_" + str(i), f'maxchi_likelihood_{exp_id}.png')))

    with open(str(Likelihood_file), 'w') as f:
        L_str = '\t'.join(map(str, Likelihood_arr[max_accuracy_n])) + '\tMaxChi\n'
        f.write(L_str)

    with open(str(Path(p_dir, exp_id, f"iter_{i}", 'maxChi_informative_sites.txt')), 'w') as f2:
        L_str = '\t'.join(map(str, info_site_code_list[max_accuracy_n].astype(int)))
        L_str += '\tmaxChi\n'
        f2.write(L_str)

    res = np.array([power_arr[max_accuracy_n], t_child_flag_arr[max_accuracy_n],
                    bp_dist_arr[max_accuracy_n], bp_CI_arr[max_accuracy_n],
                    bp_precision_arr[max_accuracy_n], bp_precision_arr2[max_accuracy_n],
                    n_informative_sites_arr[max_accuracy_n], no_weighted_info_sites_arr[max_accuracy_n],
                    d1_arr[max_accuracy_n], can_pos_arr[max_accuracy_n]])

    # print(res)
    return res


def maxchi_detection_2(treeSeq_info, exp_id, iter_i, rerun, local, left=None, right=None, max_bp=3, p_dir=""):
    genome_length = treeSeq_info['genome_length']
    sample_size = treeSeq_info['sample_size']

    if left is None:
        left = 0
    if right is None:
        right = genome_length

    seqfile = str(Path(p_dir, exp_id) / ("iter_" + str(iter_i)) / "concatenated_seq.fasta")

    n_comb = int(sample_size * (sample_size - 1) / 2)

    bp = treeSeq_info['bp']
    recomb_node = treeSeq_info['recomb_node']
    recomb_type = treeSeq_info['type']
    print("The recombination node is {0}".format(recomb_node))

    seq_dict = SeqIO.to_dict(SeqIO.parse(seqfile, 'fasta'))

    # method 1

    pairIter = combinations(range(sample_size), 2)

    n_res = 0
    seq_arange = np.arange(genome_length)
    plt.clf()
    power = 0

    maxChi_output_file = \
        Path(p_dir, exp_id, "iter_" + str(iter_i), exp_id + "_iter_" + str(iter_i) + "_maxChi.txt")
    Likelihood_file = \
        Path(p_dir, exp_id, "iter_" + str(iter_i), "maxChi_Likelihood.txt")

    ts_file = str(Path(p_dir, exp_id, "iter_" + str(iter_i), 'ts.txt'))
    ts = tskit.load(ts_file)

    window_flag = True

    if rerun:
        if Likelihood_file.is_file():
            Likelihood_file.unlink()
        if maxChi_output_file.is_file():
            maxChi_output_file.unlink()

        for item in pairIter:

            seqpair_list = [list(seq_dict["n" + str(node)]) for node in item]

            seqCode = list(map(lambda x, y: x != y, seqpair_list[0], seqpair_list[1]))
            pair_list = sorted(list(item))
            pair_name = "_".join("n" + str(node) for node in pair_list)

            bp_list = [0, genome_length]
            res_dict = {}

            for i in range(max_bp):
                p_list = np.ones(int(genome_length))
                chi_list = np.zeros(int(genome_length))
                n_bp = len(bp_list) - 2

                if i > n_bp or max_bp <= n_bp:
                    break

                for j in range(n_bp + 1):
                    curr_bp_range = [bp_list[j], bp_list[j + 1]]
                    range_length = bp_list[j + 1] - bp_list[j]
                    range_SNP = sum(seqCode[bp_list[j]:bp_list[j + 1]])
                    curr_bp_range = [bp_list[j] + range_SNP,
                                     bp_list[j + 1] - range_SNP]

                    for p in range(curr_bp_range[0], curr_bp_range[1]):
                        if window_flag:
                            left_seq = seqCode[p-range_SNP:p]
                            right_seq = seqCode[p:p+range_SNP]
                            contingency_table = np.array([[sum(left_seq), sum(right_seq)],
                                                          [range_SNP - sum(left_seq),
                                                           range_SNP - sum(right_seq)]])
                        else:
                            left_seq = seqCode[curr_bp_range[0]:p]
                            right_seq = seqCode[p:curr_bp_range[1]]
                            contingency_table = np.array([[sum(left_seq), sum(right_seq)],
                                                          [(p - curr_bp_range[0]) - sum(left_seq),
                                                           (range_length - (p - curr_bp_range[0])) - sum(right_seq)]])
                        if any(contingency_table.sum(axis=1) == 0) or any(contingency_table.sum(axis=0) == 0):
                            p_list[p] = 1
                            chi_list[p] = 0
                        else:
                            tmp_chi, tmp_chi_pval, _, _ = chi2_contingency(contingency_table)
                            p_list[p] = tmp_chi_pval
                            chi_list[p] = tmp_chi

                can_pos = np.argmax(chi_list)
                can_p = p_list[can_pos]

                segment_length = 0
                for k in range(n_bp + 1):
                    if bp_list[k] <= can_p <= bp_list[k + 1]:
                        segment_length = bp_list[k + 1] - bp_list[k]
                        break

                P_threshold = 1 / n_comb / segment_length
                # P_threshold = 1 / n_comb

                if can_p <= P_threshold:
                    bp_list.append(can_pos)
                    # can_recomb = np.append(can_recomb, can_p)
                    bp_list = sorted(bp_list)
                    likelihood_arr = (1 / p_list) / np.sum(1 / p_list)
                    power = 1
                    plt.plot(np.arange(genome_length), likelihood_arr)

                    with open(maxChi_output_file, 'a') as f:
                        chi_str = '\t'.join(list(map(str, chi_list)))
                        chi_str += f'\t{pair_name}\tchi\t{i}\n'
                        f.write(chi_str)
                        p_str = '\t'.join(list(map(str, p_list)))
                        p_str += f'\t{pair_name}\tp_value\t{i}\n'
                        f.write(p_str)
                        likelihood_str = '\t'.join(list(map(str, likelihood_arr)))
                        likelihood_str += f'\t{pair_name}\tlikelihood\t{i}\n'
                        f.write(likelihood_str)

                    with open(str(Likelihood_file), 'a') as f:
                        f.write("\t".join(map(str, likelihood_arr)) + '\tmaxChi\n')

            res_dict[pair_name] = bp_list

        plt.savefig(str(Path(p_dir, exp_id, "iter_" + str(iter_i), 'Chi.png')))
        plt.close()

    bp_CI = 0
    n_bp = len(bp)
    bp_width_avg = 1000
    bp_width_arr = []
    bp_precision_arr = []
    bp_precision_arr2 = []
    bp_precision_avg = 0
    bp_dist = []
    bp_precision_list = []
    bp_precision_list2 = []
    bp_CI_list = [0 for x in range(n_bp)]
    can_recomb = []
    can_bp_recombtype_list = []
    subalignment_length_list = []
    maxChi_window_RF_distance1 = []
    no_missed_bp_list = []

    if maxChi_output_file.is_file():
        maxChi_table = pd.read_csv(str(maxChi_output_file), header=None, sep='\t')
        chi_table = maxChi_table[(maxChi_table.iloc[:, genome_length + 1] == "chi")]
        p_table = maxChi_table[(maxChi_table.iloc[:, genome_length + 1] == "p_value")]
        likelihood_table = maxChi_table[(maxChi_table.iloc[:, genome_length + 1] == "likelihood")]

        can_recomb = np.repeat(genome_length * 1.0, n_res)
        n_res = maxChi_table.shape[0] // 3
        bp_dist = np.repeat(genome_length * 1.0, n_res)

        power = 1

        n_cnt = 0
        for index, row in likelihood_table.iterrows():
            likelihood_list = np.array(row[0:1000])
            can_recomb = np.append(can_recomb, np.argmax(likelihood_list))
            can_recomb[n_cnt] = np.argmax(likelihood_list)
            bp_width_l, bp_width_r = utility.find_confidence_interval(likelihood_list,
                                                                         np.argmax(likelihood_list))
            bp_width_arr.append([bp_width_l, bp_width_r])
            plt.plot(np.arange(genome_length), likelihood_list)
            pos_seq = np.arange(genome_length)
            pos_expectation = np.inner(pos_seq, likelihood_list)
            pos_var = np.inner(np.power(pos_seq - pos_expectation, 2), likelihood_list)
            bp_precision = 1 / pos_var
            bp_precision_arr.append(bp_precision)
            # print(pos_expectation, pos_var, bp_precision)

            likelihood_entropy = - np.sum(np.nan_to_num(np.log2(likelihood_list.astype(float))) * likelihood_list)
            bp_precision_arr2.append(likelihood_entropy)

            seq_pair = list(map(lambda x: int(x.replace("n", "")), row[1000].split("_")))
            # print(seq_pair)
            relevant_breakpoints_list = []
            for sample_node in list(ts.samples()):
                # print(sample_node)
                tmp_sample_nodes = [x for x in seq_pair]
                if sample_node not in seq_pair:
                    tmp_sample_nodes.append(sample_node)

                    simplified_ts = ts.simplify(tmp_sample_nodes)
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
                    if len(relevant_breakpoints) > 0:
                        relevant_breakpoints_list.append(set(relevant_breakpoints))

            if len(relevant_breakpoints_list) > 0:
                relevant_breakpoints = functools.reduce(lambda x, y: x.union(y), relevant_breakpoints_list)

                bp_dist[n_cnt] = genome_length + 1
                hitten_bp = -1
                for rbp in list(relevant_breakpoints):
                    if abs(can_recomb[n_cnt] - rbp) < bp_dist[n_cnt]:
                        bp_dist[n_cnt] = abs(can_recomb[n_cnt] - rbp)
                        hitten_bp = rbp
                bp_round = list(map(lambda x: round(x, 3), bp))
                can_bp_recombtype_list.append(recomb_type[bp_round.index(round(hitten_bp, 3))])

                n_cnt += 1

        can_recomb = can_recomb.astype(int).tolist()
        bp_dist = bp_dist.tolist()
        bp_precision_list = bp_precision_arr[np.argmin(bp_dist)]

        maxChi_window_RF_distance1, subalignment_length_list, no_missed_bp_list = RF_distance(exp_id, iter_i, "maxChi", genome_length, bp, sample_size,
                                                            local=local, rege_tree=True, can_bp_list=can_recomb)

    plt.savefig(str(Path(p_dir, exp_id, "iter_" + str(iter_i), 'maxChi_siteprob.png')))
    plt.close()

    res = [power, can_recomb, bp_dist,
           bp_CI_list, bp_precision_list, can_bp_recombtype_list, 0, maxChi_window_RF_distance1, subalignment_length_list, no_missed_bp_list]
    print(res)
    return res


if __name__ == "__main__":
    exp_id = sys.argv[1]
    iter_i = int(float(sys.argv[2]))
    max_bp = int(float(sys.argv[3]))
    # p_dir = sys.argv[4]

    maxChi_multiproc_task(iter_i, exp_id, False, True, max_bp)

