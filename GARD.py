# -*- coding: utf-8 -*-
"""
-------------------------------------------------
   File Name：     GARD
   Description :
   Author :       shicen
   date：          11/10/22
-------------------------------------------------
   Change Activity:
                   11/10/22:
-------------------------------------------------
"""
import os
import re
import sys
from pathlib import Path
import numpy as np
import subprocess
import json

from matplotlib import pyplot as plt

from analysis_treeSeq import analysis_treeSeq
from Bio import SeqIO
import utility
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


def GARD_multiproc_task(i, exp_id, rerun, local, max_bp, p_dir=""):
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

    print(f"GARD: {i}th iteration for {exp_id}")

    try:
        if type(treeSeq_info['bp']) != list and treeSeq_info['sample_size'] == 3:
            res_arr = gard_detection(treeSeq_info, exp_id, i, rerun=rerun, local=local)
        else:
            res_arr = gard_detection_2(treeSeq_info, exp_id, i, local=local,
                                       rerun=rerun, max_bp=max_bp)

        if res_arr[0] >= 0:
            res_file_path = Path(p_dir, exp_id, f"iter_{i}", "GARD_result.txt")
            with open(res_file_path, 'w') as f:
                f.write("\t".join(map(str, res_arr)) + f'\t{recomb_rate}\t{mut_rate}\t{i}\tGARD\n')
    except Exception as e:
        print(e)
        exc_type, exc_obj, exc_tb = sys.exc_info()
        fname = os.path.split(exc_tb.tb_frame.f_code.co_filename)[1]
        print(exc_type, fname, exc_tb.tb_lineno)
        res_arr = np.array([-1, 0, 0, 0, 0, 0, 0, 0, 0, 0])
        res_file_path = Path(p_dir, exp_id, f"iter_{i}", "GARD_result.txt")
        with open(res_file_path, 'w') as f:
            f.write("\t".join(map(str, res_arr)) + f'\t{recomb_rate}\t{mut_rate}\t{i}\tGARD\n')


def gard_siteProb(siteProb):

    positive_p = np.arange(1000)[siteProb > 0]
    positive_p_prob = siteProb[positive_p]

    step_p_interval = np.array([])
    for i in range(len(positive_p) - 1):
        step_p_interval = np.append(step_p_interval, (positive_p[i] + positive_p[i + 1]) / 2)

    site_p_step = np.zeros(1000)
    if step_p_interval[0] % 1 != 0:
        site_p_step[0:int(step_p_interval[0]) + 1] = positive_p_prob[0]
    else:
        site_p_step[0:int(step_p_interval[0])] = positive_p_prob[0]
        site_p_step[int(step_p_interval[0])] = (positive_p_prob[0] + positive_p_prob[1]) / 2

    for i in range(1, len(step_p_interval)):
        step_breakpoint = step_p_interval[i]
        previous_step_breakpoint = int(step_p_interval[i - 1])
        if step_breakpoint % 1 != 0:
            # print(i, step_breakpoint, positive_p_prob[i])
            site_p_step[(previous_step_breakpoint + 1):int(step_breakpoint) + 1] = positive_p_prob[i]
        else:
            site_p_step[(previous_step_breakpoint + 1):int(step_breakpoint)] = positive_p_prob[i]
            site_p_step[int(step_breakpoint)] = (positive_p_prob[i] + positive_p_prob[i + 1]) / 2

    if step_p_interval[-1] % 1 != 0:
        site_p_step[int(step_p_interval[-1]):] = positive_p_prob[-1]
    else:
        site_p_step[int(step_p_interval[-1]) + 1:] = positive_p_prob[-1]

    return site_p_step


def gard_detection(treeSeq_info, exp_id, i, rerun, local, p_dir=""):
    genome_length = treeSeq_info['genome_length']
    sample_size = treeSeq_info['sample_size']

    seq_dir = str(Path(p_dir, exp_id, f"iter_{i}"))

    print(f"iter_{i}")
    if local:
        cmd = f'hyphY gard --alignment {exp_id}/iter_{i}/concatenated_seq.fasta --max-breakpoints 1'
    else:
        cmd = f'hyphy gard --alignment {exp_id}/iter_{i}/concatenated_seq.fasta --max-breakpoints 1'

    if rerun:
        if Path(p_dir, seq_dir, "concatenated_seq.fasta.best-gard").is_file():
            Path(p_dir, seq_dir, "concatenated_seq.fasta.best-gard").unlink()

        if Path(p_dir, seq_dir, "concatenated_seq.fasta.best-gard.fit.bf").is_file():
            Path(p_dir, seq_dir, "concatenated_seq.fasta.best-gard.fit.bf").unlink()

        if Path(p_dir, seq_dir, "concatenated_seq.fasta.GARD.json").is_file():
            Path(p_dir, seq_dir, "concatenated_seq.fasta.GARD.json").unlink()

        run_cmd(cmd)

    json_output_file = str(Path(p_dir, exp_id, f"iter_{i}", "concatenated_seq.fasta.GARD.json"))

    with open(json_output_file, 'r') as f:
        json_parser = json.load(f)

    siteProb_raw = np.zeros(1000)
    for k, v in json_parser['siteBreakPointSupport'].items():
        siteProb_raw[int(k)] = v

    # new
    siteProb = gard_siteProb(siteProb_raw)
    siteProb = siteProb / np.sum(siteProb)

    # old
    siteProb = siteProb_raw / np.sum(siteProb_raw)

    bp = treeSeq_info['bp']

    gard_bp = json_parser['breakpointData']['0']['bps'][0][1]
    # gard_bp = np.argmax(siteProb)
    bp_dist = abs(gard_bp - 1 - bp)
    bp_width = 0
    if gard_bp != genome_length:
        power_flag = 1
    else:
        power_flag = 0

    Likelihood_file = \
        Path(p_dir, exp_id, "iter_" + str(i), "GARD_Likelihood.txt")

    with open(str(Likelihood_file), 'w') as f:
        f.write("\t".join(map(str, siteProb)) + '\tGARD\n')

    seq_file = str(Path(exp_id, f"iter_{i}", 'concatenated_seq.fasta'))
    record_dict = SeqIO.to_dict(SeqIO.parse(seq_file, 'fasta'))

    recomb_node = treeSeq_info['recomb_node']
    sample_nodes = set([i for i in range(sample_size)])
    non_recomb_node = sample_nodes - recomb_node

    c_label = "n" + str(list(recomb_node)[0])
    p_label = "n" + str(non_recomb_node.pop())
    q_label = "n" + str(non_recomb_node.pop())

    pseq_list = list(record_dict[p_label])
    qseq_list = list(record_dict[q_label])
    cseq_list = list(record_dict[c_label])

    SNP_mask = list(map(lambda x, y, z: not (x == y == z), pseq_list, qseq_list, cseq_list))

    n_SNP = sum(SNP_mask)

    seq_arange = np.arange(genome_length)

    weighted_SNP = list(map(lambda x, y, z, p: 1 - abs(p - treeSeq_info['bp']) / genome_length
    if not (x == y == z) else 0,
                            pseq_list, qseq_list, cseq_list, seq_arange))
    weighted_n_SNP = sum(weighted_SNP)

    informative_sites_code = \
        np.array(list(map(
            lambda p, q, c: 1 if (c == p and c != q) else
            -1 if (c != p and c == q) else
            2 if (c != p and c != q) else 0
            , pseq_list, qseq_list, cseq_list)))

    with open(str(Path(p_dir, exp_id, "iter_" + str(i), 'GARD_informative_sites.txt')), 'w') as f2:
        L_str = '\t'.join(map(str, informative_sites_code.astype(int)))
        L_str += '\tGARD\n'
        f2.write(L_str)

    bp_CI_l, bp_CI_r = utility.find_confidence_interval(siteProb, bp=np.argmax(siteProb_raw))
    bp_width = bp_CI_r - bp_CI_l

    if bp_CI_l <= bp <= bp_CI_r:
        bp_CI = 1
    else:
        bp_CI = 0

    t_child_flag = -1

    d_1 = np.inner(abs(seq_arange - bp), SNP_mask)
    d_2 = np.inner((seq_arange - bp) ** 2, SNP_mask)

    d_1 = d_1 / n_SNP
    d_2 = np.sqrt(d_2) / n_SNP

    pos_seq = np.arange(genome_length)
    pos_expectation = np.inner(pos_seq, siteProb)
    pos_var = np.inner(np.power(pos_seq - bp, 2), siteProb)
    bp_precision2 = 1 / pos_var

    bp_precision = 1 / bp_width

    if power_flag:
        plt.clf()
        plt.plot(np.arange(genome_length), siteProb)
        plt.savefig(str(Path(p_dir, exp_id, "iter_" + str(i), f'GARD_likelihood_{exp_id}.png')))

    res_arr = np.array([power_flag, t_child_flag, bp_dist,
                        bp_CI, bp_precision, bp_precision2, n_SNP, weighted_n_SNP, d_1, d_2])
    print(res_arr)

    return res_arr


def gard_detection_2(treeSeq_info, exp_id, iter_i, rerun, local, max_bp, p_dir=""):

    bp = treeSeq_info['bp']
    genome_length = treeSeq_info['genome_length']
    recomb_type = treeSeq_info['type']
    sample_size = treeSeq_info['sample_size']

    seq_dir = str(Path(p_dir, exp_id, "iter_" + str(iter_i)))

    print(f"iter_{iter_i}")
    if local:
        cmd = f'hyphY gard --alignment {exp_id}/iter_{iter_i}/concatenated_seq.fasta --max-breakpoints %d' % max_bp
    else:
        cmd = f'hyphy gard --alignment {exp_id}/iter_{iter_i}/concatenated_seq.fasta --max-breakpoints %d' % max_bp

    if rerun:
        if Path(p_dir, seq_dir, "concatenated_seq.fasta.best-gard").is_file():
            Path(p_dir, seq_dir, "concatenated_seq.fasta.best-gard").unlink()

        if Path(p_dir, seq_dir, "concatenated_seq.fasta.best-gard.fit.bf").is_file():
            Path(p_dir, seq_dir, "concatenated_seq.fasta.best-gard.fit.bf").unlink()

        if Path(p_dir, seq_dir, "concatenated_seq.fasta.GARD.json").is_file():
            Path(p_dir, seq_dir, "concatenated_seq.fasta.GARD.json").unlink()

        run_cmd(cmd)

    json_output_file = str(Path(p_dir, exp_id, "iter_" + str(iter_i), "concatenated_seq.fasta.GARD.json"))
    with open(json_output_file, 'r') as f:
        json_parser = json.load(f)

    siteProb_raw = np.zeros(genome_length)
    for k, v in json_parser['siteBreakPointSupport'].items():
        siteProb_raw[int(k)-1] = v
    siteProb = siteProb_raw / np.sum(siteProb_raw)

    pos_seq = np.arange(genome_length)
    pos_expectation = np.inner(pos_seq, siteProb)
    pos_var = np.inner(np.power(pos_seq - pos_expectation, 2), siteProb)
    bp_precision = 1 / pos_var
    print(pos_expectation, pos_var, bp_precision)

    likelihood_entropy = - np.sum(np.nan_to_num(np.log2(siteProb.astype(float))) * siteProb)
    bp_precision2 = likelihood_entropy

    bp_list_item = json_parser['improvements'][str(len(json_parser['improvements']) - 1)]
    can_bp_recombtype = []

    try:
        gard_bp_list = [i[0] for i in bp_list_item['breakpoints']]
        power_flag = 1
    except:
        gard_bp_list = []
        power_flag = 0

    if power_flag:

        bp_dist = np.repeat(genome_length * 1.0, len(gard_bp_list))
        for i in range(len(gard_bp_list)):
            hitten_bp = genome_length + 1
            for x in bp:
                if abs(x-gard_bp_list[i]) < bp_dist[i]:
                    bp_dist[i] = abs(x-gard_bp_list[i])
                    hitten_bp = x
            can_bp_recombtype.append(recomb_type[bp.index(hitten_bp)])

        GARD_window_RF_distance1, subalignment_length_list, no_missed_bp_list = RF_distance(exp_id, iter_i, "GARD", genome_length, bp, sample_size,
                                                            local=local, rege_tree=True, can_bp_list=gard_bp_list)
    else:
        bp_dist = np.repeat(genome_length, len(bp))
        GARD_window_RF_distance1, subalignment_length_list, no_missed_bp_list = [], [], []

    with open(str(Path(p_dir, exp_id, "iter_" + str(iter_i), "Likelihood.txt")), 'w') as f:
        f.write("\t".join(map(str, siteProb)) + '\tGARD\n')

    bp_width_arr = []
    for g_bp in gard_bp_list:
        bp_width_l, bp_width_r = utility.find_confidence_interval(siteProb, g_bp)
        bp_width_arr.append([bp_width_l, bp_width_r])

    n_bp = len(bp)
    bp_CI_list = [0 for x in range(n_bp)]

    cover_cnt = 0
    for k in range(n_bp):
        t_bp = bp[k]
        for bp_range in bp_width_arr:
            if bp_range[0] <= t_bp <= bp_range[1]:
                bp_CI_list[k] = 1
                break

    t_child_flag = -1

    res_arr = [power_flag, gard_bp_list, bp_dist.tolist(),
               bp_CI_list, bp_precision, can_bp_recombtype, 0, GARD_window_RF_distance1, subalignment_length_list, no_missed_bp_list]
    print(res_arr)
    return res_arr


if __name__ == "__main__":
    exp_id = sys.argv[1]
    iter_i = int(float(sys.argv[2]))
    max_bp = int(float(sys.argv[3]))
    # p_dir = sys.argv[4]

    GARD_multiproc_task(iter_i, exp_id, rerun=False, local=False, max_bp=max_bp)
