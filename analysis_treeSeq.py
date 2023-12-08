# -*- coding: utf-8 -*-
"""
-------------------------------------------------
   File Name：     analysis_treeSeq
   Description :
   Author :       shicen
   date：          8/1/22
-------------------------------------------------
   Change Activity:
                   8/1/22:
-------------------------------------------------
"""
import re
import subprocess
import sys

import pandas as pd
import numpy as np
from pathlib import Path

from Bio import SeqIO
from scipy.stats import ranksums

from check_topo import print_topology
from functools import reduce

import dendropy
from dendropy.calculate import treecompare

import math


def try_run_cmd(cmd):
    try:
        output = subprocess.check_output(cmd, shell=True, stderr=subprocess.STDOUT)
        # sys.stdout.write(output.decode("UTF-8"))
    except subprocess.CalledProcessError:
        print('Execution of "%s" failed!\n' % cmd)
        sys.exit(1)


def generate_ts_table(ts, exp_id, iter_n, p_dir=""):
    node_table = pd.DataFrame.from_dict({
        'flags': ts.tables.nodes.flags,
        'time': ts.tables.nodes.time
    })

    edge_table = pd.DataFrame.from_dict({
        'left': ts.tables.edges.left,
        'right': ts.tables.edges.right,
        'parent': ts.tables.edges.parent,
        'child': ts.tables.edges.child
    })

    mut_table = pd.DataFrame.from_dict({
        'site': ts.tables.mutations.site,
        'node': ts.tables.mutations.node,
        'time': ts.tables.mutations.time
    })

    site_table = pd.DataFrame.from_dict({
        'position': ts.tables.sites.position
    })

    working_path = Path(p_dir, exp_id, f"iter_{iter_n}")

    node_table.to_csv(str(Path(working_path, "node_table.txt")), sep='\t')
    edge_table.to_csv(str(Path(working_path, "edge_table.txt")), sep='\t')
    mut_table.to_csv(str(Path(working_path, "mut_table.txt")), sep='\t')
    site_table.to_csv(str(Path(working_path, "site_table.txt")), sep='\t')

    return node_table, edge_table, mut_table, site_table


def calculate_branch_length(edge_table, node_table, lca_list=[]):
    L = 0.0
    L_2 = 0.0
    L_3 = 0.0
    L_4 = 0.0

    genome_length = np.max(edge_table.right)

    for index, row in edge_table.iterrows():
        child = row['child'].astype(int)
        parent = row['parent'].astype(int)

        child_t = node_table.iloc[child].time.item()
        parent_t = node_table.iloc[parent].time.item()

        L += parent_t - child_t

        segment_length = row['right'] - row['left']

        L_2 += (parent_t - child_t) * segment_length / genome_length

        if segment_length == genome_length:
            L_3 += (parent_t - child_t)
        else:
            L_3 += (parent_t - child_t) / 2

    return np.array([L, L_2, L_3, L_4])


def write_treeSeq(exp_id, iter_n, treeSeq_info,
                  rank_stat, info_sites_global, confounding_info_sites,
                  consistent_info_sites,
                  weighted_consistent_info_sites, rf_dist, no_neu_recurr_mut, no_adv_recurr_mut,
                  info_sites_imbalance,
                  p_dir=""):
    file_path = Path(p_dir, exp_id, "treeSeq_info.txt")

    recomb_rate = float(re.findall(r'recombRate_(.*?)_mutRate_(.*)', exp_id)[0][0])
    mut_rate = float(re.findall(r'recombRate_(.*?)_mutRate_(.*)', exp_id)[0][1])

    bp_pos = 0.0
    if type(treeSeq_info['bp']) != list:
        bp_pos = abs(treeSeq_info['bp'] - 500)
    else:
        seg_fraction = np.array([])
        bp_list = np.copy(treeSeq_info['bp'])
        bp_list = np.insert(bp_list, 0, 0)
        bp_list = np.append(bp_list, treeSeq_info['genome_length'])

        for i in range(len(bp_list) - 1):
            seg_fraction = np.append(seg_fraction, (bp_list[i + 1] - bp_list[i]) / treeSeq_info['genome_length'])
        for i in seg_fraction:
            bp_pos += -1 * i * math.log(i)

    with open(str(file_path), 'a') as f:
        f.write(
            f"{iter_n}\t{recomb_rate}\t{mut_rate}\t{treeSeq_info['type']}\t{treeSeq_info['recomb_node']}"
            f"\t{treeSeq_info['branch_length'][0]}\t{treeSeq_info['branch_length'][1]}"
            f"\t{treeSeq_info['branch_length'][2]}\t{treeSeq_info['branch_length'][3]}"
            f"\t{treeSeq_info['root_height_change']}\t{treeSeq_info['bp']}\t{bp_pos}"
            f"\t{rank_stat}\t{info_sites_global}\t{confounding_info_sites}"
            f"\t{consistent_info_sites}\t{weighted_consistent_info_sites}"
            f"\t{rf_dist}\t{no_neu_recurr_mut}\t{no_adv_recurr_mut}\t{info_sites_imbalance}\n")


def cal_RF_distance(exp_id, iter_n, genome_length, sample_size, rerun_tree=False, p_dir=""):
    seq_file = Path(exp_id, f"iter_{iter_n}", "concatenated_seq.fasta")
    cmd_raxml = '~/raxml-ng/bin/raxml-ng --msa ' + str(seq_file) + ' --model GTR+G --redo'

    if rerun_tree:
        try_run_cmd(cmd_raxml)
    seq_tree_file = Path(exp_id, f"iter_{iter_n}",
                         "concatenated_seq.fasta.raxml.bestTree")
    with open(seq_tree_file) as f:
        seq_tree = f.readline().strip()

    total_RF_dist = 0
    total_RF_dist2 = 0

    tree_paths = sorted(list(Path(exp_id, f"iter_{iter_n}").glob('tree_*.tre')))

    for path in tree_paths:
        with open(str(path), 'r') as f:
            local_tree = f.readline().strip()
        local_seq_file = str(path).replace(".tre", ".fasta")
        records = list(SeqIO.parse(local_seq_file, "fasta"))
        local_seq_length = len(records[0].seq)

        tree_list = [seq_tree, local_tree]
        tns = dendropy.TaxonNamespace()
        dendropy_tree_list = list(map(lambda x: dendropy.Tree.get(data=x, schema="newick",
                                                                  taxon_namespace=tns), tree_list))
        rf_dist = dendropy.calculate.treecompare.symmetric_difference(dendropy_tree_list[0],
                                                                      dendropy_tree_list[1])
        total_RF_dist += (rf_dist/(2*(sample_size-2))) * local_seq_length / genome_length
        total_RF_dist2 += (rf_dist/(2*(sample_size-2)))

    return total_RF_dist, total_RF_dist2


def analysis_treeSeq(exp_id, iter_n, rewrite=True, ts=None, bp=None, tree_list=None, p_dir=""):
    skip_flag = False
    if tree_list is None:
        assert ts is not None and bp is not None
        print(str(iter_n) + "th iteration")

        with open(str(Path(p_dir, exp_id, "iter_" + str(iter_n), 'ts_anc_SVG.svg')), 'w') as f:
            print(ts.draw_svg(), file=f)

        node_table, edge_table, mut_table, site_table = \
            generate_ts_table(ts, exp_id, iter_n, p_dir)

        genome_length = np.max(edge_table.right).astype(int)

        tree_list = [t.as_newick(include_branch_lengths=True) for t in ts.trees()]

        ts.dump(str(Path(p_dir, exp_id, f"iter_{iter_n}", 'ts.txt')))

        root_set = set()
        for t in ts.trees():
            root_set.add(t.root)

        if len(root_set) > 1:
            skip_flag = True

        if skip_flag:
            return {}, skip_flag

    else:

        node_table = pd.read_csv(str(Path(p_dir, exp_id, f"iter_{iter_n}", 'node_table.txt')), sep='\t')
        edge_table = pd.read_csv(str(Path(p_dir, exp_id, f"iter_{iter_n}", 'edge_table.txt')), sep='\t')
        genome_length = int(np.max(edge_table.right))

    sample_nodes = set(node_table[node_table['flags'] == 1].index.to_list())
    sample_size = len(sample_nodes)

    recomb_internal_node = node_table[(node_table['flags'] != 1) & (node_table['flags'] != 0)
                                      & (node_table['time'].duplicated(keep=False))].index.to_list()

    bp_list = sorted(edge_table.left.unique()[1:])
    n_bp = len(bp_list)

    treeSeq_info = {
        'sample_size': sample_size,
        'genome_length': genome_length,
        'bp': bp_list,
        'recomb_node': set(),
        'branch_length': np.zeros(sample_size),
        'type': [],
        'root_height_change': 0
    }  # Assume only two trees

    topology_list = [set(sorted(print_topology(t)[0]))
                     for t in tree_list]
    i = 1
    while i < len(topology_list):
        if topology_list[i] != topology_list[i - 1]:
            treeSeq_info['type'].append(3)
        else:
            internal_node = node_table[~node_table.index.isin(sample_nodes)]
            duplicated_time = internal_node[internal_node.time.duplicated()].time
            for t in duplicated_time:
                recomb_node_list = internal_node[internal_node.time == t].index.to_list()
                recomb_parents = edge_table[edge_table.child.isin(recomb_node_list)]
                recomb_grandparents = edge_table[edge_table.child.isin(recomb_parents)]
                if (len(set(recomb_parents.parent)) == 1 and np.max(recomb_parents.left) == np.min(
                        recomb_parents.right) == bp_list[i - 1]) \
                        or (len(set(recomb_grandparents.parent)) == 1 and np.max(recomb_grandparents.left) == np.min(
                    recomb_grandparents.right) == bp_list[i - 1]
                            and np.max(recomb_parents.left) == np.min(recomb_parents.right) == bp_list[i - 1]):
                    treeSeq_info['type'].append(1)
                    break
                elif np.sum((recomb_parents.left <= bp_list[i - 1]) & (bp_list[i - 1] <= recomb_parents.right)) >= 2:
                    treeSeq_info['type'].append(2)
                    break
                elif np.max(recomb_parents.left) > bp_list[i - 1] or np.min(recomb_parents.right) < bp_list[i - 1]:
                    continue
                else:
                    treeSeq_info['type'].append(0)
        i += 1

    treeSeq_info['branch_length'] = calculate_branch_length(edge_table, node_table)

    if 1 in treeSeq_info['type'] and sample_size > 3:
        skip_flag = True
        return {}, skip_flag

    no_adv_recurr_mut = 0

    if n_bp == 1 and sample_size == 3:
        treeSeq_info['bp'] = treeSeq_info['bp'][0]
        for index, row in edge_table.iterrows():
            if row.parent in recomb_internal_node:
                treeSeq_info['recomb_node'].update([row.child.astype(int)])

        if not (treeSeq_info['recomb_node'] & sample_nodes):
            recomb_node_p = list(treeSeq_info['recomb_node'])[0]
            for index, row in edge_table.iterrows():
                if row.parent == recomb_node_p:
                    treeSeq_info['recomb_node'].update([row.child.astype(int)])

    print(f"{exp_id}, {iter_n}")
    print(treeSeq_info)

    if not treeSeq_info['recomb_node'] and sample_size == 3:
        print("No recombination?? {0} {1}".format(exp_id, iter_n))

    if len(treeSeq_info['recomb_node']) > 1 and sample_size == 3:
        skip_flag = True
        return treeSeq_info, skip_flag

    if rewrite and not skip_flag:
        if sample_size == 3 and type(treeSeq_info['bp']) != list:
            c_label = "n" + str(list(treeSeq_info['recomb_node'])[0])
            parent_label = sample_nodes - treeSeq_info['recomb_node']
            p_label = "n" + str(parent_label.pop())
            q_label = "n" + str(parent_label.pop())
            seq_file = str(Path(p_dir, exp_id, "iter_" + str(iter_n), 'concatenated_seq.fasta'))
            record_dict = SeqIO.to_dict(SeqIO.parse(seq_file, 'fasta'))
            p = record_dict[p_label]
            q = record_dict[q_label]
            c = record_dict[c_label]
            seqpair_list = [list(p), list(q), list(c)]
            seq_range = np.arange(genome_length)
            seqCode_p = list(map(lambda p, q, c: c == p and c != q,
                                 seqpair_list[0], seqpair_list[1], seqpair_list[2]))
            seqCode_q = list(map(lambda p, q, c: c != p and c == q,
                                 seqpair_list[0], seqpair_list[1], seqpair_list[2]))
            seqCode_pq = np.array(
                list(map(lambda p, q, c: 1 if (c == p and c != q) else -1 if (c != p and c == q) else 0,
                         seqpair_list[0], seqpair_list[1], seqpair_list[2])))
            p_seq_pos = seq_range[seqCode_p]
            q_seq_pos = seq_range[seqCode_q]
            rank_stat = abs(ranksums(p_seq_pos, q_seq_pos)[0])
            info_sites_global = sum(seqCode_p) + sum(seqCode_q)
            confounding_info_sites = abs(
                sum(seqCode_pq[0:int(treeSeq_info['bp'])]) - sum(seqCode_pq[int(treeSeq_info['bp']):genome_length])
            )
            info_sites_imbalance_left = abs(sum(
                seqCode_pq[0:int(treeSeq_info['bp'])]
            ))
            info_sites_imbalance_right = abs(sum(
                seqCode_pq[int(treeSeq_info['bp']):genome_length]
            ))
            info_sites_imbalance = abs(info_sites_imbalance_left - info_sites_imbalance_right)
            consistent_info_sites_code1 = \
                list(map(lambda p, q, c, pos: 1 if ((c == p and c != q) and pos < treeSeq_info['bp']) or
                                                   ((c != p and c == q) and pos > treeSeq_info['bp']) else 0,
                         seqpair_list[0], seqpair_list[1], seqpair_list[2], seq_range))

            consistent_info_sites_code2 = \
                list(map(lambda p, q, c, pos: 1 if ((c == p and c != q) and pos > treeSeq_info['bp']) or
                                                   ((c != p and c == q) and pos < treeSeq_info['bp']) else 0,
                         seqpair_list[0], seqpair_list[1], seqpair_list[2], seq_range))

            info_sites_imbalance = abs(sum(consistent_info_sites_code1) - sum(consistent_info_sites_code2))

            if sum(consistent_info_sites_code1) > sum(consistent_info_sites_code2):
                consistent_info_sites = sum(consistent_info_sites_code1)
                weighted_consistent_info_sites = sum(
                    list(map(lambda c, pos: (1 - abs(pos - treeSeq_info['bp']) / genome_length) * c,
                             consistent_info_sites_code1, seq_range))
                )
            else:
                consistent_info_sites = sum(consistent_info_sites_code2)
                weighted_consistent_info_sites = sum(
                    list(map(lambda c, pos: (1 - abs(pos - treeSeq_info['bp']) / genome_length) * c,
                             consistent_info_sites_code2, seq_range))
                )

            with open(str(Path(p_dir, exp_id, "iter_" + str(iter_n), 'informative_sites.txt')), 'w') as f2:
                pos = '\t'.join(map(str, seq_range.astype(int)))
                pos += '\tpos\n'
                f2.write(pos)
                L_str = '\t'.join(map(str, seqCode_pq.astype(int)))
                L_str += '\tglobal\n'
                f2.write(L_str)

        else:
            rank_stat = 0
            info_sites_global = 0
            confounding_info_sites = 0
            consistent_info_sites = 0
            weighted_consistent_info_sites = 0
            info_sites_imbalance = 0

        try:
            rf_dist1, rf_dist2 = cal_RF_distance(exp_id, iter_n, genome_length, sample_size, rerun_tree=rewrite)
        except Exception as e:
            rf_dist1, rf_dist2 = -1, -1

        write_treeSeq(exp_id, iter_n,
                      treeSeq_info,
                      rank_stat, info_sites_global, confounding_info_sites, consistent_info_sites,
                      weighted_consistent_info_sites, rf_dist1, rf_dist2, no_adv_recurr_mut,
                      info_sites_imbalance,
                      p_dir)

    return treeSeq_info, skip_flag
