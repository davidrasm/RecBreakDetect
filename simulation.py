# -*- coding: utf-8 -*-
"""
-------------------------------------------------
   File Name：     msprime_sim
   Description :
   Author :       shicen
   date：          7/3/22
-------------------------------------------------
   Change Activity:
                   7/3/22:
-------------------------------------------------
"""
import subprocess
import msprime
import numpy as np
import traceback
import pyvolve
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from pathlib import Path
import shutil
from analysis_treeSeq import analysis_treeSeq


def run_cmd(cmd):
    # try:
    #     output = subprocess.check_output(cmd, shell=True, stderr=subprocess.STDOUT)
    #     sys.stdout.write(output.decode("UTF-8"))
    # except subprocess.CalledProcessError:
    #     print('Execution of "%s" failed!\n' % cmd)
    #     sys.exit(1)
    output = subprocess.run(cmd, shell=True, text=True, input='Y', stdout=subprocess.DEVNULL, stderr=subprocess.STDOUT)


def sim_ARG(sample_size=10, Ne=100, length=1e3, recombination_rate=5e-6, mutation_rate=0, min_breakpoints=1,
            max_breakpoints=1000,
            min_breakpoint_coord=0, max_breakpoint_coord=1e3, plot=False):
    """
        Simulate ARGs with local trees in a tree sequence using msprime.

        Parameters:
           sample_size (int): Number of tips/samples in ARG
           Ne (float): Haploid effective population size
           length (int): Genome length
           recombination_rate (float): recombination rate per site per lineage
           min_breakpoints (int): minumum number of allowed breakpoints in simulated ARG
           max_breakpoints (int): maximum number of allowed breakpoints in simulated ARG
           plot (boolean): display ARG as local trees and ts tables

        Returns:
           ts (tskit.TreeSequence): TreeSequence representation of simulated ARG
    """

    breaks = 0
    breakpoint_coord = 0
    while breaks < min_breakpoints or breaks > max_breakpoints or \
            breakpoint_coord < min_breakpoint_coord or max_breakpoint_coord <= breakpoint_coord:
        ts = msprime.simulate(sample_size=sample_size, Ne=Ne, length=length, recombination_rate=recombination_rate,
                              mutation_rate=mutation_rate, record_full_arg=True)
        breakpoint_arr = ts.breakpoints(as_array=True)
        # print(breakpoint_arr)
        breaks = len(breakpoint_arr) - 2  # -2 because tskit counts ends as breakpoints
        breakpoint_coord = breakpoint_arr[1]

    for tr_num, tree in enumerate(ts.trees()):
        if plot:
            print("-" * 20)
            print("tree {}: interval = {}".format(tree.index, tree.interval))
            print(tree.draw(format="unicode"))

    '''
    # Count topo changes in tree sequences
    topo_changes = count_topo_changes(ts)
    '''

    if plot:
        print(ts.tables.nodes)
        print()
        print(ts.tables.edges)
        print()
        print("Recombination breakpoints: " + str(breaks))
        print()
        # print("Topological changes in tree sequence: " + str(topo_changes))

    return ts


def sim_seqs(tree_file, seq_file, seq_length=1000, mut_rate=1, freqs=[0.25, 0.25, 0.25, 0.25], kappa=2.75):
    """
        Simulate sequences for local trees in ARG under a HKY model in pyvolve

        Parameters:
           tree_file (str): Newick file with local tree for which sequences will be simulated
           seq_file (str): Output fasta file for simulated sequences

        Optional keyword arguements:
           mut_rate (float): Mutation rate per site
           seq_length (int): Genome length
           freqs (list[float]): Equilibrium nucleotide frequencies
           kappa (float): transition/transversion ratio for HKY model

    """

    # Read in phylogeny along which Pyvolve should simulate seqs
    # my_tree = pyvolve.read_tree(file=tree_file, scale_tree=mut_rate)  # scale_tree sets absolute mutation rate
    my_tree = pyvolve.read_tree(tree=tree_file, scale_tree=mut_rate)
    # print(tree_file)
    # pyvolve.print_tree(my_tree) # Print the parsed phylogeny

    # Parameterize HKY substitution model for sim
    nuc_model = pyvolve.Model("nucleotide", {"kappa": kappa, "state_freqs": freqs})

    # Define a Partition object which evolves set # of positions
    my_partition = pyvolve.Partition(models=nuc_model, size=seq_length)

    # Define an Evolver instance to evolve a single partition
    my_evolver = pyvolve.Evolver(partitions=my_partition, tree=my_tree)

    # Evolve sequences with custom file names
    # my_evolver(ratefile = "AMR_ratefile.txt", infofile = "AMR_infofile.txt", seqfile = "AMR-seqsim.fasta" )
    # **ratefile** is a custom name for the "site_rates.txt" file. Provide None or False to suppress file creation.
    # **infofile** is a custom name for the "site_rates_info.txt" file. Provide None or False to suppress file creation.
    my_evolver(seqfile=seq_file, ratefile=None, infofile=None, write_anc=True)


def concatenate_seq(seq_files, sample_size, output_seq):
    sim_seq_record_dict = {}
    node_name_list = ["n" + str(x) for x in range(sample_size)]
    for seq_file in seq_files:
        # print(seq_file)
        with open(seq_file) as f:
            for record in SeqIO.parse(f, "fasta"):
                node_name = record.id
                node_sim_seq = record.seq
                if node_name in node_name_list:
                    if node_name in sim_seq_record_dict.keys():
                        sim_seq_record_dict[node_name] += node_sim_seq
                    else:
                        sim_seq_record_dict[node_name] = node_sim_seq

    # ouput_seq_file = str(Path(exp_id) / 'concatenated_seq.fasta')
    seq_rec = []
    for id, seq in sim_seq_record_dict.items():
        # print(id, seq)
        c_seq_rec = SeqRecord(seq, id=id, description="")
        seq_rec.append(c_seq_rec)
    SeqIO.write(seq_rec, output_seq, 'fasta')


def msprime_sim(sample_size, Ne, genome_length, recombination_rate, mutation_rate,
                exp_id, iter_N, min_breakpoints, max_breakpoints,
                min_breakpoint_coord, max_breakpoint_coord, p_dir):
    # type_arr = np.zeros(iter_N)
    # ts_list = []
    # treeSeq_info_list = []

    # if Path(exp_id).is_dir():
    #     shutil.rmtree(exp_id)
    # Path(exp_id).mkdir()
    if (Path(p_dir, exp_id)).is_dir():
        # shutil.rmtree(Path(exp_id))
        cmd = 'rm -r {0}'.format(str(Path(p_dir, exp_id)))
        run_cmd(cmd)
    # print(str(Path(exp_id)))
    (Path(p_dir, exp_id)).mkdir()

    treeSeq_info_path = Path(p_dir, exp_id, "treeSeq_info.txt")
    # if treeSeq_info_path.is_file():
    #     run_cmd("rm " + str(treeSeq_info_path))
    with open(str(treeSeq_info_path), 'w') as f:
        f.write(
            "iter\trecomb_rate\tmut_rate\trecomb_type\trecomb_node"
            "\tARG_length_1\tARG_length_2"
            "\tARG_length_3\tARG_length_4"
            "\troot_height_change\tbp\tbp_pos_deviation\trank_stat\tinfo_sites_global"
            "\tconfounding_info_sites"
            "\tconsistent_info_sites\tweighted_consistent_info_sites"
            "\trf_dist\tno_neu_recurr_mut\tno_adv_recurr_mut\tinfo_sites_imbalance\n")

    for i in range(iter_N):
        flag = True
        while flag:
            try:
                if Path(p_dir, exp_id, "iter_" + str(i)).is_dir():
                    shutil.rmtree(Path(exp_id, "iter_" + str(i)))
                (Path(p_dir, exp_id) / ("iter_" + str(i))).mkdir()
                ts = sim_ARG(sample_size=sample_size, Ne=Ne, length=genome_length,
                             recombination_rate=recombination_rate,
                             # mutation_rate=mutation_rate,
                             min_breakpoints=min_breakpoints, max_breakpoints=max_breakpoints,
                             min_breakpoint_coord=min_breakpoint_coord, max_breakpoint_coord=max_breakpoint_coord)

                breaks = ts.breakpoints(as_array=True)
                break_point = breaks[1:-1]
                print(break_point)
                segments = len(breaks) - 1

                # type_arr[i] = treeSeq_info['type']
                print(p_dir, exp_id)

                tree_files = [str(Path(p_dir, exp_id) / ("iter_" + str(i)) / ("tree_" + str(x) + ".tre")) for x in
                              range(segments)]
                seq_files = [str(Path(p_dir, exp_id) / ("iter_" + str(i)) / ("tree_" + str(x) + ".fasta")) for x in
                             range(segments)]

                for tr_num, tree in enumerate(ts.trees()):
                    seq_length = round(tree.interval[1] - tree.interval[0])
                    with open(tree_files[tr_num], "w") as text_file:
                        print(tree.as_newick(), file=text_file)
                    # sim_seqs(tree_files[tr_num], seq_files[tr_num], mut_rate=mutation_rate, seq_length=seq_length)
                    tree_newick = tree.as_newick()
                    # sim_seqs(tree_newick, seq_files[tr_num], mut_rate=mutation_rate, seq_length=seq_length)
                    sim_seqs(tree_newick, seq_files[tr_num], seq_length=seq_length, mut_rate=mutation_rate)

                flag = False

                output_seq = str(Path(p_dir, exp_id) / ("iter_" + str(i)) / "concatenated_seq.fasta")
                concatenate_seq(seq_files, sample_size, output_seq)
                _, skip_flag = analysis_treeSeq(exp_id, i, ts=ts, bp=break_point, p_dir=p_dir)

                if skip_flag:
                    flag = True

            except Exception as e:
                print(traceback.format_exc())
                # exit(1)
                print("Failed to simulate, try again")

    # with open(str(Path(exp_id, "type_arr.txt")), 'w') as f:
    #     f.write('\t'.join(type_arr.astype(str)) + '\n')

    # return treeSeq_info_list
