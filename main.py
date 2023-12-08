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

import multiprocessing
from SeqPair_multithread import seqPair

if __name__ == '__main__':

    # "Sim params"
    sample_size = 10
    Ne = 1.0 / 2  # divide by two because msprime assumes individuals are diploid
    genome_length = 1e3
    rho = [0.1, 1]  # rate per genome
    recomb_rate = [x / genome_length for x in rho]  # rate per site
    # mu = [0.5, 1, 5, 10, 100, 500]
    # mut_rate = [x / genome_length for x in mu]
    pi = [0.01, 0.025, 0.05, 0.075, 0.1, 0.125, 0.15, 0.175, 0.2]
    mut_rate = [x / (2 * 2 * Ne) for x in pi]
    mu = [x * genome_length for x in mut_rate]
    iter_N = 5
    local = False
    regenerate = False
    rerun_list = ["threeSeq", "maxChi", "GARD"]
    # rerun_list = []

    seqpair = seqPair(n=sample_size, Ne=Ne,
                      L=genome_length,
                      r=recomb_rate, m=mut_rate, local=local)
    seqpair.simulate(iter_N, 1, 10, 0, 1000, regenerate, rerun_list)

    mult_thread = False

    if mult_thread:

        seqpair.initialization(rerun_list)

        threeSeq_thread = multiprocessing.Process(target=seqpair.threeSeq)
        maxChi_thread = multiprocessing.Process(target=seqpair.maxChi)
        GARD_thread = multiprocessing.Process(target=seqpair.gard)

        threeSeq_thread.start()
        maxChi_thread.start()
        GARD_thread.start()

        threeSeq_thread.join()
        maxChi_thread.join()
        GARD_thread.join()

    seqpair.analyze_result(rho, mu, ["threeSeq", "maxChi", "GARD"])

