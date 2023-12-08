import itertools
import json
import subprocess
from functools import reduce
from itertools import combinations
from simulation import msprime_sim
from scipy.stats import chi2_contingency
from pathlib import Path
import pandas as pd
import multiprocessing
import clean_data
from itertools import accumulate
from threeSeq import *
from maxChi import maxchi_detection_2, maxchi_detection
from GARD import gard_detection_2, gard_detection
from threeSeq import threeSeq_detection_2, threeSeq_detection


NUM_PROC = 8


def run_cmd(cmd, input=None):
    # try:
    #     output = subprocess.check_output(cmd, shell=True, stderr=subprocess.STDOUT)
    #     sys.stdout.write(output.decode("UTF-8"))
    # except subprocess.CalledProcessError:
    #     print('Execution of "%s" failed!\n' % cmd)
    #     sys.exit(1)
    output = subprocess.run(cmd, shell=True, text=True, input=input, stdout=subprocess.DEVNULL,
                            stderr=subprocess.STDOUT)


def analysis_treeSeq_multiproc_task(iter_list, exp_id, rewrite):

    for i in iter_list:
        # giving directory name to Path() function
        paths = sorted(list(Path(exp_id, f"iter_{i}").glob('*.tre')))

        tree_list = []

        # iterating over all files
        for path in paths:

            with open(str(path), 'r') as f:
                tree = f.readline()
                tree_list.append(tree)

        analysis_treeSeq(exp_id, i, rewrite, tree_list=tree_list)


class seqPair:

    def __init__(self, n, Ne, L, r, m, iter_N=1, local=True, p_dir=""):
        # self.n_informative_sites_3seq = None
        self.sample_size = n
        self.Ne = Ne
        self.genome_length = int(L)
        self.mutation_rate = m
        self.recombination_rate = r
        self.iter_N = iter_N
        self.res = {}
        self.ts = {}
        self.local = local
        self.regenerate = True
        self.rerun_list = []
        self.p_dir = p_dir
        self.min_bp = 0
        self.max_bp = 0

    @staticmethod
    def calculate_metrics(N, power_value,
                          bp_dist_value, bp_range_value):

        power = power_value / N
        bp_dist_avg = bp_dist_value / power_value if power_value != 0 else -1
        bp_range_avg = bp_range_value / power_value if power_value != 0 else -1

        return [power, bp_dist_avg, bp_range_avg]

    def simulate(self, iter_N=1, min_breakpoints=1, max_breakpoints=1,
                 min_breakpoint_coord=None, max_breakpoint_coord=None, regenerate=True, rerun_list=[]):

        if min_breakpoint_coord is None:
            min_breakpoint_coord = 0
        if max_breakpoint_coord is None:
            max_breakpoint_coord = self.genome_length

        self.min_bp = min_breakpoints
        self.max_bp = max_breakpoints

        self.regenerate = regenerate

        self.rerun_list = rerun_list

        treeSeq_dict = {}
        self.iter_N = iter_N
        for item in itertools.product(self.recombination_rate,
                                      self.mutation_rate):
            r = item[0]
            m = item[1]
            exp_id = "recombRate_" + str(r * self.genome_length) + "_mutRate_" + str(m * self.genome_length)
            treeSeq_dict[exp_id] = []

            if regenerate:

                msprime_sim(sample_size=self.sample_size,
                            Ne=self.Ne,
                            genome_length=self.genome_length,
                            recombination_rate=r,
                            mutation_rate=m,
                            iter_N=self.iter_N,
                            exp_id=exp_id,
                            min_breakpoints=min_breakpoints,
                            max_breakpoints=max_breakpoints,
                            min_breakpoint_coord=min_breakpoint_coord,
                            max_breakpoint_coord=max_breakpoint_coord,
                            p_dir=self.p_dir)
            else:
                # treeSeq_info_list = []
                treeSeq_info_path = Path(self.p_dir, exp_id, "treeSeq_info.txt")
                rewrite = False
                if rewrite:
                    if treeSeq_info_path.is_file():
                        run_cmd("rm " + str(treeSeq_info_path))
                    with open(str(treeSeq_info_path), 'w') as f:
                        f.write(
                            "iter\trecomb_rate\tmut_rate\trecomb_type\trecomb_node"
                            "\tARG_length_1\tARG_length_2"
                            "\tARG_length_3\tARG_length_4"
                            "\troot_height_change\tbp\tbp_pos_deviation"
                            "\trank_stat\tinfo_sites_global"
                            "\tconfounding_info_sites"
                            "\tconsistent_info_sites\tweighted_consistent_info_sites"
                            "\trf_dist\tno_neu_recurr_mut\tno_adv_recurr_mut\tinfo_sites_imbalance\n")

                    iter_per_proc = self.iter_N // NUM_PROC

                    jobs = []

                    for i in range(NUM_PROC + 1):
                        if (i + 1) * iter_per_proc < self.iter_N:
                            iter_list = range(i * iter_per_proc, (i + 1) * iter_per_proc)
                        else:
                            iter_list = range(i * iter_per_proc, self.iter_N)

                        process = multiprocessing.Process(
                            target=analysis_treeSeq_multiproc_task,
                            args=(iter_list, exp_id, rewrite)
                        )
                        jobs.append(process)

                    for j in jobs:
                        j.start()

                    for j in jobs:
                        j.join()

    def initialization(self, method_list=[]):

        for item in itertools.product(self.recombination_rate,
                                      self.mutation_rate):
            exp_id = "recombRate_" + str(item[0] * self.genome_length) + "_mutRate_" + str(item[1] * self.genome_length)

            for i in range(self.iter_N):

                seq_dir = Path(self.p_dir, exp_id, "iter_" + str(i))

                if "threeSeq" in self.rerun_list:
                    if Path(self.p_dir, seq_dir, exp_id + ".3s.log").is_file():
                        run_cmd("rm " + str(Path(self.p_dir, seq_dir, exp_id + ".3s.log")))

                    if Path(self.p_dir, seq_dir, exp_id + ".3s.pvalHist").is_file():
                        run_cmd("rm " + str(Path(self.p_dir, seq_dir, exp_id + ".3s.pvalHist")))

                    if Path(self.p_dir, seq_dir, exp_id + ".3s.rec").is_file():
                        run_cmd("rm " + str(Path(self.p_dir, seq_dir, exp_id + ".3s.rec")))

                    if Path(self.p_dir, seq_dir, exp_id + ".3s.longRec").is_file():
                        run_cmd("rm " + str(Path(self.p_dir, seq_dir, exp_id + ".3s.longRec")))

                else:
                    res_file_path = Path(self.p_dir, exp_id, "result_table_threeSeq.txt")

                    with open(str(res_file_path), 'w') as f:
                        f.write(
                            "power\ttrue_child\taccuracy\tprecision\tcoverage\trecomb_rate\tmut_rate\titer\tmethod\n")

                if "GARD" in self.rerun_list:
                    if Path(self.p_dir, seq_dir, "concatenated_seq.fasta.best-gard").is_file():
                        run_cmd("rm " + str(Path(self.p_dir, seq_dir, "concatenated_seq.fasta.best-gard")))

                    if Path(self.p_dir, seq_dir, "concatenated_seq.fasta.best-gard.fit.bf").is_file():
                        run_cmd("rm " + str(Path(self.p_dir, seq_dir, "concatenated_seq.fasta.best-gard.fit.bf")))

                    if Path(self.p_dir, seq_dir, "concatenated_seq.fasta.GARD.json").is_file():
                        run_cmd("rm " + str(Path(self.p_dir, seq_dir, "concatenated_seq.fasta.GARD.json")))
                else:
                    res_file_path = Path(self.p_dir, exp_id, "result_table_GARD.txt")
                    with open(str(res_file_path), 'w') as f:
                        f.write(
                            "power\ttrue_child\taccuracy\tprecision\tcoverage\trecomb_rate\tmut_rate\titer\tmethod\n")

                if "maxChi" in self.rerun_list:
                    maxChi_file = Path(self.p_dir, exp_id, "iter_" + str(i), exp_id + "_iter_" + str(i) + "_maxChi.txt")
                    if maxChi_file.is_file():
                        maxChi_file.unlink()
                    maxChi_likelihood = Path(self.p_dir, exp_id, "iter_" + str(i), "_maxChi_Likelihood.txt")
                    if maxChi_likelihood.is_file():
                        maxChi_likelihood.unlink()
                else:
                    res_file_path = Path(self.p_dir, exp_id, "result_table_maxChi.txt")
                    with open(str(res_file_path), 'w') as f:
                        f.write(
                            "power\ttrue_child\taccuracy\tprecision\tcoverage\trecomb_rate\tmut_rate\titer\tmethod\n")

    def maxChi_multiproc_task(self, iter_list, exp_id):

        recomb_rate = float(re.findall(r'recombRate_(.*?)_mutRate_(.*)', exp_id)[0][0])
        mut_rate = float(re.findall(r'recombRate_(.*?)_mutRate_(.*)', exp_id)[0][1])

        for i in iter_list:

            with open(str(Path(exp_id, f"iter_{i}", 'tree_0.tre')), 'r') as f2:
                tree_0 = f2.readline()
            with open(str(Path(exp_id, f"iter_{i}", 'tree_1.tre')), 'r') as f2:
                tree_1 = f2.readline()

            treeSeq_info, _ = analysis_treeSeq(exp_id, i, rewrite=False, tree_list=[tree_0, tree_1])

            print(f"maxChi: {i}th iteration for {exp_id}")
            n_bp = len(treeSeq_info['bp'])
            if n_bp == 1 and treeSeq_info['sample_size'] == 3:
                res_arr = maxchi_detection(treeSeq_info, exp_id, i)
            else:
                if "maxChi" in self.rerun_list:
                    res_arr = maxchi_detection_2(treeSeq_info, exp_id, i, True,
                                                 max_bp=self.max_bp)
                else:
                    res_arr = maxchi_detection_2(treeSeq_info, exp_id, i, False,
                                                 max_bp=self.max_bp)
            if res_arr[0] >= 0:
                print(res_arr)
                res_file_path = Path(self.p_dir, exp_id, f"iter_{i}", "maxChi_result.txt")
                with open(res_file_path, 'w') as f:
                    f.write("\t".join(map(str, res_arr)) + f'\t{recomb_rate}\t{mut_rate}\t{i}\tmaxChi\n')

    def maxChi(self, left=None, right=None):

        # self.res = {}
        for item in itertools.product(self.recombination_rate,
                                      self.mutation_rate):
            exp_id = "recombRate_" + str(item[0] * self.genome_length) + "_mutRate_" + str(item[1] * self.genome_length)

            iter_per_proc = self.iter_N // NUM_PROC

            jobs = []

            for i in range(NUM_PROC + 1):
                if (i + 1) * iter_per_proc < self.iter_N:
                    iter_list = range(i * iter_per_proc, (i + 1) * iter_per_proc)
                else:
                    iter_list = range(i * iter_per_proc, self.iter_N)

                process = multiprocessing.Process(
                    target=self.maxChi_multiproc_task,
                    args=(iter_list, exp_id)
                )
                jobs.append(process)

            for j in jobs:
                j.start()

            for j in jobs:
                j.join()


    def threeSeq_multiproc_task(self, iter_list, exp_id):

        recomb_rate = float(re.findall(r'recombRate_(.*?)_mutRate_(.*)', exp_id)[0][0])
        mut_rate = float(re.findall(r'recombRate_(.*?)_mutRate_(.*)', exp_id)[0][1])

        for i in iter_list:

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
            if type(treeSeq_info['bp']) != list and treeSeq_info['sample_size'] == 3:
                # res_arr = self.threeSeq_detection(treeSeq_info, exp_id, i)
                res_arr = threeSeq_detection(treeSeq_info, exp_id, i, "threeSeq" in self.rerun_list, self.local)
            else:
                res_arr = threeSeq_detection_2(treeSeq_info, exp_id, i, "threeSeq" in self.rerun_list, self.local)
            if res_arr[0] >= 0:
                print(res_arr)
                res_file_path = Path(self.p_dir, exp_id, f"iter_{i}", "threeSeq_result.txt")
                with open(res_file_path, 'w') as f:
                    f.write("\t".join(map(str, res_arr)) + f'\t{recomb_rate}\t{mut_rate}\t{i}\tthreeSeq\n')

    def threeSeq(self):
        self.res = {}

        for item in itertools.product(self.recombination_rate,
                                      self.mutation_rate):
            exp_id = "recombRate_" + str(item[0] * self.genome_length) + "_mutRate_" + str(item[1] * self.genome_length)

            print(exp_id)

            # if not res_file_path.is_file():
            #     with open(str(res_file_path), 'w') as f:
            #         f.write("power\ttrue_child\taccuracy\tprecision\tcoverage\trecomb_rate\tmut_rate\titer\tmethod\n")

            iter_per_proc = self.iter_N // NUM_PROC

            jobs = []

            for i in range(NUM_PROC + 1):
                if (i + 1) * iter_per_proc < self.iter_N:
                    iter_list = range(i * iter_per_proc, (i + 1) * iter_per_proc)
                else:
                    iter_list = range(i * iter_per_proc, self.iter_N)

                process = multiprocessing.Process(
                    target=self.threeSeq_multiproc_task,
                    args=(iter_list, exp_id)
                )
                jobs.append(process)

            for j in jobs:
                j.start()

            for j in jobs:
                j.join()

    def GARD_multiproc_task(self, iter_list, exp_id):

        recomb_rate = float(re.findall(r'recombRate_(.*?)_mutRate_(.*)', exp_id)[0][0])
        mut_rate = float(re.findall(r'recombRate_(.*?)_mutRate_(.*)', exp_id)[0][1])

        for i in iter_list:

            with open(str(Path(self.p_dir, exp_id, f"iter_{i}", 'tree_0.tre')), 'r') as f2:
                tree_0 = f2.readline()
            with open(str(Path(self.p_dir, exp_id, f"iter_{i}", 'tree_1.tre')), 'r') as f2:
                tree_1 = f2.readline()

            treeSeq_info, _ = analysis_treeSeq(exp_id, i, rewrite=False, tree_list=[tree_0, tree_1])

            print(f"GARD: {i}th iteration for {exp_id}")
            n_bp = len(treeSeq_info['bp'])
            if n_bp == 1 and treeSeq_info['sample_size'] == 3:
                res_arr = gard_detection(treeSeq_info, exp_id, i)
            else:
                if "GARD" in self.rerun_list:
                    res_arr = gard_detection_2(treeSeq_info, exp_id, i, local=self.local,
                                               rerun=True, max_bp=self.max_bp)
                else:
                    res_arr = gard_detection_2(treeSeq_info, exp_id, i, local=self.local,
                                               rerun=False, max_bp=self.max_bp)
            if res_arr[0] >= 0:
                res_file_path = Path(self.p_dir, exp_id, f"iter_{i}", "GARD_result.txt")
                with open(res_file_path, 'w') as f:
                    f.write("\t".join(map(str, res_arr)) + f'\t{recomb_rate}\t{mut_rate}\t{i}\tGARD\n')


    def gard(self):
        for item in itertools.product(self.recombination_rate,
                                      self.mutation_rate):
            exp_id = "recombRate_" + str(item[0] * self.genome_length) + "_mutRate_" + str(item[1] * self.genome_length)
            print(exp_id)

            iter_per_proc = self.iter_N // NUM_PROC

            jobs = []

            for i in range(NUM_PROC + 1):
                if (i + 1) * iter_per_proc < self.iter_N:
                    iter_list = range(i * iter_per_proc, (i + 1) * iter_per_proc)
                else:
                    iter_list = range(i * iter_per_proc, self.iter_N)

                process = multiprocessing.Process(
                    target=self.GARD_multiproc_task,
                    args=(iter_list, exp_id)
                )
                jobs.append(process)

            for j in jobs:
                j.start()

            for j in jobs:
                j.join()

    def printResult(self, method):
        if not self.res:
            print("Please run detection method!")
            return

        pd.DataFrame.from_dict(self.res, orient='index',
                               columns=["power", "breakpoint distance", "breakpoint width"]).to_csv(
            '{0}_result.csv'.format(method))

    def analyze_result(self, rho, mu, method_list):

        clean_data.get_performance_table(rho, mu, self.iter_N, method_list, local=self.local)
