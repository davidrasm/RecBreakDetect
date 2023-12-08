# -*- coding: utf-8 -*-
"""
-------------------------------------------------
   File Name：     clean_data
   Description :
   Author :       shicen
   date：          8/14/22
-------------------------------------------------
   Change Activity:
                   8/14/22:
-------------------------------------------------
"""

import numpy as np
import pandas as pd
from pathlib import Path
import itertools


def get_performance_table(rho, mu, iter_N, method_list, p_dir="", local=True):
    table_list = []
    # method_list = ["GARD", "threeSeq", "maxChi"]

    for method in method_list:
        table_list = []
        for r_m in itertools.product(rho, mu):
            exp_str = f"recombRate_{float(r_m[0])}_mutRate_{float(r_m[1])}"

            recomb_param_file = str(Path(p_dir, exp_str, 'treeSeq_info.txt'))
            recomb_param_table = pd.read_csv(recomb_param_file, sep='\t')
            iteration_table_list = []
            for i in range(iter_N):
                exp_id_dir = Path(p_dir, exp_str,
                                  f"iter_{i}")
                print(f"{method} {exp_str} {i}th iteration")
                method_performance_file = (exp_id_dir / f"{method}_result.txt")
                if method_performance_file.is_file():
                    method_performance_table = pd.read_csv(method_performance_file, sep='\t', header=None)
                    iteration_table_list.append(method_performance_table)
            if len(iteration_table_list) > 0:
                method_exp_df = pd.concat(iteration_table_list)
                method_exp_df.columns = ["power", "true_child", "accuracy", "coverage", "precision", "precision2",
                                         "no_of_info_sites", "no_of_weighted_info_sites",
                                         "dist_1", "dist_2",
                                         "recomb_rate", "mut_rate", "iter", "method"]

                recomb_param_table[['recomb_rate', 'mut_rate', 'iter']] = \
                    recomb_param_table[['recomb_rate', 'mut_rate', 'iter']].astype(float)
                method_exp_df[['recomb_rate', 'mut_rate', 'iter']] = \
                    method_exp_df[['recomb_rate', 'mut_rate', 'iter']].astype(float)
                method_table = recomb_param_table.merge(method_exp_df, how='inner',
                                                        on=['recomb_rate', 'mut_rate', 'iter'])
                table_list.append(method_table)
        if len(table_list) > 0:
            method_df = pd.concat(table_list)
            method_df['method'] = method

            method_df = method_df.sort_values(by=["recomb_rate", "mut_rate", "iter"])

            if local:
                method_df.to_csv(str(Path(p_dir, f"{method}_performance_local.txt")), sep='\t', index=False)
            else:
                method_df.to_csv(str(Path(p_dir, f"{method}_performance.txt")), sep='\t', index=False)


if __name__ == '__main__':
    rho = [0.1, 0.5, 1, 2, 4]
    mu = [0.5, 1, 5, 10, 100, 500]

    p_dir = 'C:\\tmp_recomb\\recomb_100'
    method_list = ['threeSeq', 'maxChi', 'GARD']
    get_performance_table(rho, mu, 100, method_list, p_dir)
