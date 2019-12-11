import os
from datetime import date
from itertools import product

L=1000
muvals =  [0.005, 0.01, 0.02, 0.05, 0.1, 0.15, 0.2, 0.25, 0.35, 0.5]

import argparse
parser = argparse.ArgumentParser(description = "", usage="generate_toy_data")
parser.add_argument("--nvals", nargs='+', type=int,  help="number of taxa")
parser.add_argument("--muvals", nargs='+', default=muvals, type=float,  help="number of taxa")
parser.add_argument("--submit", action='store_true', help="submit")
parser.add_argument("--date", type=str,  help="date prefix")
parser.add_argument("--rate-alpha", type=float, default=1.5, help="parameter of the rate distribution")
args=parser.parse_args()
#dsets = [(1.0, False, False)] #, (0.0, True, False)]
dsets = [(0.2, False, True), (0.5, False, True)]

date_prefix = args.date or date.today().strftime('%Y-%m-%d')

for pc, model, rates in [(10.0, False, False),(1.0, False, False), (0.1, False, False),
                         (10.0, False, True),(1.0, False, True), (0.1, False, True), (0,True, False)]:
    for alpha, JC, aa in dsets:
        prefix = '%s_simulatedData_L%d_ratealpha%1.1f_alpha%1.1f_%s'%(date_prefix, L, args.rate_alpha, alpha, "aa" if aa else "nuc") + ("_JC" if JC else "")

        for n in args.nvals:
            submit_script = "submit.sh" if n>500 else "submit_small.sh"
            for mu in args.muvals:
                for t,g,s in product([0,1], [0,1], [0,1]):
                    cmd = f"sbatch {submit_script} src/calculate_branch_length.py -L {L} -n {n} -m {mu} -t {t} -g {g} -s {s} --prefix {prefix} {'--aa' if aa else ''} --pc {pc} {'--true-model' if model else ''} {'--true-rates' if rates else ''} --true-tree"
                    print(cmd)
                    if args.submit:
                        os.system(cmd)

