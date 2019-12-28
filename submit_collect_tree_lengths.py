import os, glob
from datetime import date
from itertools import product
import argparse
parser = argparse.ArgumentParser(description = "", usage="generate_toy_data")
parser.add_argument("--nval", type=int,  help="number of taxa")
parser.add_argument("--submit", action='store_true', help="submit")
parser.add_argument("--date", type=str,  help="date prefix")
parser.add_argument("--rate-alpha", type=float, default=1.5, help="parameter of the rate distribution")
args=parser.parse_args()
L=1000
dsets = [(1.0, False, False), (0.0, True, False), (0.2, False, True), (0.5, False, True)]

date_prefix = args.date or date.today().strftime('%Y-%m-%d')

for alpha, JC, aa in dsets:
    dirname = f'{date_prefix}_simulatedData_L{L}_ratealpha{args.rate_alpha:1.1f}_freqalpha{alpha:1.1f}_{"aa" if aa else "nuc"}' + ("_JC" if JC else "")
    output = dirname + f"_results/collected_tree_lengths_n{args.nval}.tsv"
    cmd = f"sbatch submit_small.sh src/collect_tree_lengths.py --output {output} --trueTrees {dirname} --reoptimizedTrees {dirname+'_results'}  --nval {args.nval}"
    print(cmd)
    if args.submit:
        os.system(cmd)
