import os, glob
from datetime import date
from itertools import product
from submit_reconstruction import make_chunks
from math import ceil

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

dsets = [(1.0, False, False), (0.0, True, False), (0.2, False, True), (0.5, False, True)]

date_prefix = args.date or date.today().strftime('%Y-%m-%d')

for pc in [0.01, 0.1, 0.5, 1.0]:
    for alpha, JC, aa in dsets:
        dirname = f'{date_prefix}_simulatedData_L{L}_ratealpha{args.rate_alpha:1.1f}_freqalpha{alpha:1.1f}_{"aa" if aa else "nuc"}' + ("_JC" if JC else "")
        print(dirname)
        for n in args.nvals:
            submit_script = "submit.sh" if n>200 else "submit_small.sh"
            for mu in args.muvals:
                aln_mask = os.path.join(dirname, f"L{L}_n{n}_m{mu}_*.fasta.gz")
                files = glob.glob(aln_mask)
                output = dirname+f"_results"

                if len(files):
                    chunks = make_chunks(files, n=ceil(n/100))
                    for ci,chunk in enumerate(chunks):
                        if len(chunk):
                            cmd = f'sbatch {submit_script} src/calculate_branch_length.py --files {" ".join(chunk)} --output {output} {"--aa" if aa else ""} --pc {pc}'
                            print(cmd)
                            if args.submit:
                                os.system(cmd)
                                os.system(cmd + ' --true-model')
                                os.system(cmd + ' --true-rates')

