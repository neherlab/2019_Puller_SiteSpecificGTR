import os
from datetime import date
from itertools import product

L=1000

import argparse
parser = argparse.ArgumentParser(description = "", usage="generate_toy_data")
parser.add_argument("--nvals", nargs='+', type=int,  help="number of taxa")
parser.add_argument("--muvals", nargs='+', type=float,  help="mutation rate for panelB")
parser.add_argument("--submit", action='store_true', help="submit")
parser.add_argument("--date", type=str,  help="date prefix")
parser.add_argument("--rate-alpha", type=float, default=1.5, help="parameter of the rate distribution")

args=parser.parse_args()
aa = False
JC = False
alpha=1.0

date_prefix = args.date or date.today().strftime('%Y-%m-%d')

prefix = '%s_simulatedData_L%d_ratealpha%1.1f_alpha%1.1f_%s'%(date_prefix, L, args.rate_alpha, alpha, "aa" if aa else "nuc") + ("_JC" if JC else "")
for n in args.nvals:
    submit_script = "submit.sh" if n>500 else "submit_small.sh"
    for mu in args.muvals:
        cmd = f"sbatch {submit_script} src/model_deviation.py -L {L} -n {n} -m {mu} --prefix {prefix}"
        print(cmd)
        if args.submit:
            os.system(cmd)

