import os
from datetime import date

L=1000
muvals =  [0.005, 0.01, 0.02, 0.05, 0.1, 0.15, 0.2, 0.25, 0.35, 0.5]

import argparse
parser = argparse.ArgumentParser(description = "", usage="generate_toy_data")
parser.add_argument("--nvals", nargs='+', type=int,  help="number of taxa")
parser.add_argument("--muvals", nargs='+', default=muvals, type=float,  help="number of taxa")
parser.add_argument("--date", type=str,  help="date prefix")
parser.add_argument("--submit", action='store_true', help="submit")
args=parser.parse_args()
today = args.date or date.today().strftime('%Y-%m-%d')

#           alpha, rate_alpha, JC, amino-acid
datasets = [(1.0, 1.5, False, False) ]#,
#            (0.0, 1.5, True, False)],
#            (0.2, 1.5, False, True),
#            (0.5, 1.5, False, True)]

for alpha, rate_alpha, JC, aa in datasets:
    prefix = f'{today}_simulatedData_L{L}_ratealpha{rate_alpha:1.1f}_freqalpha{alpha:1.1f}_{"aa" if aa else "nuc"}' + ("_JC" if JC else "")

    for n in args.nvals:
        submit_script = "submit.sh" if n>500 else "submit_small.sh"
        for mu in args.muvals:
            cmd = 'sbatch {script} src/generate_toy_data.py -L {L} -n {n} -m {mu} --prefix {prefix} {JC} {aa} --alpha {alpha} --rate-alpha {rate_alpha}'.format(
                   script=submit_script, L=L, n=n, mu=mu, prefix=prefix, JC="--JC" if JC else "", aa="--aa" if aa else "", alpha=alpha, rate_alpha=rate_alpha)
            target = prefix+'/L{L}_n{n}_m{m}_tree00.nwk'.format(n=n, m=mu, L=L)
            print(target)
            if not os.path.isfile(target):
                print(cmd)
                if args.submit:
                    os.system(cmd)

