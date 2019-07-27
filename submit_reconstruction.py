import os
from datetime import date

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

date_prefix = args.date or date.today().strftime('%Y-%m-%d')

for pc in [0.1, 0.5, 1.0]:
    for alpha, JC, aa in [(1.0, False, False), (0.0, True, False), (0.2, False, True), (0.5, False, True)]:
        prefix = '%s_simulatedData_L%d_ratealpha%1.1f_alpha%1.1f_%s'%(date_prefix, L, args.rate_alpha, alpha, "aa" if aa else "nuc") + ("_JC" if JC else "")

        for n in args.nvals:
            submit_script = "submit.sh" if n>200 else "submit_small.sh"
            for mu in args.muvals:
                target = prefix+'_results_pc_%1.2f/L%d_n%d_m%s.pkl'%(pc, L,n, str(mu))
                print(target)
                if not os.path.isfile(target):
                    cmd = 'sbatch {script} src/reconstruct_toy_data.py -L {L} -n {n} -m {mu} --prefix {prefix} {aa} --pc {pc}'.format(
                            script=submit_script, L=L, n=n, mu=mu, prefix=prefix, aa="--aa" if aa else "", pc=pc)
                    print(cmd)
                    if args.submit:
                        os.system(cmd)

