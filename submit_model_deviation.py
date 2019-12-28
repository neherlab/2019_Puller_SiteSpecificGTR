import os, glob
from datetime import date

if __name__ == '__main__':

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

    for alpha, JC, aa in [(1.0, False, False), (0.0, True, False), (0.2, False, True), (0.5, False, True)]:
        dirname = f'{date_prefix}_simulatedData_L{L}_ratealpha{args.rate_alpha:1.1f}_freqalpha{alpha:1.1f}_{"aa" if aa else "nuc"}' + ("_JC" if JC else "")

        for n in args.nvals:
            submit_script = "submit.sh" if n>500 else "submit_small.sh"
            for mu in args.muvals:
                aln_mask = os.path.join(dirname, f"L{L}_n{n}_m{mu}_*.fasta.gz")
                files = glob.glob(aln_mask)
                if len(files):
                    output = os.path.join(dirname+f"_results", f"ModelDeviation_L{L}_n{n}_m{mu}.tsv")
                    cmd = f'sbatch {submit_script} src/model_deviation_scan_epsilon.py --files {" ".join(files)} --output {output} {"--aa" if aa else ""}'
                    print(cmd)
                    if args.submit:
                        os.system(cmd)
