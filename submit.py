import os

prefix = '2019-01-12_simulated_data_L1000_aa'

for n in [3000]: #[100,300,1000,3000]:
    for mu in [0.005, 0.01, 0.02, 0.05, 0.1, 0.15, 0.2, 0.25, 0.35, 0.5]:
        os.system('sbatch submit.sh src/reconstruct_toy_data.py -L 1000 -n {n} -m {mu} --prefix {prefix} --aa'.format(n=n, mu=mu, prefix=prefix))
