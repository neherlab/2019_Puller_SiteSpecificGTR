import os

for n in [300, 1000]:
    for mu in [0.005, 0.01, 0.02, 0.05, 0.1, 0.15, 0.2, 0.25, 0.35, 0.5]:
        os.system('sbatch submit.sh src/generate_toy_data.py -L 300 -n {n} -m {mu}'.format(n=n, mu=mu))
