import os

L=1000
n_list = [100,300,1000,3000]
mu_list = [0.005, 0.01, 0.02, 0.05, 0.1, 0.15, 0.2, 0.25, 0.35, 0.5]

for alpha, JC, aa in [(1.0, False, False), (0.0, True, False), (0.2, False, True), (0.5, False, True)]:
    prefix = '2019-04-28_simulatedData_L%d_alpha%1.1f_%s'%(L, alpha, "aa" if aa else "nuc") + ("_JC" if JC else "")

    for n in n_list:
        for mu in mu_list:
            cmd = 'sbatch submit.sh src/reconstruct_toy_data.py -L {L} -n {n} -m {mu} --prefix {prefix} {JC} {aa} --alpha {alpha}'.format(L=L, n=n, mu=mu, prefix=prefix, JC="--JC" if JC else "", aa="--aa" if aa else "", alpha=alpha)
            print(cmd)
            os.system(cmd)
            
