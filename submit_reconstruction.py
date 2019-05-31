import os

L=1000
n_list = [100, 300] #, 1000,3000]
mu_list = [0.005, 0.01, 0.02, 0.05, 0.1, 0.15, 0.2, 0.25, 0.35, 0.5]

for pc in [0.1, 0.5, 1.0]:
    for alpha, JC, aa in [(1.0, False, False), (0.0, True, False), (0.2, False, True), (0.5, False, True)]:
        prefix = '2019-05-31_simulatedData_L%d_alpha%1.1f_%s'%(L, alpha, "aa" if aa else "nuc") + ("_JC" if JC else "")

        for n in n_list:
            for mu in mu_list:
            	target = prefix+'_results_pc_%1.2f/L%d_n%d_m%s.pkl'%(pc, L,n, str(mu))
            	print(target)
            	if not os.path.isfile(target):
	                cmd = 'sbatch submit_small.sh src/reconstruct_toy_data.py -L {L} -n {n} -m {mu} --prefix {prefix} {aa} --pc {pc}'.format(L=L, n=n, mu=mu, prefix=prefix, aa="--aa" if aa else "", pc=pc)
	                print(cmd)
	                os.system(cmd)

