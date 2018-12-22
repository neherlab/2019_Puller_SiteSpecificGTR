import glob, pickle
import numpy as np
from collections import defaultdict

if __name__ == '__main__':

    mu_dist = defaultdict(lambda: defaultdict(list))
    p_dist = defaultdict(lambda: defaultdict(list))
    W_dist = defaultdict(lambda: defaultdict(list))

    mu_vals = set()
    n_vals = set()

    for fname in glob.glob('2018-12-17_simulated_data_results/*.pkl'):
        with open(fname, 'rb') as fh:
            tmp_mu_vals, tmp_n_vals, tmp_p_dist, tmp_mu_dist, tmp_W_dist = pickle.load(fh)

        n_vals.update(tmp_n_vals)
        mu_vals.update(tmp_mu_vals)

        for d,tmp_d in [(mu_dist, tmp_mu_dist), (W_dist, tmp_W_dist), (p_dist, tmp_p_dist)]:
            for x in tmp_d:
                for y in tmp_d[x]:
                    d[x][y].extend(tmp_d[x][y])


    L=300
    n_vals = sorted(n_vals)
    mu_vals = sorted(mu_vals)

    from matplotlib import pyplot as plt
    for n in n_vals:
        plt.figure()
        for label, data in p_dist.items():
            plt.errorbar(mu_vals, [np.mean(data[(L,n,mu)]) for mu in mu_vals],
                        [np.std(data[(L,n,mu)]) for mu in mu_vals], label=label)

        plt.yscale('log')
        plt.legend()
        plt.xlabel('average rate')


        plt.figure()
        mu_vals = sorted(mu_vals)
        for label, data in mu_dist.items():
            plt.errorbar(mu_vals, [np.mean(data[(L,n,mu)]) for mu in mu_vals],
                        [np.std(data[(L,n,mu)]) for mu in mu_vals], label=label)
        for label, data in W_dist.items():
            plt.errorbar(mu_vals, [np.mean(data[(L,n,mu)]) for mu in mu_vals],
                        [np.std(data[(L,n,mu)]) for mu in mu_vals], label=label)

        plt.yscale('log')
        plt.legend()
        plt.xlabel('average rate')


    plt.figure()
    for li, (label, data) in enumerate(p_dist.items()):
        for n in n_vals:
            d = []
            for mu in mu_vals:
                d.append((n*mu, np.mean(data[(L,n,mu)]), np.std(data[(L,n,mu)])))

            d = np.array(sorted(d, key=lambda x:x[0]))
            plt.errorbar(d[:,0], d[:,1], d[:,2], label=label if n==n_vals[0] else '', c='C%d'%(li+1))

    plt.yscale('log')
    plt.xscale('log')
    plt.legend()
    plt.xlabel('tree length')
