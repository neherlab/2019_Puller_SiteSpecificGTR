import glob, pickle
import numpy as np
from collections import defaultdict

if __name__ == '__main__':

    mu_dist = defaultdict(lambda: defaultdict(list))
    p_dist = defaultdict(lambda: defaultdict(list))
    W_dist = defaultdict(lambda: defaultdict(list))
    delta_LH = defaultdict(lambda: defaultdict(list))
    avg_rate = defaultdict(lambda: defaultdict(list))

    mu_vals = set()
    n_vals = set()

    for fname in glob.glob('2018-12-17_simulated_data_results/*.pkl'):
        with open(fname, 'rb') as fh:
            tmp_mu_vals, tmp_n_vals, tmp_p_dist, tmp_mu_dist, tmp_W_dist, tmp_LH, tmp_avg_mu = pickle.load(fh)

        n_vals.update(tmp_n_vals)
        mu_vals.update(tmp_mu_vals)

        for d,tmp_d in [(mu_dist, tmp_mu_dist), (W_dist, tmp_W_dist), (p_dist, tmp_p_dist),
                        (delta_LH, tmp_LH), (avg_rate, tmp_avg_mu)]:
            for x in tmp_d:
                for y in tmp_d[x]:
                    d[x][y].extend(tmp_d[x][y])


    L=300
    n_vals = np.array(sorted(n_vals))
    mu_vals = np.array(sorted(mu_vals))

    from matplotlib import pyplot as plt

    # for each data set size, plot the distance of the inferred equilibrium frequencies
    # and the substitution rates from the true frequencies and rates
    for n in n_vals:
        rtt = np.log(n)
        plt.figure()
        for label, data in p_dist.items():
            plt.errorbar(mu_vals*rtt, [np.mean(data[(L,n,mu)]) for mu in mu_vals],
                        [np.std(data[(L,n,mu)]) for mu in mu_vals], label=label)

        plt.yscale('log')
        plt.legend()
        plt.xlabel('average root-to-tip distance')
        plt.ylabel('squared deviation')


        plt.figure()
        for label, data in mu_dist.items():
            plt.errorbar(mu_vals*rtt, [np.mean(data[(L,n,mu)])/L for mu in mu_vals],
                        [np.std(data[(L,n,mu)])/L for mu in mu_vals], label=label)
        # for label, data in W_dist.items():
        #     plt.errorbar(mu_vals*rtt, [np.mean(data[(L,n,mu)]) for mu in mu_vals],
        #                 [np.std(data[(L,n,mu)]) for mu in mu_vals], label=label)

        plt.yscale('log')
        plt.legend()
        plt.xlabel('average root-to-tip distance')
        plt.ylabel('relative squared deviation of site specific rates')


        plt.figure()
        for label, data in delta_LH.items():
            plt.errorbar(mu_vals*rtt, [-np.mean(np.diff(data[(L,n,mu)], axis=1)) for mu in mu_vals],
                        [np.std(np.diff(data[(L,n,mu)], axis=1)) for mu in mu_vals], label=label)

        #plt.yscale('log')
        plt.legend()
        plt.xlabel('average root-to-tip distance')
        plt.ylabel('delta LH')

        plt.figure()
        for label, data in avg_rate.items():
            if label=='phylo': continue
            plt.plot([np.mean(data[(L,n,mu)], axis=0)[0] for mu in mu_vals],
                     [np.mean(data[(L,n,mu)], axis=0)[1] for mu in mu_vals], label=label)

        plt.plot([0, 0.35], [0, 0.35], c='k', label='correct')
        plt.legend()
        plt.xlabel('average substitution rate')
        plt.ylabel('inferred substitution rate')

    # plot the squared distance of the inferred equilibrium frequencies
    # from the true frequencies as a function of total tree length. we expect this to be 1/n
    plt.figure()
    for li, (label, data) in enumerate(p_dist.items()):
        if label not in ['naive', 'dressed', 'phylo']:
            continue
        for n in n_vals:
            d = []
            for mu in mu_vals:
                d.append((n*mu, np.mean(data[(L,n,mu)]), np.std(data[(L,n,mu)])))

            d = np.array(sorted(d, key=lambda x:x[0]))
            plt.errorbar(d[:,0], d[:,1], d[:,2], label=label if n==n_vals[0] else '', c='C%d'%(li+1))

    plt.plot([10,1000], [0.1,0.001], label=r'$\sim x^{-1}$', c='k')
    plt.yscale('log')
    plt.xscale('log')
    plt.legend()
    plt.xlabel('average number of substitutions per site')
    plt.ylabel('squared deviation')


