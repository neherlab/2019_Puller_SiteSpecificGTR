import glob, pickle
import numpy as np
from collections import defaultdict

fmts = ['.png', '.pdf']
fs = 12

labels = {'naive':"Alignment frequencies", "dressed":"True substitution counts", "regular":"Reconstructed",
          "marginal":"ancestral sum", "iterative":"Iterative, reconstructed tree", 'branch_length':"non-linear",
          "iterative_true":"iterative, true tree", 'marginal_true':"ancestral sum, true tree",
          "optimize_tree":"Optimize tree and model", "optimize_tree_true":"Optimize tree(true) and model"}
colors = {k:'C%d'%((i+1)%10) for i,k in enumerate(sorted(labels.keys()))}

def load_toy_data_results(path):
    mu_dist = defaultdict(lambda: defaultdict(list))
    p_dist = defaultdict(lambda: defaultdict(list))
    p_entropy = defaultdict(lambda: defaultdict(list))
    W_dist = defaultdict(lambda: defaultdict(list))
    delta_LH = defaultdict(lambda: defaultdict(list))
    avg_rate = defaultdict(lambda: defaultdict(list))

    mu_vals = set()
    n_vals = set()

    for fname in glob.glob(path+'*.pkl'):
        with open(fname, 'rb') as fh:
            tmp_mu_vals, tmp_n_vals, tmp_p_dist, tmp_p_entropy, tmp_mu_dist, tmp_W_dist, tmp_LH, tmp_avg_mu = pickle.load(fh)

        n_vals.update(tmp_n_vals)
        mu_vals.update(tmp_mu_vals)

        for d,tmp_d in [(mu_dist, tmp_mu_dist), (W_dist, tmp_W_dist), (p_dist, tmp_p_dist),
                        (p_entropy, tmp_p_entropy), (delta_LH, tmp_LH), (avg_rate, tmp_avg_mu)]:
            for x in tmp_d:
                for y in tmp_d[x]:
                    d[x][y].extend(tmp_d[x][y])

    n_vals = np.array(sorted(n_vals))
    mu_vals = np.array(sorted(mu_vals))

    return {'mu': mu_dist, 'p':p_dist, 'S':p_entropy, 'W':W_dist, 'LH':delta_LH, 'avg_rate':avg_rate}, n_vals, mu_vals


def plot_pdist_vs_tree_length(data, n_vals, mu_vals, methods=None, fname=None):
    # plot the squared distance of the inferred equilibrium frequencies
    # from the true frequencies as a function of total tree length. we expect this to be 1/n
    # We assume a Yule tree here, that is the total tree length can be calculated as:
    # dk/dt = -k, hence int_t k(t) = n
    plt.figure()
    for label, dset in data["p"].items():
        if methods and label not in methods:
            continue
        for n in n_vals:
            d = []
            for mu in mu_vals:
                # vector of expected number of subs, mean distance, and std_dev
                d.append((n*mu, np.mean(dset[(L,n,mu)]), np.std(dset[(L,n,mu)])))

            d = np.array(sorted(d, key=lambda x:x[0]))
            plt.errorbar(d[:,0], d[:,1], d[:,2],
                         label=labels[label] if n==n_vals[0] else '', c=colors[label])

    plt.plot([10,1000], [0.1,0.001], label=r'$\sim x^{-1}$', c='k')
    plt.yscale('log')
    plt.xscale('log')
    plt.legend(fontsize=fs)
    plt.tick_params(labelsize=0.8*fs)
    plt.xlabel('average number of substitutions per site', fontsize=fs)
    plt.ylabel('squared deviation of $p_i^a$', fontsize=fs)
    plt.tight_layout()

    if fname:
        for fmt in fmts:
            plt.savefig(fname+fmt)


def plot_pdist_vs_rtt(data, n_vals, mu_vals, methods=None, fname=None):
    # for each data set size, plot the distance of the inferred equilibrium frequencies
    # from the their true values. This is plotted vs the root-to-tip distance, which
    # determines the reconstruction accuracy. Given that we use Yule trees, the rtt is
    # roughly log(n)

    plt.figure()
    for n in n_vals:
        rtt = np.log(n)
        for label, d in data['p'].items():
            if methods and label not in methods:
                continue
            plt.errorbar(mu_vals*rtt, [np.mean(d[(L,n,mu)]) for mu in mu_vals],
                        [np.std(d[(L,n,mu)]) for mu in mu_vals],
                        c=colors[label], label=labels[label] if n==n_vals[0] else '')

    plt.yscale('log')
    plt.xscale('log')
    plt.legend(fontsize=fs)
    plt.tick_params(labelsize=0.8*fs)
    plt.xlabel('root-to-tip distance', fontsize=fs)
    plt.ylabel('squared deviation', fontsize=fs)
    plt.tight_layout()

    if fname:
        for fmt in fmts:
            plt.savefig(fname+fmt)


def plot_pentropy_vs_rtt(data, n_vals, mu_vals, pc_vals, methods=None, fname=None):
    # for each data set size, plot the difference in entropy between the inferred and true
    # equilibrium frequencies.
    line_styles = ['--','-.','-']
    plt.figure()
    for n in n_vals:
        rtt = np.log(n)
        for pi,pc in enumerate(pc_vals):
            for label, d in data[pc]['S'].items():
                if methods and label not in methods:
                    continue
                plt.errorbar(mu_vals*rtt, [-np.mean(np.diff(d[(L,n,mu)], axis=1)) for mu in mu_vals],
                            [np.std(np.diff(d[(L,n,mu)], axis=1)) for mu in mu_vals],
                            c=colors[label], ls=line_styles[pi],
                            label=labels[label] if n==n_vals[0] and pc==pc_vals[-1] else '')

    plt.xscale('log')
    plt.legend(fontsize=fs)
    plt.tick_params(labelsize=0.8*fs)
    plt.xlabel('root-to-tip distance', fontsize=fs)
    plt.ylabel('entropy difference (inferred-true)', fontsize=fs)
    plt.tight_layout()

    if fname:
        for fmt in fmts:
            plt.savefig(fname+fmt)


def plot_avg_rate(data, n_vals, mu_vals, methods=None, fname=None):
    plt.figure()
    for n in n_vals:
        for label, d in data['avg_rate'].items():
            if methods and label not in methods:
                continue
            plt.plot([np.mean(d[(L,n,mu)], axis=0)[0] for mu in mu_vals],
                     [np.mean(d[(L,n,mu)], axis=0)[1] for mu in mu_vals],
                      c=colors[label], label=labels[label] if n==n_vals[0] else '')

    plt.plot([0, 0.35], [0, 0.35], c='k', label='correct')
    plt.legend(fontsize=fs)
    plt.tick_params(labelsize=0.8*fs)
    plt.xlabel('average substitution rate', fontsize=fs)
    plt.ylabel('inferred substitution rate', fontsize=fs)
    plt.tight_layout()

    if fname:
        for fmt in fmts:
            plt.savefig(fname+fmt)


def plot_site_specific_rate_dist(data, n_vals, mu_vals, methods=None, fname=None):
    plt.figure()
    for n in n_vals:
        rtt = np.log(n)
        for label, d in data['mu'].items():
            if methods and label not in methods:
                continue
            plt.errorbar(mu_vals*rtt, [np.mean(d[(L,n,mu)])/L for mu in mu_vals],
                        [np.std(d[(L,n,mu)])/L for mu in mu_vals],
                         c=colors[label], label=labels[label] if n==n_vals[0] else '')

    plt.yscale('log')
    plt.legend(fontsize=fs)
    plt.tick_params(labelsize=0.8*fs)
    plt.xlabel('average root-to-tip distance', fontsize=fs)
    plt.ylabel('relative squared deviation of site specific rates', fontsize=fs)
    plt.tight_layout()

    if fname:
        for fmt in fmts:
            plt.savefig(fname+fmt)


if __name__ == '__main__':
    import argparse
    parser = argparse.ArgumentParser(description = "", usage="plot reconstructions")
    parser.add_argument("--prefix", type=str, help="prefix of data set")
    args=parser.parse_args()

    from matplotlib import pyplot as plt
    data = {}

    n_vals_to_plot = [100,300]
    L=1000
    for pc in [0.1, 0.5, 1.0]:
        tmp, n_vals, mu_vals = load_toy_data_results(args.prefix + '_results_pc_%1.2f/'%pc)
        data[pc]=tmp

    pc_general = 0.1
    ####
    plot_pdist_vs_tree_length(data[pc_general], n_vals, mu_vals, methods=['naive', 'dressed'],
                        fname='figures/p_dist_vs_treelength')
    plot_avg_rate(data[pc_general], n_vals_to_plot, mu_vals, methods=['dressed'],
                        fname='figures/avg_rate_dressed')

    ####
    plot_pdist_vs_rtt(data[pc_general], n_vals_to_plot, mu_vals,
                      methods=['naive', 'dressed', 'regular', 'marginal','iterative',
                               'iterative_true', 'optimize_tree'],
                      fname='figures/p_dist_vs_rtt')

    plot_site_specific_rate_dist(data[pc_general], n_vals_to_plot, mu_vals,
                      methods=['naive', 'dressed', 'regular', 'marginal', 'iterative',
                                'iterative_true', 'optimize_tree'],
                      fname='figures/mu_dist_vs_rtt')

    plot_pentropy_vs_rtt(data, n_vals_to_plot, mu_vals, pc_vals=[0.1, 0.5, 1.0],
                         methods=['naive', 'dressed', 'iterative', 'optimize_tree'],
                         fname='figures/p_entropy_vs_rtt')

    # for each data set size, plot the distance of the inferred equilibrium frequencies
    # and the substitution rates from the true frequencies and rates
    for n in []:

        plt.figure()
        for label, d in data["LH"].items():
            plt.errorbar(mu_vals*rtt, [-np.mean(np.diff(d[(L,n,mu)], axis=1)) for mu in mu_vals],
                        [np.std(np.diff(d[(L,n,mu)], axis=1)) for mu in mu_vals], label=label)

        #plt.yscale('log')
        plt.legend()
        plt.xlabel('average root-to-tip distance')
        plt.ylabel('delta LH')
