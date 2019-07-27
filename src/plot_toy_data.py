import glob, pickle
import numpy as np
from collections import defaultdict

fmts = ['.png', '.pdf']
fs = 12

labels = {'naive':"Alignment frequencies",
          "dressed":"Known ancestral sequences",
          "regular":"Reconstructed ancestral sequences",
          "marginal":"Marginalized ancestral sequences",
          "iterative":"Iterative model estimation",
          "optimize_tree":"Iterative tree and model optimization",
          "optimize_tree_true":"Optimize tree(true) and model",
          "iterative_true":"Iterative model estimation (true tree)",
          "regular_true":"Reconstructed ancestral sequences",
          "marginal_true":"Marginalized ancestral sequences",
    }
colors = {k:'C%d'%((i+1)%10) for i,k in enumerate(labels.keys())}

def load_toy_data_results(path):
    mu_dist = defaultdict(lambda: defaultdict(list))
    p_dist = defaultdict(lambda: defaultdict(list))
    p_entropy = defaultdict(lambda: defaultdict(list))
    W_dist = defaultdict(lambda: defaultdict(list))
    delta_LH = defaultdict(lambda: defaultdict(list))
    avg_rate = defaultdict(lambda: defaultdict(list))
    mucorr_dist = defaultdict(lambda: defaultdict(list))

    mu_vals = set()
    n_vals = set()

    for fname in glob.glob(path+'*.pkl'):
        with open(fname, 'rb') as fh:
            tmp_mu_vals, tmp_n_vals, tmp_p_dist, tmp_p_entropy, \
                tmp_mu_dist, tmp_W_dist, tmp_LH, tmp_avg_mu, tmp_mucorr = pickle.load(fh)

        n_vals.update(tmp_n_vals)
        mu_vals.update(tmp_mu_vals)

        for d,tmp_d in [(mu_dist, tmp_mu_dist), (W_dist, tmp_W_dist), (p_dist, tmp_p_dist),
                        (p_entropy, tmp_p_entropy), (delta_LH, tmp_LH),
                         (avg_rate, tmp_avg_mu), (mucorr_dist, tmp_mucorr)]:
            for x in tmp_d:
                for y in tmp_d[x]:
                    d[x][y].extend(tmp_d[x][y])

    n_vals = np.array(sorted(n_vals))
    mu_vals = np.array(sorted(mu_vals))

    return {'mu': mu_dist, 'p':p_dist, 'S':p_entropy, 'W':W_dist, 'LH':delta_LH,
            'avg_rate':avg_rate, 'mucorr':mucorr_dist}, n_vals, mu_vals


def plot_pdist_vs_tree_length(data, n_vals, mu_vals, methods=None, fname=None):
    # plot the squared distance of the inferred equilibrium frequencies
    # from the true frequencies as a function of total tree length. we expect this to be 1/n
    # We assume a Yule tree here, that is the total tree length can be calculated as:
    # dk/dt = -k, hence int_t k(t) = n
    plt.figure()
    ls = {n:l for n,l in zip(sorted(n_vals, reverse=True), ['-', '--', '-.', ':'])}
    cols = {n:'C%d'%i for i,n in enumerate(n_vals)}
    for label, dset in data["p"].items():
        if methods and label not in methods:
            continue
        for n in n_vals:
            d = []
            for mu in mu_vals:
                # vector of expected number of subs, mean distance, and std_dev
                d.append((n*mu, np.mean(dset[(L,n,mu)]), np.std(dset[(L,n,mu)])))

            d = np.array(sorted(d, key=lambda x:x[0]))
            plt.errorbar(d[:,0], d[:,1], d[:,2], lw=2, ls=ls[n],
                         label=labels[label] if n==n_vals[-1] else '', c=colors[label])
            if label=='naive':
                plt.text(d[0,0], d[0,1]*1.2, r'$n='+str(n)+'$', fontsize=fs*0.8)

    plt.plot([10,1000], [0.1,0.001], label=r'$\sim x^{-1}$', c='k')
    plt.ylim([0.001, 1.6])
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


def plot_pdist_vs_rtt(data, n_vals, mu_vals, methods, fname=None):
    # for each data set size, plot the distance of the inferred equilibrium frequencies
    # from the their true values. This is plotted vs the root-to-tip distance, which
    # determines the reconstruction accuracy. Given that we use Yule trees, the rtt is
    # roughly log(n)

    plt.figure()
    for n in n_vals:
        rtt = np.log(n)
        for label in methods:
            d=data['p'][label]
            plt.errorbar(mu_vals*rtt, [np.mean(d[(L,n,mu)]) for mu in mu_vals],
                        [np.std(d[(L,n,mu)]) for mu in mu_vals],
                        c=colors[label], label=labels[label] if n==n_vals[0] else '')

    plt.yscale('log')
    plt.xscale('log')
    plt.legend(fontsize=fs*0.8)
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
                            label=labels[label] + ', pc=%1.1f'%pc if n==n_vals[0] and  label!='naive' else '')

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


def plot_rate_correlation(data, n_vals, mu_vals, pc_vals, methods=None, fname=None):
    plt.figure()
    ls={1.5:'-', 3.0:'--'}
    for rate_alpha in data:
        for pi,pc in enumerate(pc_vals):
            dtmp = data[rate_alpha][pc]
            for n in n_vals:
                for label, d in dtmp['mucorr'].items():
                    if methods and label not in methods:
                        continue
                    plt.errorbar(np.array(mu_vals)*n,
                             [np.mean(d[(L,n,mu)], axis=0) for mu in mu_vals],
                             [np.std(d[(L,n,mu)], axis=0) for mu in mu_vals], ls=ls[rate_alpha],
                                 c="C%d"%pi,label='pc=%1.2f, alpha=%1.2f'%(pc, rate_alpha) if n==n_vals[0] else '')

    plt.xscale('log')
    plt.ylim(0,1.1)
    plt.legend(fontsize=fs)
    plt.tick_params(labelsize=0.8*fs)
    plt.xlabel('average number of substitution per site', fontsize=fs)
    plt.ylabel('inferred/true rate correlation', fontsize=fs)
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
    parser.add_argument("--nvals", nargs='+', type=int, default=[1000], help="n values to plot")
    parser.add_argument("--pc", type=float, default=0.1, help="pc value to use in plots")
    parser.add_argument("--rate-alpha", type=float, default=1.5, help="rate variation set to use")
    args=parser.parse_args()

    from matplotlib import pyplot as plt
    data = {}

    n_vals_to_plot = args.nvals
    L=1000
    for rate_alpha in [1.5, 3.0]:
        data[rate_alpha]={}
        for pc in [0.1, 0.5, 1.0]:
            tmp, n_vals, mu_vals = load_toy_data_results(args.prefix.replace('XXX', str(rate_alpha)) + '_results_pc_%1.2f/'%pc)
            data[rate_alpha][pc]=tmp

    aa = 'aa' if '_aa' in args.prefix else 'nuc'
    pc_general = args.pc
    rate_alpha = args.rate_alpha
    suffix = '_%s_ratealpha%1.1f'%(aa, rate_alpha)
    #### Fig1: equilibrium frequency accuracy for all n and mu's
    plot_pdist_vs_tree_length(data[rate_alpha][pc_general], n_vals, mu_vals, methods=['naive', 'dressed'],
                        fname='figures/p_dist_vs_treelength'+suffix)

    #### Fig 2: average rate vs true rate. shows the effect of first order
    #approximation when working off counts. uninformative for other models since
    #mu is set to one -- here we need to compare branch length
    plot_avg_rate(data[rate_alpha][pc_general], n_vals_to_plot, mu_vals, methods=['dressed'],
                        fname='figures/avg_rate_dressed'+suffix)

    plot_rate_correlation(data, n_vals_to_plot, mu_vals, [0.1, 0.5, 1.0], methods=['dressed'],
                        fname='figures/rate_correlation_dressed'+suffix)

    #### Fig 3: comparison of different models as a function of tree length for one pc
    plot_pdist_vs_rtt(data[rate_alpha][pc_general], n_vals_to_plot, mu_vals,
                      methods=['naive', 'regular', 'marginal',
                               'iterative','iterative_true', 'optimize_tree','dressed'],
                      fname='figures/p_dist_vs_rtt'+suffix)

    plot_site_specific_rate_dist(data[rate_alpha][pc_general], n_vals_to_plot, mu_vals,
                      methods=['naive', 'dressed', 'regular', 'marginal',
                                'iterative_true', 'optimize_tree'],
                      fname='figures/mu_dist_vs_rtt'+suffix)

    plot_pentropy_vs_rtt(data[rate_alpha], n_vals_to_plot, mu_vals, pc_vals=[0.1, 0.5, 1.0],
                         methods=['naive', 'dressed'],
                         fname='figures/p_entropy_vs_rtt'+suffix)

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
