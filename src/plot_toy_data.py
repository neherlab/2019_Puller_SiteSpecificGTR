import glob, pickle, os
import numpy as np
import pandas as pd
from collections import defaultdict
from matplotlib import pyplot as plt

def add_panel_label(ax, label, x_offset=-0.1, y_offset=0.95, fs=12):
    '''Add a label letter to a panel'''
    ax.text(x_offset, y_offset, label,
            transform=ax.transAxes,
            fontsize=fs*1.5)

fmts = ['.png', '.pdf']
fs = 12

labels = {'naive':"Alignment frequencies",
          "dressed":"Known ancestral sequences",
          "ml_reconstruction":"Reconstructed ancestral sequences",
          "marginalize":"Marginalized ancestral sequences",
          "iterative":"Iterative model estimation",
          "optimize_tree":"Iterative tree and model optimization",
          "optimize_tree_true":"Optimize tree(true) and model",
          "iterative_true":"Iterative model estimation (true tree)",
          "ml_reconstruction_true":"Reconstructed ancestral sequences",
          "marginalize_true":"Marginalized ancestral sequences",
    }
colors = {k:'C%d'%((i+1)%10) for i,k in enumerate(labels.keys())}


def load_toy_data_results(path):
    tsv_files = glob.glob(os.path.join(path, "L*tsv"))
    df = pd.concat([pd.read_csv(fname, sep='\t') for fname in tsv_files])
    return df

def make_subset(data, conditions):
    ind = np.ones(len(data), dtype=bool)
    for k, v in conditions.items():
        if k in data.columns:
            ind = ind&(data[k]==v)
        else:
            raise ValueError(f"column '{k}' not in dataframe")

    return data[ind]

def make_means(data, conditions):
    from itertools import product
    keys = list(conditions.keys())
    values = [conditions[x] for x in keys]
    res = {}
    for combination in product(*values):
        tmp_cond =  {k:v for k,v in zip(keys, combination)}
        print(tmp_cond)
        subset = make_subset(data, tmp_cond)
        res[tuple(tmp_cond.items())] = {'mean':subset.groupby('m').mean(),
                                        'std':subset.groupby('m').std()}
    return res


def plot_pdist_vs_tree_length(data, data_secondary, subsets, fname=None):
    # plot the squared distance of the inferred equilibrium frequencies
    # from the true frequencies as a function of total tree length. we expect this to be 1/n
    # We assume a Yule tree here, that is the total tree length can be calculated as:
    # dk/dt = -k, hence int_t k(t) = n
    plt.figure()
    ls = {n:l for n,l in zip(sorted(subsets['n'], reverse=True), ['-', '--', '-.', ':'])}
    cols = {n:'C%d'%i for i,n in enumerate(subsets['n'])}

    mean_vals = make_means(data, subsets)
    for params, d in mean_vals.items():
        tmp_p = {k:v for k,v in params}
        subs_per_site = d['mean'].index*tmp_p['n']
        plt.errorbar(subs_per_site, d['mean']['chisq_p'], d['std']['chisq_p'],
                     label=labels[tmp_p['method']] if tmp_p['n']==np.max(subsets['n']) else '',
                     c=colors[tmp_p['method']], ls=ls[tmp_p['n']])

        if tmp_p['method']=='naive':
            plt.text(subs_per_site[0], d['mean']['chisq_p'].iloc[0]*1.2,
                    f"n={tmp_p['n']}", fontsize=fs*0.8)

    if not (data_secondary is None):
        mean_vals_secondary = make_means(data_secondary, subsets)
        for params, d in mean_vals_secondary.items():
            tmp_p = {k:v for k,v in params}
            subs_per_site = d['mean'].index*tmp_p['n']
            plt.errorbar(subs_per_site, d['mean']['chisq_p'], d['std']['chisq_p'],
                         label='', c="#AAAAAA", ls=ls[tmp_p['n']])


    plt.plot([10,1000], [0.1,0.001], label=r'$\sim x^{-1}$', c='k')
    plt.ylim([0.001, 1.6])
    plt.yscale('log')
    plt.xscale('log')
    plt.legend(fontsize=fs)
    plt.tick_params(labelsize=0.8*fs)
    plt.xlabel('average number of substitutions per site', fontsize=fs)
    plt.ylabel('squared deviation of $p_i^a$', fontsize=fs)
    add_panel_label(plt.gca(),'A',  x_offset=-0.12, fs=fs*1.3)
    plt.tight_layout()

    if fname:
        for fmt in fmts:
            plt.savefig(fname+fmt)


def plot_pdist_vs_rtt(data, subsets, fname=None):
    # for each data set size, plot the distance of the inferred equilibrium frequencies
    # from the their true values. This is plotted vs the root-to-tip distance, which
    # determines the reconstruction accuracy. Given that we use Yule trees, the rtt is
    # roughly log(n)

    plt.figure()
    mean_vals = make_means(data, subsets)
    for params, d in mean_vals.items():
        tmp_p = {k:v for k,v in params}
        rtt = d['mean'].index*np.log(tmp_p['n'])
        plt.errorbar(rtt, d['mean']['chisq_p'], d['std']['chisq_p'],
                     label=labels[tmp_p['method']] if tmp_p['n']==subsets['n'][0] else '',
                     c=colors[tmp_p['method']])

    plt.yscale('log')
    plt.xscale('log')
    plt.legend(fontsize=fs*0.8, loc=3)
    plt.tick_params(labelsize=0.8*fs)
    plt.xlabel('root-to-tip distance', fontsize=fs)
    plt.ylabel('squared deviation', fontsize=fs)
    plt.title('amino acids' if 'aa' in fname else 'nucleotides', fontsize=fs*1.3)
    add_panel_label(plt.gca(), 'B' if 'aa' in fname else 'A', x_offset=-0.12, fs=fs*1.3)
    plt.tight_layout()

    if fname:
        for fmt in fmts:
            plt.savefig(fname+fmt)


def plot_pentropy_vs_rtt(data, subsets, fname=None):
    # for each data set size, plot the difference in entropy between the inferred and true
    # equilibrium frequencies.
    line_styles = ['--','-.','-', ':']
    mean_vals = make_means(data, subsets)

    plt.figure()
    for params, d in mean_vals.items():
        tmp_p = {k:v for k,v in params}
        pc, label, n = tmp_p['pc'], tmp_p['method'], tmp_p['n']

        rtt = d['mean'].index*np.log(n)
        plt.plot(rtt, d['mean']['model_entropy']-d['mean']['true_entropy'],
                            c=colors[tmp_p['method']], ls=line_styles[subsets['pc'].index(pc)],
                            label=f"{labels[label]}, c={pc:1.2f}" if n==subsets['n'][0] and  label!='naive' else '')

    plt.xscale('log')
    plt.legend(fontsize=fs)
    plt.tick_params(labelsize=0.8*fs)
    plt.xlabel('root-to-tip distance', fontsize=fs)
    plt.ylabel('entropy difference (inferred-true)', fontsize=fs)
    plt.tight_layout()

    if fname:
        for fmt in fmts:
            plt.savefig(fname+fmt)


def plot_avg_rate(data, subsets, fname=None):
    mean_vals = make_means(data, subsets)

    plt.figure()
    for params, d in mean_vals.items():
        tmp_p = {k:v for k,v in params}
        plt.errorbar(d['mean']['true_average_mu'], d['mean']['average_mu'], d['std']['average_mu'],
                      c=colors[tmp_p['method']],
                      label=labels[tmp_p['method']] if tmp_p['n']==np.max(subsets['n']) else '')

    plt.plot([0, 0.35], [0, 0.35], c='k', label='correct')
    plt.legend(fontsize=fs)
    plt.tick_params(labelsize=0.8*fs)
    plt.xlabel('average substitution rate', fontsize=fs)
    plt.ylabel('inferred substitution rate', fontsize=fs)
    plt.tight_layout()

    if fname:
        for fmt in fmts:
            plt.savefig(fname+fmt)


def plot_rate_correlation(data, subsets, fname=None):
    mean_vals = make_means(data, subsets)
    ls=['-', '--']
    plt.figure()
    for params, d in mean_vals.items():
        tmp_p = {k:v for k,v in params}
        pc, label, n = tmp_p['pc'], tmp_p['method'], tmp_p['n']
        plt.errorbar(d['mean'].index*n, d['mean']['r_mu'],d['std']['r_mu'],
                     c="C%d"%subsets['pc'].index(pc), label=f'c={pc:1.2f}, n={n}',
                     ls=ls[subsets['n'].index(n)])

    plt.xscale('log')
    plt.ylim(0,1.1)
    plt.legend(fontsize=fs)
    plt.tick_params(labelsize=0.8*fs)
    plt.xlabel('average number of substitution per site', fontsize=fs)
    plt.ylabel('inferred/true rate correlation', fontsize=fs)
    plt.plot([10,max(subsets['n'])/3], [1,1], lw=2, c='#CCCCCC')
    add_panel_label(plt.gca(), 'B', x_offset=-0.12, fs=fs*1.3)
    plt.tight_layout()

    if fname:
        for fmt in fmts:
            plt.savefig(fname+fmt)


def plot_site_specific_rate_dist(data, subsets, fname=None):
    mean_vals = make_means(data, subsets)

    plt.figure()
    for params, d in mean_vals.items():
        tmp_p = {k:v for k,v in params}
        pc, label, n = tmp_p['pc'], tmp_p['method'], tmp_p['n']
        rtt = d['mean'].index*np.log(tmp_p['n'])
        plt.errorbar(rtt, d['mean']['chisq_mu'], d['std']['chisq_mu'],
                     c=colors[label], label=labels[label] if n==subsets['n'][0] else '')

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
    args=parser.parse_args()

    aa = 'aa' if '_aa' in args.prefix else 'nuc'
    pc_general = args.pc
    rate_alpha = 1.5
    suffix = '_%s_ratealpha%1.1f'%(aa, rate_alpha)

    if aa=='nuc':
        data_3 = load_toy_data_results(args.prefix.replace('XXX', '3.0'))
    else:
        data_3 = None
    data = load_toy_data_results(args.prefix.replace('XXX', str(rate_alpha)))
    n_vals_to_plot = args.nvals
    n_vals = np.unique(data['n'])

    #### Fig1: equilibrium frequency accuracy for all n and mu's
    plot_pdist_vs_tree_length(data, data_3, subsets = {'pc':[pc_general], 'method':['naive', 'dressed'], 'n':n_vals},
                                fname='figures/p_dist_vs_treelength'+suffix)

    #### Fig 2: average rate vs true rate. shows the effect of first order
    #approximation when working off counts. uninformative for other models since
    #mu is set to one -- here we need to compare branch length
    plot_avg_rate(data, subsets = {'pc':[pc_general], 'method':['naive', 'dressed'], 'n':n_vals_to_plot},
                        fname='figures/avg_rate_dressed'+suffix)

    plot_rate_correlation(data,  subsets = {'pc':[0.01, 0.1, 0.5, 1.0], 'method':['dressed'], 'n':n_vals_to_plot},
                        fname='figures/rate_correlation_dressed'+suffix)

    #### comparison of different inference schemes as a function of tree length for one pc
    plot_pdist_vs_rtt(data, subsets={'pc':[pc_general], 'n':n_vals_to_plot,
                                    'method':['naive', 'ml_reconstruction', 'marginalize',
                                              'iterative','iterative_true', 'optimize_tree','dressed']},
                      fname='figures/p_dist_vs_rtt'+suffix)

    plot_site_specific_rate_dist(data, subsets={'pc':[pc_general], 'n':n_vals_to_plot,
                                  'method':['naive', 'dressed', 'ml_reconstruction', 'marginalize',
                                            'iterative_true', 'optimize_tree']},
                      fname='figures/mu_dist_vs_rtt'+suffix)

    plot_pentropy_vs_rtt(data,  subsets = {'pc':[0.01, 0.1, 0.5, 1.0], 'method':['naive', 'dressed'], 'n':n_vals_to_plot},
                         fname='figures/p_entropy_vs_rtt'+suffix)

