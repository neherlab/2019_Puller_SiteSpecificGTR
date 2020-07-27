'''
script that reads in data obtained from mixing Jukes-Cantor type models with the true
model to explore how sensitive reconstruction of branch length is to model
mis-specification.
'''
import glob
from matplotlib import pyplot as plt
from Bio import Phylo, AlignIO
import pandas as pd
import numpy as np
from treetime import TreeAnc
from treetime.gtr_site_specific import GTR_site_specific
from plot_toy_data import make_subset, add_panel_label


def load_data(files):
   df = pd.concat([pd.read_csv(fname, sep='\t') for fname in files], axis=0)
   return df


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

if __name__ == '__main__':
    import argparse
    parser = argparse.ArgumentParser(description = "", usage="analyze simulated data for site specific GTR reconstruction project")
    parser.add_argument("--files", nargs="+", type=str, help="tables to ")
    parser.add_argument("--plot-mu", type=float, default=0.2, help="mutation rate used for the plot vs epsilong")
    parser.add_argument("--output", type=str, help="file to save figure to")
    args=parser.parse_args()

    df = load_data(args.files)
    eps_vals = list(np.linspace(0.0,1,11))
    mu_vals = list(np.unique(df['m']))

    one_mu = make_subset(df, {"m":args.plot_mu})
    eps_graph = {x:[list() for e in eps_vals] for x in ['pi', 'mu', 'all']}
    mu_graph = {x:[list() for mu in mu_vals] for x in ['pi', 'mu', 'all']}
    # collect data for focal mu value
    for ri, row in one_mu.iterrows():
        for ei,eps in enumerate(eps_vals):
            for x in ['pi', 'mu', 'all']:
                eps_graph[x][ei].append(row[f'eps_{x}_{eps:1.2f}']/row['true_value'])
    # collect data across mu values for epsilon=1
    for ri, row in df.iterrows():
        mi = mu_vals.index(row['m'])
        for x in ['pi', 'mu', 'all']:
            mu_graph[x][mi].append(row[f'eps_{x}_1.00']/row['true_value'])


    # make figure
    fig, axs = plt.subplots(1,2, sharey=True, figsize=(11,5))

    fs=16
    for k, label_str in [('pi', 'frequencies'), ('mu', 'rates'), ('all', 'both')]:
        m = np.mean(eps_graph[k],axis=1)
        std = np.std(eps_graph[k],axis=1)
        axs[0].errorbar(eps_vals, m, std, label=label_str, lw=3)

    axs[0].set_xlabel(r'mixing ratio $\gamma$ true/flat, $\langle \mu \rangle=$' +str(args.plot_mu), fontsize=fs)
    axs[0].set_ylabel('rel. tree length error', fontsize=fs)
    axs[0].tick_params(labelsize=fs*0.8)
    axs[0].set_ylim([0.4,1.1])
    add_panel_label(axs[0],'A',  x_offset=+0.02,y_offset=0.9, fs=fs)

    for k, label_str in [('pi', 'frequencies'), ('mu', 'rates'), ('all', 'both')]:
        m = np.mean(mu_graph[k],axis=1)
        std = np.std(mu_graph[k],axis=1)
        axs[1].errorbar(mu_vals, m, std, label=label_str, lw=3)
    axs[1].set_xlabel(r'evolutionary rate $\langle \mu \rangle$ ($\gamma=1$)', fontsize=fs)
    axs[1].legend(fontsize=fs)
    axs[1].tick_params(labelsize=fs*0.8)
    add_panel_label(axs[1],'B',  x_offset=+0.02, y_offset=0.9, fs=fs)
    plt.tight_layout()
    plt.savefig(args.output)
