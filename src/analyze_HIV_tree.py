import os, gzip, glob, pickle
from Bio import AlignIO
import numpy as np
from scipy.stats import pearsonr, spearmanr
from collections import defaultdict
import matplotlib.pyplot as plt
import seaborn as sns

from treetime.treeanc import TreeAnc
from treetime.gtr_site_specific import GTR_site_specific
from treetime.gtr import GTR
from treetime.seq_utils import seq2prof, profile_maps, alphabets

from estimation import *
from generate_toy_data import save_model, load_model

def convert_vals(x):
    try:
        return float(x)
    except:
        if '<' in x:
            return 0.0005
        else:
            return 0.2

def load_fitness_landscape(gene, aa=False, subtype='B'):
    import pandas as pd
    if aa:
        d = pd.read_csv('HIV_fitness_landscape/data/fitness_pooled_aa/aa_%s_fitness_costs_uncensored_st_%s.tsv'%(gene, subtype),
                    sep='\t', skiprows=1)
    else:
        d = pd.read_csv('HIV_fitness_landscape/data/fitness_pooled/nuc_%s_selection_coeffcients_uncensored_%s.tsv'%(gene, subtype),
                    sep='\t', skiprows=1)
    for col in d.columns[3:6]:
        d.loc[:,col] = d.loc[:,col].apply(convert_vals)

    return d


def calculate_in_out_ratio(gtr, major_allele):
    '''
    e^f_a = sum_j pi_ia W_ij/sum W_ij pi_ja
    where i is the major allele at position i
    '''
    major_index = np.zeros_like(major_allele, dtype=int)
    for ni, n in enumerate(gtr.alphabet):
        major_index[major_allele==n]=ni
    #major_index = np.argmax(gtr.Pi, axis=0)

    r_ind = np.arange(gtr.Pi.shape[-1])
    return np.einsum('a,a,aj->a', gtr.mu, gtr.Pi[major_index, r_ind], gtr.W[major_index])/\
           np.einsum('a,aj,ja->a', gtr.mu, gtr.W[major_index], gtr.Pi)


def diversity(p):
    return 1.0-np.sum(p**2, axis=0)

if __name__ == '__main__':
    import argparse
    parser = argparse.ArgumentParser(description = "", usage="analyze simulated data for site specific GTR reconstruction project")
    hiv_out = 'HIV_inferences/'
    parser.add_argument("--prefix", type=str, help="data set")
    parser.add_argument("--gene", type=str, default='pol')
    parser.add_argument("--subtype", type=str, default='B')
    parser.add_argument("--aa", action='store_true', help="use amino acid alphabet")
    parser.add_argument("--pc", type=float, default=1.0, help="pseudocount for reconstruction")
    parser.add_argument("--redo", action='store_true', default=False, help="re-estimate")
    args=parser.parse_args()

    alphabet='aa' if args.aa else 'nuc_nogap'
    pc = args.pc
    fs = 16

    aln = args.prefix +  ('_aa-aligned.fasta' if args.aa else '_aligned.fasta')
    tree = args.prefix + '_tree.nwk'
    model_name = hiv_out + os.path.basename(aln)[:-5]+'%1.2f_inferred_model%s.npz'%(args.pc, '_aa' if args.aa else '')
    aln_freq_name = hiv_out + os.path.basename(aln[:-5])+'alignment_frequencies%s'%( '_aa' if args.aa else '')

    gtr=None
    aln_freq=None
    if not args.redo:
        try:
            aln_freq = np.loadtxt(aln_freq_name)
            gtr = load_model(model_name)
        except:
            pass
    if gtr is None:
        tt = TreeAnc(tree=tree, aln = aln, compress=False, alphabet=alphabet, verbose=3)
        # tt.optimize_tree(branch_length_mode='marginal', max_iter=10,
        #                  infer_gtr=True, site_specific_gtr=True, pc=pc)

        tt.infer_gtr_iterative(max_iter=10, site_specific=True, pc=pc)

        gtr = tt.gtr
        save_model(gtr, model_name)

    if aln_freq is None:
        aln_freq = p_from_aln(AlignIO.read(aln, 'fasta'), alphabet=alphabet)
        np.savetxt(aln_freq_name, aln_freq)

    fabio_fitness = load_fitness_landscape(args.gene, subtype=args.subtype, aa=args.aa)
    valid = np.isfinite(fabio_fitness['median'])

    fitness_estimates_gtr = calculate_in_out_ratio(gtr, fabio_fitness.loc[:,"consensus"])

    max_bin = int(50*len(fitness_estimates_gtr)/3000)


    opa = 0.3
    msize = 8
    ccoef=pearsonr
    #ccoef=spearmanr
    div_rec = diversity(gtr.Pi)
    cc_rec = ccoef(np.log(fabio_fitness['median'])[valid], np.log(1e-10+div_rec)[valid])
    plt.figure()
    # g = plt.hexbin(fabio_fitness['median'], div_rec, xscale='log', yscale='log')
#    g.ax_joint.set_yscale('log')
#    g.ax_joint.set_xscale('log')
#    g.ax_joint.set_xlim([3e-5, 3e0])
#    g.ax_joint.set_ylim([3e-5, 1.5])

    plt.scatter(fabio_fitness['median'], div_rec, label='reconstructed, r2=%1.2f'%cc_rec[0]**2, alpha=opa, s=msize)
    plt.legend()
    plt.xlim([3e-5, 3e0])
    plt.ylim([3e-5, 1.5])
    plt.yscale('log')
    plt.xscale('log')

    div_aln = diversity(aln_freq)
    cc_aln = ccoef(np.log(fabio_fitness['median'][valid]), np.log(1e-10+div_aln)[valid])
    plt.figure()
    plt.scatter(fabio_fitness['median'], div_aln, label='alignment, r2=%1.2f'%cc_aln[0]**2, alpha=opa, s=msize)
    plt.legend()
    plt.xlim([3e-5, 3e-0])
    plt.ylim([3e-5, 1.5])
    plt.yscale('log')
    plt.xscale('log')


    plt.figure()
    cc_fit = ccoef(np.log(fabio_fitness['median'][valid]), np.log(1e-10+fitness_estimates_gtr)[valid])
    plt.scatter(fabio_fitness['median'], fitness_estimates_gtr, label=r'$r^2=$'+'%1.2f'%cc_fit[0]**2, alpha=opa, s=msize)
    plt.legend(fontsize=fs)
    # plt.hexbin(fabio_fitness['median'], fitness_estimates_gtr, xscale='log', yscale='log', gridsize=16, mincnt=1, vmin=1, vmax=max_bin)
    plt.xlim([3e-5, 3e0])
    plt.ylim([3e-2, 3e3])
    # plt.text(5e-5, 1e3, r'$r^2=$'+'%1.2f'%cc_fit[0]**2, fontsize=fs)
    plt.title('B: amino acids' if args.aa else 'A: nucleotides',fontsize=fs*1.2)
    plt.xlabel('intra-host fitness cost estimates', fontsize=fs)
    plt.ylabel(r'GTR based fitness estimates $\Gamma_{in}/\Gamma_{out}$', fontsize=fs)
    plt.xscale('log')
    plt.yscale('log')
    plt.tick_params(labelsize=0.8*fs)
    plt.tight_layout()
    plt.savefig('figures/'+args.prefix.split('/')[-1]+'_fitness_pc_%1.3f'%args.pc+ ('_aa' if args.aa else '')+'.pdf')

    plt.figure()
    plt.scatter(fabio_fitness['median'], gtr.mu,  label='fitness vs rate')
    plt.legend()
    plt.xscale('log')
    plt.yscale('log')

    plt.figure()
    cc_div = ccoef(np.log(1e-10+div_aln), np.log(1e-10+fitness_estimates_gtr))
    plt.scatter(div_aln, fitness_estimates_gtr, label=r'$r^2=$'+'%1.2f'%cc_div[0]**2, alpha=opa, s=msize)
    plt.legend(fontsize=fs)
    # plt.hexbin(div_aln+1e-4, fitness_estimates_gtr, xscale='log', yscale='log', gridsize=16, mincnt=1, vmin=1, vmax=max_bin)
    # plt.text(3e-4, 1e-1, r'$r^2=$'+'%1.2f'%cc_div[0]**2, fontsize=fs)
    plt.title('amino acids' if args.aa else 'nucleotides',fontsize=fs)
    plt.xlim([2e-4, 3e0])
    plt.ylim([3e-2, 9e3])
    plt.xlabel('diversity in alignment', fontsize=fs)
    plt.ylabel(r'GTR based fitness estimates $\Gamma_{in}/\Gamma_{out}$', fontsize=fs)
    plt.xscale('log')
    plt.yscale('log')
    plt.tick_params(labelsize=0.8*fs)
    plt.tight_layout()
    plt.savefig('figures/'+args.prefix.split('/')[-1]+'_alignment_div_pc_%1.3f'%args.pc+ ('_aa' if args.aa else '')+'.pdf')
