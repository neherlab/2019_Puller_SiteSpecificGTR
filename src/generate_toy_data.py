import os, gzip, argparse, sys
import numpy as np
from Bio import Phylo, AlignIO

from treetime.treeanc import TreeAnc
from treetime.gtr_site_specific import GTR_site_specific
from treetime.gtr import GTR
from treetime.seqgen import SeqGen
from treetime.seq_utils import seq2prof, profile_maps

from betatree import betatree

from filenames import *


def save_model(gtr_model, fname):
    np.savez(fname, pi=gtr_model.Pi, mu=gtr_model.mu, W=gtr_model.W, alphabet=gtr_model.alphabet)


def save_mutation_count(T, fname):
    n_ija,T_ia,root = get_ancestral_mutation_count(T, T.gtr.alphabet)
    np.savez(fname, n_ija=n_ija, T_ia=T_ia, root_sequence=root)


def load_mutation_count(fname):
    d = np.load(fname)
    return d['n_ija'], d['T_ia'], d['root_sequence']


def load_model(fname, flatten_p=0.0, flatten_mu=0.0, flatten_W=0.0):
    d = np.load(fname)
    return GTR_site_specific.custom(alphabet=d['alphabet'],
                                    mu=d['mu']*(1-flatten_mu) + flatten_mu*np.mean(d['mu']),
                                    pi=d['pi']*(1-flatten_p)+flatten_p/d['pi'].shape[0],
                                    W=d['W']*(1-flatten_W) + 2*flatten_W/d['W'].shape[0]/(d['W'].shape[0]-1))

def simplex(params, out_prefix = None, yule=True, n_model = 5, n_seqgen=5, JC=False, alphabet='nuc_nogap', alpha=1.0, rate_alpha=1.5):
    from Bio import AlignIO
    # generate a model
    T = betatree(params['n'], alpha=2.0)
    T.yule=yule
    T.coalesce()
    # ladderize the tree and name internal nodes via loading into TreeAnc
    T.BioTree.ladderize()
    tt = TreeAnc(tree=T.BioTree)
    if out_prefix:
        Phylo.write(tt.tree, tree_name(out_prefix, params), 'newick')

    for mi in range(n_model):
        params['model']=mi
        if JC:
            myGTR = GTR_site_specific.random(L=params['L'], alphabet=alphabet,
                                             pi_dirichlet_alpha=0, W_dirichlet_alpha=0, mu_gamma_alpha=rate_alpha)
        else:
            myGTR = GTR_site_specific.random(L=params['L'], alphabet=alphabet, pi_dirichlet_alpha = alpha, mu_gamma_alpha=rate_alpha)

        myGTR.mu*=params['m']

        if out_prefix:
            save_model(myGTR, model_name(out_prefix, params))

        for si in range(n_seqgen):
            params['seqgen']=si
            # generate sequences
            mySeq = SeqGen(gtr=myGTR, tree=T.BioTree)
            mySeq.evolve()

            if out_prefix:
                save_mutation_count(mySeq, mutation_count_name(out_prefix, params))
                with open(alignment_name_raw(out_prefix, params), 'wt') as fh:
                    AlignIO.write(mySeq.get_aln(), fh, 'fasta')
                reconstruct_tree(out_prefix, params, aa='aa' in alphabet)
                os.system('gzip '+alignment_name_raw(out_prefix, params))


def reconstruct_tree(prefix, params, aa=False):
    aln_file = alignment_name_raw(prefix, params)
    out_tree = reconstructed_tree_name(prefix, params)
    if aa:
        call = ["fasttree", aln_file, ">", out_tree]
        os.system(" ".join(call))
    else:
        fast_opts = [
            "-ninit", "2",
            "-n",     "2",
            "-me",    "0.05"
        ]
        call = ["iqtree"] + fast_opts +["-nt 1", "-s", aln_file, "-m", 'GTR+R10',
                ">", "iqtree.log"]
        os.system(" ".join(call))
        os.system("mv %s.treefile %s"%(aln_file, out_tree))
        os.system("rm %s.*"%aln_file)

    rec_tree = Phylo.read(out_tree, 'newick')
    rec_tree.root_at_midpoint()
    rec_tree.ladderize()
    Phylo.write(rec_tree, out_tree, 'newick')

def get_ancestral_mutation_count(tree, alphabet):
    alphabet_to_index = {a:ai for ai,a in enumerate(alphabet)}
    L = tree.seq_len
    q = len(alphabet)
    positions = np.arange(L)
    n_ija = np.zeros((q,q,L), dtype=int)
    T_ia = np.zeros((q,L),dtype=float)
    for n in tree.tree.get_nonterminals():
        parent_profile = np.zeros(L, dtype=int)
        for ai,a in enumerate(alphabet):
            parent_profile[n.ancestral_sequence==a] = ai

        for c in n:
            child_profile = np.zeros(L, dtype=int)
            for ai,a in enumerate(alphabet):
                child_profile[c.ancestral_sequence==a] = ai

            T_ia[parent_profile,positions] += 0.5*c.branch_length
            T_ia[child_profile,positions] += 0.5*c.branch_length

            n_ija[child_profile, parent_profile, positions] += (1-(parent_profile==child_profile))


    return n_ija, T_ia, tree.tree.root.ancestral_sequence



if __name__ == '__main__':
    parser = argparse.ArgumentParser(description = "", usage="analyze simulated data for site specific GTR reconstruction project")
    parser.add_argument("-m", type=float, help="simulated mutation rate")
    parser.add_argument("-n", type=int, help="number of taxa")
    parser.add_argument("-L", type=int, default=300, help="length of sequence")
    parser.add_argument("--JC", action='store_true', help="simulate JC model")
    parser.add_argument("--aa", action='store_true', help="use amino acid alphabet")
    parser.add_argument("--alpha", default=1.0, type=float,  help="parameter of the dirichlet distribution for preferences")
    parser.add_argument("--rate-alpha", default=1.5, type=float,  help="parameter of the gamma distribution for rates")
    parser.add_argument("--prefix", type=str, help="folder to save data")
    args=parser.parse_args()

    L=args.L
    alphabet='aa_nogap' if args.aa else 'nuc_nogap'
    prefix = args.prefix
    if not os.path.isdir(prefix):
        os.mkdir(prefix)

    n = args.n
    mu = args.m
    for ti in range(2):
        params = {'L':L, 'n':n, 'm':mu, 'tree':ti}
        simplex(params, out_prefix=prefix, n_model=2, n_seqgen=2, yule=True, JC=args.JC, alphabet=alphabet, alpha=args.alpha, rate_alpha=args.rate_alpha)

