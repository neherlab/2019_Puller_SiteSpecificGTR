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
from estimation import get_mutation_count

def save_model(gtr_model, fname):
    np.savez(fname, pi=gtr_model.Pi, mu=gtr_model.mu, W=gtr_model.W, alphabet=gtr_model.alphabet)


def save_mutation_count(T, fname):
    n_ija,T_ia,root = get_mutation_count(T, T.gtr.alphabet)
    np.savez(fname, n_ija=n_ija, T_ia=T_ia, root_sequence=root)


def load_mutation_count(fname):
    d = np.load(fname)
    return d['n_ija'], d['T_ia'], d['root_sequence']


def load_model(fname):
    d = np.load(fname)
    return GTR_site_specific.custom(alphabet=d['alphabet'], mu=d['mu'], pi=d['pi'], W=d['W'])


def simplex(params, out_prefix = None, yule=True, n_model = 5, n_seqgen=5, JC=False, alphabet='nuc_nogap'):
    from Bio import AlignIO
    # generate a model
    T = betatree(params['n'], alpha=2.0)
    T.yule=yule
    T.coalesce()
    if out_prefix:
        Phylo.write(T.BioTree, tree_name(out_prefix, params), 'newick')

    for mi in range(n_model):
        params['model']=mi
        if JC:
            myGTR = GTR_site_specific.random(L=params['L'], alphabet=alphabet,
                                             pi_dirichlet_alpha=0, W_dirichlet_alpha=0)
        else:
            myGTR = GTR_site_specific.random(L=params['L'], alphabet=alphabet)
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


def reconstruct_tree(prefix, params, aa=False):
    aln_file = alignment_name_raw(prefix, params)
    fast_opts = [
        "-ninit", "2",
        "-n",     "2",
        "-me",    "0.05"
    ]
    call = ["iqtree"] + fast_opts +["-nt 1", "-s", aln_file, "-m", 'LG' if aa else 'GTR+G10',
            ">", "iqtree.log"]
    os.system(" ".join(call))
    os.system("mv %s.treefile %s"%(aln_file, reconstructed_tree_name(prefix, params)))
    os.system("rm %s.*"%aln_file)


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description = "", usage="analyze simulated data for site specific GTR reconstruction project")
    parser.add_argument("-m", type=float, help="simulated mutation rate")
    parser.add_argument("-n", type=int, help="number of taxa")
    parser.add_argument("-L", type=int, default=300, help="length of sequence")
    parser.add_argument("--JC", action='store_true', help="simulate JC model")
    parser.add_argument("--aa", action='store_true', help="use amino acid alphabet")
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
        simplex(params, out_prefix=prefix, n_model=2, n_seqgen=2, yule=True, JC=args.JC, alphabet=alphabet)

