import os, gzip, argparse, sys
import numpy as np
from Bio import Phylo, AlignIO

from treetime.treeanc import TreeAnc
from treetime.gtr_site_specific import GTR_site_specific
from treetime.gtr import GTR
from treetime.seqgen import SeqGen
from treetime.seq_utils import seq2prof, profile_maps

from betatree import betatree

def parse_alignment_name(fname):
    base = os.path.basename(fname)[:-9]
    params = {}
    for x in base.split('_'):
        try:
            params[x[0]]=int(x[1:])
        except:
            try:
                params[x[0]]=float(x[1:])
            except:
                try:
                    params[x[:-2]]=int(x[-2:])
                except:
                    pass

    return params


def model_name(prefix, params):
    return prefix+'/L{L}_n{n}_m{m}_tree{tree:02d}_model{model:02d}.npz'.format(**params)

def mutation_count_name(prefix, params):
    return prefix+'/L{L}_n{n}_m{m}_tree{tree:02d}_model{model:02d}_seqgen{seqgen:02}_mutations.npz'.format(**params)

def alignment_name(prefix, params):
    return prefix+'/L{L}_n{n}_m{m}_tree{tree:02d}_model{model:02d}_seqgen{seqgen:02}.fasta.gz'.format(**params)

def tree_name(prefix, params):
    return prefix+'/L{L}_n{n}_m{m}_tree{tree:02d}.nwk'.format(**params)

def save_model(gtr_model, fname):
    np.savez(fname, pi=gtr_model.Pi, mu=gtr_model.mu, W=gtr_model.W, alphabet=gtr_model.alphabet)

def save_mutation_count(T, fname):
    n_ija,T_ia,root = get_mutation_count(T.tree, T.gtr.alphabet)
    np.savez(fname, n_ija=n_ija, T_ia=T_ia, root_sequence=root)

def load_mutation_count(fname):
    d = np.load(fname)
    return d['n_ija'], d['T_ia'], d['root_sequence']

def load_model(fname):
    d = np.load(fname)
    return GTR_site_specific.custom(alphabet=d['alphabet'], mu=d['mu'], pi=d['pi'], W=d['W'])

def get_mutation_count(tree, alphabet):
    alphabet_to_index = {a:ai for ai,a in enumerate(alphabet)}
    L = len(tree.root.sequence)
    q=len(alphabet)
    positions = np.arange(L)
    n_ija = np.zeros((q,q,L), dtype=int)
    T_ia = np.zeros((q,L),dtype=float)
    for n in tree.get_nonterminals():
        parent_profile = np.zeros(L, dtype=int)
        for ai,a in enumerate(alphabet):
            parent_profile[n.sequence==a] = ai

        for c in n:
            child_profile = np.zeros(L, dtype=int)
            for ai,a in enumerate(alphabet):
                child_profile[c.sequence==a] = ai

            T_ia[parent_profile,positions] += 0.5*c.branch_length
            T_ia[child_profile,positions] += 0.5*c.branch_length

            n_ija[child_profile, parent_profile, positions] += (1-(parent_profile==child_profile))

    return n_ija, T_ia, tree.root.sequence


def simplex(params, out_prefix = None, yule=True, n_model = 5, n_seqgen=5):
    from Bio import AlignIO
    # generate a model
    T = betatree.betatree(params['n'], alpha=2.0)
    T.yule=yule
    T.coalesce()
    if out_prefix:
        Phylo.write(T.BioTree, tree_name(out_prefix, params), 'newick')

    for mi in range(n_model):
        params['model']=mi
        myGTR = GTR_site_specific.random(L=params['L'], alphabet='nuc_nogap')
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
                with gzip.open(alignment_name(out_prefix, params), 'wt') as fh:
                    AlignIO.write(mySeq.get_aln(), fh, 'fasta')


if __name__ == '__main__':
    L=100
    n=3000
    for mu in [0.01, 0.02, 0.05, 0.1, 0.15, 0.25, 0.5]:
        for ti in range(3):
            prefix = 'simulated_data/'
            simplex({'L':L, 'n':n, 'm':mu, 'tree':ti}, out_prefix = prefix, n_model=3, n_seqgen=3)
