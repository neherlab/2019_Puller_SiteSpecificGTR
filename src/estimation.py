import gzip
import numpy as np
from Bio import Phylo, AlignIO

from treetime.treeanc import TreeAnc
from treetime.gtr_site_specific import GTR_site_specific
from treetime.gtr import GTR
from treetime.seq_utils import seq2prof, profile_maps, alphabets
from filenames import *

def p_from_aln(in_prefix, params, alphabet='nuc_nogap'):
    with gzip.open(alignment_name(in_prefix, params), 'rt') as fh:
        aln = AlignIO.read(fh, 'fasta')

    alpha = alphabets[alphabet]
    aln_array = np.array(aln)
    af = []
    for a in alpha:
        af.append(np.mean(aln_array==a, axis=0))
    return np.array(af)


def reconstruct_counts(in_prefix, params, gtr='JC69', alphabet='nuc_nogap',
                       marginal=False, reconstructed_tree=False):
    with gzip.open(alignment_name(in_prefix, params), 'rt') as fh:
        aln = AlignIO.read(fh, 'fasta')
    tree_fname =  reconstructed_tree_name(in_prefix, params) if reconstructed_tree else tree_name(in_prefix, params)
    myTree = TreeAnc(gtr=gtr, alphabet=alphabet,
                     tree=tree_fname, aln=aln,
                     reduce_alignment=False, verbose = 0)

    if type(gtr)==str:
        myTree.infer_ancestral_sequences(marginal=True, infer_gtr=True, normalized_rate=False)
    else:
        myTree.infer_ancestral_sequences(marginal=True)

    mutation_counts = get_mutation_count(myTree, alphabet=myTree.gtr.alphabet, marginal=marginal)
    return mutation_counts, myTree.sequence_LH(), myTree


def estimate_GTR(mutation_counts, pc=0.1, single_site=False, bl=None):
    n_ija, T_ia, root_sequence = mutation_counts[:3]
    root_prof = seq2prof(root_sequence, profile_maps['nuc_nogap']).T

    if single_site:
        inferred_model = GTR.infer(n_ija.sum(axis=-1), T_ia.sum(axis=-1),
                          root_state=root_prof.sum(axis=-1), alphabet='nuc_nogap')
    else:
        inferred_model = GTR_site_specific.infer(n_ija, T_ia, pc=pc,
                          root_state=root_prof, alphabet='nuc_nogap', bl=bl)
    return inferred_model


def get_mutation_count(tree, alphabet, marginal=False):
    if marginal:
        return get_marginal_mutation_count(tree, alphabet)
    else:
        return get_ML_mutation_count(tree, alphabet)


def get_ML_mutation_count(tree, alphabet):
    alphabet_to_index = {a:ai for ai,a in enumerate(alphabet)}
    L = tree.seq_len
    q = len(alphabet)
    positions = np.arange(L)
    n_ija = np.zeros((q,q,L), dtype=int)
    T_ia = np.zeros((q,L),dtype=float)
    for n in tree.tree.get_nonterminals():
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


    return n_ija, T_ia, tree.tree.root.sequence


def get_marginal_mutation_count(tree, alphabet):
    L = tree.seq_len
    q = len(alphabet)
    n_ija = np.zeros((q,q,L), dtype=float)
    n_ija_eff = np.zeros((q,q,L), dtype=float)
    T_ia = np.zeros((q,L),dtype=float)
    for n in tree.tree.get_nonterminals():
        for c in n:
            mut_stack = np.transpose(tree.get_branch_mutation_matrix(c, full_sequence=True), (2,1,0))
            neff = np.transpose(tree.get_effective_mutation_matrix(c, full_sequence=True), (2,1,0))
            T_ia += 0.5*c.branch_length * mut_stack.sum(axis=0)
            T_ia += 0.5*c.branch_length * mut_stack.sum(axis=1)

            n_ija_eff += neff
            n_ija += mut_stack

    return n_ija, T_ia, tree.tree.root.sequence, n_ija_eff


def get_average_transition_matrix(tree, model):

    eQt = np.zeros_like(model.v)
    tQeQt = np.zeros_like(model.v)
    for n in tree.find_clades():
        if n==tree.root:
            continue
        else:
            tmp = model.expQt(n.branch_length)
            eQt+=tmp
            tQeQt += n.branch_length*model.Q*tmp

    return eQt, tQeQt
