import gzip
import numpy as np
from Bio import Phylo, AlignIO

from treetime.treeanc import TreeAnc
from treetime.gtr_site_specific import GTR_site_specific
from treetime.gtr import GTR
from treetime.seq_utils import seq2prof, profile_maps, alphabets
from filenames import *

def p_from_aln(aln, alphabet='nuc_nogap'):
    """computes frequencies of characters in alignment columns

    Parameters
    ----------
    aln : alignment
        Alignment
    alphabet : str, optional
        nucleotide or amino acid

    Returns
    -------
    np.array
        frequency matrix
    """
    alpha = alphabets[alphabet]
    aln_array = np.array(aln)
    af = []
    for a in alpha:
        af.append(np.mean(aln_array==a.decode(), axis=0))
    return np.array(af)


def estimate_GTR(mutation_counts, pc=0.1, single_site=False, tt=None, alphabet='nuc_nogap'):
    n_ija, T_ia, root_sequence = mutation_counts[:3]
    root_prof = seq2prof(root_sequence, profile_maps[alphabet]).T


    if single_site:
        inferred_model = GTR.infer(n_ija.sum(axis=-1), T_ia.sum(axis=-1), pc=pc,
                                   root_state=root_prof.sum(axis=-1), alphabet=alphabet)
    else:
        inferred_model = GTR_site_specific.infer(n_ija, T_ia, pc=pc,
                                                 root_state=root_prof, alphabet=alphabet)

    if tt:
        (m_denom, m_num), (t_denom, t_num) = exact_mu_t_update(tt)
        inferred_model.mu *= (m_num+pc)/(m_denom+pc)

    return inferred_model


def exact_mu_t_update(tree):
    L = tree.seq_len
    q = len(tree.gtr.alphabet)
    mu_diag = np.zeros((L), dtype=float)
    mu_offdiag = np.zeros((L), dtype=float)

    t_diag = []
    t_offdiag = []
    for n in tree.tree.get_nonterminals():
        for c in n:
            (m1,m2), mu, t = tree.get_rate_length_updates(c)
            mu_diag += t*m1
            mu_offdiag += t*m2

            t_diag.append(np.sum(mu*m1))
            t_offdiag.append(np.sum(mu*m2))

    return (mu_diag, mu_offdiag), (np.array(t_diag), np.array(t_offdiag))
