"""
By default, the reconstructed trees are rooted at midpoint. 
This script was only used in the early days when this was
not yet done.
"""
import glob
from Bio import Phylo
from generate_toy_data import *
from filenames import *

if __name__ == '__main__':
    import argparse
    parser = argparse.ArgumentParser(description = "", usage="reroot all reconstructed trees")
    parser.add_argument("--prefix", type=str, required=True, help="folder to retrieve data from")
    args=parser.parse_args()

    prefix = args.prefix
    mask = "/L*_n*_m*_*fasta.gz"
    files = glob.glob(prefix+mask)

    for fname in files:
        print(fname)
        params = parse_alignment_name(fname)
        iq_tree = Phylo.read(reconstructed_tree_name(prefix, params), 'newick')
        iq_tree.root_at_midpoint()
        Phylo.write(iq_tree, reconstructed_tree_name(prefix, params), 'newick')
