import glob
from collections import defaultdict
import numpy as np
from matplotlib import pyplot as plt
import pandas as pd
from Bio import Phylo, AlignIO
from filenames import *

def analyze_tree(tree, label_stub):
    res = {}
    tree.root_at_midpoint()
    res[label_stub+'_treelength'] = tree.total_branch_length()
    res[label_stub+'_depth'] = np.mean([x for c,x in tree.depths().items() if c.is_terminal()])
    res[label_stub+'_terminals'] = np.mean([x.branch_length for x in tree.get_terminals()])
    return res

if __name__ == '__main__':
    import argparse
    parser = argparse.ArgumentParser(description = "", usage="collect statistics of re-optimized trees")
    parser.add_argument("--reoptimizedTrees", type=str, help="directory with reoptimized trees")
    parser.add_argument("--trueTrees", type=str, help="directory with true trees")
    parser.add_argument("--nval", type=int, help="tree size to glob")
    parser.add_argument("--true-rates", action='store_true', help="include estimates with true rates")
    parser.add_argument("--output", type=str, help="file to save data")
    args=parser.parse_args()

    res = defaultdict(dict)
    if args.nval:
        reoptimizedTrees = glob.glob(args.reoptimizedTrees+f'/*n{args.nval}_*reoptimized.nwk')
        recTrees = glob.glob(args.trueTrees+f'/*n{args.nval}_*reconstructed.nwk')
    else:
        reoptimizedTrees = glob.glob(args.reoptimizedTrees+'/*reoptimized.nwk')
        recTrees = glob.glob(args.trueTrees+'/*reconstructed.nwk')

    for treename in reoptimizedTrees:
        dset = parse_alignment_name(treename)
        pc = 'None' if "trueModel" in treename else dset['c']
        if 'c' in dset:
            dset.pop('c')
        tree_label = tuple(dset.items())
        tree = Phylo.read(treename, 'newick')

        if "trueModel" in treename:
            label_stub = "trueModel"
        elif "trueRates" in treename:
            label_stub = f"trueRates_c{pc:1.2}"
        else:
            label_stub = f"inferred_c{pc:1.2}"

        dset.update(analyze_tree(tree, label_stub))
        res[tree_label].update(dset)


    for treename in recTrees:
        dset = parse_alignment_name(treename)
        tree_label = tuple(dset.items())
        tree = Phylo.read(treename, 'newick')
        label_stub = "reconstructed"
        dset.update(analyze_tree(tree, label_stub))
        res[tree_label].update(dset)

        tree = Phylo.read(tree_name(os.path.dirname(treename), dset), 'newick')
        for n in tree.find_clades():
            if n.branch_length:
                n.branch_length *= dset['m']
        label_stub = "trueTree"
        dset.update(analyze_tree(tree, label_stub))
        res[tree_label].update(dset)


    out_dir = os.path.dirname(args.output) or '.'
    if not os.path.isdir(out_dir):
        os.mkdir(out_dir)

    df = pd.DataFrame(res.values())
    df.to_csv(args.output, sep='\t')
