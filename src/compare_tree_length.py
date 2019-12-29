import glob
import pandas as pd
from matplotlib import pyplot as plt
from plot_toy_data import add_panel_label

if __name__ == '__main__':
    import argparse
    parser = argparse.ArgumentParser(description = "", usage="graph inferred tree lengths")
    parser.add_argument("--treelengths", type=str, help="simulated mutation rate")
    parser.add_argument("--nval", type=int, help="data set size")
    parser.add_argument("--aa", action='store_true', help="inference from amino acid sequences")
    args=parser.parse_args()

    df = pd.read_csv(args.treelengths, sep='\t')

    # make scatter plots comparing true to inferred tree length
    n = args.nval
    n_branch = 2*n-2

    fs = 16
    fig = plt.figure()
    df_slice = df.loc[df['n']==n,:]
    true_val = df_slice['trueTree_treelength']/n_branch
    plt.scatter(true_val, df_slice['trueModel_treelength']/n_branch, label='true model')
    plt.scatter(true_val, df_slice['reconstructed_treelength']/n_branch, label="FastTree" if args.aa else "IQ-tree GTR+10")
    plt.scatter(true_val, df_slice['inferred_c0.1_treelength']/n_branch, label='c=0.1')
    plt.scatter(true_val, df_slice['inferred_c1.0_treelength']/n_branch, label='c=1.0')
    plt.scatter(true_val, df_slice['trueRates_c0.1_treelength']/n_branch, label='c=0.1 (true rates)')
#    plt.scatter(true_val, df_slice['trueRates_c1.0_treelength']/n_branch, label='c=1.0 (true rates)')
    plt.plot([0,.3], [0,.3])
    plt.xlim([0,.3])
    plt.ylim([0,.3])
    plt.ylabel('Inferred avg branch length', fontsize=fs)
    plt.xlabel('True avg branch length', fontsize=fs)
    plt.tick_params(labelsize=0.8*fs)
    add_panel_label(plt.gca(),'A',  x_offset=-0.16, fs=fs)
    plt.legend(fontsize=0.8*fs)
    plt.tight_layout()
    plt.savefig(f"figures/{'aa' if args.aa else 'nuc'}_length_n{args.nval}.pdf")


    fig = plt.figure()
    df_slice = df.loc[df['n']==n,:]
    true_val = df_slice['trueTree_depth']
    plt.scatter(true_val, df_slice['trueModel_depth'], label='true model')
    plt.scatter(true_val, df_slice['reconstructed_depth'], label="FastTree" if args.aa else "IQ-tree GTR+10")
    plt.scatter(true_val, df_slice['inferred_c0.1_depth'], label='c=0.1')
    plt.scatter(true_val, df_slice['inferred_c1.0_depth'], label='c=1.0')
    plt.scatter(true_val, df_slice['trueRates_c0.1_depth'], label='c=0.1 (true rates)')
#    plt.scatter(true_val, df_slice['trueRates_c1.0_depth'], label='c=1.0 (true rates)')
    plt.ylabel('Inferred avg root-to-tip', fontsize=fs)
    plt.xlabel('True avg root-to-tip distance', fontsize=fs)
    plt.tick_params(labelsize=0.8*fs)
    plt.plot([0,3.4], [0,3.4])
    plt.xlim([0,3.4])
    plt.ylim([0,3.4])
    add_panel_label(plt.gca(),'B',  x_offset=-0.12, fs=fs)
    plt.legend(fontsize=0.8*fs)
    plt.tight_layout()
    plt.savefig(f"figures/{'aa' if args.aa else 'nuc'}_depth_n{args.nval}.pdf")


    fig = plt.figure()
    df_slice = df.loc[df['n']==n,:]
    true_val = df_slice['trueTree_terminals']
    plt.scatter(true_val, df_slice['trueModel_terminals'], label='true model')
    plt.scatter(true_val, df_slice['reconstructed_terminals'], label="FastTree" if args.aa else "IQ-tree GTR+10")
    plt.scatter(true_val, df_slice['inferred_c0.1_terminals'], label='c=0.1')
    plt.scatter(true_val, df_slice['inferred_c1.0_terminals'], label='c=1.0')
#    plt.scatter(true_val, df_slice['trueRates_c0.1_terminals'], label='c=0.1 (true rates)')
    plt.ylabel('Inferred avg terminal length', fontsize=fs)
    plt.xlabel('True avg terminal length', fontsize=fs)
    plt.tick_params(labelsize=0.8*fs)
    plt.plot([0,.3], [0,.3])
    plt.xlim([0,.3])
    plt.ylim([0,.3])
    add_panel_label(plt.gca(),'C',  x_offset=-0.16, fs=fs)
    plt.legend(fontsize=0.8*fs)
    plt.tight_layout()
    plt.savefig(f"figures/{'aa' if args.aa else 'nuc'}_terminal_n{args.nval}.pdf")

    # plt.figure(3)
    # plt.plot([0,0.25], [0,0.25])
    # plt.ylabel('inferred average terminal branch length')
    # plt.xlabel('true average terminal branch length')
    # plt.legend()
    # plt.savefig("figures/"+prefix+"_terminal.pdf")
    # # for dset in all_bl:
    # #     plt.figure()
    # #     plt.title(str(dset))
    # #     for x in all_bl[dset]:
    # #         plt.scatter(x[:,0], x[:,1])
    # #     plt.plot([0,1],[0,1])
    # #     plt.ylim([0,1.5])
    # #     plt.xlim([0,1.5])
