import glob
import pandas as pd
from matplotlib import pyplot as plt

if __name__ == '__main__':
    import argparse
    parser = argparse.ArgumentParser(description = "", usage="graph inferred tree lengths")
    parser.add_argument("--treelengths", type=str, help="simulated mutation rate")
    args=parser.parse_args()

    df = pd.read_csv(args.treelengths, sep='\t')

    # make scatter plots comparing true to inferred tree length
    n = 100
    n_branch = 2*n-2

    fig = plt.figure(1)
    df_slice = df.loc[df['n']==n,:]
    true_val = df_slice['trueTree_treelength']/n_branch
    plt.scatter(true_val, df_slice['reconstructed_treelength']/n_branch)
    plt.scatter(true_val, df_slice['inferred_c0.1_treelength']/n_branch)
    plt.scatter(true_val, df_slice['inferred_c1.0_treelength']/n_branch)
    plt.scatter(true_val, df_slice['trueRates_c0.1_treelength']/n_branch)
    plt.scatter(true_val, df_slice['trueModel_treelength']/n_branch)
    plt.plot([0,.3], [0,.3])

    fig = plt.figure(2)
    df_slice = df.loc[df['n']==n,:]
    true_val = df_slice['trueTree_depth']
    plt.scatter(true_val, df_slice['reconstructed_depth'])
    plt.scatter(true_val, df_slice['inferred_c0.1_depth'])
    plt.scatter(true_val, df_slice['inferred_c1.0_depth'])
    plt.scatter(true_val, df_slice['trueRates_c0.1_depth'])
    plt.scatter(true_val, df_slice['trueModel_depth'])
    plt.plot([0,.3], [0,.3])


    # for dset in length:
    #     length_a = np.array(length[dset])
    #     depth_a = np.array(depth[dset])
    #     terminal_a = np.array(terminal_bl[dset])
    #     if dset=='ML':
    #         label =  'FastTree JTT+CAT20' if aa else 'IQ-tree GTR+R10'
    #     elif dset=='true':
    #         label='True model'
    #     else:
    #         label=f"Inferred model {'(true rates)' if dset[1]=='true-rates' else ''}, pc={dset[0]:1.1f}"
    #     plt.figure(1)
    #     plt.scatter(length_a[:,0], length_a[:,1], label=label)
    #     plt.figure(2)
    #     plt.scatter(depth_a[:,0], depth_a[:,1], label=label)
    #     plt.figure(3)
    #     plt.scatter(terminal_a[:,0], terminal_a[:,1], label=label)

    # plt.figure(1)
    # plt.plot([0,0.25], [0,0.25])
    # plt.ylabel('inferred average branch length')
    # plt.xlabel('true average branch length')
    # plt.legend()
    # plt.savefig("figures/"+prefix+"_length.pdf")

    # plt.figure(2)
    # plt.plot([0,3.5], [0,3.5])
    # plt.ylabel('inferred average root-to-tip distance')
    # plt.xlabel('true average root-to-tip distance')
    # plt.legend()
    # plt.savefig("figures/"+prefix+"_depth.pdf")

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
