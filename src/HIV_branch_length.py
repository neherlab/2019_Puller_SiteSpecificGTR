import os, gzip, glob, pickle
from matplotlib import pyplot as plt

from treetime import TreeAnc, GTR_site_specific, GTR
from treetime.seq_utils import seq2prof, profile_maps, alphabets, prof2seq

from estimation import *
from generate_toy_data import save_model, load_model

from plot_toy_data import add_panel_label


if __name__ == '__main__':
    import argparse
    parser = argparse.ArgumentParser(description = "", usage="analyze simulated data for site specific GTR reconstruction project")
    hiv_out = 'HIV_inferences/'
    parser.add_argument("--prefix", type=str, help="data set")
    parser.add_argument("--gene", type=str, default='pol')
    parser.add_argument("--subtype", type=str, default='B')
    parser.add_argument("--aa", action='store_true', help="use amino acid alphabet")
    parser.add_argument("--pc", type=float, default=1.0, help="pseudocount for reconstruction")
    args=parser.parse_args()

    # load model computed by script "analyze_HIV_tree.py"
    alphabet='aa' if args.aa else 'nuc_nogap'
    pc = args.pc
    fs = 16
    dset_name=os.path.basename(args.prefix)
    model_name = hiv_out + dset_name+'_%saligned.%1.2f_inferred_model%s.npz'%('aa-' if args.aa else '', args.pc, '_aa' if args.aa else '')
    gtr = load_model(model_name)

    S = -np.sum(np.log(gtr.Pi)*gtr.Pi, axis=0)
    print('Entropy:', S.mean(), S.std())

    # make a replica of the model with frequencies averaged over the sequence
    gtr.mu /= gtr.average_rate().mean()
    gtr_flat = GTR_site_specific.custom(alphabet=gtr.alphabet,  mu=gtr.mu, W=gtr.W, pi=gtr.Pi.mean(axis=1))
    gtr_flat.mu /= gtr_flat.average_rate().mean()


    branch_length = []
    for r in range(20): # repeat 20 times
        # sample a random sequence from the true profile
        s1 = prof2seq(gtr.Pi.T, gtr, sample_from_prof=True)
        s1_prof = seq2prof(s1[0], profile_map=gtr.profile_map)
        tmp_branch_length = []
        for t in np.linspace(0,1,31):
            # evolve the ancestral sequence for time t
            s2_prof = gtr.evolve(s1_prof, t)
            s2 = prof2seq(s2_prof, gtr, sample_from_prof=True)
            s2_prof_q = seq2prof(s2[0], profile_map=gtr.profile_map)

            # compute ML estimate for true model, flat model, and hamming distance
            tmp_branch_length.append((t,
                gtr.optimal_t_compressed((s1_prof, s2_prof_q), multiplicity = np.ones_like(gtr.mu), profiles=True),
                gtr_flat.optimal_t_compressed((s1_prof, s2_prof_q), multiplicity = np.ones_like(gtr.mu), profiles=True),
                1-np.mean(s1[-1]==s2[-1])))
        branch_length.append(tmp_branch_length)

    branch_length = np.array(branch_length)
    branch_length_mean = branch_length.mean(axis=0)
    branch_length_std = branch_length.std(axis=0)
    max_distance_flat = 1-np.sum(gtr_flat.Pi**2, axis=0).mean()
    max_distance = 1-np.sum(gtr.Pi**2, axis=0).mean()
    print("max_distance", max_distance, max_distance_flat)

    fig, axs = plt.subplots(1,2, figsize=(12,6))
    # plot estimate using true model
    axs[0].errorbar(branch_length_mean[:,0], branch_length_mean[:,1], branch_length_std[:,1], label='generating model')
    # plot estimate using flat model
    axs[0].errorbar(branch_length_mean[:,0], branch_length_mean[:,2], branch_length_std[:,2], label='average preferences')
    # plot hamming distance
    axs[0].errorbar(branch_length_mean[:,0], branch_length_mean[:,3], branch_length_std[:,3], label='p-distance')
    # add assymptotic similarity of heavily diverged sequences
    axs[0].plot([0.7,1], np.ones(2)*max_distance, ls='--', c='C2')
    # axs[0].plot([0,1], np.ones(2)*max_distance_flat, ls='--', c='C1')
    axs[0].set_xlabel('branch length', fontsize=fs)
    axs[0].set_ylabel('inferred branch length', fontsize=fs)
    axs[0].tick_params(labelsize=0.8*fs)
    axs[0].legend(fontsize=0.8*fs)
    add_panel_label(axs[0],'A',  x_offset=-0.15,y_offset=0.95, fs=fs)

    # add a second axis that coresponds roughly to years for HIV
    sec_axx = axs[0].twinx()
    sec_axy = axs[0].twiny()
    rate = 0.002
    sec_axx.set_ylim([x/rate for x in axs[0].get_ylim()])
    sec_axy.set_xlim([x/rate for x in axs[0].get_xlim()])
    #sec_axx.set_ylabel('inferred time [y]', fontsize=fs)
    sec_axy.set_xlabel('time [y]', fontsize=fs, labelpad=fs*0.3)
    sec_axx.tick_params(labelsize=0.8*fs)
    sec_axy.tick_params(labelsize=0.8*fs)

    axs[1].plot(sorted(gtr.mu), np.linspace(0,1,gtr.seq_len))
    axs[1].plot(sorted(gtr_flat.mu), np.linspace(0,1,gtr.seq_len))
    axs[1].set_xlabel('rate', fontsize=fs)
    axs[1].set_ylabel('fraction below rate x', fontsize=fs)
    axs[1].set_xscale('log')
    axs[1].tick_params(labelsize=0.8*fs)
    add_panel_label(axs[1],'B',  x_offset=-0.15,y_offset=0.95, fs=fs)
    plt.tight_layout()

    for fmt in ['.png', '.pdf']:
        plt.savefig(f"figures/{dset_name}_branch_length_error_pc_{args.pc:1.3f}{'_aa' if args.aa else ''}{fmt}")

