import glob
import matplotlib.pyplot as plt
from generate_toy_data import *

def entropy(a):
    eps = 1e-10
    return -np.sum(np.log(a+eps)*(a+eps), axis=0)

def partip(a):
    return 1/np.sum(a**2, axis=0)

if __name__ == '__main__':

    import argparse
    parser = argparse.ArgumentParser(description = "", usage="analyze simulated data for site specific GTR reconstruction project")
    parser.add_argument('--mask', type=str, nargs='+')
    parser.add_argument('--labels', type=str, nargs='+')
    args = parser.parse_args()

    fs = 16
    plt.figure('S')
    plt.figure('N')
    for label, mask in zip(args.labels, args.mask):
        files = glob.glob(mask)
        allS = []
        allN = []
        for file in files:
            m = load_model(file)
            S = entropy(m.Pi)
            N = partip(m.Pi)
            allS.extend(S)
            allN.extend(N)

        plt.figure('S')
        plt.plot(sorted(allS), np.linspace(0,1,len(allS)), label=label)
        plt.figure('N')
        plt.plot(sorted(allN), np.linspace(0,1,len(allN)), label=label)
        print(f"{mask}, average entropy: {np.mean(allS):1.3f}, average N: {np.mean(allN):1.3f}")

    plt.figure('S')
    plt.ylabel('fraction with less entropy', fontsize=fs)
    plt.xlabel('fraction of sites', fontsize=fs)
    plt.legend()
    plt.savefig('figures/entropy_stats.pdf')

    plt.figure('N')
    plt.ylabel('fraction with fewer effective states', fontsize=fs)
    plt.xlabel('fraction of sites', fontsize=fs)
    plt.legend()
    plt.savefig('figures/partip_stats.pdf')
