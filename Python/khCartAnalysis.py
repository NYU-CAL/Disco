import sys
import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt
import discopy.util as util

def analyze(filenames):

    N = len(filenames)

    eB = np.empty(N)
    t = np.empty(N)

    fig, ax = plt.subplots(5, 3, figsize=(9, 9))

    for i,f in enumerate(filenames):
        rep = util.DiscoReport(f)

        t = rep.t
        w = rep.dist_int[0]

        wvxc = rep.dist_int[3::4, :]
        wvxs = rep.dist_int[4::4, :]
        wvyc = rep.dist_int[5::4, :]
        wvys = rep.dist_int[6::4, :]

        mvx0 = 2*np.fabs(rep.dist_int[1]/w)
        mvy0 = 2*np.fabs(rep.dist_int[2]/w)

        mvx1 = 2*np.sqrt((wvxc[0]/w)**2 + (wvxs[0]/w)**2)
        mvy1 = 2*np.sqrt((wvyc[0]/w)**2 + (wvys[0]/w)**2)

        mvx2 = 2*np.sqrt((wvxc[1]/w)**2 + (wvxs[1]/w)**2)
        mvy2 = 2*np.sqrt((wvyc[1]/w)**2 + (wvys[1]/w)**2)

        mvx3 = 2*np.sqrt((wvxc[2]/w)**2 + (wvxs[2]/w)**2)
        mvy3 = 2*np.sqrt((wvyc[2]/w)**2 + (wvys[2]/w)**2)

        mvx4 = 2*np.sqrt((wvxc[3]/w)**2 + (wvxs[3]/w)**2)
        mvy4 = 2*np.sqrt((wvyc[3]/w)**2 + (wvys[3]/w)**2)

        mask = t > t[0]

        gy0 = np.log(mvy0[1:]/mvy0[:-1]) / (t[1:]-t[:-1])
        gy1 = np.log(mvy1[1:]/mvy1[:-1]) / (t[1:]-t[:-1])
        gy2 = np.log(mvy2[1:]/mvy2[:-1]) / (t[1:]-t[:-1])
        gy3 = np.log(mvy3[1:]/mvy3[:-1]) / (t[1:]-t[:-1])
        gy4 = np.log(mvy4[1:]/mvy4[:-1]) / (t[1:]-t[:-1])

        t2 = 0.5*(t[1:] + t[:-1])

        mask2 = t2 > t2[4]

        ax[0, 0].plot(t[mask], mvx0[mask])
        ax[0, 1].plot(t[mask], mvy0[mask])
        ax[1, 0].plot(t[mask], mvx1[mask])
        ax[1, 1].plot(t[mask], mvy1[mask])
        ax[2, 0].plot(t[mask], mvx2[mask])
        ax[2, 1].plot(t[mask], mvy2[mask])
        ax[3, 0].plot(t[mask], mvx3[mask])
        ax[3, 1].plot(t[mask], mvy3[mask])
        ax[4, 0].plot(t[mask], mvx4[mask])
        ax[4, 1].plot(t[mask], mvy4[mask])
        ax[0, 2].plot(t2[mask2], gy0[mask2])
        ax[1, 2].plot(t2[mask2], gy1[mask2])
        ax[2, 2].plot(t2[mask2], gy2[mask2])
        ax[3, 2].plot(t2[mask2], gy3[mask2])
        ax[4, 2].plot(t2[mask2], gy4[mask2])

    ax[0,0].set(xlabel=r'$t$', ylabel=r'$|\tilde{v}_x^0|$',
                xscale='linear', yscale='log')
    ax[0,1].set(xlabel=r'$t$', ylabel=r'$|\tilde{v}_y^0|$',
                xscale='linear', yscale='log')
    ax[0,2].set(xlabel=r'$t$', ylabel=r'$\Gamma_y^0$',
                xscale='linear', yscale='linear')
    ax[1,0].set(xlabel=r'$t$', ylabel=r'$\tilde{v}_x^1$',
                xscale='linear', yscale='log')
    ax[1,1].set(xlabel=r'$t$', ylabel=r'$\tilde{v}_y^1$',
                xscale='linear', yscale='log')
    ax[1,2].set(xlabel=r'$t$', ylabel=r'$\Gamma_y^1$',
                xscale='linear', yscale='linear')
    ax[2,0].set(xlabel=r'$t$', ylabel=r'$\tilde{v}_x^2$',
                xscale='linear', yscale='log')
    ax[2,1].set(xlabel=r'$t$', ylabel=r'$\tilde{v}_y^2$',
                xscale='linear', yscale='log')
    ax[2,2].set(xlabel=r'$t$', ylabel=r'$\Gamma_y^2$',
                xscale='linear', yscale='linear')
    ax[3,0].set(xlabel=r'$t$', ylabel=r'$\tilde{v}_x^3$',
                xscale='linear', yscale='log')
    ax[3,1].set(xlabel=r'$t$', ylabel=r'$\tilde{v}_y^3$',
                xscale='linear', yscale='log')
    ax[3,2].set(xlabel=r'$t$', ylabel=r'$\Gamma_y^3$',
                xscale='linear', yscale='linear')
    ax[4,0].set(xlabel=r'$t$', ylabel=r'$\tilde{v}_x^4$',
                xscale='linear', yscale='log')
    ax[4,1].set(xlabel=r'$t$', ylabel=r'$\tilde{v}_y^4$',
                xscale='linear', yscale='log')
    ax[4,2].set(xlabel=r'$t$', ylabel=r'$\Gamma_y^4$',
                xscale='linear', yscale='linear')
        
    fig.tight_layout()


    figname = "khModes.png"
    print("Saving " + figname)
    fig.savefig(figname)

if __name__ == "__main__":

    if len(sys.argv) < 2:
        print("Need a report dude!")
        sys.exit()

    filenames = sys.argv[1:]
    analyze(filenames)
