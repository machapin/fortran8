import numpy as np
import matplotlib.pyplot as plt
# import japanize_matplotlib
import glob
import os

plt.rcParams['font.family'] = 'Times New Roman'
plt.rcParams['mathtext.fontset'] = 'stix'
plt.rcParams['font.size'] = 14

def UC_fig():
    data = np.loadtxt(dir+'/kabe_'+path+'/debug.d')
    x1 = data[:, 0]
    y1 = data[:, 1]
    c1 = data[:, 2]

    data = np.loadtxt(dir+'/ibm01_'+path+'/debug.d')
    x2 = data[:, 0]
    y2 = data[:, 1]
    c2 = data[:, 2]

    data = np.loadtxt(dir+'/ibm02_'+path+'/debug.d')
    x3 = data[:, 0]
    y3 = data[:, 1]
    c3 = data[:, 2]

    if os.path.exists(dir+'/kabe_'+path+'/koide_theory.d'):
        data = np.loadtxt(dir+'/kabe_'+path+'/koide_theory.d')
        x4 = data[:, 0]
        y4 = data[:, 1]

    fig = plt.figure(facecolor='w', figsize=(6, 4))
    ax = fig.add_subplot(1, 1, 1)

    if (index == 'U'):
        ax.plot(x1, y1, lw=1, c='black', label='A', marker='o', markersize=2)
        ax.plot(x2, y2, lw=1, c='red', label='B', marker='o', markersize=2)
        ax.plot(x3, y3, lw=1, c='blue', label='C', marker='o', markersize=2)
        if os.path.exists(dir+'/kabe_'+path+'/koide_theory.d'):
            ax.plot(x4, y4, lw=1, c='green', label='theory', linestyle="dashed")
    if (index == 'C'):
        ax.plot(x1, c1, lw=1, c='black', label='A', marker='o', markersize=2)
        ax.plot(x2, c2, lw=1, c='red', label='B', marker='o', markersize=2)
        ax.plot(x3, c3, lw=1, c='blue', label='C', marker='o', markersize=2)

    ax.plot([np.pi/2, np.pi/2], [-1, 10], lw=1, c='black', linestyle="dashed")
    ax.plot([np.pi*3/2, np.pi*3/2], [-1, 10], lw=1, c='black', linestyle="dashed")

    ax.set_xlim([0, 2*np.pi])
    ax.set_xticks(np.linspace(0, 2*np.pi, 5))
    ax.set_xticklabels(['$0$', '$\pi/2$', '$\pi$', '$3\pi/2$', '$2\pi$'], size=15)
    ax.set_xlabel(r'$y$', size=15)

    if (index == 'U'): ax.set_ylim([min([min(y1), min(y2), min(y3)])-0.1, max([max(y1), max(y2), max(y3)])+0.1])
    if (index == 'C'): ax.set_ylim([2-0.1, max([max(c1), max(c2), max(c3)])+0.1])
    if (index == 'U'): ax.set_ylabel(r'$u$', size=15)
    if (index == 'C'): ax.set_ylabel(r'$\rm{tr}{C}$', size=15)
    ax.legend(loc='upper right', fontsize=14)

    plt.tight_layout()
    if (index == 'U'): plt.savefig(dir+'/kabe_'+path+'/debug_U_dt'+dt+'.png')
    if (index == 'U'): plt.savefig(dir+'/kabe_'+path+'/debug_U_dt'+dt+'.svg')
    if (index == 'C'): plt.savefig(dir+'/kabe_'+path+'/debug_C_dt'+dt+'.png')
    if (index == 'C'): plt.savefig(dir+'/kabe_'+path+'/debug_C_dt'+dt+'.svg')
    plt.show()
    plt.close()

if __name__ == '__main__':
    dir = 'houkokukai'
    index = 'U'
    # index = 'C'
    # path = 'poiseuille_newton'
    path = 'poiseuille'
    # path = 'couette'
    dt = '0001'
    UC_fig()