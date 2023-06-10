import numpy as np
import matplotlib.pyplot as plt
import glob

plt.rcParams['font.family'] = 'Times New Roman'
plt.rcParams['mathtext.fontset'] = 'stix'
plt.rcParams['font.size'] = 14

def theory():
    data = np.loadtxt(path+'/koide_theory.d')
    x1 = data[:, 0]
    y1 = data[:, 1]

    data = np.loadtxt(path+'/debug.d')
    x2 = data[:, 0]
    y2 = data[:, 1]


    fig = plt.figure(facecolor='w', figsize=(6, 4))
    ax = fig.add_subplot(1, 1, 1)

    ax.plot(x1, y1, lw=1, c='black', label='theory')
    ax.plot(x2, y2, lw=1, c='red', label='fenep')

    ax.plot([np.pi/2, np.pi/2], [-1, 10], lw=1, c='black', linestyle="dashed")
    ax.plot([np.pi*3/2, np.pi*3/2], [-1, 10], lw=1, c='black', linestyle="dashed")

    ax.set_xlim([0, 2*np.pi])
    ax.set_xticks(np.linspace(0, 2*np.pi, 5))
    ax.set_xticklabels(['$0$', '$\pi/2$', '$\pi$', '$3\pi/2$', '$2\pi$'], size=15)
    ax.set_xlabel(r'$y$', size=15)

    ax.set_ylim([min([min(y1), min(y2)])-0.1, max([max(y1), max(y2)])+0.1])
    ax.set_ylabel(r'$u$', size=15)
    ax.legend(loc='upper right', fontsize=15)

    plt.tight_layout()
    # plt.savefig(path+'_kabe/debug_U.png')
    # plt.savefig(path+'_kabe/debug_U.svg')
    plt.show()
    plt.close()


def fenep_dx():
    xmax = 2*np.pi
    x = [16, 32, 64, 128]
    if (pattern == 1): y = [0.192765710958E-01, 0.481914277392E-02, 0.120478569346E-02, 0.301196423369E-03]  # Re=1, beta=1.0
    if (pattern == 2): y = [0.193725255379E-01, 0.489131784533E-02, 0.123485327133E-02, 0.313037769414E-03]  # Re=1, beta=0.9, Wi=1
    if (pattern == 3): y = [0.204397650666E-01, 0.523836050399E-02, 0.132141205045E-02, 0.000348]  # Re=1, beta=0.9, Wi=10
    
    x = [xmax/x_ for x_ in x]

    a, b = np.polyfit(np.log10(x), np.log10(y), 1)  # 1次関数でフィッティング
    print(a, b)


    fig = plt.figure(facecolor='w')
    ax = fig.add_subplot(1, 1, 1)
    # ax = fig.add_axes([0.15, 0.15, 0.75, 0.75])
    ax.plot([min(x), max(x)], [min(x)**a*10**b, max(x)**a*10**b], lw=1, c='black')
    # ax.scatter(x, y, s=10, marker='x', c='black')
    ax.scatter(x, y, s=10, c='black')

    ax.set_xlim([min(x), max(x)])
    ax.set_xscale('log')
    ax.set_yscale('log')
    # ax.legend()
    # ax.set_xticks(np.array([0.1, 0.2, 0.4]))
    # ax.set_xticklabels([r'$10^{-1}$', r'$2 \times 10^{-1}$', r'$4 \times 10^{-1}$'])
    ax.set_xlabel(r'$\Delta x$', size=15)

    # ax.set_yticks(np.linspace(0, 2*np.pi, 3))
    # ax.set_yticklabels(['0', '$\pi$', '$2\pi$'])
    ax.set_ylabel(r'$\rm{RMSE}$', size=12)

    if (pattern == 1): fig.text(0.7, 0.5, 'slope : 2.000')
    if (pattern == 2): fig.text(0.7, 0.5, 'slope : 1.984')
    if (pattern == 3): fig.text(0.7, 0.5, 'slope : 1.962')

    plt.tight_layout()
    if (pattern == 1): fig.savefig('fenep_dx/fenep_dx_1.png')
    if (pattern == 1): fig.savefig('fenep_dx/fenep_dx_1.svg')
    if (pattern == 2): fig.savefig('fenep_dx/fenep_dx_2.png')
    if (pattern == 2): fig.savefig('fenep_dx/fenep_dx_2.svg')
    if (pattern == 3): fig.savefig('fenep_dx/fenep_dx_3.png')
    if (pattern == 3): fig.savefig('fenep_dx/fenep_dx_3.svg')
    plt.show()



if __name__ == '__main__':
    # path = 'newton_poiseuille'
    # path = 'newton_couette'
    path = 'debug'
    # theory()
    pattern = 3  # どのパラメータの結果を表示するか
    fenep_dx()