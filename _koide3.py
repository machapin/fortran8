import numpy as np
import matplotlib.pyplot as plt

plt.rcParams['font.family'] = 'Times New Roman'
plt.rcParams['mathtext.fontset'] = 'stix'
plt.rcParams['font.size'] = 15


def time_trC():
    data = np.loadtxt(dir+'/koide3_cardano.d')
    x1 = data[:, 0]
    y1 = data[:, 1]

    data = np.loadtxt(dir+'/koide3_lapack.d')
    x2 = data[:, 0]
    y2 = data[:, 1]

    data = np.loadtxt(dir+'/koide3_sylve.d')
    x3 = data[:, 0]
    y3 = data[:, 1]
    
    data = np.loadtxt(dir+'/koide3_vieta.d')
    x4 = data[:, 0]
    y4 = data[:, 1]

    fig = plt.figure(facecolor='w', figsize=(6, 4))
    ax = fig.add_subplot(1, 1, 1)


    # ax.plot(x1, y1, lw=1, c='black', label='Cardano')
    ax.plot(x2, y2, lw=1, c='red', label='Lapack')
    ax.plot(x3, y3, lw=1, c='blue', label='Sylvester')
    ax.plot(x4, y4, lw=1, c='green', label='Vieta')
    ax.plot(x1, y1, lw=1, c='black', label='Cardano', linestyle="dashed")

    ax.set_xlim([0, 20000])
    ax.set_xlabel(r'step', size=15)
    ax.set_ylabel(r'$\rm{tr}{C}$', size=15)
    ax.legend(loc='upper right', fontsize=14)


    plt.tight_layout()
    plt.savefig(dir+'/koide3.png')
    plt.show()
    plt.close()



def time_dif_trC():
    data = np.loadtxt(dir+'/koide3_cardano.d')
    x1 = data[:, 0]
    y1 = data[:, 1]

    data = np.loadtxt(dir+'/koide3_lapack.d')
    x2 = data[:, 0]
    y2 = data[:, 1]

    data = np.loadtxt(dir+'/koide3_sylve.d')
    x3 = data[:, 0]
    y3 = data[:, 1]
    
    data = np.loadtxt(dir+'/koide3_vieta.d')
    x4 = data[:, 0]
    y4 = data[:, 1]

    fig = plt.figure(facecolor='w', figsize=(6, 4))
    ax = fig.add_subplot(1, 1, 1)


    ax.plot(x3, abs(y1-y2), lw=1, c='black', label='Cardano vs Lapack'); name = 1
    # ax.plot(x4, abs(y4-y1), lw=1, c='black', label='Cardano vs Vieta'); name = 2
    # ax.plot(x3, abs(y3-y1), lw=1, c='black', label='Cardano vs Sylvester'); name = 3
    # ax.plot(x3, abs(y3-y4), lw=1, c='black', label='Vieta vs Sylvester'); name = 4

    ax.set_xlim([0, 20000])
    ax.set_xlabel(r'step', size=15)
    ax.set_ylabel(r'dif $\rm{tr}{C}$', size=15)
    ax.set_yscale('log')
    ax.legend(loc='upper right', fontsize=14)


    plt.tight_layout()
    plt.savefig(dir+'/koide3_dif'+str(name)+'.png')
    plt.show()
    plt.close()

if __name__ == '__main__':
    dir = 'koide3'
    time_trC()
    # time_dif_trC()