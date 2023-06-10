import numpy as np
import matplotlib.pyplot as plt
import glob

plt.rcParams['font.family'] = 'Times New Roman'
plt.rcParams['mathtext.fontset'] = 'stix'
plt.rcParams['font.size'] = 15

def energy_spectrum(file):
    path = sorted(glob.glob(f'{file}/energy_*.d'))[0]
    data = np.loadtxt(path)
    x = data[:, 0]
    y = data[:, 1]
    p = 0.01  # E(1)の大きさ
    xmax = np.max(x[y>0])  # 値があるところの最大値

    # plt.scatter(x, y, s=0.2, c='r')
    # plt.plot(x, y, lw=2, c='r')
    # plt.plot([1, 100], [10*1**(-5.0/3.0), 10*100**(-5.0/3.0)], c='black', lw=1, linestyle="dashed", label='-5/3')
    # plt.plot([1, 100], [0.001*1**(-3), 0.001*100**(-3)], c='black', lw=1, linestyle="dashed")
    # plt.plot([1, 100], [100*1**(-4.7), 100*100**(-4.7)], c='black', lw=1, linestyle="dashed", label='-4.7')

    fig = plt.figure(facecolor='w', figsize=(5, 5))
    ax = fig.add_axes([0.18, 0.15, 0.75, 0.8])

    ax.plot(x, y, lw=1, c='black')
    ax.plot([1, xmax], [p*1**(-5.0/3.0), p*xmax**(-5.0/3.0)], c='black', lw=0.7, linestyle="dashed", label='-5/3')
    # ax.plot([1, xmax], [500*1**(-4.7), 500*xmax**(-4.7)], c='black', lw=0.7, linestyle="dashed", label='-4.7')

    ax.set_xlim([1, xmax])
    ax.set_ylim([p*1.0E-8, p*10])
    ax.set_xscale('log')
    ax.set_yscale('log')
    ax.set_xlabel(r'$k$', size=17)
    ax.set_ylabel(r'$E(k)$', size=17)
    fig.text(0.7, 0.65, r'$E(k) \propto k^{-5/3}} $', size=15)
    plt.tight_layout()
    plt.savefig(path[:-2] + '.png')
    # plt.savefig(path[:-2] + '.svg')
    plt.show()
    plt.close()
    
def energy_time(file, start=None, end=None):
    path = file+'/all_time.d'
    # path = file+'/all_time_dif.d'
    data = np.loadtxt(path)

    # 読み込むステップ数
    if (start==None) : start = 0
    if (end==None) : end = data.shape[0]
    # x = data[start:end, 0]
    x = data[start:end, 1]
    y = data[start:end, 2]; tilte = 'K'  # 6/4以降で変更あり
    y = data[start:end, 3]; tilte = '\lambda'
    # y = data[start:end, 4]; tilte = 'R_\lambda'
    # y = data[start:end, 5]; tilte = '\epsilon'
    # y = data[start:end, 6]; tilte = '\eta'
    # y = data[start:end, 7]; tilte = '\tau'
    # y = data[start:end, 8]; tilte = 'v'
    # y = data[start:end, 9]; tilte = 'Ep_Energy'
    # y = data[start:end, 10]; tilte = 'trC_mean'
    # y = data[start:end, 11]; tilte = 'trC_std'


    fig = plt.figure(facecolor='w', figsize=(8, 4))
    ax = fig.add_axes([0.18, 0.15, 0.75, 0.8])
    # ax.plot(x, y, lw=1, c='black', label='Newton')
    ax.plot(x, y, lw=1, c='black')

    ax.set_xlim(min(x), max(x))
    ax.set_xlabel(r'$t$', size=17)
    # ax.set_ylim(min(y)-0.5, max(y)+0.5)
    ax.set_ylabel(r'${:}$'.format(tilte), size=17)
    if (tilte == '\tau') : ax.set_ylabel(r'$\tau$', size=17)  # tauだけうまく行かないなぜ
    # ax.legend(loc='best', fontsize=15)

    # plt.tight_layout()
    if (tilte == 'R_\lambda'): tilte = 'Re'
    if (tilte == '\lambda'): tilte = 'lambda'
    if (tilte == '\epsilon'): tilte = 'epsilon'
    if (tilte == '\eta'): tilte = 'eta'
    if (tilte == '\tau'): tilte = 'tau'

    fig.savefig(path[:-2] + '_' + tilte + '.png')
    # fig.savefig(path[:-2] + '.svg')
    # plt.show()
    plt.close()


def Ktime_twinx(file, start=None, end=None):
    path = file+'/all_time.d'
    # path = file+'/all_time_dif.d'
    data = np.loadtxt(path)

    # 読み込むステップ数
    if (start==None) : start = 0
    if (end==None) : end = data.shape[0]
    # x = data[start:end, 0]
    x = data[start:end, 1]
    # y1 = data[start:end, 2]; tilte1 = 'K'
    y1 = data[start:end, 4]; tilte1 = 'R_\lambda'
    y2 = data[start:end, 5]; tilte2 = '\epsilon'

    
    fig = plt.figure(facecolor='w', figsize=(8, 4))
    ax1 = fig.add_axes([0.15, 0.15, 0.72, 0.8])
    ax2 = ax1.twinx()

    # ax1.scatter(x, y1, s=0.2, c='r')
    # ax2.scatter(x, y2, s=0.2, c='b')
    p1, = ax1.plot(x, y1, lw=1, c='black')
    # ax2.plot(x, y2, lw=1, c='red')
    p2, = ax2.plot(x, y2, lw=1, c='red')

    ax1.set_xlim(min(x), max(x))
    ax1.set_xlabel(r'$t$', size=15)
    # ax1.set_ylim(min(y1)-(max(y1)-min(y1))*0.1, max(y1)+(max(y1)-min(y1))*0.1)
    # ax2.set_ylim(min(y2)-(max(y2)-min(y2))*0.1, max(y2)+(max(y2)-min(y2))*0.1)
    ax1.set_ylabel(r'${:}$'.format(tilte1), size=17)
    ax2.set_ylabel(r'${:}$'.format(tilte2), size=17)

    # from matplotlib.ticker import ScalarFormatter
    # ax2.yaxis.set_major_formatter(ScalarFormatter(useMathText=True))
    # ax2.ticklabel_format(style="sci",  axis="y",scilimits=(0,0))

    # ax2.set_yticks(np.linspace(0.001, 0.0012,3))
    # ax2.set_yticklabels(['0.0010', '0.0011', '0.0012'])
    # ax2.set_ylim(0.00095, 0.00125)
    # ax2.legend(loc='best')
    ax1.legend(handles = [p1, p2], labels = [r'${:}$'.format(tilte1), r'${:}$'.format(tilte2)], loc='best', fontsize=15)

    # plt.tight_layout()
    fig.savefig(path[:-2] + '_twin.png')
    # fig.savefig(path[:-2] + '.svg')
    plt.show()
    plt.close()



if __name__ == '__main__':
    file = 'v005/data3'
    # energy_spectrum(file)
    energy_time(file, start=None, end=None)
    # Ktime_twinx(file, start=None, end=None)