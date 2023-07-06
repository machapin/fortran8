import numpy as np
import matplotlib.pyplot as plt
import glob

plt.rcParams['font.family'] = 'Times New Roman'
plt.rcParams['mathtext.fontset'] = 'stix'
plt.rcParams['font.size'] = 15

def energy_spectrum(file):
    for path in glob.glob(f'{file}/energy*50000.d'):
        # path = sorted(glob.glob(f'{file}/energy_0*.d'))[0]
        # path = sorted(glob.glob(f'{file}/energy_dif_0*.d'))[0]
        # path = sorted(glob.glob(f'{file}/energy_mask_0*.d'))[0]
        # path = sorted(glob.glob(f'{file}/energy_z_0*.d'))[0]
        data = np.loadtxt(path)
        x = data[:, 0]
        y = data[:, 1]
        p = 0.0001  # 点線の開始位置
        xmax = np.max(x[y>0])  # 値があるところの最大値

        fig = plt.figure(facecolor='w', figsize=(5, 5))
        ax = fig.add_axes([0.18, 0.15, 0.75, 0.8])

        ax.plot(x, y, lw=1, c='black')
        ax.plot([1, xmax], [p*1**(-5.0/3.0), p*xmax**(-5.0/3.0)], c='black', lw=0.7, linestyle="dashed", label='-5/3')

        ax.set_xlim([1, xmax])
        ax.set_ylim([p*1.0E-8, p*10])
        ax.set_xscale('log')
        ax.set_yscale('log')
        ax.set_xlabel(r'$k$', size=17)
        ax.set_ylabel(r'$E(k)$', size=17)
        fig.text(0.7, 0.65, r'$E(k) \propto k^{-5/3}} $', size=15)
        # plt.tight_layout()
        plt.savefig(path[:-2] + '.png')
        # plt.savefig(path[:-2] + '.svg')
        # plt.show()
        plt.close()

def energy_spectrum_sum(file):
    i = 0
    # for path in glob.glob(f'{file}/energy_same_0*.d'):
    # for path in glob.glob(f'{file}/energy_diff_0*.d'):
    for path in glob.glob(f'{file}/energy_same_mask_0*.d'):
    # for path in glob.glob(f'{file}/energy_diff_mask_0*.d'):
    # for path in glob.glob(f'{file}/energy_same_z_0*.d'):
    # for path in glob.glob(f'{file}/energy_diff_z_0*.d'):
        data = np.loadtxt(path)
        x = data[:, 0]
        if(i==0): y = data[:, 1]
        if(i!=0): y += data[:, 1]
        i += 1
    y = y/i
    p = 0.0001  # 点線の開始位置
    xmax = np.max(x[y>0])  # 値があるところの最大値

    fig = plt.figure(facecolor='w', figsize=(5, 5))
    ax = fig.add_axes([0.18, 0.15, 0.75, 0.8])

    ax.plot(x, y, lw=1, c='black')
    ax.plot([1, xmax], [p*1**(-5.0/3.0), p*xmax**(-5.0/3.0)], c='black', lw=0.7, linestyle="dashed", label='-5/3')

    ax.set_xlim([1, xmax])
    ax.set_ylim([p*1.0E-8, p*10])
    ax.set_xscale('log')
    ax.set_yscale('log')
    ax.set_xlabel(r'$k$', size=17)
    ax.set_ylabel(r'$E(k)$', size=17)
    fig.text(0.7, 0.65, r'$E(k) \propto k^{-5/3}} $', size=15)
    # plt.tight_layout()
    plt.savefig(file + '/' + path[9:-9] + '_sum.png')
    # plt.show()
    plt.close()
    
def energy_time(file, start=None, end=None):
    for path in glob.glob(f'{file}/time*.d'):
        # path = file+'/time_same.d'
        data = np.loadtxt(path)

        # 読み込むステップ数
        if (start==None) : start = 0
        if (end==None) : end = data.shape[0]
        # x = data[start:end, 0]
        x = data[start:end, 1]
        # y = data[start:end, 2]; title = 'K'
        # y = data[start:end, 3]/(256**3); title = 'K_s'
        # y = data[start:end, 4]; title = '\lambda'
        # y = data[start:end, 5]; title = 'R_\lambda'
        # y = data[start:end, 6]; title = '\epsilon'
        y = data[start:end, 7]; title = '\eta'
        # y = data[start:end, 8]; title = '\tau'
        # y = data[start:end, 9]; title = 'v'
        # y = data[start:end, 10]; title = 'Ep_Energy'
        # y = data[start:end, 11]; title = 'trC_mean'
        # y = data[start:end, 12]; title = 'trC_std'


        fig = plt.figure(facecolor='w', figsize=(8, 4))
        ax = fig.add_axes([0.18, 0.15, 0.75, 0.8])
        # ax.plot(x, y, lw=1, c='black', label='Newton')
        ax.plot(x, y, lw=1, c='black')

        ax.set_xlim(min(x), max(x))
        ax.set_xlabel(r'$t$', size=17)
        # ax.set_ylim(min(y)-0.5, max(y)+0.5)
        ax.set_ylabel(r'${:}$'.format(title), size=17)
        if (title == '\tau') : ax.set_ylabel(r'$\tau$', size=17)  # tauだけうまく行かないなぜ
        # ax.legend(loc='best', fontsize=15)

        # plt.tight_layout()
        if (title == 'R_\lambda'): title = 'Re'
        if (title == '\lambda'): title = 'lambda'
        if (title == '\epsilon'): title = 'epsilon'
        if (title == '\eta'): title = 'eta'
        if (title == '\tau'): title = 'tau'

        fig.savefig(path[:-2] + '_' + title + '.png')
        # fig.savefig(path[:-2] + '_' + title + '.svg')
        # plt.show()
        plt.close()


def Ktime_twinx(file, start=None, end=None):
    for path in glob.glob(f'{file}/all_time*.d'):
        # path = file+'/all_time.d'
        # path = file+'/all_time_dif.d'
        # path = file+'/all_time_dif_map.d'
        data = np.loadtxt(path)

        # 読み込むステップ数
        if (start==None) : start = 0
        if (end==None) : end = data.shape[0]
        # x = data[start:end, 0]
        x = data[start:end, 1]
        # y1 = data[start:end, 2]; title1 = 'K'
        # y1 = data[start:end, 3]; title1 = 'K_staggered'
        y1 = data[start:end, 5]; title1 = 'R_\lambda'
        y2 = data[start:end, 6]; title2 = '\epsilon'

        
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
        ax1.set_ylabel(r'${:}$'.format(title1), size=17)
        ax2.set_ylabel(r'${:}$'.format(title2), size=17)

        # from matplotlib.ticker import ScalarFormatter
        # ax2.yaxis.set_major_formatter(ScalarFormatter(useMathText=True))
        # ax2.ticklabel_format(style="sci",  axis="y",scilimits=(0,0))

        # ax2.set_yticks(np.linspace(0.001, 0.0012,3))
        # ax2.set_yticklabels(['0.0010', '0.0011', '0.0012'])
        # ax2.set_ylim(0.00095, 0.00125)
        # ax2.legend(loc='best')
        ax1.legend(handles = [p1, p2], labels = [r'${:}$'.format(title1), r'${:}$'.format(title2)], loc='best', fontsize=15)

        # plt.tight_layout()
        fig.savefig(path[:-2] + '_twin.png')
        # fig.savefig(path[:-2] + '.svg')
        # plt.show()
        plt.close()


def power_spectrum(file):
    for path in glob.glob(f'{file}/Spectrum_*.d'):
        # path = file+'/Spectrum_x.d'
        # path = file+'/Spectrum_y.d'
        # path = file+'/Spectrum_z.d'
        data = np.loadtxt(path)
        x = data[:, 0]
        y = data[:, 1]
        p = 1  # 点線の開始位置
        xmax = min(np.max(x), 10000000)  # xの描画範囲

        fig = plt.figure(facecolor='w', figsize=(5, 5))
        ax = fig.add_axes([0.18, 0.15, 0.75, 0.8])

        ax.plot(x, y, lw=1, c='black')
        ax.plot([1, xmax], [p*1**(-5.0/3.0), p*xmax**(-5.0/3.0)], c='black', lw=0.7, linestyle="dashed", label='-5/3')

        ax.set_xlim([1, xmax])
        ax.set_ylim([p*1.0E-8, p*10])
        ax.set_xscale('log')
        ax.set_yscale('log')
        ax.set_xlabel(r'$f$', size=17)
        ax.set_ylabel(r'$E(f)$', size=17)
        fig.text(0.7, 0.65, r'$E(f) \propto f^{-5/3}} $', size=15)
        # plt.tight_layout()
        plt.savefig(path[:-2] + '.png')
        # plt.savefig(path[:-2] + '.svg')
        # plt.show()
        plt.close()


def step_velocity(file):
    for index in range(3):
        path = file+'/step_velocity.d'
        data = np.loadtxt(path)

        # x = data[:, 0]
        x = data[:, 1]
        title_name = ['U', 'V', 'W']
        y = data[:, index+2]; title = title_name[index]

        fig = plt.figure(facecolor='w', figsize=(8, 4))
        ax = fig.add_axes([0.18, 0.15, 0.75, 0.8])
        # ax.plot(x, y, lw=1, c='black', label='Newton')
        ax.plot(x, y, lw=1, c='black')

        ax.set_xlim(min(x), max(x))
        ax.set_ylim([-0.2, 0.2])
        ax.set_xlabel(r'$t$', size=17)
        ax.set_ylabel(r'${:}$'.format(title), size=17)

        # plt.tight_layout()
        fig.savefig(path[:-2] + '_' + title + '.png')
        # fig.savefig(path[:-2] + '_' + title + '.svg')
        # plt.show()
        plt.close()

if __name__ == '__main__':
    file = 'v008/data4'
    # energy_spectrum(file)
    # energy_spectrum_sum(file)
    energy_time(file, start=None, end=None)
    # Ktime_twinx(file, start=None, end=None)
    # power_spectrum(file)
    # step_velocity(file)