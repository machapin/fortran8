import numpy as np
import matplotlib.pyplot as plt
import matplotlib.cm as cm
import matplotlib.colors as colors
from tqdm import tqdm
import glob

plt.rcParams['font.family'] = 'Times New Roman'
plt.rcParams['mathtext.fontset'] = 'stix'
plt.rcParams['font.size'] = 16

def mapping(vmax=None, vmin=None):
    if (form == 'd'):
        # for filename in tqdm(glob.glob(f'{dir}/{ziku}_*.d')):
            filename = sorted(glob.glob(f'{dir}/{ziku}_*.d'))[-1]
            data = np.loadtxt(filename)
            step = int(filename[-10:-2])
            map_fig(data, step, vmax=vmax, vmin=vmin)

    if (form == 'bin'):
        # for filename in tqdm(glob.glob(f'{dir}/{ziku}_*.bin')):
            filename = sorted(glob.glob(f'{dir}/{ziku}_*.bin'))[-1]
            with open(filename, mode='rb') as f:
                data = np.fromfile(f, dtype=np.float64)
            data = data.reshape(-1, 18)
            step = int(filename[-12:-4])
            map_fig(data, step, vmax=vmax, vmin=vmin)


def map_fig(data, step, vmax, vmin):
    # color = 'plasma'  # 特に指定ないとき
    color = 'seismic'  # 特に指定ないとき
    # Z = data[:, 0]; title = 'x'  # デバッグ用
    # Z = data[:, 1]; title = 'y'
    # Z = data[:, 2]; title = 'z'
    # Z = data[:, 3]; title = 'U'
    # Z = data[:, 4]; title = 'V'
    Z = data[:, 5]; title = 'W'
    
    # Z = data[:, 6]; title = 'Qti'  # 6/15以降はQ値
    # if (ziku=='x') : Z = data[:, 7]; title = 'omegax'; color = 'seismic'
    # if (ziku=='y') : Z = data[:, 8]; title = 'omegay'; color = 'seismic'
    # if (ziku=='z') : Z = data[:, 9]; title = 'omegaz'; color = 'seismic'
    # Z = data[:, 10]; title = 'enstrophy'
    
    # Z = data[:, 11]; title = 'trC'; color = 'plasma'
    # Z = data[:, 12]; title = 'C11'; color = 'plasma'
    # Z = data[:, 13]; title = 'C12'; color = 'plasma'
    # Z = data[:, 14]; title = 'C13'; color = 'plasma'
    # Z = data[:, 15]; title = 'C22'; color = 'plasma'
    # Z = data[:, 16]; title = 'C23'; color = 'plasma'
    # Z = data[:, 17]; title = 'C33'; color = 'plasma'

    if (ziku == 'm'):
        # Z = data[:, 2]; title = 'K'
        # Z = data[:, 3]; title = '\lambda'
        Z = data[:, 4]; title = 'R_\lambda'; color = 'plasma'
        # Z = data[:, 5]; title = '\epsilon'
        # Z = data[:, 6]; title = '\eta'
        # Z = data[:, 7]; title = '\tau'
        
        # Z = data[:, 8]; title = 'v'
        # Z = data[:, 9]; title = 'Ep_Energy'
        # Z = data[:, 10]; title = 'trC_mean'
        # Z = data[:, 11]; title = 'trC_std'
        # Z = data[:, 12]; title = 'mask'; color = 'binary'
    

    size = int(len(Z)**0.5)
    Z = Z.reshape(size, size)[::-1, :]

    fig = plt.figure(facecolor='w', figsize=(7, 5))
    ax = fig.add_subplot(1, 1, 1)

    # ax.set_title(f'step:{step}')
    xmax = 2*np.pi
    
    if (vmax == None): vmax = np.max(np.abs(Z))
    if (vmin == None):
        if (np.min(Z) <= 0) : vmin = -vmax  # 正負対称なら
        if (np.min(Z) >= 0) : vmin = np.min(Z)  # 正のみなら


    ax.set_xticks(np.linspace(0, 2*np.pi, 3))
    if (ibm <= 0): ax.set_xticklabels(['$0$', '$\pi$', '$2\pi$'], size=15)
    if (ibm >= 1): ax.set_xticklabels(['$0$', '$2D$', '$4D$'], size=15)
    ax.set_yticks(np.linspace(0, 2*np.pi, 3))
    if (ibm <= 0): ax.set_yticklabels(['$0$', '$\pi$', '$2\pi$'], size=15)
    if (ibm >= 1): ax.set_yticklabels(['$0$', '$2D$', '$4D$'], size=15)
    # ax.set_ylabel(r'$y$', size=18, rotation=0)
    if(ziku=='x'): ax.set_xlabel(r'$y$', size=18)
    if(ziku=='y'): ax.set_xlabel(r'$x$', size=18)
    if(ziku=='z'): ax.set_xlabel(r'$x$', size=18)
    if(ziku=='m'): ax.set_xlabel(r'$x$', size=18)
    if(ziku=='x'): ax.set_ylabel(r'$z$', size=18)
    if(ziku=='y'): ax.set_ylabel(r'$z$', size=18)
    if(ziku=='z'): ax.set_ylabel(r'$y$', size=18)
    if(ziku=='m'): ax.set_ylabel(r'$y$', size=18)

    # if 1:  # マスクしたいとき
    #     mask = data[:, 12]; title = 'Re_mask'
    #     mask = mask.reshape(size, size)[::-1, :]
    #     Z = Z * mask
    #     cmap = cm.plasma
    #     cmap_data = cmap(np.arange(cmap.N))
    #     print(cmap_data)
    #     cmap_data[0, 3] = 0  # 0 のときのα値を0(透明)にする
    #     print(cmap_data)
    #     customized_cool = colors.ListedColormap(cmap_data)
    #     im = ax.imshow(Z, interpolation='nearest', extent=[0, xmax, 0, xmax], cmap=customized_cool, vmin=vmin, vmax=vmax)

    im = ax.imshow(Z, interpolation='nearest', extent=[0, xmax, 0, xmax], cmap=color, vmin=vmin, vmax=vmax)
    fig.colorbar(im)
    

    # if (ziku=='z' and ibm==2):
    #     for i in range(1, ibm+1):
    #         data = np.loadtxt('ibm'+str(ibm)+'/{:08d}.d'.format(i))
    #         xc = data[:, 0]
    #         yc = data[:, 1]
    #         ax.plot(xc, yc, lw=1.0, color='black', linestyle="dashed")
    #         # ax.scatter(xc, yc, s=1.0, color='black')

    if (ziku=='z' and ibm==4):
        for i in range(1, ibm+1):
            data = np.loadtxt('ibm/ibm21/{:08d}.d'.format(i))
            xc = data[:, 0]
            yc = data[:, 1]
            ax.plot(xc, yc, lw=1.0, color='black', linestyle="dashed")
            # ax.scatter(xc, yc, s=1.0, color='black')
    
    if (title == 'R_\lambda'): title = 'Re'
    if (title == '\lambda'): title = 'lambda'
    if (title == '\epsilon'): title = 'epsilon'
    if (title == '\eta'): title = 'eta'
    if (title == '\tau'): title = 'tau'

    fig.savefig(f'{dir}/{ziku}_{title}_' + '{:08d}.png'.format(step))
    fig.savefig(f'{dir}/{ziku}_{title}_' + '{:08d}.svg'.format(step))
    # plt.show()
    plt.close()



if __name__ == '__main__':
    dir = 'v008/data3'
    form = 'bin'  # 'd' or 'bin'
    ziku = 'm'
    ibm = 1  # 1以上：軸メモリD、4：外力点表示(できない)
    mapping(vmax=None, vmin=None)