import numpy as np
import matplotlib.pyplot as plt
import matplotlib.patches as patches
from tqdm import tqdm
import glob

# export XDG_RUNTIME_DIR=/tmp/runtime-suzuki
plt.rcParams['font.family'] = 'Times New Roman'
plt.rcParams['mathtext.fontset'] = 'stix'
plt.rcParams['font.size'] = 15

def vector():
    if (form == 'd'):
        # for filename in tqdm(glob.glob(f'{dir}/{ziku}_*.d')):
            filename = sorted(glob.glob(f'{dir}/{ziku}_*.d'))[-1]
            data = np.loadtxt(filename)
            step = int(filename[-10:-2])
            vector_fig(data, step)

    if (form == 'bin'):
        # for filename in tqdm(glob.glob(f'{dir}/{ziku}_*.bin')):
            filename = sorted(glob.glob(f'{dir}/{ziku}_*.bin'))[-1]
            with open(filename, mode='rb') as f:
                data = np.fromfile(f, dtype=np.float64)
            data = data.reshape(-1, 18)
            step = int(filename[-12:-4])
            vector_fig(data, step)


def vector_fig(data, step):
        N = len(set(data[:, 0]))
        
        if ziku=='z': x = data[:, 0]; y = data[:, 1]; U = data[:, 3]; V = data[:, 4]; xlabel = '$x$'; ylabel = '$y$'  # テイラー渦
        if ziku=='y': x = data[:, 0]; y = data[:, 2]; U = data[:, 3]; V = data[:, 5]; xlabel = '$x$'; ylabel = '$z$'  # クエット流れ
            
        s = 2  # 矢印のサイズ
        l = 8  # 間引き
        U = U * s
        V = V * s
        fig = plt.figure(figsize=(5, 5))
        ax = fig.add_subplot(1, 1, 1)
        # ax = fig.add_axes([0.13, 0.15, 0.4, 0.8])
        
        if (l == 1):
            ax.quiver(x, y, U, V, color='black', angles='xy', scale_units='xy', scale=5.0)
        else:
            x_quiver = []
            y_quiver = []
            U_quiver = []
            V_quiver = []
            for i in range(int(N/l)):
                for j in range(int(N/l)):
                    index = i * N * l + j * l + int(N*l/2) + int(l/2)
                    # index = i * N * l + j * l
                    x_quiver.append(x[index])
                    y_quiver.append(y[index])
                    U_quiver.append(U[index])
                    V_quiver.append(V[index])

            ax.quiver(x_quiver, y_quiver, U_quiver, V_quiver, color='black', angles='xy', scale_units='xy', scale=5.0)

        xmax = 2*np.pi
        if int(max(x)) == int(xmax):
            ax.set_xticks(np.linspace(0, 2*np.pi, 3))
            ax.set_xticklabels(['$0$', '$\pi$', '$2\pi$'], size=15)
            ax.set_xlabel(r'$x$', size=18)
            ax.set_yticks(np.linspace(0, 2*np.pi, 3))
            ax.set_yticklabels(['$0$', '$\pi$', '$2\pi$'], size=15)
            # ax.set_ylabel(r'$y$', size=18, rotation=0)
            if(ziku=='z'): ax.set_ylabel(r'$y$', size=18)
            if(ziku=='y'): ax.set_ylabel(r'$z$', size=18)
            ax.set_xlim([0.0, xmax])
            ax.set_ylim([0.0, xmax])

        else:
            xmax = max(x) + min(x)
            ax.set_xticks(np.linspace(0, xmax, 3))
            ax.set_xticklabels(['$0$', '$2D$', '$4D$'], size=15)
            ax.set_xlabel(r'$x$', size=18)
            ax.set_yticks(np.linspace(0, xmax, 3))
            ax.set_yticklabels(['$0$', '$2D$', '$4D$'], size=15)
            # ax.set_ylabel(r'$y$', size=18, rotation=0)
            if(ziku=='z'): ax.set_ylabel(r'$y$', size=18)
            if(ziku=='y'): ax.set_ylabel(r'$z$', size=18)
            ax.set_xlim([0.0, xmax])
            ax.set_ylim([0.0, xmax])


        # xmax = max(x) + min(x)
        # ax.set_xticks(np.linspace(0, xmax, 4))
        # ax.set_xticklabels(['$0$', '$D$', '$2D$', '$3D$'], size=15)
        # ax.set_xlabel(r'$x$', size=18)
        # ax.set_yticks(np.linspace(0, xmax, 4))
        # ax.set_yticklabels(['$0$', '$D$', '$2D$', '$3D$'], size=15)
        # # ax.set_ylabel(r'$y$', size=18, rotation=0)
        # if(ziku=='z'): ax.set_ylabel(r'$y$', size=18)
        # if(ziku=='y'): ax.set_ylabel(r'$z$', size=18)
        # ax.set_xlim([0.0, xmax])
        # ax.set_ylim([0.0, xmax])


        # if 0:
        #     ax.set_xticks(np.linspace(xmax*3/8, xmax*5/8, 3))
        #     # ax.set_xticklabels(['$\\frac{3\pi}{4}$', '$\pi$', '$\\frac{5\pi}{4}$'], size=15)
        #     ax.set_xticklabels(['$0$', '$D$', '$2D$'], size=15)
        #     ax.set_xlabel(r'$x$', size=18)
        #     ax.set_yticks(np.linspace(xmax*2/8, xmax*6/8, 3))
        #     ax.set_yticklabels(['$0$', '$2D$', '$4D$'], size=15)
        #     # ax.set_ylabel(r'$y$', size=18, rotation=0)
        #     ax.set_ylabel(r'$z$', size=18)
        #     ax.set_xlim([xmax*3/8, xmax*5/8])
        #     ax.set_ylim([xmax*2/8, xmax*6/8])
        #     ziku = 'tyouhoukei_'+ziku


        # ax.set_title(f'step:{step}')

        # if (ibm == 1):
        #     r = patches.Rectangle(xy=(0, 0), width=xmax, height=xmax/4, fc='black', ec='black', alpha=0.5)
        #     ax.add_patch(r)
        #     r = patches.Rectangle(xy=(0, xmax*3/4), width=xmax, height=xmax/4, fc='black', ec='black', alpha=0.5)
        #     ax.add_patch(r)

        # if (ibm == 2):
        #     for i in range(1, ibm+1):
        #         data = np.loadtxt('ibm'+str(ibm)+'/{:08d}.d'.format(i))
        #         xc = data[:, 0]
        #         yc = data[:, 1]
        #         # ax.plot(xc, yc, lw=1.0, color='black', linestyle="dashed")
        #         ax.plot(xc, yc, lw=1.0, color='black')
        #         # ax.scatter(xc, yc, s=1.0, color='black')
        
        # if (ibm == 5):
        #     for i in range(1, 66):
        #         data = np.loadtxt('ibm'+str(ibm)+'/{:08d}.d'.format(i))
        #         xc = data[:, 0]
        #         yc = data[:, 1]
        #         # ax.plot(xc, yc, lw=1.0, color='black', linestyle="dashed")
        #         ax.plot(xc, yc, lw=1.0, color='black')
        #         # ax.scatter(xc, yc, s=1.0, color='black')
        
        # if (ibm == 6):
        #     for i in range(1, 33):
        #         data = np.loadtxt('ibm'+str(ibm)+'/{:08d}.d'.format(i))
        #         xc = data[:, 0]
        #         yc = data[:, 1]
        #         # ax.plot(xc, yc, lw=1.0, color='black', linestyle="dashed")
        #         ax.plot(xc, yc, lw=1.0, color='black')
        #         # ax.scatter(xc, yc, s=1.0, color='black')

        if (ziku=='z' and ibm==4):
            for i in range(1, 5):
                data = np.loadtxt('ibm/ibm13/{:08d}.d'.format(i))
                xc = data[:, 0]
                yc = data[:, 1]
                # ax.plot(xc, yc, lw=1.0, color='black', linestyle="dashed")
                # ax.plot(xc, yc, lw=1.0, color='black')
                ax.scatter(xc, yc, s=0.3, lw=0.3, color='red')
                
        # if (ziku=='z' and ibm==22):  # 円柱の外力点表示
        #     data = np.loadtxt('./ibm/ibm22/00000001.d')
        #     xc = data[:, 0]
        #     yc = data[:, 1]
        #     # ax.plot(xc, yc, lw=1.0, color='black', linestyle="dashed")
        #     # ax.plot(xc, yc, lw=1.0, color='black')
        #     ax.scatter(xc, yc, s=10.0, lw=1.0, color='red')

        # if (ziku=='z' and ibm==11):  # 円柱の外力点表示
        #     data = np.loadtxt('./ibm/ibm11/00000001.d')
        #     xc = data[:, 0]
        #     yc = data[:, 1]
        #     # ax.plot(xc, yc, lw=1.0, color='black', linestyle="dashed")
        #     # ax.plot(xc, yc, lw=1.0, color='black')
        #     ax.scatter(xc, yc, s=10.0, lw=1.0, color='red')

        fig.savefig(f'{dir}/{ziku}_f_'+'{:08d}.png'.format(step))
        # fig.savefig(f'{dir}/{ziku}_f_'+'{:08d}.svg'.format(step))
        # fig.savefig(f'{dir}/{ziku}_method_{N}.png'.format(step))
        # fig.savefig(f'{dir}/{ziku}_method_{N}.svg'.format(step))
        # fig.savefig(f'{dir}/method.png')
        # fig.savefig(f'{dir}/method.svg')
        plt.show()
        plt.close()


if __name__ == '__main__':
    dir = 'v008/data3'
    form = 'bin'  # 'd' or 'bin'
    ziku = 'y'  # 'z' or 'y'
    ibm = 0
    vector()