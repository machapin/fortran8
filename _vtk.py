import numpy as np
from pyevtk.hl import gridToVTK
import glob


def vtk(file):
    path = sorted(glob.glob(f'{file}/all_*.d'))[-1]
    data = np.loadtxt(path)

    N = len(set(data[:, 0]))

    x = data[:N, 0]
    y = data[::N, 1]
    y = y[:N]
    z = data[::N*N, 2]

    U = data[:, 3]
    V = data[:, 4]
    W = data[:, 5]
    omega = data[:, 6]
    enstrophy = data[:, 7]
    sendan = data[:, 8]
    energy = data[:, 9]  # 適当
    Q = data[:, 10]

    # print(np.max(omega), np.median(omega), np.min(omega), np.mean(omega), np.std(omega))
    print(np.max(enstrophy), np.median(enstrophy), np.min(enstrophy), np.mean(enstrophy), np.std(enstrophy))
    # print(np.max(Q), np.median(Q), np.min(Q), np.mean(Q), np.std(Q))

    # trc = data[:, 11]
    # C11 = data[:, 12]
    # C12 = data[:, 13]
    # C13 = data[:, 14]
    # C22 = data[:, 15]
    # C23 = data[:, 16]
    # C33 = data[:, 17]

    # gridToVTK("./v10/omega", np.array(x), np.array(y), np.array(z), pointData={"omega": np.array(omega)})
    gridToVTK(file+"/enstrophy", np.array(x), np.array(y), np.array(z), pointData={"enstrophy": np.array(enstrophy)})
    # gridToVTK("./v10/Q", np.array(x), np.array(y), np.array(z), pointData={"Q": np.array(Q)})


if __name__ == '__main__':
    file = 'v170'
    vtk(file)