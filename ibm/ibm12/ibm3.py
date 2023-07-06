import numpy as np
import matplotlib.pyplot as plt
plt.rcParams['font.family'] = 'Times New Roman'
plt.rcParams['mathtext.fontset'] = 'stix'
plt.rcParams['font.size'] = 15


fig = plt.figure(figsize=(5, 5))
ax = fig.add_axes([0.15, 0.15, 0.8, 0.8])
# L_C = 8*0.03/2.0/np.pi*1000  # [mm]
# U_C = 0.06*np.pi*1000  # [mm/s]

for i in range(1, 5):
    data = np.loadtxt('./ibm/ibm12/{:08d}.d'.format(i))
    xc = data[:, 0]
    yc = data[:, 1]
    ax.plot(xc, yc, lw=1.0, color='blue')
    # ax.scatter(xc, yc, s=1.0, color='blue')

xmax = 2*np.pi

ax.set_xlim([0, xmax])
ax.set_ylim([0, xmax])

ax.set_xticks(np.linspace(0, 2*np.pi, 3))
ax.set_xticklabels(['0', '$2D$', '$4D$'])
ax.set_xlabel('$x$', size=20)
ax.set_yticks(np.linspace(0, 2*np.pi, 3))
ax.set_yticklabels(['0', '$2D$', '$4D$'])
ax.set_ylabel('$y$', size=20)

ax.annotate('', xy=[0, xmax*3/4], xytext=[xmax/8, xmax*3/4], arrowprops=dict(arrowstyle='<|-|>', connectionstyle='arc3', facecolor='blue', edgecolor='blue'))
ax.annotate('', xy=[xmax/8, xmax*3/4], xytext=[xmax*3/8, xmax*3/4], arrowprops=dict(arrowstyle='<|-|>', connectionstyle='arc3', facecolor='blue', edgecolor='blue'))
ax.annotate('', xy=[xmax*3/8, xmax*3/4], xytext=[xmax*5/8, xmax*3/4], arrowprops=dict(arrowstyle='<|-|>', connectionstyle='arc3', facecolor='blue', edgecolor='blue'))
ax.annotate('', xy=[xmax*5/8, xmax*3/4], xytext=[xmax*7/8, xmax*3/4], arrowprops=dict(arrowstyle='<|-|>', connectionstyle='arc3', facecolor='blue', edgecolor='blue'))
ax.annotate('', xy=[xmax*7/8, xmax*3/4], xytext=[xmax, xmax*3/4], arrowprops=dict(arrowstyle='<|-|>', connectionstyle='arc3', facecolor='blue', edgecolor='blue'))

ax.annotate('', xy=[xmax/4, 0], xytext=[xmax/4, xmax/8], arrowprops=dict(arrowstyle='<|-|>', connectionstyle='arc3', facecolor='blue', edgecolor='blue'))
ax.annotate('', xy=[xmax/4, xmax/8], xytext=[xmax/4, xmax*3/8], arrowprops=dict(arrowstyle='<|-|>', connectionstyle='arc3', facecolor='blue', edgecolor='blue'))
ax.annotate('', xy=[xmax/4, xmax*3/8], xytext=[xmax/4, xmax*5/8], arrowprops=dict(arrowstyle='<|-|>', connectionstyle='arc3', facecolor='blue', edgecolor='blue'))
ax.annotate('', xy=[xmax/4, xmax*5/8], xytext=[xmax/4, xmax*7/8], arrowprops=dict(arrowstyle='<|-|>', connectionstyle='arc3', facecolor='blue', edgecolor='blue'))
ax.annotate('', xy=[xmax/4, xmax*7/8], xytext=[xmax/4, xmax], arrowprops=dict(arrowstyle='<|-|>', connectionstyle='arc3', facecolor='blue', edgecolor='blue'))


fig.text(0.18, 0.77, '$D/2$', size=15)
fig.text(0.37, 0.77, '$D$', size=15)
fig.text(0.53, 0.77, '$D$', size=15)
fig.text(0.74, 0.77, '$D$', size=15)
fig.text(0.87, 0.77, '$D/2$', size=15)

fig.text(0.37, 0.20, '$D/2$', size=15)
fig.text(0.37, 0.35, '$D$', size=15)
fig.text(0.37, 0.54, '$D$', size=15)
fig.text(0.37, 0.88, '$D/2$', size=15)


# ax.set_aspect('equal')
# plt.tight_layout()
fig.savefig('./ibm/ibm12/ibm.png', bbox_inches='tight')
# fig.savefig('./ibm/ibm12/ibm.pdf', bbox_inches='tight')
fig.savefig('./ibm/ibm12/ibm.svg', bbox_inches='tight')
plt.show()