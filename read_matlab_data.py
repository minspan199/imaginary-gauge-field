# -*- coding: utf-8 -*-
"""
Created on Wed May  9 16:15:52 2018

@author: Mingsen
"""
import numpy as np
import scipy.io
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D

mat = scipy.io.loadmat('data1.mat')
Y = mat['Y']
T = mat['T']

Nd = 16
N_span = len(T)
fig_1 = plt.figure()
ax1 = fig_1.gca(projection='3d')
L = np.linspace(0, Nd - 1, Nd).astype(int)

for i in reversed(L):
    ax1.plot(T, i * np.ones(N_span), np.angle(Y[:, i]), linewidth=1.5)

ax1.ticklabel_format(style='sci', axis='x', scilimits=(0, 0))
ax1.set_facecolor((1, 1, 1))
# make the panes transparent
ax1.xaxis.set_pane_color((1.0, 1.0, 1.0, 0.0))
ax1.yaxis.set_pane_color((1.0, 1.0, 1.0, 0.0))
ax1.zaxis.set_pane_color((1.0, 1.0, 1.0, 0.0))
# make the grid lines transparent
# ax1.xaxis._axinfo["grid"]['color'] =  (1,1,1,0)
# ax1.yaxis._axinfo["grid"]['color'] =  (1,1,1,0)
# ax1.zaxis._axinfo["grid"]['color'] =  (1,1,1,0)
ax1.set_zlim(-3.5, 3.5)
ax1.set_title('Eigenstates evolution of the laser array')
ax1.set_xlabel('Time')
ax1.set_ylabel('Site number')
ax1.set_zlabel('Phase')

fig_2 = plt.figure()
ax2 = fig_2.gca(projection='3d')
Y[abs(Y) > 0.5] = 0.5

for i in L:
    ax2.plot(T, i * np.ones(N_span), np.abs(Y[:, i]), linewidth=1.5)

ax2.ticklabel_format(style='sci', axis='x', scilimits=(0, 0))
ax2.xaxis.set_pane_color((1.0, 1.0, 1.0, 0.0))
ax2.yaxis.set_pane_color((1.0, 1.0, 1.0, 0.0))
ax2.zaxis.set_pane_color((1.0, 1.0, 1.0, 0.0))
ax2.set_zlim(0, 0.5)
ax2.set_title('Eigenstates evolution of the laser array')
ax2.set_xlabel('Time')
ax2.set_ylabel('Site number')
ax2.set_zlabel('Amplitudes')

fig_3 = plt.figure()
Y[abs(Y) > 0.5] = 0.5

for i in L:
    plt.plot(T, np.abs(Y[:, i]), linewidth=1.5)

plt.ticklabel_format(style='sci', axis='x', scilimits=(0, 0), useMathText=True)
plt.ylim(0, 0.5)
plt.title('Eigenstates evolution of the laser array')
plt.xlabel('Normalized time')
plt.ylabel('Amplitudes')
plt.tight_layout()
plt.show();
