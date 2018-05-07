import h5py
import matplotlib.pyplot as plt
import numpy as np
from matplotlib import animation


f_mc = h5py.File('output.h5', 'r');
f = h5py.File('output_TRMM.h5', 'r');

#===============================================================================
# alpha eigenvalues
#===============================================================================

alpha = np.array(f['alpha']).transpose()[0]
alpha_adj = np.array(f['alpha_adj']).transpose()[0]
N = len(alpha)
J = 6;
G = N - J;

idx = alpha.argsort()
idx_adj = alpha_adj.argsort()

alpha = alpha[idx]
alpha_adj = alpha_adj[idx_adj]

plt.plot(alpha.real,alpha.imag,'o');
plt.plot(alpha_adj.real,alpha_adj.imag,'x');
plt.xlabel("Re");
plt.ylabel("Im");
plt.grid();
plt.show(); 

#===============================================================================
# eigenvectors
#===============================================================================

phi_mode = np.array(f['phi_mode'])
phi_mode_adj = np.array(f['phi_mode_adj'])

phi_mode = phi_mode[:,idx]
phi_mode_adj = phi_mode_adj[:,idx_adj]

# Inverse speed
v_inv = np.array(f_mc['inverse_speed'])

# Energy bin and lethargy step
energy_grid = np.array(f_mc['TRM_simple/energy'])
energy = energy_grid
du     = np.log(energy[-1] / energy[-2])
energy = np.array(energy)*1E-6;
energy = (energy[1:] + energy[:-1])/2;

#===============================================================================
# Verification with TDMC
#===============================================================================

# Initial condition
phi_initial = np.zeros(N)
phi_initial[-1-J] = 13831.5926439 * 14.1**0.5 * 1000 * 100.0
for g in range(G):
    phi_initial[g] = phi_initial[g] * v_inv[g]

# Expansion coefficients
A = np.zeros(N,dtype=complex)
for i in range(N):
    num = complex(0,0)
    gamma = complex(0,0)
    for g in range(G):
        num = num + phi_mode_adj[g][i] * phi_initial[g]
        gamma = gamma + phi_mode_adj[g][i] * v_inv[g] * phi_mode[g][i]
    for g in range(G,G+J):
        num = num + phi_mode_adj[g][i] * phi_initial[g]
        gamma = gamma + phi_mode_adj[g][i] * phi_mode[g][i]
    A[i] = num / gamma

# Time grid
time = [ 3e-8, 15e-8, 4e-6, 1e-4 ]

# Verification solution
for t in time:
    phi_ver = np.zeros(N,dtype=complex)
    for g in range(N):
        for n in range(N):
            phi_ver[g] = phi_ver[g] + A[n] * phi_mode[g][n] \
                                      * np.e**(alpha[n] * t)
    phi_ver = phi_ver / du
    plt.loglog(energy,phi_ver[:G]);

# TDMC solution
f_td = h5py.File('TDMC.h5', 'r')
flux_td = np.array(f_td['spectrum/flux/mean']);
energy_td = np.array(f_td['spectrum/energy']);
du_td     = np.log(energy_td[-1] / energy_td[-2])

for i in range(len(flux_td)):
    flux_td[i][0] = flux_td[i][0] / du_td
energy_td = energy_td*1E-6;
energy_td = (energy_td[1:] + energy_td[:-1])/2;
plt.scatter(energy_td,flux_td[0][0], s=80, facecolors='none', edgecolors='blue');
plt.scatter(energy_td,flux_td[1][0], s=80, facecolors='none', edgecolors='orange');
plt.scatter(energy_td,flux_td[2][0], s=80, facecolors='none', edgecolors='g');
plt.scatter(energy_td,flux_td[3][0], s=80, facecolors='none', edgecolors='r');
plt.xlabel("Energy, MeV");
plt.ylabel("Scalar flux");
plt.grid();
plt.xlim(1E-9,20)
plt.ylim(1,1E10)
plt.show();

#===============================================================================
# animation
#===============================================================================

energy = np.array(f_mc['TRM_simple/energy'])
energy = np.array(energy)*1E-6;
energy = (energy[1:] + energy[:-1])/2;

time = np.logspace(-10,-2,500)

fig = plt.figure()
ax = plt.axes(xlim=(1E-9, 20), ylim=(1, 2E10))
ax.set_xscale('log')
ax.set_yscale('log')
line, = ax.plot([], [], '-', lw=2)
time_text = ax.text(0.02, 0.95, '', transform=ax.transAxes)
plt.xlabel("Energy, MeV");
plt.ylabel("Scalar flux");

def init():
    line.set_data([], [])
    time_text.set_text('')
    return time_text, line

def animate(i):
    phi = np.zeros(N,dtype=complex)
    for g in range(N):
        for n in range(N):
            phi[g] = phi[g] + A[n] * phi_mode[g][n] * np.e**(alpha[n] * time[i])
    phi = phi / du
    line.set_data(energy, phi[:G])
    time_text.set_text('time = %.9f s' %time[i])
    return time_text, line

inter = 5000 / len(time)
anim = animation.FuncAnimation(fig, animate, init_func=init,
                               frames=len(time), interval=inter, blit=True)

plt.show()
