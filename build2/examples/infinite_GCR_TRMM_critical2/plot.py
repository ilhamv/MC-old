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
psi_initial = np.array(f_mc['psi_initial'])
C_initial = np.array(f_mc['C_initial'])
for g in range(G):
    phi_initial[g] = psi_initial[g]
for j in range(J):
    phi_initial[G+j] = C_initial[j]

# Expansion coefficients
A = np.zeros(N,dtype=complex)
for i in range(N):
    num = complex(0,0)
    gamma = complex(0,0)
    for g in range(G):
        num = num + phi_mode_adj[g][i] * phi_initial[g] * v_inv[g] 
        gamma = gamma + phi_mode_adj[g][i] * v_inv[g] * phi_mode[g][i]
    for g in range(G,G+J):
        num = num + phi_mode_adj[g][i] * phi_initial[g]
        gamma = gamma + phi_mode_adj[g][i] * phi_mode[g][i]
    A[i] = num / gamma

val=[]
v_i=[]
v_j=[]
for i in range(N):
    for j in range(N):
        if i != j:
            gamma = complex(0,0)
            for g in range(G):
                gamma = gamma + phi_mode_adj[g][i] * v_inv[g] * phi_mode[g][j]
            for g in range(G,G+J):
                gamma = gamma + phi_mode_adj[g][i] * phi_mode[g][j]
            val.append(gamma)
            v_i.append(i)
            v_j.append(j)
val=np.array(val)
v_i=np.array(v_i)
v_j=np.array(v_j)
idx = val.argsort()
val = val[idx]
v_i = v_i[idx]
v_j = v_j[idx]

for i in range(len(val)):
    print(val[i],v_i[i],v_j[i])

# Initial expansion
phi = np.zeros(N,dtype=complex)
for g in range(N):
    for n in range(N):
        phi[g] = phi[g] + A[n] * phi_mode[g][n]
plt.loglog(energy, phi[:G])
plt.loglog(energy, phi_initial[:G])
plt.show()

#===============================================================================
# animation
#===============================================================================

energy = np.array(f_mc['TRM_simple/energy'])
energy = np.array(energy)*1E-6;
energy = (energy[1:] + energy[:-1])/2;

time = np.logspace(-10,-2,500)

fig = plt.figure()
ax = plt.axes(xlim=(1E-9, 20), ylim=(1, 2E12))
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
