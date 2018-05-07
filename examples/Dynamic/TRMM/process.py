import h5py
import matplotlib.pyplot as plt
import numpy as np
from matplotlib import animation


f_mc = h5py.File('output.h5', 'r');
f = []
f.append(h5py.File('output_TRMM1.h5', 'r'))
f.append(h5py.File('output_TRMM2.h5', 'r'))
f.append(h5py.File('output_TRMM1.h5', 'r'))

t_switch = [0., 10., 50.]
#t_switch = [0.]

NC = len(t_switch)

# Inverse speed
v_inv = np.array(f_mc['inverse_speed'])

# Energy bin and lethargy step
energy_grid = np.array(f_mc['TRM_simple/energy'])
energy = energy_grid
du     = np.log(energy[-1] / energy[-2])
energy = np.array(energy)*1E-6;
energy = (energy[1:] + energy[:-1])/2;


#===============================================================================
# eigenpairs
#===============================================================================

alpha_case = []
alpha_adj_case = []

phi_mode_case = []
phi_mode_adj_case = []

idx_case = []
idx_adj_case = []

for i in range(NC):
    alpha = np.array(f[i]['alpha']).transpose()[0]
    alpha_adj = np.array(f[i]['alpha_adj']).transpose()[0]
    N = len(alpha)
    J = 6;
    G = N - J;

    idx = alpha.argsort()
    idx_adj = alpha_adj.argsort()

    alpha = alpha[idx]
    alpha_adj = alpha_adj[idx_adj]

    phi_mode = np.array(f[i]['phi_mode'])
    phi_mode_adj = np.array(f[i]['phi_mode_adj'])

    phi_mode = phi_mode[:,idx]
    phi_mode_adj = phi_mode_adj[:,idx_adj]

    plt.plot(alpha.real,alpha.imag,'o');
    plt.plot(alpha_adj.real,alpha_adj.imag,'x');
    plt.xlabel("Re");
    plt.ylabel("Im");
    plt.grid();
    plt.show(); 
    
    alpha_case.append(alpha)
    alpha_adj_case.append(alpha_adj)
    phi_mode_case.append(phi_mode)
    phi_mode_adj_case.append(phi_mode_adj)

#===============================================================================
# ICs and Expansion Coeffisients
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
A_case = []
A = np.zeros(N,dtype=complex)
for i in range(N):
    num = complex(0,0)
    gamma = complex(0,0)
    for g in range(G):
        num = num + phi_mode_adj_case[0][g][i] * phi_initial[g] * v_inv[g]
        gamma = gamma + phi_mode_adj_case[0][g][i] * v_inv[g] * phi_mode_case[0][g][i]
    for g in range(G,G+J):
        num = num + phi_mode_adj_case[0][g][i] * phi_initial[g]
        gamma = gamma + phi_mode_adj_case[0][g][i] * phi_mode_case[0][g][i]
    A[i] = num / gamma
A_case.append(A)

phi = np.zeros(N,dtype=complex)
for g in range(N):
    for n in range(N):
        phi[g] = phi[g] + A_case[0][n] * phi_mode_case[0][g][n]
'''
plt.loglog(energy,phi_initial[:G])
plt.loglog(energy,phi[:G],'o')
plt.show()
'''
for j in range(1,NC):
    phi_initial = np.zeros(N,dtype=complex)
    for g in range(N):
        for n in range(N):
            phi_initial[g] = phi_initial[g] + A_case[j-1][n] * phi_mode_case[j-1][g][n] * np.e**(alpha_case[j-1][n] * (t_switch[j] - t_switch[j-1]))

    A = np.zeros(N,dtype=complex)
    for i in range(N):
        num = complex(0,0)
        gamma = complex(0,0)
        for g in range(G):
            num = num + phi_mode_adj_case[j][g][i] * phi_initial[g] * v_inv[g]
            gamma = gamma + phi_mode_adj_case[j][g][i] * v_inv[g] * phi_mode_case[j][g][i]
        for g in range(G,G+J):
            num = num + phi_mode_adj_case[j][g][i] * phi_initial[g]
            gamma = gamma + phi_mode_adj_case[j][g][i] * phi_mode_case[j][g][i]
        A[i] = num / gamma
    A_case.append(A)
    phi = np.zeros(N,dtype=complex)
    for g in range(N):
        for n in range(N):
            phi[g] = phi[g] + A_case[j][n] * phi_mode_case[j][g][n]
    '''
    plt.loglog(energy,phi_initial[:G])
    plt.loglog(energy,phi[:G],'o')
    plt.show()
    '''


#===============================================================================
# animation
#===============================================================================
Nd = 1000
time = np.linspace(0.,100.,Nd)

fig = plt.figure()
ax = plt.axes(xlim=(1E-9, 20), ylim=(1E-4, 1E6))
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
P = []
t = []
def animate(i):
    phi = np.zeros(N,dtype=complex)
    for j in reversed(range(NC)):
        if time[i] >= t_switch[j]:
            A = A_case[j]
            alpha = alpha_case[j]
            alpha_adj = alpha_adj_case[j]
            phi_mode = phi_mode_case[j]
            phi_mode_adj = phi_mode_adj_case[j]
            case = j
            break
    for g in range(N):
        for n in range(N):
            phi[g] = phi[g] + A[n] * phi_mode[g][n] * np.e**(alpha[n] * (time[i]-t_switch[case]))
    phi = phi / du
    line.set_data(energy, phi[:G])
    time_text.set_text('time = %.2f s' %time[i])

    Pw = 0
    for g in range(G):
        Pw = Pw + phi[g]
    P.append(Pw)
    t.append(time[i])
    return time_text, line

inter = 1000 / len(time)
anim = animation.FuncAnimation(fig, animate, init_func=init,
                               frames=len(time), interval=inter, blit=True)

plt.show()
plt.plot(t[:Nd],P[:Nd],"-*")
plt.xlabel('time, s')
plt.ylabel('total flux')
plt.show()
