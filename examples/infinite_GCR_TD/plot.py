import h5py
import matplotlib.pyplot as plt
import numpy as np
from matplotlib import animation

f = h5py.File('output.h5', 'r');

flux = f['spectrum/flux/mean'];
energy = f['spectrum/energy'];
time = np.array(f['spectrum/time']);
print(time)

energy = np.array(energy)*1E-6;
energy = (energy[1:] + energy[:-1])/2;


for i in range(0, len(flux)):
    plt.scatter(energy,flux[i][0], s=80, facecolors='none', edgecolors='k');
    ax = plt.gca()
    ax.set_xscale('log')
    ax.set_yscale('log')
    plt.ylim(1,1E9);
    plt.xlim(1E-9,2E1);
    plt.xlabel("Energy, MeV");
    plt.ylabel("Scalar flux");
    plt.show(); 

fig = plt.figure()
ax = plt.axes(xlim=(1E-9, 20), ylim=(1, 1E9))
ax.set_xscale('log')
ax.set_yscale('log')
line, = ax.plot([], [], 'o', lw=2)
time_text = ax.text(0.02, 0.95, '', transform=ax.transAxes)
plt.xlabel("Energy, MeV");
plt.ylabel("Scalar flux");

def init():
    line.set_data([], [])
    time_text.set_text('')
    return time_text, line

def animate(i):
    line.set_data(energy, flux[i][0])
    time_text.set_text('time = %.9f s' %time[i])
    return time_text, line

inter = 1000 / len(flux)
anim = animation.FuncAnimation(fig, animate, init_func=init,
                               frames=len(flux), interval=inter, blit=True)

plt.show()
