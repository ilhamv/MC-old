import h5py
import matplotlib.pyplot as plt
import numpy as np

f = h5py.File('output.h5', 'r');

k = f['ksearch/k_cycle'];
N = np.array(f['summary/Nsample']);
NG = int(N**0.5);
plt.plot(k,'x');
plt.xlabel("cycle#");
plt.ylabel("k");
plt.ylim(0)
plt.title('G=%s, sample=%s'%(NG,N))
plt.show(); 

flux = f['MG/flux/mean'];
energy = f['MG/energy'];

energy = np.array(energy)*1E-6;
energy = (energy[1:] + energy[:-1])/2;

plt.loglog(energy,flux[0]);
#plt.ylim(1,1E9);
#plt.xlim(1E-9,2E1);
plt.xlabel("Energy, MeV");
plt.ylabel("Scalar flux");
plt.grid();
plt.title('G=%s, sample=%s'%(NG,N))
plt.show(); 

x = f['TRMM/alpha/real'];
y = f['TRMM/alpha/imag'];

plt.plot(x,y,'o');
plt.xlabel("Re");
plt.ylabel("Im");
plt.grid();
plt.title('G=%s, sample=%s'%(NG,N))
plt.show(); 

flux = f['TRMM/flux'];
energy = f['MG/energy'];

f_td = h5py.File('TDMC.h5', 'r');
flux_td = f_td['spectrum/flux/mean'];
energy_td = f_td['spectrum/energy'];

energy_td = np.array(energy_td)*1E-6;
energy_td = (energy_td[1:] + energy_td[:-1])/2;

plt.scatter(energy_td,flux_td[0][0], s=80, facecolors='none', edgecolors='blue');
plt.scatter(energy_td,flux_td[1][0], s=80, facecolors='none', edgecolors='orange');
plt.scatter(energy_td,flux_td[2][0], s=80, facecolors='none', edgecolors='g');
plt.scatter(energy_td,flux_td[3][0], s=80, facecolors='none', edgecolors='r');



energy = np.array(energy)*1E-6;
energy = (energy[1:] + energy[:-1])/2;

#plt.loglog(energy,flux[0]);
plt.loglog(energy,flux[1]);
plt.loglog(energy,flux[2]);
plt.loglog(energy,flux[3]);
plt.loglog(energy,flux[4]);
plt.ylim(1,1E12);
plt.xlim(1E-9,2E1);
plt.xlabel("Energy, MeV");
plt.ylabel("Scalar flux");
plt.grid();
plt.title('G=%s, sample=%s'%(NG,N))
plt.show();
