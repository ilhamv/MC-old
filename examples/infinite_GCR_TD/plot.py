import h5py
import matplotlib.pyplot as plt
import numpy as np

f = h5py.File('output.h5', 'r');

flux = f['spectrum/flux/mean'];
energy = f['spectrum/energy'];

energy = np.array(energy)*1E-6;
energy = (energy[1:] + energy[:-1])/2;


plt.loglog(energy,flux[0][0]);
plt.loglog(energy,flux[1][0]);
plt.loglog(energy,flux[2][0]);
plt.loglog(energy,flux[3][0]);
plt.ylim(1,1E9);
plt.xlim(1E-9,2E1);
plt.xlabel("Energy, MeV");
plt.ylabel("Scalar flux");
plt.show(); 
