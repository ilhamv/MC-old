import h5py
import matplotlib.pyplot as plt
import numpy as np

f = h5py.File('output.h5', 'r');

#===============================================================================
# k search and label
#===============================================================================

k = f['ksearch/k_cycle'];
N = np.array(f['summary/Nsample']);
plt.plot(k,'x');
plt.xlabel("cycle#");
plt.ylabel("k");
plt.ylim(0)
plt.show(); 

H = f['ksearch/H_cycle'];
N = np.array(f['summary/Nsample']);
plt.plot(H,'x');
plt.xlabel("cycle#");
plt.ylabel("H");
plt.ylim(0)
plt.show(); 
