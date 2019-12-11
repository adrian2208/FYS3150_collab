import numpy as np
import matplotlib.pyplot as plt
import matplotlib
import sys
infile=np.loadtxt("../results/energy_distribution.csv",delimiter=",",skiprows=1)
index=infile[:,0].astype(float)[:-1]
energy=infile[:,1].astype(float)[:-1]
alpha=infile[0,3]
omega=infile[0,2]
fig=plt.figure(figsize=(20,10))
colors=["blue","orange","green","red"]
font = {'family' : 'normal',
        'weight' : 'normal',
        'size'   : 30}
matplotlib.rc('font', **font)
plt.plot(index[:1000],energy[:1000],label=r"$\Psi_{T1}$")
plt.xlabel("Simulation steps")
plt.ylabel(r"<E> [au] at $\omega$=1.0, $\alpha=0.88$")
plt.legend()
plt.tight_layout()
plt.savefig("../plots/Energy_variation.pdf")
plt.show()
