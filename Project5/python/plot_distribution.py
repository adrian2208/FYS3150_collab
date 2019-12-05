import numpy as np
import matplotlib.pyplot as plt
import matplotlib
import sys
infile=np.loadtxt("../results/energy_distribution.csv",delimiter=",",skiprows=1)
index=infile[:,0].astype(float)[:-1]
energy=infile[:,1].astype(float)[:-1]
alpha=infile[0,3]
omega=infile[0,2]
plt.plot(index[:1000],energy[:1000])
plt.show()
