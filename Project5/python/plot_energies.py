import numpy as np
import matplotlib.pyplot as plt
infile=np.loadtxt("../results/function1.csv",delimiter=",",skiprows=1)
alpha=infile[:int(len(infile)/3),1]
energy_001=infile[:int(len(infile)/3),2].astype(float); sigma_001=infile[:int(len(infile)/3),3].astype(float)

energy_05=infile[int(len(infile)/3):int(2*len(infile)/3),2].astype(float);sigma_05=infile[int(len(infile)/3):int(2*len(infile)/3),3].astype(float)

energy_1=infile[int(2*len(infile)/3):,2].astype(float);sigma_1=infile[int(2*len(infile)/3):,3].astype(float)
plt.plot(alpha,energy_1,label="1")
plt.plot(alpha,energy_05,label="0.5")
plt.plot(alpha,energy_001,label="0.01")
plt.legend()
plt.xlabel(r"$\alpha$")
plt.ylabel(r"$E_0$")
plt.yscale("log")
plt.show()
