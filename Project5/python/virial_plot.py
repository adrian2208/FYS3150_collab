import numpy as np
import matplotlib.pyplot as plt
import matplotlib
import sys
infile=np.loadtxt("../results/virial.csv",delimiter=",",skiprows=1)
omegas=infile[:int(len(infile)/3),0]
print(omegas)
kinetic=[]
potential=[]
for i in range(3):
    print(int((i)*len(omegas)),int((i+1)*len(omegas)))
    kinetic.append(infile[int((i)*len(omegas)):int((i+1)*len(omegas)),1].astype(float))
    potential.append(infile[int((i)*len(omegas)):int((i+1)*len(omegas)),2].astype(float))
print(kinetic)
labels=["No repulsion",r"$\Psi_{T1}$",r"$\Psi_{T2}$"]
for i in range(3):
    plt.plot(omegas[:-1],(kinetic[i]/potential[i])[:-1],label=labels[i])
plt.legend()
plt.xlabel(r"$\omega$")
plt.ylabel(r"$<T>/<V>$")

plt.show()
