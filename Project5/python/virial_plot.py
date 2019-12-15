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
fig=plt.figure(figsize=(20,20))
colors=["blue","orange","green","red"]
font = {'family' : 'normal',
        'weight' : 'normal',
        'size'   : 40}
matplotlib.rc('font', **font)
for i in range(3):
    plt.plot(omegas[:-1],(kinetic[i]/potential[i])[:-1],label=labels[i])
plt.legend()
plt.xlabel(r"$\omega$")
plt.ylabel(r"$<T>/<V>$")
plt.tight_layout()
plt.xlim(-0.05,1.05)
plt.ylim(0, 1.05)
plt.savefig("../plots/Virial.pdf")
plt.show()
