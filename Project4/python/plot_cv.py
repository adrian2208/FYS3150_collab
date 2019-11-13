import numpy as np
import matplotlib.pyplot as plt
import matplotlib
data=np.loadtxt("../results/results_calculations_general.csv",skiprows=1,dtype=float,delimiter=",")

temperatures=data[:,0][0:int(len(data[:,0])/4)]
lens=[40,60,80,100]
print(temperatures)
cv=[[] for lenny in lens]
energy=[[] for lenny in lens]
xi=[[] for lenny in lens]
magnetisation=[[] for lenny in lens]

print(data[:,0])
print(data[0,:])
for i in range(len(lens)):
    for temperature in temperatures:
        relevant=data[np.logical_and(abs(data[:,0]-temperature)<1e-10,data[:,1]==lens[i])][0][-2]
        cv[i].append(relevant)
        relevant=data[np.logical_and(abs(data[:,0]-temperature)<1e-10,data[:,1]==lens[i])][0][-1]
        xi[i].append(relevant)
        relevant=data[np.logical_and(abs(data[:,0]-temperature)<1e-10,data[:,1]==lens[i])][0][3]
        energy[i].append(relevant)
        relevant=data[np.logical_and(abs(data[:,0]-temperature)<1e-10,data[:,1]==lens[i])][0][7]
        magnetisation[i].append(relevant)
font = {'family' : 'normal',
        'weight' : 'normal',
        'size'   : 22}

matplotlib.rc('font', **font)
plt.figure(figsize=(20,20))

plt.subplot(221)
for i in range(len(cv)):
    plt.plot(temperatures,cv[i],label="L=%s"%str(lens[i]))
plt.xlabel("Temperature [kT/J]")
plt.ylabel(r"Heat capacity $C_V$ []")
plt.legend()
plt.subplot(222)

for i in range(len(xi)):
    plt.plot(temperatures,xi[i],label="L=%s"%str(lens[i]))
plt.xlabel("Temperature [kT/J]")
plt.ylabel(r"Susceptibility $\chi$")
plt.legend()
plt.subplot(223)

for i in range(len(energy)):
    plt.plot(temperatures,energy[i],label="L=%s"%str(lens[i]))
plt.xlabel("Temperature [kT/J]")
plt.ylabel("Energy [J]")
plt.legend()
plt.subplot(224)

for i in range(len(magnetisation)):
    plt.plot(temperatures,magnetisation[i],label="L=%s"%str(lens[i]))
plt.xlabel("Temperature [kT/J]")
plt.ylabel("Absolute Magnetisation")
plt.legend()
plt.savefig("")
plt.show()
