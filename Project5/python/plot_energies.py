import numpy as np
import matplotlib.pyplot as plt
import matplotlib
import sys
infile=np.loadtxt("../results/function1.csv",delimiter=",",skiprows=1)
alpha=infile[:int(len(infile)/3),1]
energies=[]
sigmas=[]
distances=[]
omegas=[0.01,0.5,1]
for i in range(3):
    energies.append(infile[int((i)*len(infile)/3):int((i+1)*len(infile)/3),2].astype(float))
    sigmas.append(infile[int((i)*len(infile)/3):int((i+1)*len(infile)/3),3].astype(float))
    distances.append(infile[int((i)*len(infile)/3):int((i+1)*len(infile)/3),4].astype(float))
fig=plt.figure(figsize=(30,20))
colors=["blue","orange","green","red"]
font = {'family' : 'normal',
        'weight' : 'normal',
        'size'   : 22}
matplotlib.rc('font', **font)
for i in range(3):
    plt.subplot(3,2,2*i+1)
    plt.plot(alpha,energies[i],label=r"$\omega=%.2f$"%omegas[i],color=colors[i])
    mini=np.argmin(energies[i])
    plt.plot(alpha[mini],energies[i][mini],"o",color=colors[i])
    plt.text(alpha[mini],energies[i][mini],"Minimum at %.2f\n with E=%.2f"%(alpha[mini],energies[i][mini]))
    plt.xlabel(r"$\alpha$")
    plt.ylabel(r"$E_0$")
    plt.legend()
    plt.subplot(3,2,2*i+2)
    plt.plot(alpha,sigmas[i],label=r"$\omega=%.2f$"%omegas[i],color=colors[i])
    mini=np.argmin(sigmas[i])
    plt.plot(alpha[mini],sigmas[i][mini],"o",color=colors[i])
    plt.text(alpha[mini],sigmas[i][mini],"Minimum at %.2f\n with $\sigma$=%.2f"%(alpha[mini],sigmas[i][mini]))
    plt.xlabel(r"$\alpha$")
    plt.ylabel(r"$\sigma$")
    plt.legend()
plt.tight_layout()
plt.savefig("../plots/function1_plot.pdf")

plt.show()
