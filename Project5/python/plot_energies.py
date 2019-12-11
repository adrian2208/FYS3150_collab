import numpy as np
import matplotlib.pyplot as plt
import matplotlib
import sys
infile=np.loadtxt("../results/function1.csv",delimiter=",",skiprows=1)
omegas=[0.01,0.5,1.0,5.0]
alpha=infile[:int(len(infile)/len(omegas)),1]
energies=[]
sigmas=[]
distances=[]
for i in range(len(omegas)):
    energies.append(infile[int((i)*len(infile)/len(omegas)):int((i+1)*len(infile)/len(omegas)),2].astype(float))
    sigmas.append(infile[int((i)*len(infile)/len(omegas)):int((i+1)*len(infile)/len(omegas)),3].astype(float))
    distances.append(infile[int((i)*len(infile)/len(omegas)):int((i+1)*len(infile)/len(omegas)),4].astype(float))
fig=plt.figure(figsize=(20,20))
colors=["blue","orange","green","red"]
font = {'family' : 'normal',
        'weight' : 'normal',
        'size'   : 25}
matplotlib.rc('font', **font)
for i in range(len(omegas)):
    plt.subplot(len(omegas),2,2*i+1)
    plt.plot(alpha,energies[i],label=r"$\omega=%.2f$"%omegas[i],color=colors[i])
    mini=np.argmin(energies[i])
    plt.plot(alpha[mini],energies[i][mini],"o",color=colors[i],markersize=15)
    #plt.text(alpha[mini],energies[i][mini],"Minimum at %.3f\n with E=%.3f"%(alpha[mini],energies[i][mini]))
    print("Omega: %.4f Alpha: %.4f Energy: %.4f Sigma: %.4f Distance: %.4f"%(omegas[i],alpha[mini],energies[i][mini],sigmas[i][mini],distances[i][mini]))
    plt.xlabel(r"$\alpha$")
    plt.ylabel(r"$E_0$")
    plt.legend()
    plt.subplot(len(omegas),2,2*i+2)
    plt.plot(alpha,sigmas[i],label=r"$\omega=%.2f$"%omegas[i],color=colors[i])
    mini=np.argmin(sigmas[i])
    plt.plot(alpha[mini],sigmas[i][mini],"o",color=colors[i],markersize=15)
    #plt.text(alpha[mini],sigmas[i][mini],"Minimum at %.3f\n with $\sigma$=%.3f"%(alpha[mini],sigmas[i][mini]))
    plt.xlabel(r"$\alpha$")
    plt.ylabel(r"$\sigma$")
    plt.legend()
plt.tight_layout()
plt.savefig("../plots/function1_plot.pdf")

plt.show()
