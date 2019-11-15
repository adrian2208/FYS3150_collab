import numpy as np
import matplotlib.pyplot as plt
import matplotlib
import sys
infilename="results_calculations_general.csv"

if len(sys.argv) >=2:
    infilename=sys.argv[1]

##There's more infilenames
data=np.loadtxt("../results/%s"%infilename,skiprows=1,dtype=float,delimiter=",")

temperatures=data[:,0][0:int(len(data[:,0])/4)]
lens=[40,60,80,100]
inv_lens=np.asarray(lens,dtype=float)**(-1)
crit_temp_cv=[]
crit_temp_chi=[]
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

colors=["blue","orange","green","red"]
font = {'family' : 'normal',
        'weight' : 'normal',
        'size'   : 22}

matplotlib.rc('font', **font)
plt.figure(figsize=(20,20))

plt.subplot(221)
for i in range(len(cv)):
    plt.plot(temperatures,cv[i],color=colors[i],label="L=%s"%str(lens[i]))
    maxindex=np.argmax(cv[i])
    crit_temp_cv.append(temperatures[maxindex])
    plt.plot(temperatures[maxindex],cv[i][maxindex], "o",color=colors[i])
plt.xlabel("Temperature [kT/J]")
plt.ylabel(r"Heat capacity $C_V$ per particle [$k_B$]")
plt.legend()
plt.subplot(222)

for i in range(len(xi)):
    plt.plot(temperatures,xi[i],color=colors[i],label="L=%s"%str(lens[i]))
    maxindex=np.argmax(xi[i])
    crit_temp_chi.append(temperatures[maxindex])
    plt.plot(temperatures[maxindex],xi[i][maxindex], "o",color=colors[i])
plt.xlabel("Temperature [kT/J]")
plt.ylabel(r"Susceptibility $\chi$ per particle [$\frac{1}{J}$]")
plt.legend()
plt.subplot(223)

for i in range(len(energy)):
    plt.plot(temperatures,energy[i],color=colors[i],label="L=%s"%str(lens[i]))
plt.xlabel("Temperature [kT/J]")
plt.ylabel("Avg. energy per particle [J]")
plt.legend()
plt.subplot(224)

for i in range(len(magnetisation)):
    plt.plot(temperatures,magnetisation[i],color=colors[i],label="L=%s"%str(lens[i]))
plt.xlabel("Temperature [kT/J]")
plt.ylabel("Avg. abs. Magnetisation per particle [unitless]")
plt.legend()
plt.tight_layout()

plt.savefig("../plots/%s.pdf"%infilename[:-4])
plt.show()

regression_cv=np.poly1d(np.polyfit(inv_lens,crit_temp_chi,deg=1))
print(inv_lens,crit_temp_chi)
print(regression_cv)
regression_chi=np.poly1d(np.polyfit(inv_lens,crit_temp_cv,deg=1))
x=np.linspace(0,1,1000)

plt.plot(x,regression_cv(x))
plt.plot(x,regression_chi(x))
plt.ylim(2.2,2.3)
plt.show()
