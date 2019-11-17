import numpy as np
import matplotlib.pyplot as plt
import matplotlib
import sys
def readfromtwo():
    infilename="results_calculations_40.csv"
    data=np.loadtxt("../results/%s"%infilename,skiprows=1,dtype=float,delimiter=",")
    temperatures2=data[:,0]
    print(temperatures)
    xi=data[:,-1]
    return temperatures2,xi
def calculate_regression(x_val,y_val):
    regression=np.poly1d(np.polyfit(x_val,y_val,deg=1))
    mse=0;
    for i in range(len(x_val)):
        mse+=(y_val[i]-regression(x_val[i]))**2
    return regression(0),np.sqrt(mse/(len(x_val)-2))
infilename="results_calculations_superaccurate.csv"
if len(sys.argv) >=2:
    infilename=sys.argv[1]

##There's more infilenames
data=np.loadtxt("../results/%s"%infilename,skiprows=1,dtype=float,delimiter=",")

temperatures=data[:,0][0:int(len(data[:,0])/4)]
lens=[40,60,80,100]
inv_lens=np.asarray(lens,dtype=float)**(-1)
crit_temp_cv=[]
crit_temp_chi=[]
cv=[[] for lenny in lens]
energy=[[] for lenny in lens]
xi=[[] for lenny in lens]
magnetisation=[[] for lenny in lens]

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

temperatures2=temperatures
if infilename=="results_calculations_superaccurate.csv":
    temperatures2,xi[0] = readfromtwo()
    print(temperatures2)
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
    if i==0:
        plt.plot(temperatures2,xi[i],color=colors[i],label="L=%s"%str(lens[i]))
        maxindex=np.argmax(xi[i])
        crit_temp_chi.append(temperatures2[maxindex])
        plt.plot(temperatures2[maxindex],xi[i][maxindex], "o",color=colors[i])
    else:
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

krit_cv,rmse_cv=calculate_regression(inv_lens,crit_temp_cv)
krit_chi,rmse_chi=calculate_regression(inv_lens,crit_temp_chi)
print("The critical temperature using Cv is found to be %f, with an RMSD value of %f "%(krit_cv,rmse_cv))
print("The critical temperature using chi is found to be %f, with an RMSD value of %f "%(krit_chi,rmse_chi))
print("\n\n%8s%8s%8s"%("L", "Tc(Cv)", "Tc(chi)"))

for i in range(len(lens)):
    print("%8d%8.4f%8.4f"%(lens[i],crit_temp_cv[i],crit_temp_chi[i]))
