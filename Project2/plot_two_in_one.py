import numpy as np
import matplotlib.pyplot as plt
import sys
def density(rho,func):
    return func**2
omegas=["1.000000"]
omega_val=[float(omega) for omega in omegas]
filenames=["solutions_one_electron_","solutions_two_electrons_"]
position_peak_two=[] #Two electrons
position_peak_one=[] #one electron'
for filename in filenames:
    omega=omegas[0]
    infile=open(filename+omega+".txt","r")
    n=int(infile.readline().split()[-1])
    rhomax=float(infile.readline().split()[-1])
    rhomin=float(infile.readline().split()[-1])
    omega=""
    try:
        omega=float(infile.readline().split()[-1])
    except:
        print("fail")
        pass
    eigenval=[float(i) for i in infile.readline().split()];
    values=infile.readlines();
    for i in range(len(values)):
        values[i]=list(map(float, values[i].split()))
    eigenvec=np.array(values)
    print(eigenvec)
    i=np.argsort(eigenval)
    eigenvec=eigenvec[:,i]
    eigenvec=np.transpose(eigenvec)
    rho=np.linspace(rhomin,rhomax,n)
    label="Ground state with one electron"
    if filename=="solutions_two_electrons_":
        label="Ground state with two electrons"
    densityp=density(rho,eigenvec[0])
    densityp=densityp/sum(densityp)
    plt.plot(rho,densityp,label=label)
plt.title("Comparision of wave functions for omega=%s"%omegas[0])
plt.xlabel(r"$\rho$")
plt.ylabel(r"$(u(\rho))^2$")
plt.legend()
plt.savefig("ground_state_comparision_omega=%s.png"%omegas[0])
plt.show()
