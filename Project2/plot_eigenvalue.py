import numpy as np
import matplotlib.pyplot as plt
omegas=["0.010000","0.500000","1.000000","5.000000"]
omega_val=[float(omega) for omega in omegas]
filenames=["solutions_one_electron_","solutions_two_electrons_"]
lowest_eigenval_two=[] #Two electrons
lowest_eigenval_one=[] #one electron
for filename in filenames:
    for omega in omegas:
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
        if filename=="solutions_one_electron_":
            lowest_eigenval_one.append(np.sort(eigenval)[0])
        else:
            lowest_eigenval_two.append(np.sort(eigenval)[0])
plt.title("Lowest eigenvalue")
plt.xlabel(r"$\rho$")
plt.ylabel(r"$(R(r)*r)^2$")
print(omega_val)
print(lowest_eigenval_one)
print(lowest_eigenval_two)
outfile=open("eigenvalues.txt","w")
outfile.write("%15s %15s %15s\n"%("Eigenvalue_one", "Eigenvalue_two","omega"))
for i in range(len(omega_val)):
    outfile.write("%15.3f %15.3f %15.3f\n"%(lowest_eigenval_one[i],lowest_eigenval_two[i],omega_val[i]))
plt.plot(omega_val,lowest_eigenval_one,"o",label="Lowest eigenvalue as a fucntion of omega for one electron")
plt.plot(omega_val,lowest_eigenval_two,"o",label="Lowest eigenvalue as a fucntion of omega for two electrons")
plt.legend()
plt.savefig("eigenvalues.png")
plt.show()
