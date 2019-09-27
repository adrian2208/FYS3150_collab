"""For an existing file, finds the maximum value of a given wave function and writes it to file. Here, many files.
"""

import numpy as np
import matplotlib.pyplot as plt
import sys
def density(rho,func): #The function squared. Rho for rhoference
    return func**2
omegas=["0.010000","0.500000","1.000000","5.000000"] #for file name purposes
omega_val=[float(omega) for omega in omegas]
filenames=["solutions_one_electron_","solutions_two_electrons_"]
position_peak_two=[] #Two electrons
position_peak_one=[] #one electron
outfile=open("maximum_peaks_wave_functions.txt","w") #written to this file
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
            pass
        eigenval=[float(i) for i in infile.readline().split()];
        values=infile.readlines();
        for i in range(len(values)):
            values[i]=list(map(float, values[i].split()))
        eigenvec=np.array(values)
        i=np.argsort(eigenval)
        eigenvec=eigenvec[:,i]
        eigenvec=np.transpose(eigenvec)
        rho=np.linspace(rhomin,rhomax,n)
        groundstate=eigenvec[0]
        maxpeak=np.argmax(density(rho,groundstate))
        if filename=="solutions_one_electron_":
            position_peak_one.append(rho[maxpeak])
        else:
            position_peak_two.append(rho[maxpeak])
outfile.write("%15s %15s %15s\n"%("Maxpeak_one", "Maxpeak_two","omega"))
for i in range(len(omega_val)):
    outfile.write("%15.3f %15.3f %15.3f\n"%(position_peak_one[i],position_peak_two[i],omega_val[i]))
outfile.close()
plt.show()
