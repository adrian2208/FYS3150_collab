import numpy as np
import matplotlib.pyplot as plt
A=np.zeros((2,2),dtype=int)
E_avg=0
E_squared=0
M_avg=0
M_abs_avg=0
M_squared=0

T=10
beta=1/T
energies=[]
magnetics=[]
def calculate_energy(A):
	return -(2*A[0,0]*A[0,1]+2*A[0,0]*A[1,0]+2*A[1,0]*A[1,1]+2*A[0,1]*A[1,1])
def calculate_magnetic(A):
	return A[0][0]+A[0][1]+A[1][0]+A[1][1]

for i in range(-1,2,2):
	for j in range(-1,2,2):
		for k in range(-1,2,2):
			for l in range(-1,2,2):
				A[0,0]=i
				A[0,1]=j
				A[1,0]=k
				A[1,1]=l
				energies.append(calculate_energy(A))
				magnetics.append(calculate_magnetic(A))
Z=0
for i in range(len(energies)):
	Z=Z+np.exp(-beta*energies[i])
for i in range(len(energies)):
	E_avg=E_avg+np.exp(-beta*energies[i])*energies[i]
	E_squared=E_squared+np.exp(-beta*energies[i])*energies[i]**2
	M_avg=M_avg+np.exp(-beta*energies[i])*magnetics[i]
	M_abs_avg=M_abs_avg+np.exp(-beta*energies[i])*np.abs(magnetics[i])
	M_squared=M_squared+np.exp(-beta*energies[i])*magnetics[i]**2
E_avg/=Z;E_squared/=Z;M_avg/=Z;M_abs_avg/=Z;M_squared/=Z;
print("Average Energy: %.5f\nAverage Magnetisation: %.5f\nAverage absolute magnetisation: %.5f\nE**2: %.5f\nM**2: %.5f"%(E_avg,M_avg,M_abs_avg,E_squared,M_squared))
