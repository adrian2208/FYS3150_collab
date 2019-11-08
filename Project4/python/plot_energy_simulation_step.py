#This...somewhat...broke down
import numpy as np
import matplotlib.pyplot as plt
infile_random=np.loadtxt("../results/all_results_random.csv",skiprows=1,dtype=float,delimiter=",")
infile_up=np.loadtxt("../results/all_results_up.csv",skiprows=1,dtype=float,delimiter=",")
temp=infile_random[0][0]
infile_random=infile_random.transpose()
infile_up=infile_up.transpose()

matrix_size=infile_random[1][0]
indeces=(infile_random[2][:])[1:]
energies_random=(infile_random[:][3])[1:]/(matrix_size**2)
energies_up=(infile_up[:][3])[1:]/(matrix_size**2)

accepted_configurations_random=(infile_random[:][5])[1:]/(matrix_size**2)
magnetics_random=(infile_random[:][4])[1:]/(matrix_size**2)
accepted_configurations_up=(infile_up[:][5])[1:]/(matrix_size**2)
magnetics_up=(infile_up[:][4])[1:]/(matrix_size**2)
plt.figure(figsize=(20,10))
plt.subplot(131)
plt.title("Avg. Energy, %d spins, T=%.1f"%(matrix_size**2,temp))
plt.xlabel("Simulation step")
plt.ylabel("Energy [J] per part.")
plt.plot(indeces,energies_random,label="Random")
plt.plot(indeces,energies_up,label="Up")
plt.grid(True)

plt.legend()

plt.subplot(132)
plt.title("Avg. Magnetisation, %d spins, T=%.1f"%(matrix_size**2,temp))
plt.xlabel("Simulation step")
plt.ylabel("Magnetisation per part.")
plt.plot(indeces,magnetics_random,label="Random")
plt.plot(indeces,magnetics_up,label="Up")
plt.grid(True)

plt.legend()

plt.subplot(133)
plt.title("Acc. configurations, %d spins, T=%.1f"%(matrix_size**2,temp))
plt.xlabel("Simulation step")
plt.ylabel("Percentage of accepted configurations")
plt.plot(indeces,100*accepted_configurations_random/indeces,label="Random")
plt.plot(indeces,100*accepted_configurations_up/indeces,label="Up")
plt.legend()
plt.grid(True)
#plt.tight_layout()

plt.savefig("../plots/configurations_%s_size_%.1f_temp.pdf"%(matrix_size,temp))
plt.show()
