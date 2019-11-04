#This...somewhat...broke down
import numpy as np
import matplotlib.pyplot as plt
infile=np.loadtxt("../results/all_results.csv",skiprows=1,dtype=float,delimiter=",")
temp=infile[0][0]
infile=infile.transpose()
matrix_size=infile[1][0]
indeces=(infile[2][:])[1:]
print(indeces)
energies=(infile[:][3])[1:]/(matrix_size**2)
accepted_configurations=(infile[:][5])[1:]
print(energies)
magnetics=-(infile[:][4])[1:]/(matrix_size**2)
plt.plot(indeces,energies,label="Energy")
plt.plot(indeces,magnetics,label="Magnetisation")
plt.title("Random,Temp: %5.2f N: %2.0f"%(temp,matrix_size**2))
plt.xlabel("Simulation step")
plt.ylabel("Energy [J]/Magnetisation p.p.")
plt.legend()
plt.savefig("../plots/TIMEDEV_%s_size_%s_temp_rand.pdf"%(matrix_size,temp))
plt.show()
plt.plot(indeces,100*accepted_configurations/indeces)
plt.title("Random,Temp: %5.2f N: %2.0f"%(temp,matrix_size**2))
plt.xlabel("Simulation step")
plt.xlim(indeces[0],indeces[-1])
plt.ylabel("Percentage of accepted configurations")
plt.savefig("../plots/configurations_%s_size_%s_temp_rand.pdf"%(matrix_size,temp))
plt.show()
