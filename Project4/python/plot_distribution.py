import numpy as np
import matplotlib.pyplot as plt
infile_1=open("../results/individual_energies_1.000000.txt").readlines()
infile_24=open("../results/individual_energies_2.400000.txt").readlines()
energies=[]
temp_1_dist=[]
temp_24_dist=[]
for line in infile_1[2:]:
    energy,dist=line.split(",")
    energies.append(int(energy))
    temp_1_dist.append(int(dist))
for line in infile_24[2:]:
    energy,dist=line.split(",")
    temp_24_dist.append(int(dist))
energies_1=[]
energies=np.asarray(energies)
temp_1_dist=np.asarray(temp_1_dist)/sum(temp_1_dist)
temp_24_dist=np.asarray(temp_24_dist)/sum(temp_24_dist)
"""
for i in range(len(temp_1_dist)):
    while  temp_1_dist[i]>0:
        energies_1.append(energies[i])
        temp_1_dist[i]-=1;
plt.hist(energies_1)
"""
print(temp_1_dist)
plt.bar(energies[:10],temp_1_dist[:10],width=4)
plt.xticks(energies[:10])
plt.title("Distribution at T=1.0")
plt.xlabel("Energy in units of J")
plt.ylabel("Relative occurence")
plt.savefig("../plots/Dist_1.0.pdf")
plt.show()
plt.bar(energies,temp_24_dist,width=4)
#plt.xticks(energies)
plt.title("Distribution at T=2.4")
plt.xlabel("Energy in units of J")
plt.ylabel("Relative occurence")
plt.savefig("../plots/Dist_2.4.pdf")
plt.show()
