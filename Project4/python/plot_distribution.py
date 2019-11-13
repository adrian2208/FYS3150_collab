import numpy as np
import matplotlib.pyplot as plt
import matplotlib
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
font = {'family' : 'normal',
        'weight' : 'normal',
        'size'   : 22}

matplotlib.rc('font', **font)

plt.figure(figsize=(20,10))

plt.subplot(121)
plt.bar(energies[:10],temp_1_dist[:10],width=4)
plt.xticks(energies[:10:2])
plt.title("T=1.0, spins 400, runs:1e6")
plt.xlabel("Energy in units of J")
plt.ylabel("Relative occurence")
plt.subplot(122)
plt.bar(energies[:int((len(energies)/2))],temp_24_dist[:int((len(energies)/2))],width=4)
plt.title("T=2.4, spins 400, runs:1e6")
plt.xlabel("Energy in units of J")
plt.ylabel("Relative occurence")
plt.savefig("../plots/Dist.pdf",bbox_inches='tight')
plt.show()
