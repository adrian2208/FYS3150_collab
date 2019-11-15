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
temp_1_dist=np.asarray(temp_1_dist)
temp_24_dist=np.asarray(temp_24_dist)
average24=average10=var24=var10=0;
for i in range(len(energies)):
    average24+=energies[i]*temp_24_dist[i]
    average10+=energies[i]*temp_1_dist[i]
average24=average24/sum(temp_24_dist)
average10=average10/sum(temp_1_dist)
for i in range(len(energies)):
    var24+=temp_24_dist[i]*(energies[i]-average24)**2
    var10+=temp_1_dist[i]*(energies[i]-average10)**2
var24=var24/sum(temp_24_dist)
var10=var10/sum(temp_1_dist)
print(sum(temp_24_dist))
print(sum(temp_1_dist))
print("Average: %f , Standard deviation for T=2.4: %f"%(average24,np.sqrt(var24)))
print("Average: %f , Standard deviation for T=1.0: %f"%(average10,np.sqrt(var10)))

font = {'family' : 'normal',
        'weight' : 'normal',
        'size'   : 22}
temp_1_dist=np.asarray(temp_1_dist)/sum(temp_1_dist)
temp_24_dist=np.asarray(temp_24_dist)/sum(temp_24_dist)
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
