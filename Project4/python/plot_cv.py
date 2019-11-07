import numpy as np
import matplotlib.pyplot as plt
data=np.loadtxt("../results/results_calculations.csv",skiprows=1,dtype=float,delimiter=",")

temperatures=data[:,0][0:int(len(data[:,0])/4)]
lens=[40,60,80,100]
print(temperatures)
cv=[[] for lenny in lens]
print(data[:,0])
print(data[0,:])
for i in range(len(lens)):
    for temperature in temperatures:
        relevant=data[np.logical_and(abs(data[:,0]-temperature)<1e-10,data[:,1]==lens[i])][0][-2]
        cv[i].append(relevant)
        print(temperature,lens[i],relevant)
for i in range(len(cv)):
    plt.plot(temperatures,cv[i],label="Size=%s"%str(lens[i]))
plt.xlabel("Temperature")
plt.ylabel("Cv")
plt.legend()
plt.show()
