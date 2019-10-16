import numpy as np
import matplotlib.pyplot as plt
infile=np.loadtxt("results/time_info.csv",delimiter=",",dtype="str")
outfile1=open("results/accuracy_time.csv","w")
gauleg=infile[infile[:,0]=="gauleg"]
gaulag=infile[infile[:,0]=="gauleg"]
monteBad=infile[infile[:,0]=="badMonteCarlo"]
monteGood=infile[infile[:,0]=="goodMonteCarlo"]
monteGood_n=np.unique(monteGood[:,3]).astype(np.int) #N values
monteGood_vectorized_parallel=monteGood[np.logical_and(monteGood[:,1]=="yes",monteGood[:,2]=="4")]
monteGood_vectorized_noparallel=monteGood[np.logical_and(monteGood[:,1]=="yes",monteGood[:,2]=="1")]
monteGood_novectorized_parallel=monteGood[np.logical_and(monteGood[:,1]=="no",monteGood[:,2]=="4")]
monteGood_novectorized_noparallel=monteGood[np.logical_and(monteGood[:,1]=="no",monteGood[:,2]=="1")]
#time_vectorized_parallel=np.avg()
time_vectorized_parallel=[]
time_vectorized_noparallel=[]
time_novectorized_parallel=[]
time_novectorized_noparallel=[]
monteGood_n_short=np.unique(monteGood_vectorized_noparallel[:,3]).astype(np.int)
print(monteGood_n_short)
for n in monteGood_n_short:
    time_vectorized_noparallel.append(np.mean(monteGood_vectorized_noparallel[monteGood_vectorized_noparallel[:,3]==str(n)][:,4].astype(np.float)))
    time_novectorized_noparallel.append(np.mean(monteGood_novectorized_noparallel[monteGood_novectorized_noparallel[:,3]==str(n)][:,4].astype(np.float)))
for n in monteGood_n:
    time_vectorized_parallel.append(np.mean(monteGood_vectorized_parallel[monteGood_vectorized_parallel[:,3]==str(n)][:,4].astype(np.float)))
    time_novectorized_parallel.append(np.mean(monteGood_novectorized_parallel[monteGood_novectorized_parallel[:,3]==str(n)][:,4].astype(np.float)))
plt.plot(monteGood_n_short,time_vectorized_noparallel,label="time_vectorized_noparallel")
plt.plot(monteGood_n_short,time_novectorized_noparallel,label="time_novectorized_noparallel")
plt.plot(monteGood_n[:6],time_vectorized_parallel[:6],label="time_vectorized_parallel")
plt.plot(monteGood_n[:6],time_novectorized_parallel[:6],label="time_novectorized_parallel")
plt.xscale("log")
#plt.yscale("log")
plt.legend()
plt.show()
plt.title("Relative difference in time use as function of n")
plt.xlabel("N")
plt.ylabel(r"Relative difference $\frac{nonparallelized}{parallelized}$")
plt.xscale("log")
print(time_vectorized_noparallel)
rel=np.asarray(time_novectorized_noparallel)/np.asarray(time_novectorized_parallel[:5])

plt.plot(monteGood_n[:5],rel)
plt.savefig("Relative_difference_time.png")
plt.show()
