import numpy as np
import matplotlib.pyplot as plt
infile=np.loadtxt("results/time_info.csv",delimiter=",",dtype="str")
outfile1=open("results/accuracy_time.csv","w")
gauleg=infile[infile[:,0]=="gauleg"]
gaulag=infile[infile[:,0]=="gaulag"]
monteBad=infile[infile[:,0]=="badMonteCarlo"]
monteGood=infile[infile[:,0]=="goodMonteCarlo"]
monteGood_n=np.unique(monteGood[:,3]).astype(np.int) #N values
monteGood_vectorized_parallel=monteGood[np.logical_and(monteGood[:,1]=="yes",monteGood[:,2]=="4")]
monteGood_vectorized_noparallel=monteGood[np.logical_and(monteGood[:,1]=="yes",monteGood[:,2]=="1")]
monteGood_novectorized_parallel=monteGood[np.logical_and(monteGood[:,1]=="no",monteGood[:,2]=="4")]
monteGood_novectorized_noparallel=monteGood[np.logical_and(monteGood[:,1]=="no",monteGood[:,2]=="1")]
gauleg_novectorized_parallel=gauleg[np.logical_and(gauleg[:,1]=="no",gauleg[:,2]=="4")]
gaulag_novectorized_parallel=gaulag[np.logical_and(gaulag[:,1]=="no",gaulag[:,2]=="4")]
monteBad_novectorized_parallel=monteBad[np.logical_and(monteBad[:,1]=="no",monteBad[:,2]=="4")]
n_gauleg=[10,15,20,25,30,40,45]
n_monte=[1000,10000,100000,1000000,10000000,100000000,1000000000]
time_gauleg=[]
error_gauleg=[]
time_gaulag=[]
error_gaulag=[]
time_monteGood=[]
error_monteGood=[]
time_monteBad=[]
error_monteBad=[]
for n in n_gauleg:
    time_gauleg.append(np.mean(gauleg_novectorized_parallel[gauleg_novectorized_parallel[:,3]==str(n)][:,4].astype(np.float)))
    time_gaulag.append(np.mean(gaulag_novectorized_parallel[gaulag_novectorized_parallel[:,3]==str(n)][:,4].astype(np.float)))
    error_gauleg.append(np.mean(gauleg_novectorized_parallel[gauleg_novectorized_parallel[:,3]==str(n)][:,6].astype(np.float)))
    error_gaulag.append(np.mean(gaulag_novectorized_parallel[gaulag_novectorized_parallel[:,3]==str(n)][:,6].astype(np.float)))
for n in n_monte:
    time_monteGood.append(np.mean(monteGood_novectorized_parallel[monteGood_novectorized_parallel[:,3]==str(n)][:,4].astype(np.float)))
    time_monteBad.append(np.mean(monteBad_novectorized_parallel[monteBad_novectorized_parallel[:,3]==str(n)][:,4].astype(np.float)))
    error_monteGood.append(np.max(monteGood_novectorized_parallel[monteGood_novectorized_parallel[:,3]==str(n)][:,6].astype(np.float)))
    error_monteBad.append(np.max(monteBad_novectorized_parallel[monteBad_novectorized_parallel[:,3]==str(n)][:,6].astype(np.float)))
plt.plot(time_monteGood,error_monteGood,label="MC with importance sampling")
plt.plot(time_monteBad,error_monteBad,label="MC without importance sampling")
plt.plot(time_gaulag,error_gaulag,label="Gauss Laguerre")
print(error_monteBad)
print(error_monteGood)
plt.plot(time_gauleg,error_gauleg,label="Gauss Legendre")
plt.legend()
plt.xlabel("Time [s]")
plt.ylabel("Relative error")
plt.xscale("log")
plt.yscale("log")
plt.savefig("Error_time_use.png")
plt.show()
#time_vectorized_parallel=np.avg()
time_vectorized_parallel=[]
time_vectorized_noparallel=[]
time_novectorized_parallel=[]
time_novectorized_noparallel=[]
monteGood_n_short=np.unique(monteGood_vectorized_noparallel[:,3]).astype(np.int)
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
rel=np.asarray(time_novectorized_noparallel)/np.asarray(time_novectorized_parallel[:6])

plt.plot(monteGood_n[:6],rel)
plt.savefig("Relative_difference_time.png")
plt.show()
