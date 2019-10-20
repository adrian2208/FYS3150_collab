import numpy as np
import matplotlib.pyplot as plt
infile=np.loadtxt("../results/accuracy_time.csv",delimiter=",",dtype="str")
infile2=np.loadtxt("../results/time_info.csv",delimiter=",",dtype="str")
gauleg=infile[infile[:,0]=="gauleg"].transpose()
gaulag=infile[infile[:,0]=="gaulag"].transpose()
monteBad=infile[infile[:,0]=="no_importance_sampling"].transpose()
monteGood=infile[infile[:,0]=="importance_sampling"].transpose()

plt.plot(monteBad[4,:].astype(np.float),monteBad[3,:].astype(np.float),label="MC without importance sampling")
plt.plot(monteGood[4,:].astype(np.float),monteGood[3,:].astype(np.float),label="MC with importance sampling")
plt.plot(gaulag[4,:].astype(np.float),gaulag[3,:].astype(np.float),label="Gauss Laguerre")
plt.plot(gauleg[4,:].astype(np.float),gauleg[3,:].astype(np.float),label="Gauss Legendre")
plt.legend()
plt.xlabel("Time [s]")
plt.ylabel("Relative error")
plt.xscale("log")
plt.yscale("log")
plt.savefig("../plots/Error_time_use.pdf")
plt.show()
plt.plot(monteGood[1,:].astype(np.float),monteGood[5,:].astype(np.float)/monteGood[4,:].astype(np.float))
plt.title("Relative difference in time use as function of n")
plt.xlabel("N")
plt.ylabel(r"Relative difference $\frac{nonparallelized}{parallelized}$")
plt.xscale("log")
plt.savefig("../plots/Relative_difference_time.pdf")
plt.show()
infile2_monte=infile2[infile2[:,0]=="goodMonteCarlo"]
monteGood_vectorized_parallel=infile2_monte[np.logical_and(infile2_monte[:,1]=="yes",infile2_monte[:,2]=="4")]
monteGood_vectorized_noparallel=infile2_monte[np.logical_and(infile2_monte[:,1]=="yes",infile2_monte[:,2]=="1")]
time_vectorized_parallel=[]
time_vectorized_noparallel=[]
for n in monteGood[1,:].astype(np.int):
    time_vectorized_noparallel.append(np.mean(monteGood_vectorized_noparallel[monteGood_vectorized_noparallel[:,3]==str(n)][:,4].astype(np.float)))
    time_vectorized_parallel.append(np.mean(monteGood_vectorized_parallel[monteGood_vectorized_parallel[:,3]==str(n)][:,4].astype(np.float)))

#plt.plot(monteGood[1,:].astype(np.float),monteGood[5,:].astype(np.float).astype(np.float),label="non-parallel")
#plt.plot(monteGood[1,:].astype(np.float),monteGood[4,:].astype(np.float).astype(np.float),label="non-vectorized, non-parallel")
plt.plot(monteGood[1,:].astype(np.float),monteGood[5,:].astype(np.float).astype(np.float)/np.asarray(time_vectorized_noparallel),label="vectorized, non-parallel")
plt.plot(monteGood[1,:].astype(np.float),monteGood[4,:].astype(np.float).astype(np.float)/np.asarray(time_vectorized_parallel),label="vectorized, parallel")
plt.xscale("log")
plt.xlabel("N")
plt.ylabel(r"Relative difference $\frac{nonvectorized}{vectorized}$")

#plt.yscale("log")
plt.legend()
plt.savefig("../plots/Parrallelisation_flag.pdf")
plt.show()

"""
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
gauleg_novectorized_noparallel=gauleg[np.logical_and(gauleg[:,1]=="no",gauleg[:,2]=="1")]
gaulag_novectorized_noparallel=gaulag[np.logical_and(gaulag[:,1]=="no",gaulag[:,2]=="1")]
monteBad_novectorized_parallel=monteBad[np.logical_and(monteBad[:,1]=="no",monteBad[:,2]=="4")]
monteBad_novectorized_noparallel=monteBad[np.logical_and(monteBad[:,1]=="no",monteBad[:,2]=="1")]

n_gauleg=[10,15,20,25,30,40,45,55,65]
n_monte=[1000,10000,100000,1000000,10000000,100000000,1000000000]







time_gauleg=[]
error_gauleg=[]
time_gaulag=[]
error_gaulag=[]
time_monteGood=[]
error_monteGood=[]
time_monteBad=[]
error_monteBad=[]
result_monteGood=[]
result_monteBad=[]
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
plt.savefig("../plots/Error_time_use.pdf")
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
plt.savefig("../plots/Relative_difference_time.pdf")
plt.show()
"""
