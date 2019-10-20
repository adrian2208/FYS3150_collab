import numpy as np
import matplotlib.pyplot as plt
infile=np.loadtxt("../results/time_info.csv",delimiter=",",dtype="str")
gauleg=infile[infile[:,0]=="gauleg"]
gaulag=infile[infile[:,0]=="gaulag"]
monteBad=infile[infile[:,0]=="badMonteCarlo"]
monteGood=infile[infile[:,0]=="goodMonteCarlo"]
outfile=open("../results/accuracy_time.csv","w")
outfile.write("Type,n,Result,Relative_error,time_use_parallel,time_use_nonparallel,StandardError\n")
n_gauleg=[10,15,20,25,30,40,45,55,65]
n_monte=[1000,10000,100000,1000000,10000000,100000000,1000000000]
monteGood_novectorized_parallel=monteGood[np.logical_and(monteGood[:,1]=="no",monteGood[:,2]=="4")]
monteGood_novectorized_noparallel=monteGood[np.logical_and(monteGood[:,1]=="no",monteGood[:,2]=="1")]
gauleg_novectorized_parallel=gauleg[np.logical_and(gauleg[:,1]=="no",gauleg[:,2]=="4")]
gaulag_novectorized_parallel=gaulag[np.logical_and(gaulag[:,1]=="no",gaulag[:,2]=="4")]
gauleg_novectorized_noparallel=gauleg[np.logical_and(gauleg[:,1]=="no",gauleg[:,2]=="1")]
gaulag_novectorized_noparallel=gaulag[np.logical_and(gaulag[:,1]=="no",gaulag[:,2]=="1")]
monteBad_novectorized_parallel=monteBad[np.logical_and(monteBad[:,1]=="no",monteBad[:,2]=="4")]
monteBad_novectorized_noparallel=monteBad[np.logical_and(monteBad[:,1]=="no",monteBad[:,2]=="1")]
for n in n_gauleg:
    outfile.write("gauleg,")
    outfile.write(str(n)+",")
    outfile.write(str(np.mean(gauleg_novectorized_parallel[gauleg_novectorized_parallel[:,3]==str(n)][:,5].astype(np.float)))) #Result
    outfile.write(",")
    outfile.write(str(np.mean(gauleg_novectorized_parallel[gauleg_novectorized_parallel[:,3]==str(n)][:,6].astype(np.float)))) #Relative Error
    outfile.write(",")
    outfile.write(str(np.mean(gauleg_novectorized_parallel[gauleg_novectorized_parallel[:,3]==str(n)][:,4].astype(np.float)))) #time_parallel
    outfile.write(",")
    try:
        outfile.write(str(np.mean(gauleg_novectorized_noparallel[gauleg_novectorized_noparallel[:,3]==str(n)][:,4].astype(np.float)))) #time_noparallel
    except:
        outfile.write("")
    outfile.write(",0\n")
for n in n_gauleg:
    outfile.write("gaulag,")
    outfile.write(str(n)+",")
    outfile.write(str(np.mean(gaulag_novectorized_parallel[gaulag_novectorized_parallel[:,3]==str(n)][:,5].astype(np.float)))) #Result
    outfile.write(",")
    outfile.write(str(np.mean(gaulag_novectorized_parallel[gaulag_novectorized_parallel[:,3]==str(n)][:,6].astype(np.float)))) #Relative Error
    outfile.write(",")
    outfile.write(str(np.mean(gaulag_novectorized_parallel[gaulag_novectorized_parallel[:,3]==str(n)][:,4].astype(np.float)))) #time_parallel
    outfile.write(",")
    try:
        outfile.write(str(np.mean(gaulag_novectorized_noparallel[gaulag_novectorized_noparallel[:,3]==str(n)][:,4].astype(np.float)))) #time_noparallel
    except:
        outfile.write("")
    outfile.write(",0\n")
for n in n_monte:
    outfile.write("no_importance_sampling,")
    outfile.write(str(n)+",")
    errors=(monteBad_novectorized_parallel[monteBad_novectorized_parallel[:,3]==str(n)][:,6].astype(np.float))
    median_error= np.where(errors==np.percentile(errors,50,interpolation='nearest'))
    biggest_error=np.argmax(monteBad_novectorized_parallel[monteBad_novectorized_parallel[:,3]==str(n)][:,6].astype(np.float))
    outfile.write(str((monteBad_novectorized_parallel[monteBad_novectorized_parallel[:,3]==str(n)][:,5].astype(np.float)[median_error])[0])) #Result
    outfile.write(",")
    outfile.write(str((monteBad_novectorized_parallel[monteBad_novectorized_parallel[:,3]==str(n)][:,6].astype(np.float)[median_error])[0])) #Error
    outfile.write(",")
    outfile.write(str(np.mean(monteBad_novectorized_parallel[monteBad_novectorized_parallel[:,3]==str(n)][:,4].astype(np.float)))) #time_parallel
    outfile.write(",")
    std_arr=monteBad_novectorized_parallel[monteBad_novectorized_parallel[:,3]==str(n)][:,7].astype(np.float)
    variance_arr=(monteBad_novectorized_parallel[monteBad_novectorized_parallel[:,3]==str(n)][:,7].astype(np.float)*np.sqrt(n))**2
    std=np.sqrt(np.mean(variance_arr))/np.sqrt(n)
    outfile.write(str(np.mean(monteBad_novectorized_noparallel[monteBad_novectorized_noparallel[:,3]==str(n)][:,4].astype(np.float)))) #time_noparallel
    outfile.write(",")
    outfile.write(str(std))
    outfile.write("\n")
for n in n_monte:
    outfile.write("importance_sampling,")
    outfile.write(str(n)+",")
    errors=(monteGood_novectorized_parallel[monteGood_novectorized_parallel[:,3]==str(n)][:,6].astype(np.float))
    median_error= np.where(errors==np.percentile(errors,50,interpolation='nearest'))
    biggest_error=np.argmax(monteGood_novectorized_parallel[monteGood_novectorized_parallel[:,3]==str(n)][:,6].astype(np.float))
    outfile.write(str((monteGood_novectorized_parallel[monteGood_novectorized_parallel[:,3]==str(n)][:,5].astype(np.float)[median_error])[0])) #Result
    outfile.write(",")
    outfile.write(str((monteGood_novectorized_parallel[monteGood_novectorized_parallel[:,3]==str(n)][:,6].astype(np.float)[median_error])[0])) #Error
    outfile.write(",")
    outfile.write(str(np.mean(monteGood_novectorized_parallel[monteGood_novectorized_parallel[:,3]==str(n)][:,4].astype(np.float)))) #time_parallel
    outfile.write(",")
    variance_arr=(monteGood_novectorized_parallel[monteGood_novectorized_parallel[:,3]==str(n)][:,7].astype(np.float)*np.sqrt(n))**2
    std=np.sqrt(np.mean(variance_arr))/np.sqrt(n)
    outfile.write(str(np.mean(monteGood_novectorized_noparallel[monteGood_novectorized_noparallel[:,3]==str(n)][:,4].astype(np.float)))) #time_noparallel
    outfile.write(",")
    outfile.write(str(std))
    outfile.write("\n")
outfile.close()
