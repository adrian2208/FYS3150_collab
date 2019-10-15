import subprocess, sys
import numpy as np
filenames=["a","b","c","d","a_parallel","b_parallel","c_parallel","d_parallel"];
helpfiles=["lib","functions","tests_main","testings"]
n_gaussisk=[10,15,20]
n_gaussisk_parallel=[10,15,20,25,30,40]
n_monte=[1000,10000,100000]
n_monte_parallel=[1000,10000,100000,100000,1000000,10000000,100000000,1000000000]
for i in range(int(len(filenames)/2)):
    subprocess.call("c++ -c "+filenames[i]+".cpp",shell=True)
for i in range(int(len(filenames)/2),len(filenames)):
    subprocess.call("mpic++ -c "+filenames[i]+".cpp",shell=True)

for filename in helpfiles:
    subprocess.call("c++ -c "+filename+".cpp",shell=True)
helpstring="c++ -o testings.exe "
for file in helpfiles:
    helpstring+=file+".o "
subprocess.call(helpstring,shell=True)
subprocess.call(["./testings.exe"])
outfile=open("results/time_info.csv","w")
outfile.write("type,parrallelisation_flag,Amount_threads,N,time_use,result,relative_error,standard_deviation")
outfile.close()
def compile(fast):
    for i in range(int(len(filenames)/4)):
        subprocess.call("c++ "+fast+" -o "+filenames[i]+".exe "+filenames[i]+".o lib.o functions.o",shell=True)
    for i in range(int(len(filenames)/4),int(len(filenames)/2)):
        subprocess.call("c++ "+fast+" -o "+filenames[i]+".exe "+filenames[i]+".o functions.o",shell=True)
    for i in range(int(len(filenames)/2),int(3*len(filenames)/4)):
        subprocess.call("mpic++ "+fast+" -o "+filenames[i]+".exe "+filenames[i]+".o lib.o functions.o",shell=True)
    for i in range(int(3*len(filenames)/4),len(filenames)):
        subprocess.call("mpic++ "+fast+" -o "+filenames[i]+".exe "+filenames[i]+".o functions.o",shell=True)
compile("-O3")
def runAll(mengde):
    for amount in range(mengde):
        for n in n_gaussisk:
            for i in range(int(len(filenames)/4)):
                subprocess.call("./"+filenames[i]+".exe "+str(n),shell=True)
        for n in n_gaussisk_parallel:
            for i in range(int(len(filenames)/2),int(3*len(filenames)/4)):
                subprocess.call("mpirun -n 4 "+filenames[i]+".exe "+str(n),shell=True)
        for n in n_monte:
            for i in range(int(len(filenames)/4),int(len(filenames)/2)):
                subprocess.call("./"+filenames[i]+".exe "+str(n),shell=True)
        for n in n_monte_parallel:
            for i in range(int(len(filenames)*3/4),len(filenames)):
                subprocess.call("mpirun -n 4 "+filenames[i]+".exe "+str(n),shell=True)
compile("-O3")
runAll(5)
infile=np.loadtxt("results/time_info.csv",delimiter=",",dtype="str")
for i in range(1,len(infile)):
    infile[i][1]='yes'
np.savetxt("results/time_info.csv",infile,delimiter=",",fmt='%s')
compile("")
runAll(5)
