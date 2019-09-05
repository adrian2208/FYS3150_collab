"""
Reads from each file that has the name "oppc_NR.txt" (where NR is 10**n), and plots the log of the maximum error against the log of the steplength
"""
import matplotlib.pyplot as plt
import numpy as np
i=10
iss=[]
maxErr=[]
while(i<=1e7):
    iss.append(np.log10(1/(i+1)))
    infile=open("oppc_"+str(i)+".txt") #given that these files already exist
    infile.readline();
    infile.readline();
    infile.readline();
    lines=infile.readlines();
    error=np.zeros(len(lines)-1)
    for rk in range(len(lines)-1):
        error[rk]=(float(lines[rk].split()[-1]))
    i*=10
    maxErr.append(np.log10(np.amax(error)))
    infile.close()
outfile=open("log_err_c.txt","w")
outfile.write("log10(h) log10(relative_error)\n")
for i in range(len(maxErr)):
    outfile.write("%7.3f %7.3f\n"%(iss[i],maxErr[i]))
outfile.close()
