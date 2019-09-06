"""
Plots the log(time use) that is given in time_info.txt for the general algorithm, the specialized algorithm
and the LU decomposition solution against log(n)
"""
import matplotlib.pyplot as plt
import numpy as np
def makedict(listy):
    dicty={}
    for element in listy:
        nval=int(element[0])
        if nval in dicty.keys():
            dicty[nval].append(float(element[1]))
        else:
            dicty[nval]=[]
            dicty[nval].append(float(element[1]))
    return dicty;
def makeplots(dicty):
    n=[]
    val=[]
    for key in dicty:
        n.append(key)
        val.append(sum(dicty[key])/len(dicty[key]))
    return np.array(n), np.array(val), len(dicty[key])
infile=open("time_info.txt","r")
lines=infile.readlines();
infile.close();
b=[]
c=[]
e=[]
for line in lines:
    line=line.split();
    if line[0]=="b":
        b.append([line[2],line[-1]]);

    elif line[0]=="c":
        c.append([line[2],line[-1]])
    elif line[0]=="e":
        e.append([line[2],line[-1]])
b_dict=makedict(b)
c_dict=makedict(c)
e_dict=makedict(e)
n,eval,amount=makeplots(e_dict)
plt.plot(np.log10(n),np.log10(eval),label="LU decomposition of general matrix (%d runs)"%amount)
n,bval,amount=makeplots(b_dict)
plt.plot(np.log10(n),np.log10(bval),label="general triangular algorithm (%d runs)"%amount)
n,cval,amount=makeplots(c_dict)
plt.plot(np.log10(n),np.log10(cval),label="specialized algorithm (%d runs)"%amount)
plt.xlabel("$log_{10}(n)$")
plt.ylabel("$log_{10}$(ellapsed time [s])")
plt.legend()
plt.savefig("time_use.png")
plt.show()
plt.plot(np.log10(n),cval/bval)
plt.xlabel("$log_{10}(n)$")
plt.ylabel("Percentual time use of specialized algorithm comp. to general")
plt.savefig("time_use2.png")
plt.show()
