import matplotlib.pyplot as plt
import numpy as np
from numpy import exp
n=[10,100,1000]
data=[]
def f(x):
    return 1-(1-exp(-10))*x-exp(-10*x)
for i in n:
    infile=open("oppb_%d.txt"%i);
    infile.readline();
    allines=infile.readlines();
    x=[]
    sol=[]
    for line in allines:
        nv,xv, accv,appv=line.split();
        x.append(float(xv));
        sol.append(float(appv));
    data.append([x,sol,i])
    infile.close();
for set in data:
    plt.plot(set[0],set[1],"--",label="Approximation with n=%d"%set[2])
x=np.linspace(0,1,1000)
plt.plot(x,f(x),label="Analytical solution")
plt.legend()
plt.xlabel("x-value of function")
plt.ylabel("y-value of function")
plt.savefig("Plottings.png")
plt.show()
