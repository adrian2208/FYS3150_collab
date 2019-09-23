import numpy as np
import matplotlib.pyplot as plt
import sys
def density(rho,func):
    return func**2
infilename=sys.argv[1]
infile=open(infilename)
n=int(infile.readline().split()[-1])
print(n)
rhomax=float(infile.readline().split()[-1])
rhomin=float(infile.readline().split()[-1])
eigenval=[float(i) for i in infile.readline().split()];
values=infile.readlines();
for i in range(len(values)):
    values[i]=list(map(float, values[i].split()))

eigenvec=np.array(values)
print(eigenvec)
i=np.argsort(eigenval)
eigenvec=eigenvec[:,i]
eigenvec=np.transpose(eigenvec)
print(eigenvec)
rho=np.linspace(rhomin,rhomax,n)
for i in range(3):
    densityp=density(rho,eigenvec[i])
    densityp=densityp/sum(densityp)
    plt.plot(rho,densityp,label="%d"%i)
plt.legend()
plt.show()
