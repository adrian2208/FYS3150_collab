import numpy as np
import matplotlib.pyplot as plt
infile=open("log_err_c.txt")
infile.readline()
lines=infile.readlines()
a=np.zeros(len(lines))
b=np.zeros(len(lines))
for i in range(len(lines)):
    a[i],b[i]=lines[i].split();
    a[i],b[i]=float(a[i]),float(b[i])
plt.plot(a,b)
plt.xlabel("$log_{10}$(h)")
plt.ylabel("$log_{10}$(relative_error)")
plt.savefig("error.png")
plt.show()
