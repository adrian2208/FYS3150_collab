import subprocess, sys
nvals=[10**i for i in range(1,8)] #10 to 10 millions
howmany=int(sys.argv[1]);
print(nvals)
subprocess.call(["rm", "time_info.txt"])
for n in nvals:
    if n<=1000:
        for i in range(howmany):
            subprocess.call(["./oppe", str(n)])
            subprocess.call(["./oppc", str(n)])
            subprocess.call(["./oppb", str(n)])
        subprocess.call(["./oppe", str(n)+" yes"])
        subprocess.call(["./oppc", str(n)+" yes"])
        subprocess.call(["./oppb", str(n)+" yes"])
    else:
        for i in range(howmany):
            subprocess.call(["./oppc", str(n)])
            subprocess.call(["./oppb", str(n)])
        subprocess.call(["./oppc", str(n)+" yes"])
        subprocess.call(["./oppb", str(n)+" yes"])
