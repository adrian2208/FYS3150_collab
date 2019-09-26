import subprocess, sys
nvals=[800,500,400,400]
omegas=[0.01,0.5,1,5]
rhomax=[40,8,6,4]
howmany=int(1);
fast=True
test=True
try:
    nvals[0]=int(sys.argv[1]);
    fast=bool(int(sys.argv[2]));
    test=bool(int(sys.argv[3]));
except:
    pass
filenames=["vecop","solve_one_electron","solve_two_electrons","solve_harmonic_oscillator"]
for filename in filenames:
    #subprocess.call(["c++", str(" -c "+filename+".cpp")])
    subprocess.call("c++ -c "+filename+".cpp",shell=True)
done_filenames=["one_electron","two_electrons","harmonic_oscillator"]
operator="-o"
if test:
    subprocess.call("c++ -c tests.cpp",shell=True)
    subprocess.call("c++ -c tests_main.cpp",shell=True)
    subprocess.call("c++ -o testing tests.o tests_main.o vecop.o",shell=True)
if fast:
    operator="-O3 -o"
for i in range(len(done_filenames)):
    subprocess.call("c++ "+operator+ " "+done_filenames[i]+" "+filenames[0]+".o "+filenames[i+1]+".o",shell=True)
if test:
    subprocess.call(["./testing",""])

for k in range(howmany):
    for file in done_filenames:
        if file != "harmonic_oscillator":
            for ir in range(len(omegas)):
                subprocess.call(("./"+file+" "+str(nvals[ir])+" "+str(rhomax[ir])+" "+str(omegas[ir])),shell=True)
        else:
            subprocess.call(["./"+file, str(nvals[-1])])
