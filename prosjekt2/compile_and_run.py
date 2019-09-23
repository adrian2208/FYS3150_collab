import subprocess, sys
nvals=[400]
omegas=[0.01,0,5,1,5]
howmany=1;
fast=True
test=True
try:
    howmany=int(sys.argv[1]));
    fast=bool(sys.argv[2]);
    test=bool(sys.argv[3]);
except:
    pass
filenames=["vecop","solve_one_electron","solve_two_electrons","solve_harmonic_oscillator"]
for filename in filenames:
    subprocess.call(["c++", "-c "+filename+".cpp"])
done_filenames=["one_electron","two_electrons","harmonic_oscillator"]
operator="-o"
if fast:
    operator="-O3 -o"
for i in range(len(done_filenames)):
    subprocess.call(["c++", operator+" "+done_filenames[i]+" "+filenames[0]+" "+filenames[i+1]])
for k in range(howmany):
    for n in nvals:
        for file in done_filenames:
            if file.equals("two_electrons"):
                for omega in omegas:
                    subprocess.call(["./"+file, str(n)+" "+str(omega)])
            else:
                subprocess.call(["./"+file, str(n)])
txtfilenames=["solutions_"+filename+".txt" for filename in done_filenames]
