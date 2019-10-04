from numpy import *
import matplotlib.pyplot as plt
def int_func(x1,x2,y1,y2,z1,z2):
    exponential_part=exp(-4*(sqrt(x1*x1+y1*y1+z1*z1)+sqrt(x2*x2+y2*y2+z2*z2)));
    denom=sqrt((x1-x2)*(x1-x2)+(y1-y2)*(y1-y2)+(z1-z2)*(z1-z2));
    returnval=0;
    if (abs(denom)<= 1e-15):
        pass
    else:
        returnval= exponential_part/denom;
    return returnval;
x=linspace(-5,5,11)
y=[int_func(val,0,0,0,0,0) for val in x]
print(y)
plt.plot(x,y)
plt.show()
