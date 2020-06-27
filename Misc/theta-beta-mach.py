import matplotlib.pyplot as mat
import numpy as np
import math
#mat.xkcd()
mach = [1.1,1.2,1.3,1.4,1.6,1.8,2,2.25,2.5,3,4,6,8,15,20,50]
gam = 1.4
beta = np.arange(0,math.radians(90),0.001)
for M in mach:
    theta = np.arctan(2*(1/np.tan(beta))*(((M**2)*(np.sin(beta)**2)-1)/((M**2)*(gam + np.cos(2*beta))+2)))
    mat.plot(np.degrees(theta),np.degrees(beta))
mat.xlim(0)
mat.rc('grid', linestyle="-", color='black')
mat.xlabel("Theta")
mat.ylabel("Beta")
mat.grid(1)
mat.title("Theta-Beta-Mach relation")
mat.show()


