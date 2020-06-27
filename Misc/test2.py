import numpy as np
import matplotlib.pyplot as mat
#import math
#mat.xkcd()
rosh = 1.4
mach = np.arange(0,5,0.1)

def areaRatio(M,gam):
    ARatio = (1/M)*(((gam+1)/2)**(-((gam+1)/(2*gam-2))))*((1+(((gam-1)/2)*M**2))**((gam+1)/(2*gam-2)))
    return ARatio

# mat.plot(areaRatio(mach,rosh),mach)
# mat.show()
while 1:
    ar = areaRatio(float(input("Enter M: ")), rosh)
    print(ar)

# x = ((aR1-aR2)/2*np.tan(theta2))*th_rad + x2
# y = self.aR*th_rad
# aR = [1.356,1.483,1.606,1.760,1.935,2.154,2.403]
# contx = [0.238,0.429,0.6585,1.021,1.8535,3.943,