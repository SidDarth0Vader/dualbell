import numpy as np
import matplotlib.pyplot as mat
#import math
#mat.xkcd()
gam = [1.33, 1.4, 1.67]
ExMa = 24
rosh = 1.4
#noc = int(input("Enter number of characteristic lines: "))
mach = np.arange(0.1,ExMa , 0.1)
#ExMa = float(input("Enter Nozzle Exit Mach number: "))
print(np.tan(np.pi/4))
mac = float(input("Mach number: "))
def Nu(m, rosh):
    #rosh = ratio of specific heats
    vM = np.sqrt((rosh+1)/(rosh-1))*np.arctan(np.sqrt(((rosh-1)*((m**2)-1))/(rosh+1)))-np.arctan(np.sqrt((m**2)-1))
    #return np.degrees(vM)
    return vM
pran = Nu(mac,rosh)
print("nu(M) = ", np.degrees(pran))
def M(nu):
    nu_inf = (np.pi * (np.sqrt(6) - 1)) / 2  # the maximum turning angle... apparently
    y = (nu / nu_inf) ** (2 / 3)
    A = 1.3604
    B = 0.0962
    C = -0.5127
    D = -0.6722
    E = -0.3278
    M = ((1 + A * y + B * (y ** 2) + C * (y ** 3)) / (1 + D * y + E * (y ** 2)))
    return M
print("M(nu(M)) = ", M(pran))
'''
for i in Nu(mach,rosh):
    mat.plot(M(Nu(mach,rosh)),Nu(mach,rosh))
'''
mat.xlabel('Mach')
mat.ylabel('v(M)')
mat.title('Prandtl-Meyer function relation b/w v(M) and Mach')
mat.plot(mach,Nu(mach,rosh))
#mat.show()
# while 1:
#     print()
