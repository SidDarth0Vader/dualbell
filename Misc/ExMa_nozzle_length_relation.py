import numpy as np
import matplotlib.pyplot as mat
from scipy.optimize import curve_fit
import xlwt
import time

#import math
#mat.xkcd()

def M(nu):  # This method of finding inverse Prandtl-Meyer Function was found at http://www.pdas.com/pm.pdf
    nu_inf = (np.pi * (np.sqrt(6) - 1)) / 2  # the maximum turning angle... apparently
    y = (np.radians(nu) / nu_inf) ** (2 / 3)
    A = 1.3604
    B = 0.0962
    C = -0.5127
    D = -0.6722
    E = -0.3278
    #M = ((1 + A * y + B * (y ** 2) + C * (y ** 3)) / (1 + D * y + E * (y ** 2)))
    return ((1 + A * y + B * (y ** 2) + C * (y ** 3)) / (1 + D * y + E * (y ** 2)))

def Nu(m, rosh):
    #rosh = ratio of specific heats
    #vM = np.sqrt((rosh+1)/(rosh-1))*np.arctan(np.sqrt(((rosh-1)*((m**2)-1))/(rosh+1)))-np.arctan(np.sqrt((m**2)-1))
    return np.degrees(np.sqrt((rosh+1)/(rosh-1))*np.arctan(np.sqrt(((rosh-1)*((m**2)-1))/(rosh+1)))-np.arctan(np.sqrt((m**2)-1)))

def areaRatio(M,gam):
    return (1/M)*(((gam+1)/2)**(-((gam+1)/(2*gam-2))))*((1+(((gam-1)/2)*M**2))**((gam+1)/(2*gam-2)))

ExMa = float(3.5)
#ExMa = float(input("Enter Nozzle Exit Mach number: ")) #ExMa = Exit Mach Number
th_rad = 0.03
#th_rad = float(input("Enter radius of cross-section of throat(in metres): ")) #th_rad = Radius at throat(m)
nozzLength = 1.5
#nozzLength = float(input("Enter length of nozzle: "))
rosh = 1.4
#rosh = float(input("Enter ratio of specific heats: ")) #rosh = Ratio of Specific Heats
noc = 25
#noc = int(input("Enter number of characteristic lines: ")) #noc = number of characteristic lines
thetamax = Nu(ExMa,rosh)/2
rem = thetamax%(noc-1)
dtheta = (thetamax - rem)/(noc-1)
i = noc + 1
nop = 0

while i>=2:
    nop = nop + i
    i = i - 1

ptloc_array = []
tmp = noc+1

while len(ptloc_array)<nop:
    i=1
    while i<=tmp and tmp>1:
        if i==tmp:
            ptloc_array.append('wall')
        elif i==1:
            ptloc_array.append('axis')
        else:
            ptloc_array.append('interior')
        i+=1
    tmp-=1

infPtLocX = 0.5 * th_rad * np.sin(np.radians(thetamax+Nu(ExMa,rosh)))
infPtLocY = 1.5 * th_rad - 0.5 * th_rad * np.cos(np.radians(thetamax+Nu(ExMa,rosh)))
print(infPtLocX, infPtLocY, Nu(ExMa, rosh))


class Points:
    def __init__(self, aloc_pno, char1, char2):
        self.pno = aloc_pno+1
        #print(aloc_pno, self.pno)
        self.pointloc = ptloc_array[self.pno - 1]
        self.Kmin = 0
        self.Kplus = 0
        self.theta = 0
        self.p_nu = 0
        if ((self.pointloc == 'axis' or self.pointloc == 'interior') and self.pno <= noc):
            self.theta = self.p_nu = rem + (self.pno-1)*dtheta
            self.Kmin = self.theta + self.p_nu
            self.Kplus = self.theta - self.p_nu
            #print(self.pno, self.pointloc, self.Kmin, self.Kplus)
        elif self.pointloc == 'axis' and self.pno != 1:
            self.Kmin = char1.Kmin
            self.Kplus = -char1.Kmin
            self.theta = 0
            self.p_nu = 0.5 * (self.Kmin - self.Kplus)
            #print(self.pno,self.pointloc,self.Kmin,self.Kplus)
        elif self.pointloc == 'wall':
            self.Kmin = char2.Kmin
            self.Kplus = char2.Kplus
            self.theta = 0.5 * (self.Kmin + self.Kplus)
            self.p_nu = 0.5 * (self.Kmin - self.Kplus)
            #print(self.pno, self.pointloc, self.Kmin, self.Kplus)
        else:
            self.Kmin = char1.Kmin
            self.Kplus = char2.Kplus
            self.theta = 0.5 * (self.Kmin + self.Kplus)
            #print(self.pno, self.pointloc, self.Kmin, self.Kplus)
        self.p_nu = 0.5 * (self.Kmin - self.Kplus)
        self.Ma = M(self.p_nu)
        self.mu = np.degrees(np.arcsin(1/self.Ma))
        self.aR = areaRatio(self.Ma,rosh)
        self.x = 0
        self.y = 0

#Calculating wall point coords using method from Gas Dynamics by Ethirajan Rathakrishnan
        if self.pointloc == 'wall' and self.pno == noc+1:
            self.x = ((self.aR-1)/(2*np.tan(np.radians(thetamax))))*th_rad
            self.y = self.aR*th_rad #- infPtLocY +th_rad
            #print(self.pno,self.aR, self.x, self.y)
        elif self.pointloc == 'wall':
            self.x = ((self.aR-char1.aR)/(2*np.tan(np.radians(char1.theta))))*th_rad + char1.x
            self.y = self.aR*th_rad
            #print(self.pno,self.aR, self.x, self.y)

flag = 0
OG_ExMa = ExMa
truncFlag = 0
count = 1
px = []
py = []
while ExMa<=3*OG_ExMa:
    #print("try ", count)
    thetamax = Nu(ExMa, rosh) / 2
    rem = thetamax % (noc - 1)
    dtheta = (thetamax - rem) / (noc - 1)
    pts=[]
    tmp = noc+1
    for x in range(nop):
        #print(x, tmp, x - tmp, x - 1)
        if ptloc_array[x] == 'axis' and x>noc and x>0:
            tmp-=1
            #print("------------tmp reduction, tmp = ",tmp)

        if x == 0:
            pts.append(Points(x,x,x))
            #print("--first cond pno={}".format(x+1))
        elif x != 0 and x < noc+1:
            pts.append(Points(x, x, pts[x-1]))
            #print("--second cond pno = {}    pno under = {}".format(x+1,x))
        elif x >= noc+1:
            pts.append(Points(x,pts[x-tmp],pts[x-1]))
            #print("--third cond pno = {}  pno left = {} pno under = {} Kmin = {:>6.2f}  Kplus = {:>6.2f}    ".format(x+1,x+1-tmp,x,pts[x-tmp].Kmin, pts[x-tmp].Kplus))
    px.append(pts[nop-1].x)
    py.append(ExMa)
    print(round(ExMa,5))
    ExMa+=0.1

mat.plot(px,py)

mat.axis('scaled')
mat.xlabel("Nozzle length (m)")
mat.ylabel("Exit Mach")
#mat.title("Truncated Ideal Contour Nozzle")
mat.grid()
mat.show()