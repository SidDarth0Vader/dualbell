import numpy as np
import matplotlib.pyplot as mat
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

ExMa = 2.4
#ExMa = float(input("Enter Nozzle Exit Mach number: ")) #ExMa = Exit Mach Number
th_rad = 1
#th_rad = float(input("Enter radius of cross-section of throat(in metres): ")) #th_rad = Radius at throat(m)
rosh = 1.4
#rosh = float(input("Enter ratio of specific heats: ")) #rosh = Ratio of Specific Heats
noc = 7
noc = int(input("Enter number of characteristic lines: ")) #noc = number of characteristic lines
mach = np.arange(0,ExMa , 0.1)
thetamax = Nu(ExMa,rosh)/2
print(thetamax)
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

class Points:
    def __init__(self, pno, char1, char2):
        self.pno = pno+1
        self.pointloc = ptloc_array[self.pno - 1]
        self.Kmin = 0
        self.Kplus = 0
        self.theta = 0
        self.p_nu = 0
        if ((self.pointloc == 'axis' or self.pointloc == 'interior') and self.pno <= noc):
            self.theta = self.p_nu =  rem + (self.pno-1)*dtheta
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
            self.y = self.aR*th_rad
            #print(self.pno,self.aR, self.x, self.y)
        elif self.pointloc == 'wall':
            self.x = ((self.aR-char1.aR)/(2*np.tan(np.radians(char1.theta))))*th_rad + char1.x
            self.y = self.aR*th_rad
            #print(self.pno,self.aR, self.x, self.y)

#Calculating coors of grid points by a wonky method derived by me
        # if self.pno == 1:
        #     self.x = th_rad*np.tan(np.radians(thetamax))
        #     self.y = 0
        # elif self.pointloc == 'interior' and self.pno > noc+1:
        #     a = np.tan(np.radians(char1.theta-char1.mu))
        #     b = np.tan(np.radians(char2.theta+char2.mu))
        #     self.x = -((((char2.x*a-char1.y+char2.y)/b)+char2.x)/(1-a/b))
        #     self.y = -(((a*((char2.y/b)-char2.x+char1.x))+char1.y)/(1-a/b))
        # elif self.pointloc == 'interior' and self.pno < noc+1:
        #     a = np.tan(np.radians(90-thetamax+dtheta*self.pno))
        #     b = np.tan(np.radians(char2.theta+char2.mu))
        #     self.x = -((((char2.x*a-th_rad+char2.y)/b)+char2.x)/(1-a/b))
        #     self.y = -(((a*((char2.y/b)-char2.x))+th_rad)/(1-a/b))
        # elif self.pointloc == 'axis':
        #     b = np.tan(np.radians(char2.theta + char2.mu))
        #     self.x = -(char1.y/b)+char1.x
        #     self.y = 0

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

    #print("Everyloop pno = {}   pno left = {}  pno under = {}".format(x+1,x+1-tmp,x))

def displayTable(): #Literally does what it looks like it's gonna do
    print("-----------|------------|------------|-----------|-----------|----------|------------")
    print("  Point no.|       Kmin |      Kplus |     Theta |        Nu |        M |        mu ")
    print("-----------|------------|------------|-----------|-----------|----------|------------")
    for x in range(nop):
        print(" {0:>9} |   {1:>8.3f} |   {2:>8.3f} |  {3:>8.3f} |  {4:>8.3f} | {5:>8.3f} |  {6:>8.3f}".format(pts[x].pno,pts[x].Kmin,pts[x].Kplus,pts[x].theta,pts[x].p_nu,pts[x].Ma,pts[x].mu))
    return 0

def minLength(): #This function draws the minimum length nozzle
    conx = []
    cony = []
    i = 0
    while i < nop:
        if pts[i].pno == 1:
            conx.append(0)
            cony.append(th_rad)

        elif pts[i].pointloc == 'wall':
            conx.append(pts[i].x)
            cony.append(pts[i].y)
        mat.plot(conx,cony, color = 'black')
        i+=1
    # print(conx)
    # print(cony)
    return 0

def drawChar(): #Draws the characteristic lines #WIP
    tmp = noc + 1
    i = 0
    while i<nop:

        if ptloc_array[x] == 'axis' and x > noc and x > 0:
            tmp -= 1

        if pts[i].pno == 1:
            mat.plot([0,pts[i].x],[th_rad,pts[i].y], color = 'Black', linewidth = '0.5')
        elif pts[i].pointloc == 'interior' and i < noc+1:
            mat.plot([0,pts[i].x],[th_rad,pts[i].y], color = 'Black', linewidth = '0.5')
            mat.plot([pts[i-1].x, pts[i].x], [pts[i-1].y, pts[i].y], color='Black', linewidth='0.5')
        elif pts[i].pointloc == 'interior' and i > noc+1:
            mat.plot([pts[i - tmp].x, pts[i].x], [pts[i - tmp].y, pts[i].y], color='Black', linewidth='0.5')
            mat.plot([pts[i - 1].x, pts[i].x], [pts[i - 1].y, pts[i].y], color='Black', linewidth='0.5')
        elif pts[i].pointloc == 'wall':
            mat.plot([pts[i - 1].x, pts[i].x], [pts[i - 1].y, pts[i].y], color='Black', linewidth='0.5')
        i+=1
    return 0

def dispLocation(): #Just for Debugging
    print("-----------|------------|------------|------------|--")
    print("  Point no.|          X |          Y |      Theta | ")
    print("-----------|------------|------------|------------|--")
    for x in range(nop):
        #if pts[x].pointloc != 'wall':
           # print(" {0:>9} |   {1:>8.3f} |   {2:>8.3f}  |".format(pts[x].pno,pts[x].x,pts[x].y))
        if pts[x].pointloc == 'wall':
            print(" {0:>9} |   {1:>8.3f} |   {2:>8.3f} |   {3:>8.3f} | Wallpoint".format(pts[x].pno, pts[x].x, pts[x].y,pts[x].theta))
    return 0

def curveGen(): #Truncated nozzle generator
    x = []
    y = []
    xy = []
    x2y = []
    x2 = []
    x3 = []
    x4 = []
    i = s_x = s_y = s_x2 = s_x3 = s_x4 = s_x5 = s_x6 = s_x7 = s_x8 = s_xy = s_x2y = s_x3y = s_x4y = 0
    j = -1
    while i < nop:
        # if pts[i].pno == 1:
        #     x.append(0)
        #     y.append(1)
        #     j+=1


        if pts[i].pointloc == 'wall':
            x.append(round(pts[i].x,4))
            y.append(round(pts[i].y,4))
            s_x += x[j]
            s_y += y[j]
            s_x2 += round(x[j] ** 2, 4)
            s_x3 += round(x[j] ** 3, 4)
            s_x4 += round(x[j] ** 4, 4)
            s_x5 += round(x[j] ** 5, 4)
            s_x6 += round(x[j] ** 6, 4)
            s_x7 += round(x[j] ** 7, 4)
            s_x8 += round(x[j] ** 8, 4)
            s_xy += round((x[j] * y[j]), 4)
            s_x2y += round(((x[j] ** 2) * y[j]), 4)
            s_x3y += round(((x[j] ** 3) * y[j]), 4)
            s_x4y += round(((x[j] ** 4) * y[j]), 4)
            j+=1
            # print(x)
            # print(y)

        # s_x += x[j]
        # s_y += y[j]
        # s_x2 += round(x[j] ** 2, 4)
        # s_x3 += round(x[j] ** 3, 4)
        # s_x4 += round(x[j] ** 4, 4)
        # s_x5 += round(x[j] ** 5, 4)
        # s_x6 += round(x[j] ** 6, 4)
        # s_x7 += round(x[j] ** 7, 4)
        # s_x8 += round(x[j] ** 8, 4)
        # s_xy += round((x[j] * y[j]), 4)
        # s_x2y += round(((x[j] ** 2) * y[j]), 4)
        # s_x3y += round(((x[j] ** 3) * y[j]), 4)
        # s_x4y += round(((x[j] ** 4) * y[j]), 4)
        i += 1

    n = len(x)

    a2 = (np.linalg.det(np.array([[s_y,s_x,s_x2],[s_xy,s_x2,s_x3],[s_x2y,s_x3,s_x4]])))/(np.linalg.det(np.array([[n,s_x,s_x2],[s_x,s_x2,s_x3],[s_x2,s_x3,s_x4]])))
    b2 = (np.linalg.det(np.array([[n,s_y,s_x2],[s_x,s_xy,s_x3],[s_x2,s_x2y,s_x4]])))/(np.linalg.det(np.array([[n,s_x,s_x2],[s_x,s_x2,s_x3],[s_x2,s_x3,s_x4]])))
    c2 = (np.linalg.det(np.array([[n,s_x,s_y],[s_x,s_x2,s_xy],[s_x2,s_x3,s_x2y]])))/(np.linalg.det(np.array([[n,s_x,s_x2],[s_x,s_x2,s_x3],[s_x2,s_x3,s_x4]])))

    a3 = (np.linalg.det(np.array([[s_y,s_x,s_x2,s_x3],[s_xy,s_x2,s_x3,s_x4],[s_x2y,s_x3,s_x4,s_x5],[s_x3y,s_x4,s_x5,s_x6]])))/(np.linalg.det(np.array([[n,s_x,s_x2,s_x3],[s_x,s_x2,s_x3,s_x4],[s_x2,s_x3,s_x4,s_x5],[s_x3,s_x4,s_x5,s_x6]])))
    b3 = (np.linalg.det(np.array([[n,s_y,s_x2,s_x3],[s_x,s_xy,s_x3,s_x4],[s_x2,s_x2y,s_x4,s_x5],[s_x3,s_x3y,s_x5,s_x6]])))/(np.linalg.det(np.array([[n,s_x,s_x2,s_x3],[s_x,s_x2,s_x3,s_x4],[s_x2,s_x3,s_x4,s_x5],[s_x3,s_x4,s_x5,s_x6]])))
    c3 = (np.linalg.det(np.array([[n,s_x,s_y,s_x3],[s_x,s_x2,s_xy,s_x4],[s_x2,s_x3,s_x2y,s_x5],[s_x3,s_x4,s_x3y,s_x6]])))/(np.linalg.det(np.array([[n,s_x,s_x2,s_x3],[s_x,s_x2,s_x3,s_x4],[s_x2,s_x3,s_x4,s_x5],[s_x3,s_x4,s_x5,s_x6]])))
    d3 = (np.linalg.det(np.array([[n,s_x,s_x2,s_y],[s_x,s_x2,s_x3,s_xy],[s_x2,s_x3,s_x4,s_x2y],[s_x3,s_x4,s_x5,s_x3y]])))/(np.linalg.det(np.array([[n,s_x,s_x2,s_x3],[s_x,s_x2,s_x3,s_x4],[s_x2,s_x3,s_x4,s_x5],[s_x3,s_x4,s_x5,s_x6]])))

    a4 = (np.linalg.det(np.array([[s_y,s_x,s_x2,s_x3,s_x4],[s_xy,s_x2,s_x3,s_x4,s_x5],[s_x2y,s_x3,s_x4,s_x5,s_x6],[s_x3y,s_x4,s_x5,s_x6,s_x7],[s_x4y,s_x5,s_x6,s_x7,s_x8]])))/(np.linalg.det(np.array([[n,s_x,s_x2,s_x3,s_x4],[s_x,s_x2,s_x3,s_x4,s_x5],[s_x2,s_x3,s_x4,s_x5,s_x6],[s_x3,s_x4,s_x5,s_x6,s_x7],[s_x4,s_x5,s_x6,s_x7,s_x8]])))
    b4 = (np.linalg.det(np.array([[n,s_y,s_x2,s_x3,s_x4],[s_x,s_xy,s_x3,s_x4,s_x5],[s_x2,s_x2y,s_x4,s_x5,s_x6],[s_x3,s_x3y,s_x5,s_x6,s_x7],[s_x4,s_x4y,s_x6,s_x7,s_x8]])))/(np.linalg.det(np.array([[n,s_x,s_x2,s_x3,s_x4],[s_x,s_x2,s_x3,s_x4,s_x5],[s_x2,s_x3,s_x4,s_x5,s_x6],[s_x3,s_x4,s_x5,s_x6,s_x7],[s_x4,s_x5,s_x6,s_x7,s_x8]])))
    c4 = (np.linalg.det(np.array([[n,s_x,s_y,s_x3,s_x4],[s_x,s_x2,s_xy,s_x4,s_x5],[s_x2,s_x3,s_x2y,s_x5,s_x6],[s_x3,s_x4,s_x3y,s_x6,s_x7],[s_x4,s_x5,s_x4y,s_x7,s_x8]])))/(np.linalg.det(np.array([[n,s_x,s_x2,s_x3,s_x4],[s_x,s_x2,s_x3,s_x4,s_x5],[s_x2,s_x3,s_x4,s_x5,s_x6],[s_x3,s_x4,s_x5,s_x6,s_x7],[s_x4,s_x5,s_x6,s_x7,s_x8]])))
    d4 = (np.linalg.det(np.array([[n,s_x,s_x2,s_y,s_x4],[s_x,s_x2,s_x3,s_xy,s_x5],[s_x2,s_x3,s_x4,s_x2y,s_x6],[s_x3,s_x4,s_x5,s_x3y,s_x7],[s_x4,s_x5,s_x6,s_x4y,s_x8]])))/(np.linalg.det(np.array([[n,s_x,s_x2,s_x3,s_x4],[s_x,s_x2,s_x3,s_x4,s_x5],[s_x2,s_x3,s_x4,s_x5,s_x6],[s_x3,s_x4,s_x5,s_x6,s_x7],[s_x4,s_x5,s_x6,s_x7,s_x8]])))
    e4 = (np.linalg.det(np.array([[n,s_x,s_x2,s_x3,s_y],[s_x,s_x2,s_x3,s_x4,s_xy],[s_x2,s_x3,s_x4,s_x5,s_x2y],[s_x3,s_x4,s_x5,s_x6,s_x3y],[s_x4,s_x5,s_x6,s_x7,s_x4y]])))/(np.linalg.det(np.array([[n,s_x,s_x2,s_x3,s_x4],[s_x,s_x2,s_x3,s_x4,s_x5],[s_x2,s_x3,s_x4,s_x5,s_x6],[s_x3,s_x4,s_x5,s_x6,s_x7],[s_x4,s_x5,s_x6,s_x7,s_x8]])))

    print(a2,b2,c2)
    print(a3,b3,c3,d3)
    print(a4,b4,c4,d4,e4)
    xplot = np.arange(0, pts[nop-1].x, 0.01)
    y2 = a2 + b2*xplot + c2*xplot**2
    y3 = a3 + b3*xplot + c3*xplot**2 + d3*xplot**3
    y4 = a4 + b4*xplot + c4*xplot**2 + d4*xplot**3 + e4*xplot**4
    mat.plot(xplot,y2,'r')
    mat.plot(xplot,y3,'b')
    mat.plot(xplot,y4,'g')


# actAR = (np.pi*pts[nop-1].y**2)/(np.pi*(th_rad*2)**2)
# print(np.pi, pts[nop-1].y)
# print("Actual Area Ratio = {:.4f}".format(actAR))
# print("Target Area Ratio = {:.4f}".format(areaRatio(ExMa,rosh)))

#Enable or disable the functions down below by removing the '#' or putting the '#' in front of the functions respectively

#minLength()      #This function draws the minimum length nozzle
#displayTable()   #This table shows flow properties along the points in the flow
#drawChar()       #This function is still in development. Probably won't ever be developed
dispLocation()   #Shows coordinates of all wallpoints as well as characteristic intersection points in Cartesian coords sys
curveGen()       #This function shows the shape of the nozzle

tempx = []
tempy = []
i=0
for i in range(nop):
    tempx.append(pts[i].x)
    tempy.append(pts[i].y)
mat.plot(tempx,tempy,'ro')
print(pts[noc].x)
x = np.arange(0,pts[noc].x,0.01)
y = 1 + 0.3097*x**2
mat.plot(x,y,'black')
mat.axis('scaled')
#mat.plot([0,pts[nop-1].x],[pts[nop-1].y,pts[nop-1].y], label = "{}m".format(pts[nop-1].x), color = 'Blue', linewidth = 0.7)
#mat.plot([pts[nop-1].x,pts[nop-1].x],[0,pts[nop-1].y], label = "{}m".format(pts[nop-1].y), color = 'Blue', linewidth = 0.7)
mat.xlabel("X-axis")
mat.ylabel("Y-axis")
mat.title("Minimum length Nozzle")
# x = np.arange(0,pts[noc+1].x,0.01)
y = 1 + 0.3097 * x ** 2
#mat.plot(x,y,'orange')
# mat.ylim(0,(pts[nop-1].y+0.1*pts[nop-1].y))
# mat.xlim(0,(pts[nop-1].x+0.1*pts[nop-1].x))
mat.grid()
mat.show()