import numpy as np
import matplotlib.pyplot as mat
from scipy.optimize import curve_fit
import xlwt

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


ExMa = 4.3
# ExMa = float(input("Enter Nozzle Exit Mach number: ")) #ExMa = Exit Mach Number
th_rad = 0.03354
#th_rad = float(input("Enter radius of cross-section of throat(in metres): ")) #th_rad = Radius at throat(m)
downWallRad = 0.5*th_rad
rosh = 1.4
#rosh = float(input("Enter ratio of specific heats: ")) #rosh = Ratio of Specific Heats
noc = 30
#noc = int(input("Enter number of characteristic lines: ")) #noc = number of characteristic lines
thetamax = Nu(ExMa,rosh)/2
print(thetamax)
#rem = thetamax%(noc-1)
dtheta = (thetamax)/(noc-1)
i = noc + 1
nop = 0

print("dtheta = ", dtheta)

while i>=2:
    nop = nop + i
    i = i - 1

ptloc_array = []
ptzo_array = []
ptzo = 1
tmp = noc+1

while len(ptloc_array)<nop:
    i=1
    while i<=tmp and tmp>1:
        if i==tmp:
            ptloc_array.append('wall')
            ptzo_array.append(ptzo)
            ptzo+=1
        elif i==1:
            ptloc_array.append('axis')
            ptzo_array.append(ptzo)
        else:
            ptloc_array.append('interior')
            ptzo_array.append(ptzo)
        i+=1
    tmp-=1

infPtLocX = 0.5 * th_rad * np.sin(np.radians(thetamax))
infPtLocY = 1.5 * th_rad - 0.5 * th_rad * np.cos(np.radians(thetamax))
print(infPtLocX, infPtLocY, Nu(ExMa, rosh))


class Points:
    def __init__(self, aloc_pno, char1, char2):
        self.pno = aloc_pno+1
        self.pzo = ptzo_array[aloc_pno]
        #print(aloc_pno, self.pno)
        self.pointloc = ptloc_array[aloc_pno]
        self.Kmin = 0
        self.Kplus = 0
        self.theta = 0
        self.p_nu = 0
        if ((self.pointloc == 'axis' or self.pointloc == 'interior') and self.pno <= noc):
            self.theta = self.p_nu = (self.pno-1)*dtheta
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
        self.aR = areaRatio(self.Ma, rosh)
        self.x = 0
        self.y = 0

#Calculating wall point coords using method from Gas Dynamics by Ethirajan Rathakrishnan
        # if self.pointloc == 'wall' and self.pno == noc+1:
        #     self.x = ((self.aR-1)/(2*np.tan(np.radians(thetamax))))*th_rad
        #     self.y = self.aR*th_rad #- infPtLocY +th_rad
        #     #print(self.pno,self.aR, self.x, self.y)
        # elif self.pointloc == 'wall':
        #     self.x = ((self.aR-char1.aR)/(2*np.tan(np.radians(char1.theta))))*th_rad + char1.x
        #     self.y = self.aR*th_rad
        if self.pointloc == 'wall' and self.pno == noc+1:
            self.y = np.sqrt(self.aR * (th_rad ** 2))  # - infPtLocY +th_rad
            self.x = (self.y-infPtLocY)/np.tan(np.radians(thetamax))+infPtLocX
            # print(self.pno,self.aR, self.x, self.y)
        elif self.pointloc == 'wall':
            self.y = np.sqrt(self.aR*(th_rad**2))
            self.x = (self.y-char1.y)/np.tan(np.radians(char1.theta)) + char1.x
            # print(self.pno,self.aR, self.x, self.y)
# Characteristic drawing
        if self.pointloc == 'axis':
            self.x = (th_rad+downWallRad)*np.tan(np.radians(self.pzo*dtheta))

        elif self.pointloc == 'interior' and self.pno <= noc:
            tanth1 = np.tan(np.radians((-(90-self.pno*dtheta)+(self.theta-self.mu))/2))
            tanth2 = np.tan(np.radians((char2.theta+char2.mu+self.theta+self.mu)/2))
            denom = 1 - tanth1/tanth2
            self.x = (char2.x + (1/tanth2)*(th_rad+downWallRad - char2.y))/denom
            self.y = (th_rad+downWallRad + tanth1*(char2.x - char2.y/tanth2))/denom

        elif self.pointloc == 'interior' and self.pno > noc:
            tanth1 = np.tan(np.radians((char1.theta+self.theta)*0.5-(char1.mu+self.mu)*0.5))
            tanth2 = np.tan(np.radians((char1.theta+self.theta)*0.5+(char1.mu+self.mu)*0.5))
            denom = 1 - tanth1 / tanth2
            self.x = (char2.x + (1/tanth2)*(char1.y - char2.y - char1.x/tanth1))/denom
            self.y = (char1.y + tanth1*(char2.x - char1.x - char2.y/tanth2))/denom


pts=[]
tmp = noc+1
ptzo = 1
for x in range(nop):
    #print(x, tmp, x - tmp, x - 1)
    if ptloc_array[x] == 'axis' and x>noc and x>0:
        tmp-=1
        ptzo += 1
        #print("------------tmp reduction, tmp = ",tmp)

    if x == 0:
        pts.append(Points(x,x,x))
        pts[x].pzo = ptzo
        #print("--first cond pno={}".format(x+1))
    elif x != 0 and x < noc+1:
        pts.append(Points(x, x, pts[x-1]))
        pts[x].pzo = ptzo
        #print("--second cond pno = {}    pno under = {}".format(x+1,x))
    elif x >= noc+1:
        pts.append(Points(x,pts[x-tmp],pts[x-1]))
        pts[x].pzo = ptzo
        #print("--third cond pno = {}  pno left = {} pno under = {} Kmin = {:>6.2f}  Kplus = {:>6.2f}    ".format(x+1,x+1-tmp,x,pts[x-tmp].Kmin, pts[x-tmp].Kplus))



    #print("Everyloop pno = {}   pno left = {}  pno under = {}".format(x+1,x+1-tmp,x))

def displayTable(): #Literally does what it looks like it's gonna do
    print("-----------|------------|------------|-----------|-----------|----------|------------|----------")
    print("  Point no.|       Kmin |      Kplus |  θ(Theta) |     ν(Nu) |        M |      μ(mu) |      pzo ")
    print("-----------|------------|------------|-----------|-----------|----------|------------|----------")
    for x in range(nop):
        print(" {0:>9} |   {1:>8.3f} |   {2:>8.3f} |  {3:>8.3f} |  {4:>8.3f} | {5:>8.3f} |  {6:>8.3f} |   {7:>9}".format(pts[x].pno,pts[x].Kmin,pts[x].Kplus,pts[x].theta,pts[x].p_nu,pts[x].Ma,pts[x].mu,pts[x].pzo))
    return 0

def idealConNozz(): #This function draws the ideal contour nozzle
    conx = []
    cony = []
    xplot = np.arange(0,infPtLocX, th_rad/100)
    for a in range(len(xplot)):
        conx.append(xplot[a])
        cony.append(round(-np.sqrt((0.5*th_rad)**2-xplot[a]**2) + 1.5*th_rad,4))

    i = 0
    while i < nop:
        if pts[i].pointloc == 'wall':
            conx.append(pts[i].x)
            cony.append(pts[i].y)
        i+=1
    mat.plot(conx, cony, color='black')
    print(conx)
    print(cony)
    return 0

def drawChar():
    colour = 'red'
    for h in range(nop):
        if pts[h].pointloc == 'axis':
            mat.plot([0, pts[h].x], [th_rad+downWallRad, pts[h].y], colour)
        elif pts[h].pointloc == 'interior' or pts[h].pointloc == 'wall':
            mat.plot([pts[h-1].x, pts[h].x], [pts[h-1].y, pts[h].y], colour)

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

def curveGen(): #Attempts to generate smooth regression curve of nozzle points
    x = []
    y = []
    xexpo = []
    yexpo = []
    i = 0
    xplot = np.arange(-50, pts[nop-1].x, 0.01)
    inExpReg_X = np.arange(0, infPtLocX, th_rad/100)
    for a in range(len(inExpReg_X)):
        x.append(inExpReg_X[a])
        xexpo.append(inExpReg_X[a])
        yexpo.append(round(-np.sqrt((0.5*th_rad)**2-inExpReg_X[a]**2) + 1.5*th_rad,4))
        y.append(round(-np.sqrt((0.5*th_rad)**2-inExpReg_X[a]**2) + 1.5*th_rad,4))

    xexpo.append(pts[noc+1].x)
    yexpo.append(pts[noc+1].y)
    #inExpReg_Y = - np.sqrt((0.5 * th_rad) ** 2 - inExpReg_X ** 2) + 1.5 * th_rad
    #print(len(xplot))

    def sigmoid(x,a,b):
        return ((pts[nop-1].y/(1+np.e**(infPtLocX+a*x)))+b)

    def exponential(x,a,b):
        return (a*np.e**(b*x))

    def parabolic(x,a,b,c):
        return (a+(b*x)+(c*(x**2)))

    def cubic(x, a, b, c, d):
        return (a+(b*x)+(c*(x**2))+(d*(x**3)))

    def quadratic(x, a, b, c, d, e):
        return (a+(b*x)+(c*(x**2))+(d*(x**3))+(e*(x**4)))

    while i < nop:
        # if pts[i].pno == 1:
        #     x.append(0)
        #     y.append(th_rad)
        #
        # if len(x) == 1:
        #     x.append(round(infPtLocX,4))
        #     y.append(round(infPtLocY,4))

        if pts[i].pointloc == 'wall':
            x.append(round(pts[i].x,4))
            y.append(round(pts[i].y,4))

        i += 1

    print(x)
    print(y)

    print(len(x))

    expDAT,pcov = curve_fit(exponential,xexpo,yexpo)
    sigDAT,pcov = curve_fit(sigmoid,x,y)
    paraDAT,pcov = curve_fit(parabolic,x,y)
    cubicDAT,pcov = curve_fit(cubic,x,y)
    quadDAT,pcov = curve_fit(quadratic,x,y)

    asig = sigDAT[0]
    bsig = sigDAT[1]
    ae2 = expDAT[0]
    be2 = expDAT[1]
    a2 = paraDAT[0]
    b2 = paraDAT[1]
    c2 = paraDAT[2]
    a3 = cubicDAT[0]
    b3 = cubicDAT[1]
    c3 = cubicDAT[2]
    d3 = cubicDAT[3]
    a4 = quadDAT[0]
    b4 = quadDAT[1]
    c4 = quadDAT[2]
    d4 = quadDAT[3]
    e4 = quadDAT[4]
    y_sig = sigmoid(xplot, asig, bsig)
    y_ex = exponential(xplot, ae2, be2)
    y2 = parabolic(xplot, a2, b2, c2)
    y3 = cubic(xplot, a3, b3, c3, d3)
    y4 = quadratic(xplot, a4, b4, c4, d4, e4)
    #mat.plot(xexpo, y_ex, 'yellow')
    #mat.plot(xplot,y_sig, 'black')
    mat.plot(xplot, y2, 'orange')
    mat.plot(xplot, y3, 'blue')
    #mat.plot(xplot, y4, 'red')

    # print(a2,b2,c2)
    # print(a3,b3,c3,d3)
    # print(a4,b4,c4,d4,e4)
    # inExpReg_X = np.arange(0, infPtLocX, 0.01)
    # inExpReg_Y = - np.sqrt((0.5*th_rad)**2 - inExpReg_X**2)+1.5*th_rad
    # mat.plot(inExpReg_X,inExpReg_Y,'orange')
    mat.plot(xplot,y2,'r')
    mat.plot(xplot,y3,'b')
    mat.plot(xplot,y4,'g')

def XLimporter():
    dir = "D:\Study Files\Aerospace Engg\Project\Output\\"
    #dir = input("Enter save directory: ")
    wb = xlwt.Workbook()
    xplot = np.arange(0, infPtLocX, th_rad / 100)
    ws = wb.add_sheet("ICN coordinates")
    ws.write(0, 0, "Point Number")
    ws.write(0, 1, "X")
    ws.write(0, 2, "Y")
    j = 0
    while j<len(xplot):
        ws.write(j+1, 0, j+1)
        ws.write(j+1, 1, xplot[j])
        ws.write(j+1, 2, round(-np.sqrt((0.5*th_rad)**2-xplot[j]**2) + 1.5*th_rad,4))
        j+=1

    i = 0
    k = 0
    while i < nop:
        if pts[i].pointloc == 'wall':
            ws.write(j+k+1, 0, j+k+1)
            ws.write(j+k+1, 1, round(pts[i].x,4))
            ws.write(j+k+1, 2, round(pts[i].y,4))
            k+=1
        i+=1

    wb.save("{0:}Ideal Con Nozzle, ExMa={1:}, Throat rad={2:}, No. of char={3:}.xls".format(dir,ExMa,th_rad,noc))
    return 0

def everyLoc():
    print("-----------|------------|------------|------------|----------")
    print("  Point no.|          X |          Y |      Theta |       mu ")
    print("-----------|------------|------------|------------|----------")
    for x in range(nop):
        print(" {0:>9} |   {1:>8.3f} |   {2:>8.3f} |   {3:>8.3f} |   {4:>8.3f} | {5:>3}".format(pts[x].pno, pts[x].x, pts[x].y, pts[x].theta, pts[x].mu, pts[x].pointloc))
    return 0

# actAR = (np.pi*pts[nop-1].y**2)/(np.pi*(th_rad*2)**2)
# print(np.pi, pts[nop-1].y)
# print("Actual Area Ratio = {:.4f}".format(actAR))
# print("Target Area Ratio = {:.4f}".format(areaRatio(ExMa,rosh)))

#Enable or disable the functions down below by removing the '#' or putting the '#' in front of the functions respectively

#drawChar()
idealConNozz()      #This function draws the minimum length nozzle
#displayTable()      #This table shows flow properties along the points in the flow
#dispLocation()      #Shows coordinates of all wallpoints as well as characteristic intersection points in Cartesian coords sys
curveGen()          #Attempts to generate smooth regression curve of nozzle points but fails miserably
#XLimporter()        #Imports coords to an excel file
#everyLoc()

tempx = []
tempy = []
# for i in range(nop):
#     if pts[i].pointloc == 'wall':
#         tempx.append(pts[i].x)
#         tempy.append(pts[i].y)
# mat.plot(tempx,tempy,'ro')
# print(pts[noc].x)

mat.axis('scaled')
#mat.plot([0,pts[nop-1].x],[pts[nop-1].y,pts[nop-1].y], label = "{}m".format(pts[nop-1].x), color = 'Blue', linewidth = 0.7)
#mat.plot([pts[nop-1].x,pts[nop-1].x],[0,pts[nop-1].y], label = "{}m".format(pts[nop-1].y), color = 'Blue', linewidth = 0.7)
mat.xlabel("X-axis")
mat.ylabel("Y-axis")
mat.title("Ideal Contour Nozzle")
mat.ylim(0,(pts[nop-1].y+0.1*pts[nop-1].y))
mat.xlim(0,(pts[nop-1].x+0.1*pts[nop-1].x))
mat.grid()

mat.show()
