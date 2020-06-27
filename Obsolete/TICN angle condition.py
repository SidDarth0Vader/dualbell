import numpy as np
import matplotlib.pyplot as mat
from scipy.optimize import curve_fit
import xlwt
import time

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
#ExMa = float(input("Enter Nozzle Exit Mach number: ")) #ExMa = Exit Mach Number
th_rad = 0.03354
#th_rad = float(input("Enter radius of cross-section of throat(in metres): ")) #th_rad = Radius at throat(m)
nozzExAng = 10
#nozzExAng = float(input("Enter nozzle exit angle(in degrees): ")) #nozzExAng = Nozzle Exit Angle(degrees)
rosh = 1.4
#rosh = float(input("Enter ratio of specific heats: ")) #rosh = Ratio of Specific Heats
noc = 30
#noc = int(input("Enter number of characteristic lines: ")) #noc = number of characteristic lines
downWallRad = 0.5*th_rad
thetamax = Nu(ExMa,rosh)/2
dtheta = round((thetamax)/(noc-1),3)
i = noc + 1
nop = 0

while i>=2:
    nop = nop + i
    i = i - 1

ptloc_array = []
ptzo_array = []
tmp = noc+1
ptzo = 1
while len(ptloc_array)<nop:
    i=1
    while i<=tmp and tmp>1:
        if i==tmp:
            ptloc_array.append('wall')
            ptzo_array.append(ptzo)
            ptzo += 1
        elif i==1:
            ptloc_array.append('axis')
            ptzo_array.append(ptzo)
        else:
            ptloc_array.append('interior')
            ptzo_array.append(ptzo)
        i+=1

    tmp-=1

# infPtLocX = 0.5 * th_rad * np.sin(np.radians(thetamax+Nu(ExMa,rosh)))
# infPtLocY = 1.5 * th_rad - 0.5 * th_rad * np.cos(np.radians(thetamax+Nu(ExMa,rosh)))
infPtLocX = downWallRad * np.sin(np.radians(thetamax))
infPtLocY = th_rad+downWallRad - downWallRad * np.cos(np.radians(thetamax))
#print(infPtLocX, infPtLocY, Nu(ExMa, rosh))


class Points:
    def __init__(self, aloc_pno, char1, char2):
        self.pno = aloc_pno+1
        #print(aloc_pno, self.pno)
        self.pzo = ptzo_array[aloc_pno]
        self.pointloc = ptloc_array[aloc_pno]
        self.Kmin = 0
        self.Kplus = 0
        self.theta = 0
        self.p_nu = 0
        self.x = 0
        self.y = 0
        if ((self.pointloc == 'axis' or self.pointloc == 'interior') and self.pno <= noc):
            self.theta = self.p_nu = (self.pno-1)*dtheta
            #rem = thetamax%(noc-1) and dtheta = (thetamax - rem)/(noc-1)
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

# Alternate method
        if self.pointloc == 'wall' and self.pno == noc+1:
            self.y = np.sqrt(self.aR * (th_rad ** 2))  # - infPtLocY +th_rad
            self.x = (self.y-infPtLocY)/np.tan(np.radians((thetamax+self.theta)*0.5))+infPtLocX

        elif self.pointloc == 'wall':
            self.y = np.sqrt(self.aR*(th_rad**2))
            self.x = (self.y-char1.y)/np.tan(np.radians((char1.theta+self.theta)*0.5)) + char1.x

# Characteristic drawing
        if self.pointloc == 'axis':
            self.x = (th_rad+downWallRad)*np.tan(np.radians(self.pzo*dtheta))

        elif self.pointloc == 'interior' and self.pno <= noc:
            x2 = char2.x
            y2 = char2.y
            m1 = np.tan(np.radians((-(90-self.pno*dtheta)+(self.theta-self.mu))/2))
            m2 = np.tan(np.radians((char2.theta+char2.mu+self.theta+self.mu)/2))
            denom = 1 - m1/m2
            self.x = (x2 + (1/m2)*(th_rad+downWallRad - y2))/denom
            self.y = (th_rad+downWallRad + m1*(x2 - y2/m2))/denom
            # print(self.x)
            # print(self.y)

        elif self.pointloc == 'interior' and self.pno > noc:
            x1 = char1.x
            y1 = char1.y
            x2 = char2.x
            y2 = char2.y
            m1 = np.tan(np.radians((char1.theta+self.theta)*0.5-(char1.mu+self.mu)*0.5))
            m2 = np.tan(np.radians((char2.theta+self.theta)*0.5+(char2.mu+self.mu)*0.5))
            denom = 1 - m1/m2
            self.x = (x2 + (1/m2)*(y1 - y2 - x1*m1))/denom
            self.y = (y1 + m1*(x2 - x1 - y2/m2))/denom
            # print(self.x)
            # print(self.y)


flag = 0
OG_ExMa = ExMa
truncFlag = 0
count = 0
limx = 0
limy = 0
trunc_pno = 0
pts = []
while flag == 0:
    #print("try ", count)
    thetamax = Nu(ExMa, rosh) / 2
    dtheta = round((thetamax)/(noc-1),4)
    pts=[]
    tmp = noc+1
    infPtLocX = downWallRad * np.sin(np.radians(thetamax))
    infPtLocY = th_rad + downWallRad - downWallRad * np.cos(np.radians(thetamax))
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
            # print("--third cond pno = {}  pno left = {} pno under = {} Kmin = {:>6.2f}  Kplus = {:>6.2f}    ".format(x+1,x+1-tmp,x,pts[x-tmp].Kmin, pts[x-tmp].Kplus))
    if ExMa <= 5*OG_ExMa:
        for x in range(nop):
            if pts[x].pointloc == 'wall':
                #print(abs(nozzExAng-pts[x].theta))
                # time.sleep(0.5)
                if (abs(nozzExAng-pts[x].theta) <= dtheta/2): # and (nozzExAng-pts[x].theta >= 0):
                    #print("I'm in")
                    if (pts[x].Ma-OG_ExMa <= 0.05) and (pts[x].Ma-OG_ExMa >= 0):
                        flag = 1
                        truncFlag = 1
                        limx = pts[x].x
                        limy = pts[x].y
                        trunc_pno = pts[x].pno
                        #print(ExMa)
                        break
                    else:
                        ExMa += 0.01
                        break
    else:
        print("Unable to generate TICN")
        flag = 1
        truncFlag = 0
    # print("ExMa = ", round(ExMa,5), x, pts[nop-1].x)
    # print("Tried ExMa = ", ExMa, " at ", pts[x].x, "pt Mach = ", pts[x].Ma," diff = ", abs(pts[x].Ma-OG_ExMa))
    count += 1
    #print(count)

if truncFlag == 0:
    ExMa = OG_ExMa
    thetamax = Nu(ExMa, rosh) / 2
    rem = thetamax % (noc - 1)
    dtheta = round(thetamax/(noc-1),3)
    pts = []
    tmp = noc + 1
    infPtLocX = downWallRad * np.sin(np.radians(thetamax))
    infPtLocY = th_rad + downWallRad - downWallRad * np.cos(np.radians(thetamax))
    for x in range(nop):
        #print(x, tmp, x - tmp, x - 1)
        if ptloc_array[x] == 'axis' and x>noc and x>0:
            tmp -= 1
            #print("------------tmp reduction, tmp = ",tmp)

        if x == 0:
            pts.append(Points(x, x, x))
            #print("--first cond pno={}".format(x+1))
        elif x != 0 and x < noc+1:
            pts.append(Points(x, x, pts[x-1]))
            #print("--second cond pno = {}    pno under = {}".format(x+1,x))
        elif x >= noc+1:
            pts.append(Points(x, pts[x-tmp], pts[x-1]))
            #print("--third cond pno = {}  pno left = {} pno under = {} Kmin = {:>6.2f}  Kplus = {:>6.2f}    ".format(x+1,x+1-tmp,x,pts[x-tmp].Kmin, pts[x-tmp].Kplus))

print(thetamax)
print(ExMa)

# if truncFlag == 1:
#     purgeVal = nop-trunc_pno
#     i = 0
#     while i<purgeVal:
#         pts.pop(nop-1-i)
#         i+=1
#     nop = len(pts)

#print("Everyloop pno = {}   pno left = {}  pno under = {}".format(x+1,x+1-tmp,x))

def displayTable(): #Literally does what it looks like it's gonna do
    print("-----------|------------|------------|-----------|-----------|----------|------------|----------")
    print("  Point no.|       Kmin |      Kplus |  θ(Theta) |     ν(Nu) |        M |      μ(mu) |      pzo ")
    print("-----------|------------|------------|-----------|-----------|----------|------------|----------")
    for x in range(nop):
        print(" {0:>9} |   {1:>8.3f} |   {2:>8.3f} |  {3:>8.3f} |  {4:>8.3f} | {5:>8.3f} |   {6:>8.3f} |   {7:>6}".format(pts[x].pno,pts[x].Kmin,pts[x].Kplus,pts[x].theta,pts[x].p_nu,pts[x].Ma,pts[x].mu,pts[x].pzo))
    return 0

def drawChar():
    colour = 'red'
    for h in range(nop):
        if pts[h].pointloc == 'axis':
            mat.plot([0, pts[h].x], [th_rad+downWallRad, pts[h].y], colour)

        elif pts[h].pointloc == 'interior' or pts[h].pointloc == 'wall':
            mat.plot([pts[h-1].x, pts[h].x], [pts[h-1].y, pts[h].y], colour)

def altDrawChar():
    colour = 'red'
    back = noc
    for h in range(nop):
        if pts[h].pointloc == 'axis':
            mat.plot([0, pts[h].x], [th_rad+downWallRad, pts[h].y], colour)

        elif pts[h].pointloc == 'wall':
            mat.plot([pts[h-back].x, pts[h].x], [pts[h-back].y, pts[h].y], colour)
            back -= 1


def TICNang(): #This function draws the ideal contour nozzle
    conx = []
    cony = []
    xplot = np.arange(0,infPtLocX, downWallRad/100)
    if truncFlag == 1:
        for a in range(len(xplot)):
            conx.append(xplot[a])
            cony.append(round(-np.sqrt(downWallRad**2 - xplot[a]**2) + downWallRad + th_rad, 4))

        i = 0
        while i < nop:
            if pts[i].pno == pts[nop-1].pno:
                conx.append(limx)
                cony.append(limy)

            elif pts[i].pointloc == 'wall':
                conx.append(pts[i].x)
                cony.append(pts[i].y)
            i += 1
    elif truncFlag == 0:
        for a in range(len(xplot)):
            conx.append(xplot[a])
            cony.append(round(-np.sqrt((downWallRad ** 2) - xplot[a] ** 2) + downWallRad + th_rad,4))

        i = 0
        while i < nop:
            if pts[i].pointloc == 'wall':
                conx.append(pts[i].x)
                cony.append(pts[i].y)
            i+=1

    mat.plot(conx, cony, color='black')
    # print(conx)
    # print(cony)
    return 0

def dispLocation(): #Just for Debugging
    print("-----------|------------|------------|------------|--")
    print("  Point no.|          X |          Y |      Theta | ")
    print("-----------|------------|------------|------------|--")
    for x in range(nop):
        #if pts[x].pointloc != 'wall':
           # print(" {0:>9} |   {1:>8.3f} |   {2:>8.3f}  |".format(pts[x].pno,pts[x].x,pts[x].y))
        if pts[x].pointloc == 'wall':
            print(" {0:>9} |   {1:>8.3f} |   {2:>8.3f} |   {3:>8.3f} | Wallpoint".format(pts[x].pno, pts[x].x, pts[x].y, pts[x].theta))
    return 0

def curveGen(): #Attempts to generate smooth regression curve of nozzle points
    xx = []
    yy = []
    def parabolic(x, a, b, c):
        return (a+(b*x)+(c*(x**2)))

    def cubic(x, a, b, c, d):
        return (a+(b*x)+(c*(x**2))+(d*(x**3)))

    def quadratic(x, a, b, c, d, e):
        return (a+(b*x)+(c*(x**2))+(d*(x**3))+(e*(x**4)))

    i = 0
    while i < nop:
        if len(xx) == 0:
            xx.append(round(infPtLocX,4))
            yy.append(round(infPtLocY,4))

        elif pts[i].pno == pts[nop-1].pno:
            xx.append(round(nozzLength,4))
            yy.append(round(limy,4))

        elif pts[i].pointloc == 'wall':
            xx.append(round(pts[i].x,4))
            yy.append(round(pts[i].y,4))

        i += 1

    # print(x)
    # print(y)
    print(len(xx))
    print(len(yy))

    paraDAT,pcov = curve_fit(parabolic,xx,yy)
    cubicDAT,pcov = curve_fit(cubic,xx,yy)
    quadDAT,pcov = curve_fit(quadratic,xx,yy)

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
    y2 = parabolic(xx, a2, b2, c2)
    y3 = cubic(xx, a3, b3, c3, d3)
    y4 = quadratic(xx, a4, b4, c4, d4, e4)
    mat.plot(xx, y2, 'orange')
    mat.plot(xx, y3, 'blue')
    mat.plot(xx, y4, 'red')

    # print(a2,b2,c2)
    # print(a3,b3,c3,d3)
    # print(a4,b4,c4,d4,e4)
    # inExpReg_X = np.arange(0, infPtLocX, 0.01)
    # inExpReg_Y = - np.sqrt((0.5*th_rad)**2 - inExpReg_X**2)+1.5*th_rad
    # mat.plot(inExpReg_X,inExpReg_Y,'orange')
    # mat.plot(x,y2,'red')
    # mat.plot(x,y3,'blue')
    # mat.plot(x,y4,'green')

def XLimporter():
    dir = "D:\Study Files\Aerospace Engg\Project\Output\\"
    #dir = input("Enter save directory: ")
    wb = xlwt.Workbook()
    xplot = np.arange(0, infPtLocX, downWallRad / 100)
    ws = wb.add_sheet("ICN coordinates")
    ws.write(0, 0, "Point Number")
    ws.write(0, 1, "X")
    ws.write(0, 2, "Y")
    j = 0
    while j<len(xplot):
        ws.write(j+1, 0, j+1)
        ws.write(j+1, 1, xplot[j])
        ws.write(j+1, 2, round(-np.sqrt((downWallRad) ** 2 - xplot[a] ** 2) + downWallRad + th_rad))
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

    wb.save("{0:}Ideal Con Nozzle, ExMa={1:}, Throat rad={2:}, No. of char={3:}.xls".format(dir,OG_ExMa,th_rad,noc))
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

print('dtheta = ', dtheta)

#drawChar()
altDrawChar()
TICNang()        #This function draws the minimum length nozzle
#altDrawChar()       #draws characterisitic lines
displayTable()   #This table shows flow properties along the points in the flow
dispLocation()   #Shows coordinates of all wallpoints as well as characteristic intersection points in Cartesian coords sys
#everyLoc()
#curveGen()       #Attempts to generate smooth regression curve of nozzle points but fails miserably
#XLimporter()     #Imports coords to an excel file

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
mat.xlabel("X-axis (m)")
mat.ylabel("Y-axis (m)")
mat.title("Truncated Ideal Contour Nozzle (angle limit)")
if truncFlag == 0:
    mat.ylim(0,(pts[nop-1].y+0.1*pts[nop-1].y))
    mat.xlim(0,(pts[nop-1].x+0.1*pts[nop-1].x))
elif truncFlag == 1:
    mat.ylim(0, (limy + 0.1*limy))
    mat.xlim(0, (limx + 0.1*limx))
mat.grid()
mat.show()