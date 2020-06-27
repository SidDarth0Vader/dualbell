import numpy as np
import matplotlib.pyplot as mat
import xlwt
import time


def M(nu):  # This method of finding inverse Prandtl-Meyer Function was found at http://www.pdas.com/pm.pdf
    nu_inf = (np.pi * (np.sqrt(6) - 1)) / 2  # the maximum turning angle
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
    return ((1/M)*(((gam+1)/2)**(-((gam+1)/(2*gam-2))))*((1+(((gam-1)/2)*M**2))**((gam+1)/(2*gam-2))))


ExMa = 4.3
#ExMa = float(input("Enter Nozzle Exit Mach number: ")) #ExMa = Exit Mach Number
th_rad = 0.03354
#th_rad = float(input("Enter radius of cross-section of throat(in metres): ")) #th_rad = Radius at throat(m)
nozzExAng = 10
#nozzExAng = float(input("Enter nozzle exit angle(in degrees): ")) #nozzExAng = Nozzle Exit Angle(degrees)
throatAng = 23.9
rosh = 1.4
#rosh = float(input("Enter ratio of specific heats: ")) #rosh = Ratio of Specific Heats
noc = 30
#noc = int(input("Enter number of characteristic lines: ")) #noc = number of characteristic lines
downWallRad = 0.5*th_rad
thetamax = Nu(ExMa,rosh)/2
dtheta = round((thetamax)/(noc),3)
altdtheta = round(throatAng/noc, 3)
i = noc + 1
nop = 0

while i>=2:
    nop = nop + i
    i = i - 1

ptloc_array = []
ptzo_array = []
ptcount = []
tmp = noc+1
ptzo = 1
while len(ptloc_array)<nop:
    i=1
    while i<=tmp and tmp>1:
        if i==tmp:
            ptloc_array.append('wall')
            ptzo_array.append(ptzo)
            ptcount.append(i)
            ptzo += 1
        elif i==1:
            ptloc_array.append('axis')
            ptzo_array.append(ptzo)
            ptcount.append(i)
        else:
            ptloc_array.append('interior')
            ptzo_array.append(ptzo)
            ptcount.append(i)
        i+=1

    tmp-=1

# infPtLocX = 0.5 * th_rad * np.sin(np.radians(thetamax+Nu(ExMa,rosh)))
# infPtLocY = 1.5 * th_rad - 0.5 * th_rad * np.cos(np.radians(thetamax+Nu(ExMa,rosh)))
infPtLocX = downWallRad * np.sin(np.radians(throatAng))
infPtLocY = th_rad+downWallRad - downWallRad * np.cos(np.radians(throatAng))
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

# # Alternate method wall point plot
#         if self.pointloc == 'wall' and self.pno == noc+1:
#             self.y = np.sqrt(self.aR * (th_rad ** 2))  # - infPtLocY +th_rad
#             self.x = (self.y-infPtLocY)/np.tan(np.radians((thetamax+self.theta)*0.5))+infPtLocX
#
#         elif self.pointloc == 'wall':
#             self.y = np.sqrt(self.aR*(th_rad**2))
#             self.x = (self.y-char1.y)/np.tan(np.radians((char1.theta+self.theta)*0.5)) + char1.x

# Alternate method wall point plot
        if self.pointloc == 'wall' and self.pno == noc+1:
            x1 = infPtLocX
            y1 = infPtLocY
            x2 = char2.x
            y2 = char2.y
            m1 = np.tan(np.radians(thetamax))
            m2 = np.tan(np.radians((char2.theta + self.theta)*0.5 + (char2.mu + self.mu)*0.5))
            denom = 1 - (m1/m2)
            self.x = (x2 + (1/m2)*(y1 - y2 - x1*m1))/denom
            self.y = (y1 + m1*(x2 - x1 - (y2/m2)))/denom

        elif self.pointloc == 'wall':
            x1 = char1.x
            y1 = char1.y
            x2 = char2.x
            y2 = char2.y
            m1 = np.tan(np.radians((char1.theta+self.theta-self.mu/1.05)/2))
            m2 = np.tan(np.radians((char2.theta+self.theta)*0.5+(char2.mu+self.mu)*0.5))
            denom = 1 - (m1/m2)
            self.x = (x2 + (1/m2)*(y1 - y2 - x1*m1))/denom
            #self.y = np.sqrt(self.aR * (th_rad ** 2))
            self.y = (y1 + m1*(x2 - x1 - y2/m2))/denom

# Characteristic drawing
        if self.pointloc == 'axis':
            self.x = (th_rad+downWallRad)*np.tan(np.radians(self.pzo*altdtheta))

        elif self.pointloc == 'interior' and self.pno <= noc:
            x2 = char2.x
            y2 = char2.y
            m1 = np.tan(np.radians(-90 + self.pno*altdtheta))
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
    dtheta = round(thetamax/noc,4)
    altdtheta = round(throatAng / noc, 3)
    pts=[]
    tmp = noc+1
    #infPtLocX = downWallRad * np.sin(np.radians(thetamax))
    #infPtLocY = th_rad + downWallRad - downWallRad * np.cos(np.radians(thetamax))
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
                if (abs(nozzExAng-pts[x].theta) <= dtheta/2): # and (pts[x].theta-nozzExAng >= 0):
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
    dtheta = round(thetamax/(noc),3)
    altdtheta = round(throatAng / noc, 3)
    pts = []
    tmp = noc + 1
    #infPtLocX = downWallRad * np.sin(np.radians(thetamax))
    #infPtLocY = th_rad + downWallRad - downWallRad * np.cos(np.radians(thetamax))
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

trunc_nop = nop
if truncFlag == 1:
    purgeVal = nop-trunc_pno
    i = 0
    while i<purgeVal:
        ptcount.pop(nop-1-i)
        i+=1
    trunc_nop = len(ptcount)

#print("Everyloop pno = {}   pno left = {}  pno under = {}".format(x+1,x+1-tmp,x))


def displayTable(): #Literally does what it looks like it's gonna do
    print("-----------|------------|------------|-----------|-----------|----------|------------")
    print("  Point no.|       Kmin |      Kplus |  θ(Theta) |     ν(Nu) |        M |      μ(mu) ")
    print("-----------|------------|------------|-----------|-----------|----------|------------")
    for x in range(trunc_nop):
        print(" {0:>9} |   {1:>8.3f} |   {2:>8.3f} |  {3:>8.3f} |  {4:>8.3f} | {5:>8.3f} |   {6:>8.3f}".format(pts[x].pno,pts[x].Kmin,pts[x].Kplus,pts[x].theta,pts[x].p_nu,pts[x].Ma,pts[x].mu))
    return 0


def displayFullTable(): #Literally does what it looks like it's gonna do
    print("-----------|------------|------------|-----------|-----------|----------|------------")
    print("  Point no.|       Kmin |      Kplus |  θ(Theta) |     ν(Nu) |        M |      μ(mu) ")
    print("-----------|------------|------------|-----------|-----------|----------|------------")
    for x in range(nop):
        print(" {0:>9} |   {1:>8.3f} |   {2:>8.3f} |  {3:>8.3f} |  {4:>8.3f} | {5:>8.3f} |   {6:>8.3f}".format(pts[x].pno,pts[x].Kmin,pts[x].Kplus,pts[x].theta,pts[x].p_nu,pts[x].Ma,pts[x].mu))
    return 0


def drawChar():
    colour = 'red'
    for h in range(nop):
        if pts[h].pointloc == 'axis':
            mat.plot([0, pts[h].x], [th_rad+downWallRad, pts[h].y], colour, linewidth = 0.5)
            mat.plot([0, pts[h].x], [th_rad+downWallRad, pts[h].y], 'o', linewidth = 0.5)

        elif pts[h].pointloc == 'interior' or pts[h].pointloc == 'wall':
            mat.plot([pts[h-1].x, pts[h].x], [pts[h-1].y, pts[h].y], colour, linewidth = 0.5)
            mat.plot([pts[h-1].x, pts[h].x], [pts[h-1].y, pts[h].y], 'o', linewidth = 0.5)


conx = []
cony = []
conPtType = []
fullConx = []
fullCony = []


def TICNang(): #This function draws the ideal contour nozzle
    infconX = []
    infconY = []
    stConX = []
    stConY = []
    xplot = np.arange(0,infPtLocX+(downWallRad/100), downWallRad/100)
    if truncFlag == 1:
        for a in range(len(xplot)):
            conx.append(xplot[a])
            cony.append(-np.sqrt(downWallRad**2 - xplot[a]**2) + downWallRad + th_rad)
            infconX.append(xplot[a])
            infconY.append(-np.sqrt(downWallRad**2 - xplot[a]**2) + downWallRad + th_rad)
            conPtType.append('Pre-det')

        i = 0
        while i < trunc_nop:
            if pts[i].pno == pts[trunc_nop-1].pno:
                conx.append(limx)
                cony.append(limy)
                stConX.append(limx)
                stConY.append(limy)
                conPtType.append('Calc')

            elif pts[i].pointloc == 'wall':
                conx.append(pts[i].x)
                cony.append(pts[i].y)
                stConX.append(pts[i].x)
                stConY.append(pts[i].y)
                conPtType.append('Calc')
            i += 1

        i = 0
        while i < nop:
            if pts[i].pointloc == 'wall':
                fullConx.append(pts[i].x)
                fullCony.append(pts[i].y)
            i+=1

    elif truncFlag == 0:
        for a in range(len(xplot)):
            conx.append(xplot[a])
            cony.append(-np.sqrt((downWallRad ** 2) - xplot[a] ** 2) + downWallRad + th_rad)
            infconX.append(xplot[a])
            infconY.append(-np.sqrt(downWallRad**2 - xplot[a]**2) + downWallRad + th_rad)
            conPtType.append('Pre-det')

        i = 0
        while i < nop:
            if pts[i].pointloc == 'wall':
                conx.append(pts[i].x)
                cony.append(pts[i].y)
                stConX.append(pts[i].x)
                stConY.append(pts[i].y)
                conPtType.append('Calc')
            i+=1

    mat.plot(fullConx, fullCony, ':', color='blue')
    # mat.plot(infconX, infconY, color='black')
    # mat.plot(stConX, stConY, color='black')
    # print(conx)
    # print(cony)
    return 0


def dispLocation(): #Just for Debugging
    print("-----------|------------|------------|------------|--")
    print("  Point no.|          X |          Y |      Theta | ")
    print("-----------|------------|------------|------------|--")

    for a in range(len(conx)):
        if conPtType[a] == 'Pre-det':
            print(" {0:>9} |   {1:>8.6f} |   {2:>8.6f} |            | {3:>7}".format(a, conx[a], cony[a], conPtType[a]))
        elif conPtType[a] == 'curve fit':
            print(" {0:>9} |   {1:>8.6f} |   {2:>8.6f} |            | {3:>8}".format(a, conx[a], cony[a], conPtType[a]))

    for x in range(trunc_nop):
        if pts[x].pointloc == 'wall':
            print(" {0:>9} |   {1:>8.3f} |   {2:>8.3f} |   {3:>8.3f} | Calc".format(pts[x].pno, pts[x].x, pts[x].y, pts[x].theta))
    return 0


def curveGen():
    xvar = []
    yvar = []
    cutoff = 0

    for xd in range(len(conx)):
        if conPtType[xd] == 'Calc':
            cutoff = xd
            xvar.append(conx[0])
            xvar.append(conx[xd-1])
            xvar.append(conx[xd])
            xvar.append(conx[xd+1])
            #xvar.append(conx[xd+2])
            #xvar.append(conx[xd+3])
            yvar.append(cony[0])
            yvar.append(cony[xd-1])
            yvar.append(cony[xd])
            yvar.append(cony[xd+1])
            #yvar.append(cony[xd+2])
            #yvar.append(cony[xd + 3])
            break

    #poly_coeff = np.polyfit(conx, cony, 3)
    poly_coeff = np.polyfit(xvar, yvar, 3)
    print(xvar, yvar)
    print(poly_coeff)
    #print(a, b, c, d)
    print(poly_coeff[0])

    xplot = np.arange(conx[cutoff-1], conx[cutoff], 0.0005)
    yplot3 = poly_coeff[3] + poly_coeff[2]*xplot + poly_coeff[1]*xplot**2 + poly_coeff[0]*xplot**3
    for lol in range(len(conx)):
        if conPtType[lol] == 'Calc':
            val = len(xplot)-1
            while val >= 0:
                conx.insert(lol, xplot[val])
                cony.insert(lol, poly_coeff[3] + poly_coeff[2]*xplot[val] + poly_coeff[1]*xplot[val]**2 + poly_coeff[0]*xplot[val]**3)
                conPtType.insert(lol, 'curve fit')
                val-=1
            break

    #mat.plot(xplot, yplot3, 'blue')
    mat.plot(conx, cony, 'black')


def XLimporter():
    dir = "D:\Study Files\Aerospace Engg\Project\Output\\"
    wb = xlwt.Workbook()
    props = wb.add_sheet("Flow properties")
    coords = wb.add_sheet("Minimum coordinates")
    props.write(0, 0, "Point no.")
    props.write(0, 1, "Kmin")
    props.write(0, 2, "Kplus")
    props.write(0, 3, "θ(Theta)")
    props.write(0, 4, "ν(Nu)")
    props.write(0, 5, "M")
    props.write(0, 6, "μ(mu)")
    xplot = np.arange(0, infPtLocX, downWallRad/100)

    for x in range(trunc_nop):
        props.write(x+1, 0, pts[x].pno)
        props.write(x+1, 1, round(pts[x].Kmin, 4))
        props.write(x+1, 2, round(pts[x].Kplus, 4))
        props.write(x+1, 3, round(pts[x].theta, 4))
        props.write(x+1, 4, round(pts[x].p_nu, 4))
        props.write(x+1, 5, round(pts[x].Ma, 4))
        props.write(x+1, 6, round(pts[x].mu, 4))

    coords.write(0, 0, "Point Number")
    coords.write(0, 1, "X")
    coords.write(0, 2, "Y")
    for val in range(len(conx)):
        coords.write(val+1, 0, val+1)
        coords.write(val+1, 1, round(conx[val],6))
        coords.write(val+1, 2, round(cony[val],6))

    wb.save("{0:}TICN(ang) modded, ExMa={1:}, Throat rad={2:}, No. of char={3:}.xls".format(dir, OG_ExMa, th_rad, noc))
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

# Enable or disable the functions down below by removing the '#' or putting the '#' in
# front of the functions respectively

print('dtheta = ', dtheta)
print('altdtheta = ', altdtheta)

drawChar()       # This function draws characteristic lines
TICNang()        # This function draws the minimum length nozzle
displayTable()   # This table shows flow properties along the points in the flow
#dispLocation()   # Shows coordinates of all wallpoints as well as characteristic intersection points in Cartesian coords sys
# everyLoc()
curveGen()       # Attempts to generate smooth regression curve of nozzle points but fails miserably
#XLimporter()     # Imports coords to an excel file
dispLocation()




mat.axis('scaled')
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

# while 1:
#     importquery = input("Import data to an excel file?(Y/N):")
#     if importquery == 'Y' or importquery == 'y':
#         XLimporter()
#     elif importquery == 'N' or importquery == 'n':
#         break
#     else:
#         print("Invalid option, please try again")