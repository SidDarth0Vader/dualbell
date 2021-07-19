import numpy as np
import matplotlib.pyplot as mat
import xlwt

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
#ExMa = float(input("Enter design Exit Mach number: ")) #ExMa = Exit Mach Number
th_rad = 1
#th_rad = float(input("Enter radius of cross-section of throat(in metres): ")) #th_rad = Radius at throat(m)
rosh = 1.4
#rosh = float(input("Enter ratio of specific heats: ")) #rosh = Ratio of Specific Heats
noc = 7
#noc = int(input("Enter number of characteristic lines: ")) #noc = number of characteristic lines
mach = np.arange(0,ExMa , 0.1)
thetamax = Nu(ExMa,rosh)/2
rem = thetamax%(noc-1)
dtheta = (thetamax - rem)/(noc-1)
i = noc + 1
nop = 0

while i>=2:
    nop = nop + i
    i -= 1

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
            self.theta = self.p_nu = rem+(self.pno-1)*dtheta
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

# Alternate method wall point plot
        if self.pointloc == 'wall' and self.pno == noc+1:
            self.y = np.sqrt(self.aR * (th_rad ** 2))  # - infPtLocY +th_rad
            self.x = (self.y-th_rad)/np.tan(np.radians(thetamax))

        elif self.pointloc == 'wall':
            self.y = np.sqrt(self.aR*(th_rad**2))
            self.x = (self.y-char1.y)/np.tan(np.radians(char1.theta)) + char1.x

# Characteristic drawing
        if self.pointloc == 'axis':
            self.x = (th_rad)*np.tan(np.radians(self.pzo*dtheta))

        elif self.pointloc == 'interior' and self.pno <= noc:
            x2 = char2.x
            y2 = char2.y
            m1 = np.tan(np.radians(-90 + self.pno*dtheta))
            m2 = np.tan(np.radians((char2.theta+char2.mu+self.theta+self.mu)/2))
            denom = 1 - m1/m2
            self.x = (x2 + (1/m2)*(th_rad - y2))/denom
            self.y = (th_rad + m1*(x2 - y2/m2))/denom
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

def displayTable(): #Displays flow properties
    print("-----------|------------|------------|-----------|-----------|----------|------------")
    print("  Point no.|       Kmin |      Kplus |  θ(Theta) |     ν(Nu) |        M |      μ(mu) ")
    print("-----------|------------|------------|-----------|-----------|----------|------------")
    for x in range(nop):
        print(" {0:>9} |   {1:>8.3f} |   {2:>8.3f} |  {3:>8.3f} |  {4:>8.3f} | {5:>8.3f} |  {6:>8.3f}".format(pts[x].pno,pts[x].Kmin,pts[x].Kplus,pts[x].theta,pts[x].p_nu,pts[x].Ma,pts[x].mu))
    return 0

def minLength(): #Draws the shape of the nozzle
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
        mat.plot(conx,cony, color = 'Red')
        i+=1
    # print(conx)
    # print(cony)
    return 0

def dispLocation(): #Displays location wall points
    print("-----------|------------|------------|------------|--")
    print("  Point no.|          X |          Y |      Theta | ")
    print("-----------|------------|------------|------------|--")
    for x in range(nop):
        #if pts[x].pointloc != 'wall':
           # print(" {0:>9} |   {1:>8.3f} |   {2:>8.3f}  |".format(pts[x].pno,pts[x].x,pts[x].y))
        if pts[x].pointloc == 'wall':
            print(" {0:>9} |   {1:>8.3f} |   {2:>8.3f}  |   {3:>8.3f} | Wallpoint".format(pts[x].pno, pts[x].x, pts[x].y,pts[x].theta))
    return 0

def XLimporter():
    dir = "D:\Study Files\Aerospace Engg\Project\Output\\"
    #dir = input("Enter save directory: ")
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

    for x in range(nop):
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
    i = 1

    for y in range(nop):
        if pts[y].pointloc == 'wall':
            coords.write(i, 0, pts[y].pno)
            coords.write(i, 1, round(pts[y].x,4))
            coords.write(i, 2, round(pts[y].y,4))
            i+=1

    wb.save("{0:}Minimum length nozzle, ExMa={1:}, Throat rad={2:}, No. of char={3:}.xls".format(dir,ExMa,th_rad,noc))
    return 0

# actAR = (np.pi*pts[nop-1].y**2)/(np.pi*(th_rad*2)**2)
# print(np.pi, pts[nop-1].y)
# print("Actual Area Ratio = {:.4f}".format(actAR))
# print("Target Area Ratio = {:.4f}".format(areaRatio(ExMa,rosh)))

#Enable or disable the functions down below by removing the '#' or putting the '#' in front of the functions respectively

minLength()      #This function shows the shape of the nozzle
displayTable()   #This table shows flow properties along the points in the flow
#dispLocation()   #Shows coordinates of all wallpoints as well as characteristic intersection points in Cartesian coords sys
XLimporter()     #Imports data to an excel file


mat.axis('scaled')
#mat.plot([0,pts[nop-1].x],[pts[nop-1].y,pts[nop-1].y], label = "{}m".format(pts[nop-1].x), color = 'Blue', linewidth = 0.7)
#mat.plot([pts[nop-1].x,pts[nop-1].x],[0,pts[nop-1].y], label = "{}m".format(pts[nop-1].y), color = 'Blue', linewidth = 0.7)
mat.xlabel("X-axis (m)")
mat.ylabel("Y-axis (m)")
mat.title("Minimum length Nozzle")
mat.ylim(0,(pts[nop-1].y+0.1*pts[nop-1].y))
mat.xlim(0,(pts[nop-1].x+0.1*pts[nop-1].x))
mat.grid()
mat.show()