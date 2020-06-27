import numpy as np
import matplotlib.pyplot as mat
import xlwt
#import time


def M(nu):  # This method of finding inverse Prandtl-Meyer Function was found at http://www.pdas.com/pm.pdf
    nu_inf = (np.pi * (np.sqrt(6) - 1)) / 2  # the maximum turning angle
    y = (np.radians(nu) / nu_inf) ** (2 / 3)
    A = 1.3604
    B = 0.0962
    C = -0.5127
    D = -0.6722
    E = -0.3278
    # M = ((1 + A * y + B * (y ** 2) + C * (y ** 3)) / (1 + D * y + E * (y ** 2)))
    return (1 + A * y + B * (y ** 2) + C * (y ** 3)) / (1 + D * y + E * (y ** 2))


def Nu(m, rosh):
    #rosh = ratio of specific heats
    #vM = np.sqrt((rosh+1)/(rosh-1))*np.arctan(np.sqrt(((rosh-1)*((m**2)-1))/(rosh+1)))-np.arctan(np.sqrt((m**2)-1))
    return np.degrees(np.sqrt((rosh+1)/(rosh-1))*np.arctan(np.sqrt(((rosh-1)*((m**2)-1))/(rosh+1)))-np.arctan(np.sqrt((m**2)-1)))


def areaRatio(M, gam):
    return (1/M)*(((gam+1)/2)**(-((gam+1)/(2*gam-2))))*((1+(((gam-1)/2)*M**2))**((gam+1)/(2*gam-2)))


def tempRatio(M, gam):
    return 1/(1+((rosh-1)/2)*M**2)


ExMa = 4.3
# ExMa = float(input("Enter Nozzle Exit Mach number: ")) #ExMa = Exit Mach Number
th_rad = 0.03354
# th_rad = float(input("Enter radius of cross-section of throat(in metres): ")) #th_rad = Radius at throat(m)
nozzExAng = 10
# nozzExAng = float(input("Enter nozzle exit angle(in degrees): ")) #nozzExAng = Nozzle Exit Angle(degrees)
throatAng = 23.9
rosh = 1.4
gas_T = 400
# gas_T = float(input("Enter temperature of gas: "))
R = 287
ivp_count = 11
# ivp_count = int(input("Enter number of initial value line points: ")) #noc = number of characteristic lines
nop_ivprob = ivp_count*(ivp_count-1)
if ivp_count == 2:
    nop_ivprob = 3

delta = 1
upWallRad = 5.5*th_rad
downWallRad = 0.5*th_rad
a_star = np.sqrt((2*rosh*R*gas_T)/(rosh+1))
a0sq = rosh*R*gas_T
a_speed = np.sqrt(rosh*R*gas_T)
alpha = np.sqrt((1+delta)/((rosh+1)*upWallRad*th_rad))
epsilon = -((th_rad/(2*(3+delta)))*np.sqrt(((rosh+1)*(1+delta))/(upWallRad/th_rad)))
thetamax = Nu(ExMa,rosh)/2
dy = th_rad/(ivp_count-1)


def xsonic(yy):
    return -((rosh+1)*alpha*(yy**2))/(2*(1+delta))


def xcorr(yy):
    return -((rosh+1)*alpha*(yy**2))/(2*(3+delta))


def ubar(xx, yy):
    return a_star*(1+alpha*xx+(((rosh+1)*(alpha**2)*(yy**2))/(2*(1+delta))))


def aspd(u):
    return np.sqrt(a0sq - ((rosh-1)/2)*(u**2))


def invWall(x3, y3, u3, v3, x1, y1, u1, v1, x4, y4, theta4):
    u2 = u3
    v2 = v3
    theta2 = V2 = a2 = mu2 = x2 = y2 = u2temp = v2temp = u4 = v4 = u4temp = flag = count = 0
    # print('x3 = {0:>5.4f}, y3 = {1:>5.4f}, u3 = {2:>5.4f}, v3 = {3:>5.4f}'.format(x3, y3, u3, v3))
    # print('x1 = {0:>5.4f}, y1 = {1:>5.4f}, u1 = {2:>5.4f}, v1 = {3:>5.4f}'.format(x1, y1, u1, v1))
    # print('x4 = {0:>5.4f}, y4= {1:>5.4f}, theta4 = {2:>5.4f}\n'.format(x4, y4, theta4))
    while flag == 0 and count <= 50:
        #print("Loop 1")
        u2temp = u2
        theta2 = np.degrees(np.arctan(v2/u2))
        Vel2 = np.sqrt(u2**2+v2**2)
        a2 = np.sqrt(rosh*R*gas_T - ((rosh-1)/2)*Vel2**2)
        #print('M = ', Vel2/a2)
        if (Vel2/a2) < 1:
            mu2 = 90
        else:
            mu2 = np.degrees(np.arcsin(a2/Vel2))
        #print('mu2 = ', mu2)
        lamPlus = np.tan(np.radians(theta2+mu2))
        s = ((y3-y1)/(x3-x1))
        A = [[-((y3-y1)/(x3-x1)), 1], [-lamPlus, 1]]
        B = [(y1-s*x1),(y4-lamPlus*x4)]
        x2 = np.linalg.solve(A, B)[0]
        y2 = np.linalg.solve(A, B)[1]
        #print("u2 = {0:6.3f}, v2 = {1:6.4f}".format(u2, v2))
        u2 = u1 + ((x2-x1)/(x3-x1))*(u3-u1)
        v2 = v1 + ((x2-x1)/(x3-x1))*(v3-v1)
        count += 1
        #print(count,'\n')
        if abs(u2-u2temp) <= 0.001:
            flag = 1

    flag = count = 0
    uplus = u2
    vplus = v2
    yplus = (y2+y4)/2
    while flag == 0 and count <= 50:
        #print("Loop 2\n", count)
        u4temp = u4
        theta2 = np.degrees(np.arctan(vplus/uplus))
        Vel2 = np.sqrt(uplus**2+vplus**2)
        a2 = np.sqrt(rosh*R*gas_T - ((rosh-1)/2)*(Vel2**2))
        if (Vel2/a2) < 1:
            muplus = 90
        else:
            muplus = np.degrees(np.arcsin(a2/Vel2))
        # print('M = ', Vel2/a2)
        # print('mu2 = ', np.degrees(np.arcsin(a2/Vel2)))
        # print('Vel2 = ', Vel2)
        # print('a2 = ', a2)
        lamPlus = np.tan(np.radians(theta2+muplus))
        Qplus = uplus**2 - a2**2
        Rplus = 2*uplus*vplus - (uplus**2 - a2**2)*lamPlus
        Splus = delta*(((a2**2)*vplus)/yplus)
        Tplus = Splus*(x4-x2)+Qplus*uplus+Rplus*vplus
        A = [[np.tan(np.radians(theta4)), -1], [Qplus, Rplus]]
        B = [0, Tplus]
        u4 = np.linalg.solve(A, B)[0]
        v4 = np.linalg.solve(A, B)[1]
        # print("u2 = {0:6.3f}, v2 = {1:6.4f}".format(u2, v2))
        # print("u4 = {0:>6.3f}, v4 = {1:>6.3f}, uplus = {2:>6.3f}, vplus = {3:>6.3f}".format(u4, v4, uplus, vplus))
        if abs(u4-u4temp) <= 0.001:
            flag = 1
        else:
            uplus = (u2 + u4)/2
            vplus = (v2 + v4)/2
        count += 1
    if flag != 0:
        return [u4, v4, y2]


def intersect_lineNcircle(x1,y1,m1):
    r = downWallRad
    b = th_rad+r
    intC = -m1*x1+y1
    A = 1 + m1**2
    B = 2*m1*intC - 2*m1*b
    C = b**2 + intC**2 - 2*b*intC - r**2
    xplus = (-B + np.sqrt(B**2 - 4*A*C))/(2*A)
    xminus = (-B - np.sqrt(B**2 - 4*A*C))/(2*A)
    yplus = m1*xplus + intC
    yminus = m1*xminus + intC
    return [xminus, yminus]

arcX = downWallRad*np.sin(np.radians(throatAng))
arcY = th_rad+downWallRad - downWallRad*np.cos(np.radians(throatAng))
arcEnd = 0

conx = []
cony = []
stConX = []
stConY = []
fullConx = []
fullCony = []
xplot = np.arange(0, arcX+(downWallRad/100), 0.0001)
for a in range(len(xplot)):
    conx.append(xplot[a])
    cony.append(-np.sqrt(downWallRad**2 - xplot[a]**2) + downWallRad + th_rad)


class InitValLine:
    def __init__(self, aloc_pno):
        self.pno = aloc_pno
        if self.pno == 0:
            self.pointloc = 'contour'
        elif self.pno == ivp_count-1:
            self.pointloc = 'axis'
        else:
            self.pointloc = 'interior'

        if self.pno == 0:
            self.y = th_rad
        elif self.pno != 0:
            self.y = th_rad - self.pno*dy
        self.xson = xsonic(self.y)
        self.v0_x = xcorr(self.y)
        self.x = self.v0_xbar = xcorr(self.y) - epsilon
        self.u = self.ucap = ubar(self.v0_x, self.y)
        self.v = 0
        self.V = self.u
        self.a_spd = aspd(self.ucap)
        self.Ma = self.ucap/self.a_spd
        self.nu = Nu(self.Ma, rosh)
        self.mu = np.degrees(np.arcsin(1/self.Ma))
        self.theta = 0
        self.Kplus = self.theta-self.nu
        self.Kmin = self.theta+self.nu
        self.corr = 1


i = 0
ivl = []
ivlX = []
ivlY = []
while i < ivp_count:
    ivl.append(InitValLine(i))
    ivlX.append(ivl[i].v0_xbar)
    ivlY.append(ivl[i].y)
    i += 1

ivpLocArray = []
ivpLA_lite = []
tmp = 1
val = 1
for cn in range(nop_ivprob):
    if cn < nop_ivprob/2:
        if tmp == 1:
            ivpLocArray.append('first char')
            ivpLA_lite.append('f')
            val += 1
            tmp = val
        else:
            ivpLocArray.append('interior')
            ivpLA_lite.append('i')
            tmp -= 1
    elif cn >= nop_ivprob/2:
        if cn == (nop_ivprob/2):
            val = ivp_count-1
        if (ivpLocArray[cn-1] == 'first char') and (cn != nop_ivprob-1):
            ivpLocArray.append('axis')
            ivpLA_lite.append('a')
            val -= 1
            tmp = val
        elif cn == nop_ivprob-1:
            ivpLocArray.append('axis')
            ivpLA_lite.append('a')
        elif tmp == 1:
            ivpLocArray.append('first char')
            ivpLA_lite.append('f')
        else:
            ivpLocArray.append('interior')
            ivpLA_lite.append('i')
            tmp -= 1

ivpLA_lite = ''.join(ivpLA_lite)

# for x in range(len(ivpLocArray)):
#     if ivpLocArray[x] == 'first char':
#         print(ivpLocArray[x])
#     else:
#         print(ivpLocArray[x], end=' ')
# print('\n')


class IniValProbReg:
    def __init__(self, aloc_pno, right, left):
        self.pno = aloc_pno+1
        self.pointloc = ivpLocArray[aloc_pno]
        self.Kmin = right.Kmin
        self.Kplus = left.Kplus
        self.x = 0
        self.y = 0
        self.u = 0
        self.v = 0
        if self.pointloc == 'axis':
            self.Kplus = -self.Kmin
        self.theta = 0.5*(self.Kmin + self.Kplus)
        self.nu = 0.5*(self.Kmin - self.Kplus)
        self.Ma = M(self.nu)
        self.T = gas_T*tempRatio(self.Ma, rosh)
        self.a = np.sqrt(rosh*R*self.T)
        self.V = self.a*self.Ma
        self.u = self.V*np.cos(np.radians(self.theta))
        self.v = self.V*np.sin(np.radians(self.theta))
        self.mu = np.degrees(np.arcsin(1/self.Ma))

        if (((ivpLocArray[aloc_pno-1] == 'first char') or (aloc_pno == 0))) and self.pointloc != 'axis':
            x1 = right.v0_xbar
            y1 = right.y
            x2 = left.v0_xbar
            y2 = left.y
            m1 = np.tan(np.radians((right.theta+self.theta)*0.5-(right.mu+self.mu)*0.5))
            m2 = np.tan(np.radians((left.theta+left.mu+self.theta+self.mu)/2))
            denom = 1 - m1/m2
            self.x = (x2 + (1/m2)*(y1 - y2 - x1*m1))/denom
            self.y = (y1 + m1*(x2 - x1 - y2/m2))/denom
            # print(self.x)
            # print(self.y)

        elif self.pointloc == 'axis':
            self.x = right.x-right.y/np.tan(np.radians(right.theta-right.mu))
            # print(right.pno)
            # print(right.theta)
            # print(right.mu)

        elif self.pointloc == 'interior' or self.pointloc == 'first char':
            x1 = right.x
            y1 = right.y
            x2 = left.x
            y2 = left.y
            m1 = np.tan(np.radians((right.theta+self.theta)*0.5-(right.mu+self.mu)*0.5))
            m2 = np.tan(np.radians((left.theta+left.mu+self.theta+self.mu)/2))
            denom = 1 - m1/m2
            self.x = (x2 + (1/m2)*(y1 - y2 - x1*m1))/denom
            self.y = (y1 + m1*(x2 - x1 - y2/m2))/denom
            # print(self.x)
            # print(self.y)


ivpReg = []
for x in range(nop_ivprob):
    adj = 0
    if x < nop_ivprob/2:
        adj = ivpLA_lite.count('f', 0, x)+1
    elif x >= nop_ivprob/2:
        adj = ivpLA_lite.count('f', x, len(ivpLA_lite))+1

    if ((ivpLocArray[x-1] == 'first char') or (x == 0)) and ivpLocArray[x] != 'axis':
        ivpReg.append(IniValProbReg(x, ivl[adj-1], ivl[adj]))
    elif (ivpLocArray[x] == 'interior') or (ivpLocArray[x] == 'first char') or (ivpLocArray[x] == 'axis'):
        ivpReg.append(IniValProbReg(x, ivpReg[x-adj], ivpReg[x-1]))


# print(" Point No. |         X |         Y |      M |     theta |        mu |    pt loc |")
# for x in range(len(ivpReg)):
#     print("{0:>11}|{1:>11.8f}|{2:>11.8f}|{3:>8.3f}|{4:>11.3f}|{5:>11.3f}| {6:<11}".format(ivpReg[x].pno, ivpReg[x].x, ivpReg[x].y, ivpReg[x].Ma, ivpReg[x].theta, ivpReg[x].mu, ivpReg[x].pointloc))

everychar = []
everycharType = [['contour']]
fc = [ivl[0]]
for x in range(len(ivpReg)):
    if ivpReg[x].pointloc == 'first char' or ivpReg[x].pno == len(ivpReg):
        fc.append(ivpReg[x])
        everycharType[0].append(ivpReg[x].pointloc)

everychar.append(fc)


class PredetPoints:
    def __init__(self, aloc_pno, o):
        self.pno = aloc_pno+1
        self.pointloc = 'contour'
        #print("o = ", o)
        if o <= int(throatAng):
            self.x = downWallRad*np.sin(np.radians(o))
            self.y = th_rad+downWallRad-downWallRad*np.cos(np.radians(o))
            self.theta = o
        else:
            self.x = downWallRad*np.sin(np.radians(throatAng))
            self.y = th_rad+downWallRad-downWallRad*np.cos(np.radians(throatAng))
            self.theta = throatAng
        #print('self.theta = ', self.theta)
        flag = iW = 0
        i = 1
        #mat.plot(x1, y1, 'o')
        while flag == 0:
            #print('i =', i)
            x3 = everychar[o-1][0].x
            y3 = everychar[o-1][0].y
            u3 = everychar[o-1][0].u
            v3 = everychar[o-1][0].v
            x1 = everychar[o-1][i].x
            y1 = everychar[o-1][i].y
            u1 = everychar[o-1][i].u
            v1 = everychar[o-1][i].v
            iW = invWall(x3,y3,u3,v3,x1,y1,u1,v1,self.x,self.y,self.theta)
            #print(iW)
            # if iW[2] > y1:
            #     mat.plot(x1, y1, 'o')
            if (y1 < iW[2]) and (iW[2] < self.theta):
                #print(y1 * 1000, iW[2] * 1000, self.y * 1000)
                flag = 1
            i += 1

        self.u = iW[0]
        self.v = iW[1]
        self.ytest = iW[2]
        self.V = np.sqrt(self.u**2+self.v**2)
        self.a = aspd(self.V)
        self.Ma = self.V/self.a
        # print("Ma = ", self.Ma)
        self.nu = Nu(self.Ma, rosh)
        self.mu = np.degrees(np.arcsin(1/self.Ma))
        self.Kmin = self.theta+self.nu
        self.corr = i-1
        self.switch = 'same'
        self.Kplus = 0
        if everychar[o-1][0].corr-self.corr >= 0:
            self.switch = 'increase'
        # print(everychar[o-1][0].corr, self.corr, self.switch)
        # print('spacer\n')


class Points:
    def __init__(self, aloc_pno, datum, right, left):
        self.pno = aloc_pno+1
        #print(aloc_pno, self.pno)
        #self.pzo = ptzo_array[aloc_pno]
        self.pointloc = everycharType[datum][aloc_pno]
        self.Kmin = 0
        self.Kplus = 0
        self.theta = 0
        self.nu = 0
        self.x = 0
        self.y = 0
        if self.pointloc == 'axis':
            self.Kmin = right.Kmin
            self.Kplus = -self.Kmin


        elif self.pointloc == 'interior':
            self.Kmin = right.Kmin
            self.Kplus = left.Kplus

        else:
            self.Kmin = right.Kmin
            self.Kplus = left.Kplus
            #print(self.pno, self.pointloc, self.Kmin, self.Kplus)

        self.theta = 0.5 * (self.Kmin + self.Kplus)
        self.nu = 0.5 * (self.Kmin - self.Kplus)
        self.Ma = M(self.nu)
        self.T = gas_T*tempRatio(self.Ma, rosh)
        self.a = np.sqrt(rosh*R*self.T)
        self.V = self.a*self.Ma
        self.u = self.V*np.cos(np.radians(self.theta))
        self.v = self.V*np.sin(np.radians(self.theta))
        # print("Ma = ", self.Ma)
        # print("Theta = ", self.theta)
        # print("V = ", self.V)
        # print("u = ", self.u)
        # print("v = ", self.v)
        self.mu = np.degrees(np.arcsin(1/self.Ma))

# Alternate method wall point plot
        if self.pointloc == 'contour':
            x1 = everychar[datum-1][0].x
            y1 = everychar[datum-1][0].y
            x2 = left.x
            y2 = left.y
            m1 = np.tan(np.radians(self.theta)) #everychar[datum-1][0].theta))
            m2 = np.tan(np.radians((left.theta+self.theta)*0.5+(left.mu+self.mu)*0.5))
            denom = 1 - (m1/m2)
            self.x = (x2 + (1/m2)*(y1 - y2 - x1*m1))/denom
            self.y = (y1 + m1*(x2 - x1 - (y2/m2)))/denom


# Characteristic drawing
#         if self.pno == 1:
#             self.x = right.x-right.y/np.tan(np.radians(right.theta-right.mu))

        if self.pointloc == 'axis':
            self.x = right.x-right.y/np.tan(np.radians(right.theta-right.mu))
            # print(right.pno)
            # print(right.theta)
            # print(right.mu)

        elif self.pointloc == 'interior':
            x1 = right.x
            y1 = right.y
            x2 = left.x
            y2 = left.y
            m1 = np.tan(np.radians((right.theta+self.theta)*0.5-(right.mu+self.mu)*0.5))
            m2 = np.tan(np.radians((left.theta+self.theta)*0.5+(left.mu+self.mu)*0.5))
            denom = 1 - m1/m2
            self.x = (x2 + (1/m2)*(y1 - y2 - x1*m1))/denom
            self.y = (y1 + m1*(x2 - x1 - y2/m2))/denom
        #print(datum, aloc_pno)
        #if self.pointloc == 'interior' and ptOnChar[aloc_pno-1].pointloc == 'contour':
            #mat.plot(self.x, self.y,'o')


class ContourPts:
    def __init__(self, aloc_pno, right, left):
        self.pno = aloc_pno+1
        self.corr = 1
        self.Kplus = left.Kplus
        self.Kmin = left.Kmin
        if self.pno == 1:
            self.theta = thetamax #(everychar[datum-1][0].theta+left.theta+left.mu)/2 #everychar[datum-1][0].theta
        else:
            self.theta = (wallpts[aloc_pno-1].theta+left.theta)/2

        # print((left.theta+left.mu)/2)
        # print(self.theta)
        self.nu = self.theta-self.Kplus
        #print(self.pno, self.pointloc, self.Kmin, self.Kplus)
        #print(self.pno,self.pointloc,self.Kmin,self.Kplus)
        self.Ma = M(self.nu)
        self.T = gas_T*tempRatio(self.Ma, rosh)
        self.a = np.sqrt(rosh*R*self.T)
        self.V = self.a*self.Ma
        self.u = self.V*np.cos(np.radians(self.theta))
        self.v = self.V*np.sin(np.radians(self.theta))
        self.mu = np.degrees(np.arcsin(1/self.Ma))
        self.aR = areaRatio(self.Ma, rosh)
        # print("Ma = ", self.Ma)
        # print("Theta = ", self.theta)
        # print("V = ", self.V)
        # print("u = ", self.u)
        # print("v = ", self.v)

# Alternate method wall point plot
        if self.pno == 1:
            x1 = right.x
            y1 = right.y
            x2 = left.x
            y2 = left.y
            m1 = np.tan(np.radians(self.theta)) #everychar[datum-1][0].theta))
            m2 = np.tan(np.radians(left.theta+left.mu))
            denom = 1 - (m1/m2)
            self.x = (x2 + (1/m2)*(y1 - y2 - x1*m1))/denom
            self.y = (y1 + m1*(x2 - x1 - (y2/m2)))/denom
        else:
            x1 = right.x
            y1 = right.y
            x2 = left.x
            y2 = left.y
            m1 = np.tan(np.radians(right.theta))  # everychar[datum-1][0].theta))
            m2 = np.tan(np.radians(left.theta+left.mu))
            denom = 1 - (m1 / m2)
            self.x = (x2 + (1 / m2) * (y1 - y2 - x1 * m1)) / denom
            self.y = (y1 + m1 * (x2 - x1 - (y2 / m2))) / denom


i = 1
nop_preDetReg = len(fc)-1
flag = 0
while flag == 0:
    j = 0
    ptOnChar = []
    ptOnChar.append(PredetPoints(j, i))
    everycharType.append([ptOnChar[j].pointloc])

    if ptOnChar[0].switch == 'increase':
        jmax = len(everychar[i-1])
        stp = 0
        #print('increase to ', jmax)
    else:
        jmax = len(everychar[i-1])
        stp = 1
        #print('same as ', jmax)

    j = 1
    #print('i = {0}, j ='.format(i),j)

    while j<=jmax:
        if j<jmax-stp:
            everycharType[i].append('interior')
        elif j == jmax-stp:
            everycharType[i].append('axis')
        j+=1
    # print('length of prev char = ', len(everycharType[i-1]), everycharType[i-1])
    # print('length of curr char = ', len(everycharType[i]), everycharType[i])
    j = 1
    while (j+stp)<=jmax:
        if everycharType[i][j] == 'interior':
            #print('des j int =', j)
            #print(everychar[i-1][j].pointloc)
            if ptOnChar[j-1].pointloc == 'contour':
                ptOnChar.append(Points(j, i, ptOnChar[0], everychar[i-1][j+stp]))
            else:
                ptOnChar.append(Points(j, i, ptOnChar[j-1], everychar[i-1][j+stp]))
        elif everycharType[i][j] == 'axis':
            #print('des j axi =', j)
            ptOnChar.append(Points(j, i, ptOnChar[j-1], ptOnChar[j-1]))
        j += 1
    nop_preDetReg = nop_preDetReg+jmax+1
    everychar.append(ptOnChar)
    if ptOnChar[0].theta == throatAng:
        flag = 1
    i += 1

# for x in range(len(everychar)):
#     print(x,'len = ', len(everychar[x]), end=' ')
#     for y in range(len(everychar[x])):
#         print(everychar[x][y].pointloc, end=' ')
#     print('')

arcEnd = i = len(everychar)-1
#print('arcEnd = ', arcEnd)
wallpts = []
j = 0
ptOnChar = []
wallpts.append(ContourPts(j, everychar[arcEnd][0], everychar[i-1][j+1]))
j = 1
jmax = len(everychar[arcEnd])-1
#print(jmax)
#print('y = ', everychar[i][35].y)
while j<=jmax:
    wallpts.append(ContourPts(j, wallpts[j-1], everychar[i][j]))
    # print("Ma = ", ptOnChar[j].Ma)
    #print('theta = ', ptOnChar[j].theta)
    # print('mu = ', ptOnChar[j].mu)
    j += 1
#everychar.append(ptOnChar)

nop_wallpts = len(wallpts)
totpts = nop_preDetReg+nop_wallpts

# for x in range(len(wallpts)):
#     print(wallpts[x].theta)
#     mat.plot(wallpts[x].x, wallpts[x].y, 'o')



# for x in range(len(everychar)):
#     for y in range(len(everychar[x])):
#         print('{0:>8.4f}'.format(everychar[x][y].x), end=' ')
#     print('')


'''

flag = 0
OG_ExMa = ExMa
truncFlag = 0
count = 0
limx = 0
limy = 0
trunc_pno = 0
pts = []
# while flag == 0:
#     #print("try ", count)
#     thetamax = Nu(ExMa, rosh) / 2
#     dtheta = round(thetamax/noc,4)
#     altdtheta = round(throatAng / noc, 3)
#     pts=[]
#     tmp = noc+1
#     M_array = np.arange(1, ExMa, 0.0001)
#     aRtab = []
#     for assign in range(len(M_array)):
#         aRtab.append(areaRatio(M_array[assign], rosh))
#     #arcX = downWallRad * np.sin(np.radians(thetamax))
#     #arcY = th_rad + downWallRad - downWallRad * np.cos(np.radians(thetamax))
#     for x in range(nop):
#         #print(x, tmp, x - tmp, x - 1)
#         if ptloc_array[x] == 'axis' and x>noc and x>0:
#             tmp-=1
#             #print("------------tmp reduction, tmp = ",tmp)
#
#         if x == 0:
#             pts.append(Points(x,x,x))
#             #print("--first cond pno={}".format(x+1))
#         elif x != 0 and x < noc+1:
#             pts.append(Points(x, x, pts[x-1]))
#             #print("--second cond pno = {}    pno under = {}".format(x+1,x))
#         elif x >= noc+1:
#             pts.append(Points(x,pts[x-tmp],pts[x-1]))
#             # print("--third cond pno = {}  pno left = {} pno under = {} Kmin = {:>6.2f}  Kplus = {:>6.2f}    ".format(x+1,x+1-tmp,x,pts[x-tmp].Kmin, pts[x-tmp].Kplus))
#     if ExMa <= 5*OG_ExMa:
#         for x in range(nop):
#             if pts[x].pointloc == 'wall':
#                 #print(abs(nozzExAng-pts[x].theta))
#                 # time.sleep(0.5)
#                 if (abs(nozzExAng-pts[x].theta) <= dtheta/2): # and (pts[x].theta-nozzExAng >= 0):
#                     #print("I'm in")
#                     if (pts[x].Ma-OG_ExMa <= 0.05) and (pts[x].Ma-OG_ExMa >= 0):
#                         flag = 1
#                         truncFlag = 1
#                         limx = pts[x].x
#                         limy = pts[x].y
#                         trunc_pno = pts[x].pno
#                         #print(ExMa)
#                         break
#                     else:
#                         ExMa += 0.01
#                         break
#     else:
#         print("Unable to generate TICN")
#         flag = 1
#         truncFlag = 0
#     # print("ExMa = ", round(ExMa,5), x, pts[nop-1].x)
#     # print("Tried ExMa = ", ExMa, " at ", pts[x].x, "pt Mach = ", pts[x].Ma," diff = ", abs(pts[x].Ma-OG_ExMa))
#     count += 1
#     #print(count)

if truncFlag == 0:
    ExMa = OG_ExMa
    thetamax = Nu(ExMa, rosh)/2
    dtheta = round(thetamax/noc, 3)
    altdtheta = round(throatAng/noc, 3)
    pts = []
    tmp = noc + 1
    M_array = np.arange(1, OG_ExMa, 0.0001)
    aRtab = []
    for assign in range(len(M_array)):
        aRtab.append(areaRatio(M_array[assign], rosh))
    #arcX = downWallRad * np.sin(np.radians(thetamax))
    #arcY = th_rad + downWallRad - downWallRad * np.cos(np.radians(thetamax))
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
'''

def displayTable(): #Literally does what it looks like it's gonna do
    print("-----------|------------|------------|-----------|-----------|----------|------------")
    print("  Point no.|       Kmin |      Kplus |  θ(Theta) |     ν(Nu) |        M |          V ")
    print("-----------|------------|------------|-----------|-----------|----------|------------")
    tot_pts = 1
    for a in range(len(everychar)):
        for b in range(len(everychar[a])):
            print(" {0:>9} |   {1:>8.3f} |   {2:>8.3f} |  {3:>8.3f} |  {4:>8.3f} | {5:>8.3f} |   {6:>8.3f}".format(tot_pts, everychar[a][b].Kmin, everychar[a][b].Kplus, everychar[a][b].theta, everychar[a][b].nu, everychar[a][b].Ma, everychar[a][b].V))
            tot_pts+=1
    for a in range(len(wallpts)):
        print(" {0:>9} |   {1:>8.3f} |   {2:>8.3f} |  {3:>8.3f} |  {4:>8.3f} | {5:>8.3f} |   {6:>8.3f}".format(tot_pts, wallpts[a].Kmin, wallpts[a].Kplus, wallpts[a].theta, wallpts[a].nu, wallpts[a].Ma, wallpts[a].V))
        tot_pts +=1
    return 0


def displayFullTable(): #Literally does what it looks like it's gonna do
    print("-----------|------------|------------|-----------|-----------|----------|------------")
    print("  Point no.|       Kmin |      Kplus |  θ(Theta) |     ν(Nu) |        M |      μ(mu) ")
    print("-----------|------------|------------|-----------|-----------|----------|------------")
    for x in range(nop_preDetReg):
        print(" {0:>9} |   {1:>8.3f} |   {2:>8.3f} |  {3:>8.3f} |  {4:>8.3f} | {5:>8.3f} |   {6:>8.3f}".format(pts[x].pno,pts[x].Kmin,pts[x].Kplus,pts[x].theta,pts[x].nu,pts[x].Ma,pts[x].mu))
    return 0


def drawChar():
    colour = 'red'
    mat.plot(ivlX, ivlY, 'green')
    mat.plot(conx, cony, 'black')
    for h in range(round(nop_ivprob)):
        k = 0
        if h < nop_ivprob/2:
            k = ivpLA_lite.count('f', 0, h)+1
            #print("before halfway point")
        elif h >= nop_ivprob/2:
            k = ivpLA_lite.count('f', h, len(ivpLA_lite))+1
            #print("after halfway point")

        #print('k = ', k)

        if ((ivpLocArray[h-1] == 'first char') or (h == 0)) and ivpLocArray[h] != 'axis':
            mat.plot([ivl[k-1].v0_xbar, ivpReg[h].x], [ivl[k-1].y, ivpReg[h].y], colour, linewidth = 0.5)
            mat.plot([ivl[k].v0_xbar, ivpReg[h].x], [ivl[k].y, ivpReg[h].y], colour, linewidth = 0.5)
            #print(h, ivl[k-1].pno, ivl[k-1].v0_xbar, ivpReg[h].x, ivl[k-1].y, ivpReg[h].y)
        elif ivpLocArray[h] == 'axis':
            mat.plot([ivpReg[h-k].x, ivpReg[h].x], [ivpReg[h-k].y, ivpReg[h].y], colour, linewidth = 0.5)

        elif (ivpLocArray[h] == 'interior') or (ivpLocArray[h] == 'first char'):
            mat.plot([ivpReg[h-k].x, ivpReg[h].x], [ivpReg[h-k].y, ivpReg[h].y], colour, linewidth = 0.5)
            mat.plot([ivpReg[h-1].x, ivpReg[h].x], [ivpReg[h-1].y, ivpReg[h].y], colour, linewidth = 0.5)
            #mat.plot([prepts[h].x, pts[h].x], [prepts[h].y, pts[h].y], 'o', linewidth = 0.5)

    i = 0
    while i < len(everychar):
        #print('\n')
        k = everychar[i][0].corr-1
        for j in range(len(everychar[i])):
            #print('i =', i, ' j =', j)
            if everychar[i][j].pointloc == 'contour' or everychar[i][j].pointloc == 'interior':
                if i != 0 or j != 0:
                    mat.plot([everychar[i][j+1].x, everychar[i][j].x],[everychar[i][j+1].y,everychar[i][j].y], colour, linewidth=0.5)
                    mat.plot([everychar[i-1][j+k].x, everychar[i][j].x],[everychar[i-1][j+k].y,everychar[i][j].y], colour, linewidth=0.5)
                else:
                    mat.plot([everychar[i][j+1].x, everychar[i][j].x],[everychar[i][j+1].y,everychar[i][j].y], colour, linewidth=0.5)

        i+=1

    i = 1
    while i < nop_wallpts:
        mat.plot([everychar[arcEnd][i].x, wallpts[i].x], [everychar[arcEnd][i].y, wallpts[i].y], colour, linewidth=0.5)
        i+=1


def ICN(): #This function draws the ideal contour nozzle
    infconX = conx.copy()
    infconY = cony.copy()
    i = 0
    while i < len(wallpts):
        conx.append(wallpts[i].x)
        cony.append(wallpts[i].y)
        stConX.append(wallpts[i].x)
        stConY.append(wallpts[i].y)
        i+=1


    #mat.plot(infconX, infconY, color='black')
    mat.plot(conx,cony, color = 'black')
    #mat.plot(stConX, stConY, color='black')
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
        print(everychar[arcEnd][0].x, conx[xd])
        if abs(conx[xd]-everychar[arcEnd][0].x)<=0.0001:
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
    poly_coeff = np.polyfit(xvar, yvar, 2)
    print(xvar, yvar)
    print(poly_coeff)
    #print(a, b, c, d)
    print(poly_coeff[0])

    xplot = np.arange(conx[cutoff], conx[cutoff+1], 0.00001)
    #yplot3 = poly_coeff[3] + poly_coeff[2]*xplot + poly_coeff[1]*xplot**2 + poly_coeff[0]*xplot**3
    for lol in range(len(conx)):
        if conx[lol] == everychar[arcEnd][0].x:
            val = len(xplot)-1
            while val >= 0:
                conx.insert(lol, xplot[val])
                cony.insert(lol, poly_coeff[3] + poly_coeff[2]*xplot[val] + poly_coeff[1]*xplot[val]**2) #+ poly_coeff[0]*xplot[val]**3)
                val-=1
            break

    #mat.plot(xplot, yplot3, 'blue')
    mat.plot(conx, cony, 'black')


def XLimporter():
    dir = "D:\Study Files\Aerospace Engg\Project\Output\\"
    wb = xlwt.Workbook()
    props = wb.add_sheet("Flow properties") #, cell_overwrite_ok=True)
    coords = wb.add_sheet("Coordinates")
    props.write(0, 0, "Point no.")
    props.write(0, 1, "Kmin")
    props.write(0, 2, "Kplus")
    props.write(0, 3, "θ(Theta)")
    props.write(0, 4, "ν(Nu)")
    props.write(0, 5, "M")
    props.write(0, 6, "μ(mu)")
    tot_pts = 1
    for a in range(len(everychar)):
        for b in range(len(everychar[a])):
            print('pt.no = ', tot_pts)
            props.write(tot_pts, 0, tot_pts)
            props.write(tot_pts, 1, round(everychar[a][b].Kmin, 4))
            props.write(tot_pts, 2, round(everychar[a][b].Kplus, 4))
            props.write(tot_pts, 3, round(everychar[a][b].theta, 4))
            props.write(tot_pts, 4, round(everychar[a][b].nu, 4))
            props.write(tot_pts, 5, round(everychar[a][b].Ma, 4))
            props.write(tot_pts, 6, round(everychar[a][b].mu, 4))
            tot_pts += 1

    for a in range(len(wallpts)):
        props.write(tot_pts, 0, tot_pts)
        props.write(tot_pts, 1, round(wallpts[a].Kmin, 4))
        props.write(tot_pts, 2, round(wallpts[a].Kplus, 4))
        props.write(tot_pts, 3, round(wallpts[a].theta, 4))
        props.write(tot_pts, 4, round(wallpts[a].nu, 4))
        props.write(tot_pts, 5, round(wallpts[a].Ma, 4))
        props.write(tot_pts, 6, round(wallpts[a].mu, 4))
        tot_pts += 1

    coords.write(0, 0, "Point Number")
    coords.write(0, 1, "X")
    coords.write(0, 2, "Y")
    for val in range(len(conx)):
        coords.write(val+1, 0, val+1)
        coords.write(val+1, 1, round(conx[val],6))
        coords.write(val+1, 2, round(cony[val],6))

    wb.save("{0:}ICN(sonic line), ExMa={1:}, Throat rad={2:}, No. of pts on sonic line={3:}.xls".format(dir, ExMa, th_rad, ivp_count))
    return 0


def everyLoc():
    print("-----------|------------|------------|------------|----------")
    print("  Point no.|          X |          Y |      Theta |       mu ")
    print("-----------|------------|------------|------------|----------")
    for x in range(noc):
        print(" {0:>9} |   {1:>8.3f} |   {2:>8.3f} |   {3:>8.3f} |   {4:>8.3f} | {5:>3}".format(pts[x].pno, pts[x].x, pts[x].y, pts[x].theta, pts[x].mu, pts[x].pointloc))
    return 0


# actAR = (np.pi*pts[nop-1].y**2)/(np.pi*(th_rad*2)**2)
# print(np.pi, pts[nop-1].y)
# print("Actual Area Ratio = {:.4f}".format(actAR))
# print("Target Area Ratio = {:.4f}".format(areaRatio(ExMa,rosh)))

# Enable or disable the functions down below by removing the '#' or putting the '#' in
# front of the functions respectively

#print('dtheta = ', dtheta)

drawChar()       # This function draws characteristic lines
ICN()        # This function draws the minimum length nozzle
displayTable()   # This table shows flow properties along the points in the flow
#dispLocation()   # Shows coordinates of all wallpoints as well as characteristic intersection points in Cartesian coords sys
#everyLoc()
#curveGen()       # Attempts to generate smooth regression curve of nozzle points but fails miserably
#XLimporter()     # Imports coords to an excel file
#dispLocation()

# print(everychar[-1][-2].Ma)
# print(everychar[-1][-2].y)


# mat.ylim(0, 0.2)
# mat.xlim(0, 1)


mat.axis('scaled')
mat.xlabel("X-axis (m)")
mat.ylabel("Y-axis (m)")
mat.title("Ideal Contour Nozzle (considering Curved Sonic Line)")
#if truncFlag == 0:
mat.ylim(0,(wallpts[nop_wallpts-1].y+0.1*wallpts[nop_wallpts-1].y))
mat.xlim(0,(wallpts[nop_wallpts-1].x+0.1*wallpts[nop_wallpts-1].x))
# elif truncFlag == 1:
#     mat.ylim(0, (limy + 0.1*limy))
#     mat.xlim(0, (limx + 0.1*limx))
mat.grid()
mat.show()
#mat.savefig("hellooooo.png", dpi=1000)

# while 1:
#     importquery = input("Import data to an excel file?(Y/N):")
#     if importquery == 'Y' or importquery == 'y':
#         XLimporter()
#         break
#     elif importquery == 'N' or importquery == 'n':
#         break
#     else:
#         print("Invalid option, please try again")
