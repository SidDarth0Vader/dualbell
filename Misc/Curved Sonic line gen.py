import matplotlib.pyplot as mat
import numpy as np

throatAng = 23.9
th_rad = 0.03354
delta = 1
rosh = 1.4
gas_T = 400
R = 287
a_star = np.sqrt((2*rosh*R*gas_T)/(rosh+1))
a0sq = rosh*R*gas_T
print("a* = ", a_star)
upWallRad = 5.5*th_rad
downWallRad = 0.5*th_rad
alpha = np.sqrt((1+delta)/((rosh+1)*upWallRad*th_rad))
print("alpha = ", alpha)
epsilon = -((th_rad/(2*(3+delta)))*np.sqrt(((rosh+1)*(1+delta))/(upWallRad/th_rad)))
print("Epsilon = ", epsilon)


conx = []
cony = []
infPtLocX = downWallRad*np.sin(np.radians(throatAng))
xplot = np.arange(0, infPtLocX+(downWallRad/100), 0.0001)
for a in range(len(xplot)):
    conx.append(xplot[a])
    cony.append(-np.sqrt(downWallRad**2 - xplot[a]**2) + downWallRad + th_rad)
    # infconX.append(xplot[a])
    # infconY.append(-np.sqrt(downWallRad**2 - xplot[a]**2) + downWallRad + th_rad)


def xsonic(yy):
    return -((rosh+1)*alpha*(yy**2))/(2*(1+delta))


def xcorr(yy):
    return -((rosh+1)*alpha*(yy**2))/(2*(3+delta))


def ubar(xx,yy):
    return a_star*(1+alpha*xx+(((rosh+1)*(alpha**2)*(yy**2))/(2*(1+delta))))


def a(xx,yy):
    return np.sqrt(a0sq - ((rosh-1)/2)*(ubar(xx,yy)**2))


y = []
xson = []
v0_x = []
v0_xbar = []
M = []
T_rat = [] # 1 + ((gamma - 1)/2)*M**2
dy = th_rad/10
ucap = []
i = 0
while i <= 10:
    if i == 0:
        y.append(0)
    elif i != 0:
        y.append(y[i-1]+dy)
    xson.append(xsonic(y[i]))
    v0_x.append((xcorr(y[i])))
    v0_xbar.append(xcorr(y[i])-epsilon)
    ucap.append(ubar(v0_x[i], y[i]))
    M.append(ubar(v0_x[i], y[i])/a(v0_x[i], y[i]))
    T_rat.append(1/(1 + ((rosh - 1)/2)*M[i]**2))
    i += 1
print(epsilon)
# x = - ((rosh+1)*alpha*(y**2))/(2*(3+delta))

print("     M |          u |      alt a |")
for x in range(len(y)):
    print("{0:>7.3f}|{1:>12.6f}|{2:>12.6f}|".format(M[x], ucap[x], (M[x]*np.sqrt(rosh*R*gas_T*T_rat[x])))) #a(v0_x[x], y[x]), np.sqrt(rosh*R*gas_T*T_rat[x])))

mat.plot(conx, cony, 'black')
# mat.plot(xson,y, ':', color = 'blue', linewidth = 0.7)
mat.plot(v0_xbar, y, 'o')
mat.plot(v0_xbar, y, color='red', linewidth=0.7)
mat.axis('scaled')
# mat.xlim(-0.01, 0.03)
mat.grid()
mat.show()
