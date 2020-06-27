import matplotlib.pyplot as mat
import numpy as np
x = [0,1,2,3,4]
y = [0,1,4,9,16]
s_x = s_y = s_x2 = s_x3 = s_x4 = s_x5 = s_x6 = s_xy = s_x2y = s_x3y = 0

n=5
print(len(x))
for i in range(len(x)):
    s_x += x[i]
    s_y += y[i]
    s_x2 += round(x[i]**2,4)
    s_x3 += round(x[i] ** 3, 4)
    s_x4 += round(x[i] ** 4, 4)
    s_x5 += round(x[i] ** 5, 4)
    s_x6 += round(x[i] ** 6, 4)
    s_xy += round((x[i]*y[i]), 4)
    s_x2y += round(((x[i]**2)*y[i]), 4)
    s_x3y += round(((x[i] ** 3) * y[i]), 4)


s_x3y = round(s_x3y,3)
s_x2y = round(s_x2y,3)
s_xy = round(s_xy,3)
s_y = round(s_y,3)
print(x)
print(y)
print(s_x, s_y, s_xy, s_x2y, s_x3y, s_x2, s_x3, s_x4,s_x5, s_x6)

a = (np.linalg.det(np.array([[s_y,s_x,s_x2,s_x3],[s_xy,s_x2,s_x3,s_x4],[s_x2y,s_x3,s_x4,s_x5],[s_x3y,s_x4,s_x5,s_x6]])))/(np.linalg.det(np.array([[n,s_x,s_x2,s_x3],[s_x,s_x2,s_x3,s_x4],[s_x2,s_x3,s_x4,s_x5],[s_x3,s_x4,s_x5,s_x6]])))
b = (np.linalg.det(np.array([[n,s_y,s_x2,s_x3],[s_x,s_xy,s_x3,s_x4],[s_x2,s_x2y,s_x4,s_x5],[s_x3,s_x3y,s_x5,s_x6]])))/(np.linalg.det(np.array([[n,s_x,s_x2,s_x3],[s_x,s_x2,s_x3,s_x4],[s_x2,s_x3,s_x4,s_x5],[s_x3,s_x4,s_x5,s_x6]])))
c = (np.linalg.det(np.array([[n,s_x,s_y,s_x3],[s_x,s_x2,s_xy,s_x4],[s_x2,s_x3,s_x2y,s_x5],[s_x3,s_x4,s_x3y,s_x6]])))/(np.linalg.det(np.array([[n,s_x,s_x2,s_x3],[s_x,s_x2,s_x3,s_x4],[s_x2,s_x3,s_x4,s_x5],[s_x3,s_x4,s_x5,s_x6]])))
d = (np.linalg.det(np.array([[n,s_x,s_x2,s_y],[s_x,s_x2,s_x3,s_xy],[s_x2,s_x3,s_x4,s_x2y],[s_x3,s_x4,s_x5,s_x3y]])))/(np.linalg.det(np.array([[n,s_x,s_x2,s_x3],[s_x,s_x2,s_x3,s_x4],[s_x2,s_x3,s_x4,s_x5],[s_x3,s_x4,s_x5,s_x6]])))

x = (s_y*((s_x2*s_x4)-(s_x3*s_x3)) - s_x*((s_xy*s_x4)-(s_x3*s_x2y)) + s_x2*((s_xy*s_x3)-(s_x2*s_x2y)))/(n*((s_x2*s_x4)-(s_x3*s_x3)) - s_x*((s_x*s_x4)-(s_x3*s_x2)) + s_x2*((s_x*s_x3)-(s_x2*s_x2)))
y = (n*((s_xy*s_x4)-(s_x3*s_x2y)) - s_y*((s_x*s_x4)-(s_x3*s_x2)) + s_x2*((s_x*s_x2y)-(s_xy*s_x2)))/(n*((s_x2*s_x4)-(s_x3*s_x3)) - s_x*((s_x*s_x4)-(s_x3*s_x2)) + s_x2*((s_x*s_x3)-(s_x2*s_x2)))
z = (n*((s_x2*s_x2y)-(s_xy*s_x3)) - s_x*((s_x*s_x2y)-(s_xy*s_x2)) + s_y*((s_x*s_x3)-(s_x2*s_x2)))/(n*((s_x2*s_x4)-(s_x3*s_x3)) - s_x*((s_x*s_x4)-(s_x3*s_x2)) + s_x2*((s_x*s_x3)-(s_x2*s_x2)))

a4 = (np.linalg.det(np.array([[s_y,s_x,s_x2,s_x3,s_x4],[s_xy,s_x2,s_x3,s_x4,s_x5],[s_x2y,s_x3,s_x4,s_x5,s_x6],[s_x3y,s_x4,s_x5,s_x6,s_x7],[s_x4y,s_x5,s_x6,s_x7,s_x8]])))/(np.linalg.det(np.array([[n,s_x,s_x2,s_x3,s_x4],[s_x,s_x2,s_x3,s_x4,s_x5],[s_x2,s_x3,s_x4,s_x5,s_x6],[s_x3,s_x4,s_x5,s_x6,s_x7],[s_x4,s_x5,s_x6,s_x7,s_x8]])))
b4 = (np.linalg.det(np.array([[n,s_y,s_x2,s_x3,s_x4],[s_x,s_xy,s_x3,s_x4,s_x5],[s_x2,s_x2y,s_x4,s_x5,s_x6],[s_x3,s_x3y,s_x5,s_x6,s_x7],[s_x4,s_x4y,s_x6,s_x7,s_x8]])))/(np.linalg.det(np.array([[n,s_x,s_x2,s_x3,s_x4],[s_x,s_x2,s_x3,s_x4,s_x5],[s_x2,s_x3,s_x4,s_x5,s_x6],[s_x3,s_x4,s_x5,s_x6,s_x7],[s_x4,s_x5,s_x6,s_x7,s_x8]])))
c4 = (np.linalg.det(np.array([[n,s_x,s_y,s_x3,s_x4],[s_x,s_x2,s_xy,s_x4,s_x5],[s_x2,s_x3,s_x2y,s_x5,s_x6],[s_x3,s_x4,s_x3y,s_x6,s_x7],[s_x4,s_x5,s_x4y,s_x7,s_x8]])))/(np.linalg.det(np.array([[n,s_x,s_x2,s_x3,s_x4],[s_x,s_x2,s_x3,s_x4,s_x5],[s_x2,s_x3,s_x4,s_x5,s_x6],[s_x3,s_x4,s_x5,s_x6,s_x7],[s_x4,s_x5,s_x6,s_x7,s_x8]])))
d4 = (np.linalg.det(np.array([[n,s_x,s_x2,s_y,s_x4],[s_x,s_x2,s_x3,s_xy,s_x5],[s_x2,s_x3,s_x4,s_x2y,s_x6],[s_x3,s_x4,s_x5,s_x3y,s_x7],[s_x4,s_x5,s_x6,s_x4y,s_x8]])))/(np.linalg.det(np.array([[n,s_x,s_x2,s_x3,s_x4],[s_x,s_x2,s_x3,s_x4,s_x5],[s_x2,s_x3,s_x4,s_x5,s_x6],[s_x3,s_x4,s_x5,s_x6,s_x7],[s_x4,s_x5,s_x6,s_x7,s_x8]])))
e4 = (np.linalg.det(np.array([[n,s_x,s_x2,s_x3,s_y],[s_x,s_x2,s_x3,s_x4,s_xy],[s_x2,s_x3,s_x4,s_x5,s_x2y],[s_x3,s_x4,s_x5,s_x6,s_x3y],[s_x4,s_x5,s_x6,s_x7,s_x4y]])))/(np.linalg.det(np.array([[n,s_x,s_x2,s_x3,s_x4],[s_x,s_x2,s_x3,s_x4,s_x5],[s_x2,s_x3,s_x4,s_x5,s_x6],[s_x3,s_x4,s_x5,s_x6,s_x7],[s_x4,s_x5,s_x6,s_x7,s_x8]])))

print(a,b,c,d)
print(x,y,z)
xplot = np.arange(-100,100,0.01)
yplot2 = a + b*xplot + c*xplot**2
yplot3 = a + b*xplot + c*xplot**2 + d*xplot**3
mat.plot(xplot,yplot2,'r')
mat.plot(x,y)
mat.plot(xplot,yplot3,'b')
mat.grid(1)
mat.show()