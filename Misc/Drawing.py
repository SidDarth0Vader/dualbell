import matplotlib.pyplot as mat
# import matplotlib.lines as lin
# import numpy as np
# # import matplotlib.pyplot as plt
# # import matplotlib.lines as mlines
# #
# # con = [(0,8.8),(7.76,10.62),(9.01,10.89),(10.25,11.13),(11.69,11.38),(13.38,11.65),(15.40,11.93)]
# #
# # plt.contour(con,0)
# #
# # x = np.arange(1, 10)
# # y = x.reshape(-1, 1)
# # h = x * y
# #
# # cs = plt.contourf(con, levels=[10,50])
# #
# #
# # def newline(p1, p2):
# #     ax = plt.gca()
# #     xmin, xmax = ax.get_xbound()
# #
# #     if(p2[0] == p1[0]):
# #         xmin = xmax = p1[0]
# #         ymin, ymax = ax.get_ybound()
# #     else:
# #         ymax = p1[1]+(p2[1]-p1[1])/(p2[0]-p1[0])*(xmax-p1[0])
# #         ymin = p1[1]+(p2[1]-p1[1])/(p2[0]-p1[0])*(xmin-p1[0])
# #
# #     l = mlines.Line2D([xmin,xmax], [ymin,ymax])
# #     ax.add_line(l)
# #     return l
# #
# # # x = np.linspace(0,100)
# # # y = x**2
# #
# # p1 = [1,20]
# # p2 = [6,70]
# #
# # # plt.plot(x, y)
# # # newline(p1,p2)
# # # newline(p2,[5,20])
# # plt.show()

from matplotlib.figure import Figure
from matplotlib.axes import Axes
from matplotlib.lines import Line2D
from matplotlib.backends.backend_agg import FigureCanvasAgg

# fig = Figure(figsize=[4, 4])
# ax = Axes(fig, [.1, .1, .8, .8])
# fig.add_axes(ax)
# l = Line2D([0, 1], [0, 1])
# ax.add_line(l)
#
# canvas = FigureCanvasAgg(fig)
# canvas.print_figure("line_ex.png")

conx = [0,7.76,9.01,10.25,11.69,13.38,15.40,17.4]
cony = [8.8,10.62,10.89,11.13,11.38,11.65,11.93,12]
mat.ylim(0,20)
mat.xlim(0,20)
mat.plot(conx,cony)
mat.show()