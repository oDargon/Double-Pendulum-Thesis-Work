import matplotlib.pyplot as plt
import numpy as np
import time
import csv
import pandas as pd
from matplotlib.colors import LinearSegmentedColormap, colorConverter, ListedColormap
from matplotlib import colors
pi = 3.14159265359 

my_colormap = LinearSegmentedColormap.from_list('my_colormap',[(0,'#FFFFFF'),(0,'#FFFFFF'),(0,'#ccff33'),(0.1,'#134611'),(0.1,'#ff002b'),(1,'#660708')],N = 50000)
my_colormapBlue = LinearSegmentedColormap.from_list('my_colormap',[(0,'#FFFFFF'),(0,'#FFFFFF'),(0,'#ccff33'),(0.01,'#134611'),(0.01,'#ff002b'),(0.1,'#660708'),(0.1,'#48cae4'),(1,'#03045e')],N = 500000)
my_colormapPurple = LinearSegmentedColormap.from_list('my_colormap',[(0,'#FFFFFF'),(0,'#FFFFFF'),(0,'#ccff33'),(0.001,'#134611'),(0.001,'#ff002b'),(0.01,'#660708'),(0.01,'#48cae4'),(0.1,'#03045e'),(0.1,'#bf1eea'),(1,'#701089')],N = 500000)

colr = LinearSegmentedColormap.from_list('my_colormap',[(0,'#FFFFFF'),(0,'#FFFFFF'),(0,'#ccff33'),(0.3333,'#134611'),(0.3333,'#ff002b'),(0.6666,'#660708'),(0.6666,'#48cae4'),(1,'#03045e')],N = 500000)


data = pd.read_csv("FullFracLyap1000.csv",header=None)
data2 = pd.read_csv("MaxLyaparm.csv",header=None)
# data2 = np.ma.masked_where(data2 == 0, data2)

plt.figure(figsize=(7, 6))
pl2 = plt.imshow(data2.iloc[::-1],extent=[100,130,70,100], cmap = "inferno")
# pl2 = plt.imshow(data2,extent=[-180,180,-180,180], cmap = "inferno")
# pl = plt.imshow(data,extent=[-102,-99.5,-117,-112.5], cmap = colr, norm=colors.LogNorm(vmin= 1, vmax=1000))
plt.colorbar(pl2,pad = 0.01)
# plt.imshow(data2,extent=[-102,-99.5,-117,-112.5],interpolation = "nearest", cmap = "inferno")
# plt.colorbar(label="Color Scale", orientation="vertical")

axes=plt.gca()
axes.set_aspect(1)
plt.xlabel(r"$\theta_1$", fontsize=16)
plt.ylabel(r"$\theta_2$", fontsize=16,rotation=0)
plt.savefig("ArmFracLyap.pdf",dpi = 2500)
plt.show()

