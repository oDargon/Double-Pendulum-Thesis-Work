from utility import gF
from numpy import pi, linspace
import time
# from Animator import animator2
from toolbox import LyapDat, LyapPlot2, FLyapDat
import matplotlib.pyplot as plt


y = gF("hamil.jl")

B1 = [1,1,1,1]
C1 = [0,60,0.01]
t1 = (0,60,6001)
u0 = (90,180,0,0)

Aa = ((pi/180)*(120), (pi/180)*230, 0, 0)
Ab = ((pi/180)*(120), (pi/180)*230.00001, 0, 0)
FA = y(Aa,B1,C1)
FB = y(Ab,B1,C1)
# Hi = rk45(u0,t1)

LyapPlot2(FA,FB,t1,20,0)
# print(FLyapDat(FA,FB,t1,20))
# plt.figure()  
# plt.plot(linspace(0,60,6001),abs(FA[0]-Hi[0]))
# plt.ylabel(r"$|\theta_{RK45}-\theta_{Vern9}|_1$ (rad)")
# plt.xlabel("time (s)")
# plt.show()

