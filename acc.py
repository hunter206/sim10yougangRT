# -*- coding: utf-8 -*-
"""
Spyder Editor

This is a temporary script file.

"""
"用来仿真双护盾张总提出的十油缸算法"

from math import sin, cos
from scipy import optimize
from numpy import *
import numpy as np
import os
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
os.system('clear') 

'杆长'
l1 = 10
l2 = 10
l3 = 12
l4 = 12
l5 = 10
l6 = 10

'前盾铰接点在前盾坐标系下坐标'
Pd = np.matrix([[0,0,0],
                [7,0,0],
                [7,1,0],
                [4,4,0],
                [3,4,0],
                [0,1,0]
                ])

'支撑盾铰接点在前盾坐标系下坐标'
Qe = np.matrix([[0,0,0],
                [1,0,0],
                [5,2,0],
                [5,3,0],
                [2,5,0],
                [1,1,0]
            ])
'''
'求质心点'
dataAp = np.matrix([sum(dataA[:,0])/6,sum(dataA[:,1])/6,0])
dataBp = np.matrix([sum(dataB[:,0])/6,sum(dataB[:,1])/6,0])
'质心化后的坐标'
dataA = dataA - dataAp
dataB = dataB - dataBp
 
H = np.matrix([[0,0,0],
               [0,0,0],
               [0,0,0]
             ])

for i in range(0,5):
    H = H + dataA[i,:].T*dataB[i,:]


#print(dataA),print(dataB)
#print(H)
u,sigma,vt=linalg.svd(H)
#print(u),print(vt)
print(vt*u.T)

'''
ax = plt.subplot(111, projection='3d')  #创建一个三维的绘图工程
#  将数据点分成三部分画，在颜色上有区分度
ax.scatter(Pd[:,0], Pd[:,1], 0, c='b')  #绘制前盾数据点
#ax.scatter(dataB[:,0], dataB[:,1], 0, c='r')  #绘制支撑盾数据点

ax.set_zlabel('Z')  #坐标轴
ax.set_ylabel('Y')
ax.set_xlabel('X')

plt.show()

              
'方程定义'
def f(x):
    '''x1,x2,x3,x4,x5,x6,x7,x8,x9,x10,x11,x12,x13,x14,x15 = x.tolist()
    '''
    t1,t2,t3,yaw,roll,pitch = x.tolist()
    '''
    r11 = int(x1)
    r12 = int(x2)
    r13 = int(x3)
    r21 = int(x4)
    r22 = int(x5)
    r23 = int(x6)
    r31 = int(x7)
    r32 = int(x8)
    r33 = int(x9)
    '''
       
    
    R3 = np.matrix([[np.cos(yaw),-np.sin(yaw),0],
            		[np.sin(yaw),np.cos(yaw),0],
            		[0,0,1]])
    R2 = np.matrix([[np.cos(roll),0,np.sin(roll)],
            		[0,1,0],
            		[-np.sin(roll),0,np.cos(roll)]])
    R1 = np.matrix([[1,0,0],
            		[0,np.cos(pitch),-np.sin(pitch)],
            		[0,np.sin(pitch),np.cos(pitch)]])
    Cbn = (R2*R1*R3).I
    #print(Cbn)
    
    r11 = Cbn[0,0]
    r12 = Cbn[0,1]
    r13 = Cbn[0,2]
    r21 = Cbn[1,0]
    r22 = Cbn[1,1]
    r23 = Cbn[1,2]
    r31 = Cbn[2,0]
    r32 = Cbn[2,1]
    r33 = Cbn[2,2]
            
    return [                   
            r11*r11 + r21*r21 + r31*r31 -1,
            r12*r12 + r22*r22 + r32*r32 -1,
            r13*r13 + r23*r23 + r33*r33 -1,
            r11*r12 + r21*r22 + r31*r32,
            r12*r13 + r22*r23 + r32*r33,
            r11*r13 + r21*r23 + r31*r33,
            pow((r11*Qe[0,0]+r12*Qe[0,1]+t1-Pd[0,1]),2)+pow((r21*Qe[0,0]+r22*Qe[0,1]+t2-Pd[0,1]),2)+pow((r31*Qe[0,0]+r32*Qe[0,1]+t3),2)-pow(l1,2),
            pow((r11*Qe[1,0]+r12*Qe[1,1]+t1-Pd[1,1]),2)+pow((r21*Qe[1,0]+r22*Qe[1,1]+t2-Pd[1,1]),2)+pow((r31*Qe[1,0]+r32*Qe[1,1]+t3),2)-pow(l2,2),
            pow((r11*Qe[2,0]+r12*Qe[2,1]+t1-Pd[2,1]),2)+pow((r21*Qe[2,0]+r22*Qe[2,1]+t2-Pd[2,1]),2)+pow((r31*Qe[2,0]+r32*Qe[2,1]+t3),2)-pow(l3,2),
            pow((r11*Qe[3,0]+r12*Qe[3,1]+t1-Pd[3,1]),2)+pow((r21*Qe[3,0]+r22*Qe[3,1]+t2-Pd[3,1]),2)+pow((r31*Qe[3,0]+r32*Qe[3,1]+t3),2)-pow(l4,2),
            pow((r11*Qe[4,0]+r12*Qe[4,1]+t1-Pd[4,1]),2)+pow((r21*Qe[4,0]+r22*Qe[4,1]+t2-Pd[4,1]),2)+pow((r31*Qe[4,0]+r32*Qe[4,1]+t3),2)-pow(l5,2),
            pow((r11*Qe[5,0]+r12*Qe[5,1]+t1-Pd[5,1]),2)+pow((r21*Qe[5,0]+r22*Qe[5,1]+t2-Pd[5,1]),2)+pow((r31*Qe[5,0]+r32*Qe[5,1]+t3),2)-pow(l6,2)
            ]
# f 计算方程组的误差
result = optimize.fsolve(f,[1,6,5,20,20,20])

print(result)
print(f(result))
#print f(result)