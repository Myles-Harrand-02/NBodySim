import matplotlib.pyplot as plt
import pandas as pd

from matplotlib import cm
from mpl_toolkits.mplot3d import axes3d
import mpl_toolkits.mplot3d.art3d as art3d

import matplotlib.animation as anim

import numpy as np

sim = pd.read_csv("D:\\LOLayto\\Programming\\N Body\\n_body_1\\output.csv", sep=",")

sun = sim[ sim["Body"] == "Sun" ]

earth = sim[ sim["Body"] == "Earth" ]
moon = sim[ sim["Body"] == "Moon" ]

#jupiter = sim[ sim["Body"] == "Jupiter" ]
#saturn = sim[ sim["Body"] == "Saturn" ]

#uranus = sim[ sim["Body"] == "Uranus" ]
#neptune = sim[ sim["Body"] == "Neptune" ]

#plot 3D
ax = plt.figure().add_subplot(projection='3d')

ax.plot([0,0],[0,0],[0,0],'ro')

ax.plot3D(earth['Disp (x)'], earth['Disp (y)'], earth['Disp (z)'])

ax.set_aspect('equal', adjustable='box')
plt.show()

#plot 2D from above
""" fig, ax = plt.subplots(figsize=(10.2,14))

ax.plot([0,0],[0,0],'ro')

ax.plot(earth['Disp (x)'],earth['Disp (y)'],'b')
#ax.plot(moon['Disp (x)'],moon['Disp (y)'],'orange')

#ax.plot(jupiter['Disp (x)'],jupiter['Disp (y)'],'y')
#ax.plot(saturn['Disp (x)'],saturn['Disp (y)'],'g')

#ax.plot(uranus['Disp (x)'],uranus['Disp (y)'],'teal')
#ax.plot(neptune['Disp (x)'],neptune['Disp (y)'],'purple')

ax.set_aspect('equal', adjustable='box')
plt.show() """

#plot 3D animation

fig = plt.figure()
ax = fig.add_subplot(projection='3d')

earthX = moon['Disp (x)'].tolist()
earthY = moon['Disp (y)'].tolist()
earthZ = moon['Disp (z)'].tolist()

ax.plot([0,0],[0,0],[0,0],'ro')

earthPos, = ax.plot(earthX[0], earthY[0], earthZ[0], color = 'blue', linewidth = 2)
earthPosTran, = ax.plot(earthX[0], earthY[0], earthZ[0], color = 'blue', linewidth = 2, alpha = 0.25)

n1, n2 = 150,300

def update(frame):
    if frame - n1 <= 0:
        xDat = earthX[:frame]
        yDat = earthY[:frame]
        zDat = earthZ[:frame]

        xDatTran = earthX[0:1]
        yDatTran = earthY[0:1]
        zDatTran = earthZ[0:1]
    else:
        xDat = earthX[frame - n1:frame]
        yDat = earthY[frame - n1:frame]
        zDat = earthZ[frame - n1:frame]

        xDatTran = earthX[:frame - n1]
        yDatTran = earthY[:frame - n1]
        zDatTran = earthZ[:frame - n1]

    earthPos.set_data_3d(np.asarray(xDat), np.asarray(yDat), np.asarray(zDat))
    earthPosTran.set_data_3d(np.asarray(xDatTran), np.asarray(yDatTran), np.asarray(zDatTran))

circ = plt.Circle((0, 0), 149.6e9, color='g', fill=False)
ax.add_patch(circ)
art3d.pathpatch_2d_to_3d(circ, z=0, zdir="z")

ani = anim.FuncAnimation(fig=fig, func=update, frames=len(earthX), interval=600)

ax.axes.set_xlim3d(left=-150e9, right=150e9) 
ax.axes.set_ylim3d(bottom=-150e9, top=150e9) 
ax.axes.set_zlim3d(bottom=-150e9, top=150e9)

ax.set_aspect('equal', adjustable='box')

#writervideo = anim.FFMpegWriter(fps=60) 
#ani.save('bodies.mp4', writer=writervideo)
          
plt.show()

""" fig, ax = plt.subplots(figsize=(10.2,14))

ax.plot(earth['Time'],earth['Disp (x)'])

plt.show()

fig, ax = plt.subplots(figsize=(10.2,14))

ax.plot(earth['Time'],earth['E'])

plt.show()

fig, ax = plt.subplots(figsize=(10.2,14))

ax.plot(earth['Time'],earth['L'])

plt.show()

fig, ax = plt.subplots(figsize=(10.2,14))

ax.plot(earth['Time'],earth['L Total'])

plt.show()

fig, ax = plt.subplots(figsize=(10.2,14))

ax.plot(earth['Time'],earth['Bound'])

plt.show()

fig, ax = plt.subplots(figsize=(10.2,14))

ax.plot(earth['Time'],earth['Ecce'])

plt.show() """