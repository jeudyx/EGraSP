#!/usr/bin/python

# -*- coding: utf-8 -*-

import numpy as np
import scipy as sp
import matplotlib.pyplot as plt
from matplotlib import animation
from mpl_toolkits.mplot3d import Axes3D
import sys
import glob
import os

if len(sys.argv) != 2:
    print "Error. You need to provide a valid folder containing csv data for a cloud."
    exit()

path = sys.argv[1]

files = sorted(glob.glob("%s/*.csv" % (path)), key=os.path.getmtime)

data_holder = []

for path in files:
    raw_data = np.loadtxt(path, delimiter=',', skiprows=1)
    data_holder.append(raw_data)

data_holder = np.array(data_holder)


fig = plt.figure("Data animation")
ax = plt.axes(projection='3d')

first_row = data_holder[0]

x = first_row[:,0]
y = first_row[:,1]
z = first_row[:,2]

sp, = ax.plot(x, y, z, 'r.')

def update(i):
    i_row = data_holder[i]
    x = i_row[:,0]
    y = i_row[:,1]
    z = i_row[:,2]
    sp.set_data(x, y)
    sp.set_3d_properties(z)
    return sp,

ani = animation.FuncAnimation(fig, update, frames=len(files), interval=50, repeat=True)

plt.show()