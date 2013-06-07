__author__ = 'jeudy'

import numpy as np
import scipy as sp
import matplotlib.pyplot as plt
from matplotlib import animation
from mpl_toolkits.mplot3d import Axes3D
import sys

if len(sys.argv) != 2:
    print "Error. You need to provide a valid file path containing csv data for a cloud."
    exit()

path = sys.argv[1]

raw_data = np.loadtxt(path, delimiter=',', skiprows=1)

x = raw_data[:,0]
y = raw_data[:,1]
z = raw_data[:,2]

# Assumes first 3 positions are X, Y, Z coordinates

fig = plt.figure("N Body")

ax = plt.axes(projection='3d')

ax.scatter(x, y, z)

plt.show()