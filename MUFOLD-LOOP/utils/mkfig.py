#!/usr/bin/python
import glob2
lpfiles = glob2.glob('./*.crd')

import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
fig = plt.figure()
ax = fig.add_subplot(111, projection='3d')

import numpy as np
crds = np.loadtxt('lpint0.crd')
plt.show()

