# -*- coding: utf-8 -*-
"""
Created on Tue Jun 23 14:23:38 2020

@author: Blaine Bateman
"""
#
# 3D visualization and animation
#
#%% libraries
#
import pandas as pd
import numpy as np
import matplotlib as mpl
from mpl_toolkits.mplot3d import Axes3D
mpl.rcParams['figure.dpi'] = 300
import matplotlib.pyplot as plt
import moviepy.video.io.ImageSequenceClip
import os
#
#%% code
#
#%% generate an ellipsoid offset from 0
#
# generate some X-Y-Z data
#
# ellipsoid:
#
# (x - A)^2 + (y - B)^2 + (z - C)^2
#   ---         ---         ---     =  1
#   a^2         b^2         c^2
#
# a longish flattish ellipsoid
#
a = 2
A = -3
b = 6
B = 2
c = 10
C = 7
#
x = np.random.uniform(-a + A, a + A, 1000)
y = np.random.uniform(-b + B, b + B, 1000)
z = np.random.uniform(-c + C, c + C, 1000)
ellipsoid = pd.DataFrame({'x' : x, 'y' : y, 'z' : z})
#
# we keep only the x-y-z that
# are inside the targeted ellipsoid
#
drop_rows = []
for i in range(ellipsoid.shape[0]):
    if ((x[i] - A)**2/a**2 + (y[i] - B)**2/b**2 + (z[i] - C)**2/c**2) >= 1:
        drop_rows.append(i)
ellipsoid.drop(drop_rows, axis = 0, inplace = True)
#
#%% define a color mapping
#
# calculate the radius of every point
#
ellipsoid['radius'] = [(ellipsoid['x'][i] - A)**2/a**2 +
                       (ellipsoid['y'][i] - B)**2/b**2 +
                       (ellipsoid['z'][i] - C)**2/c**2
                       for i in ellipsoid.index]
ellipsoid.sort_values('radius', inplace = True, ascending = False)
ellipsoid.reset_index(drop = True, inplace = True)
radius_range = list(np.arange(0, 1, 0.2))
#
# now we color the points according to radius
# we work backwards for simplicity
# following is more complex than needed here
# so it can work for general variables that might be
# associated with a set of x-y-z data, likea UMAP embedding
#
my_colors = plt.cm.jet(np.linspace(0, 1, len(radius_range)))
cbar_values = pd.Series(radius_range)
radius_range.reverse()
colors = []
color = len(radius_range) - 1
count = 0
for temp in radius_range:
    while (ellipsoid['radius'][count] >= temp):
        colors.append(color)
        count += 1
        if count >= ellipsoid.shape[0]:
            break
    color -= 1
    if count == ellipsoid.shape[0]:
        break
#
ellipsoid['color'] = colors
#
#%% rotate the ellipsoid around y axis
#
psi = 45
#
# following defines the transform for rotation about y
# this is just expanding a 3 x 3 rotation matrix
# note that y does not change since tha axis of rotation
#
ellipsoid['xp'] =  ellipsoid['x'] * np.cos(psi / (2 * np.pi)) + ellipsoid['z'] * np.sin(psi / (2 * np.pi))
ellipsoid['zp'] = -ellipsoid['x'] * np.sin(psi / (2 * np.pi)) + ellipsoid['z'] * np.cos(psi / (2 * np.pi))
ellipsoid['x']  =  ellipsoid['xp']
ellipsoid['z']  =  ellipsoid['zp']
#
#%% visualize
#
# rotation about z axis
# careful!  this can make things really confusing!
#
phi = 5
#
# range of angles to spin around x axis
#
theta_range = range(15, 365, 10)
frames = len(list(theta_range))
#
# get the extreme corner of the space
#
x_max = abs(ellipsoid['x'].max() - A)
y_max = abs(ellipsoid['y'].max() - B)
z_max = abs(ellipsoid['z'].max() - C)
#
# counter to name files and count total frames
#
plot_number = 0
for theta in theta_range:
    plt.clf()
    plt.close()
    plt.style.use('dark_background')
    fig = plt.figure(figsize = (11, 8))
    ax = fig.gca(projection = '3d')
    ax.set_axis_off()
#
# we translate to the origin to keep things simple
#
    for i in range(ellipsoid.shape[0]):
        p = ax.scatter(ellipsoid.loc[i, 'x'] - A,
                       ellipsoid.loc[i, 'y'] - B,
                       ellipsoid.loc[i, 'z'] - C,
                       c = my_colors[ellipsoid.loc[i, 'color']][0:3].reshape(1, -1),
                       alpha = 1,
                       s = 12,
                       antialiased = False)
#
# set rotation of plot
# about z and x axes
#
    ax.view_init(phi, theta)
    ax.set_xlabel(ellipsoid.columns[0], fontsize = 12)
    ax.set_ylabel(ellipsoid.columns[1], fontsize = 12)
    ax.set_zlabel(ellipsoid.columns[2], fontsize = 12)
    ax.text2D(0.1, 0.9,
              'theta: ' + str(theta) + '\nphi: ' + str(phi) + '\npsi: ' + str(psi),
              transform = ax.transAxes,
              fontsize = 14)
#
# construct a color key
# this is a bit of hack but works
# we colro by radius
#
    for color in range(0, len(cbar_values)):
        ax.scatter(x_max, y_max, z_max,
                   c = my_colors[color][0:3].reshape(1, -1),
                   label = str(round(cbar_values[color], 2)),
                   s = 13)
#
# draw a point to over-write the faux points we used
# to create the legend entries
#
    ax.scatter(x_max, y_max, z_max, c = "black", s = 15)
    ax.legend(loc = 'center right')
    plt.savefig('demo_ellipsoid_3d_' + str("{:02d}".format(plot_number)) + '.png')
    plt.show()
    plot_number += 1
#
#%% construct animations
#
frames = plot_number
#
# speed of the playback
#
fps = 10
image_files = [img for img in os.listdir() if img.endswith('.png')]
image_files.sort()
#
# moviepy does all the hard work!
#
my_clip = moviepy.video.io.ImageSequenceClip.ImageSequenceClip(image_files, fps = fps)
#
# depending on your system etc. the codec and file format may
# need to change--followng works on Windows 10 on a laptop
#
my_clip.write_videofile('3d_plot_animation_demo.avi', codec = 'mpeg4')

