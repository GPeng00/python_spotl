#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Nov  1 10:17:29 2023

@author: amt
"""

import spotl
import numpy as np
import os

leng=10
lat=39
dep=10
lon=-123
mypoly='poly.hawaii'
    
bsig = np.zeros((3, 3, leng))
lsig = np.zeros((3, 3, leng))
loadTides = np.zeros((8,))  # Assuming loadTides is an array of length 8
bodyTide = np.zeros((8,))  # Assuming bodyTide is an array of length 8

# Make a Station object (implement makeStation function)
station = spotl.make_station('station', lat, lon, 0)

# Define darwinSymbols
darwinSymbols = ['k1'] #, 'k2', 'm2', 'n2', 'o1', 'p1', 'q1', 's2']

# Implement readPolygonFile
polywc = spotl.read_polygon_file('../polys/'+mypoly)

# Implement readGreensFunctionFile
gr = '../green/green.bou.'+str(dep)  # Adjust file path
gr = '../green/green.gbaver.wef.p01.ce'
file=gr
green=spotl.read_greens_function_file(file)

#tideModelWorld = spotl.readTideModelFile(f'../tidmod/m2.osu.tpxo72.2010.txt')
tide_file='../tidmod/m2.osu.tpxo72.2010.txt'

for iDarwin in range(len(darwinSymbols) - 1, -1, -1):
    # Implement readTideModelFile
    tideModelWorld = spotl.read_tide_model_file(f'../tidmod/{darwinSymbols[iDarwin]}.osu.tpxo72.2010.txt')
    tideModel = spotl.read_tide_model_file(f'../tidmod/{darwinSymbols[iDarwin]}.osu.hawaii.2010.txt')

    loadTides[iDarwin * 2] = spotl.nloadf(station, tideModelWorld, green, 'l' , 'polygonExclude', polywc)
    # loadTides[iDarwin * 2 - 1] = nloadf(station, tideModel, green, 'l', 'polygonInclude', polywc)

'''
# Implement combineTideStruct
loadTides = combineTideStruct(loadTides)

# Load tides for a specific date
# Implement hartid function to obtain time and load_strain
time, load_strain = hartid(loadTides, 'strain', leng, dt, 'ymdhms', yr, month, day, hr, minute, sec)

v = lame['nu']  # Assuming lame is a dictionary with keys like 'nu', 'G', 'E', 'K', 'lambda'

# Calculate body tides
# Implement computeEarthTide
bodyTide = computeEarthTide(loadTides)

# Calculate body strain and adjust it
time, body_strain = hartid(bodyTide, 'strain', leng, dt, 'ymdhms', yr, month, day, hr, minute, sec)
body_strain = np.hstack((body_strain, -v / (1 - v) * (body_strain[:, 0] + body_strain[:, 1])))

G = lame['G']
E = lame['E']
K = lame['K']
lam = lame['lambda']

for i in range(leng):
    e = np.array([[body_strain[i, 0], body_strain[i, 2], 0],
                  [body_strain[i, 2], body_strain[i, 1], 0],
                  [0, 0, body_strain[i, 3]]]) * 1e-9
    bsig[:, :, i] = lam * np.trace(e) * np.eye(3) + 2 * G * e

for i in range(leng):
    e = np.array([[load_strain[i, 0], load_strain[i, 2], load_strain[i, 3]],
                  [load_strain[i, 2], load_strain[i, 1], load_strain[i, 4]],
                  [load_strain[i, 3], load_strain[i, 4], load_strain[i, 5]]]) * 1e-9
    lsig[:, :, i] = lam * np.trace(e) * np.eye(3) + 2 * G * e
    '''