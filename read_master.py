#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Author: shar_sp
Created on Wed Dec 20 12:10:42 2023
"""

import os
import numpy as np
import pickle
import csv

# User input
nFaces = 65130
startFace = 58012835
M = np.array([0.02, -0.005, 0.00])
path_base = '/scratch/ws/m0/shar_sp-LES/LEN_CNode'
file_name = 'p_AIRFOIL.raw'
patch = 'AIRFOIL'

# Paths
path = f'{path_base}/postProcessingM/surfaceSampling/'
path_points = f'{path_base}/constant/polyMesh/points'
path_faces = f'{path_base}/constant/polyMesh/faces'
path_bound = f'{path_base}/constant/polyMesh/boundary'

os.chdir(path)

# Read directories
dir = sorted(os.listdir(path))

# Read the first time step
test = open(f'{path}{dir[0]}/{file_name}').read()

# Initialize geometrical information
a = -1
with open(f'{path}{dir[0]}/{file_name}', 'r') as csvfile:
    reader = csv.reader(csvfile, delimiter=' ')
    nxyz = int(next(reader, None)[-1])
    xyz = np.empty([nxyz, 3], dtype=float)
    next(reader, None)
    for row in reader:
        a += 1
        xyz[a, 0] = float(row[0])
        xyz[a, 1] = float(row[1])
        xyz[a, 2] = float(row[2])
print('Finished Initializing geometrical information')

# Read pressure information for every time step
p_raw = np.empty([nxyz, len(dir)], dtype=float)
b = -1
for x in range(len(dir)):
    a = -1
    b += 1
    with open(f'{path}{dir[x]}/{file_name}', 'r') as csvfile:
        reader = csv.reader(csvfile, delimiter=' ')
        next(reader, None)
        next(reader, None)
        for row in reader:
            a += 1
            p_raw[a, b] = float(row[3])
    print(f'Finished timestep {dir[x]}')

print('Finished Reading out the pressure time signals')

# Read points from polyMesh/points
points = open(path_points).readlines()
xyz_mesh = points[21:-4]
xyz_anders = np.ones([len(xyz_mesh), 3], dtype=float)
for bla in range(len(xyz_mesh)):
   xyz_anders[bla, :] = np.float64(xyz_mesh[bla][1:-2].split())
close(path_points)

# Read faces from polyMesh/faces
faces = open(path_faces).readlines()
faces_mesh = faces[startFace + 20:startFace + 20 + nFaces]
close(path_faces)

# Function to calculate polygon area
def poly_area(poly, n):
    if len(poly) < 3:
        return 0
    total = [0, 0, 0]
    N = len(poly)
    for i in range(N):
        vi1 = poly[i]
        vi2 = poly[(i + 1) % N]
        prod = np.cross(vi1, vi2)
        total[0] += prod[0]
        total[1] += prod[1]
        total[2] += prod[2]
    result = np.dot(total, n)
    return abs(result / 2)

# Initialize variables
p_cxyz = np.empty([len(faces_mesh), p_raw.shape[1]], dtype=float)
c_xyz = np.empty([len(faces_mesh), 3], dtype=float)
n = np.empty([len(faces_mesh), 3], dtype=float)
A = np.empty([len(faces_mesh)], dtype=float)

xyz = np.round(xyz, decimals=6)
xyz_anders = np.round(xyz_anders, decimals=6)

# For every face: calculate centroid c_xyz, normal vector n, Area A, pressure spatial mean value p_cxyz of all points from face
for face in range(len(faces_mesh)):
    c_xyz[face, :] = np.mean(xyz_anders[np.int32(faces_mesh[face][2:-2].split())], axis=0)
    X1 = xyz_anders[np.int32(faces_mesh[face][2:-2].split())[1]] - xyz_anders[np.int32(faces_mesh[face][2:-2].split())[0]]
    X2 = xyz_anders[np.int32(faces_mesh[face][2:-2].split())[2]] - xyz_anders[np.int32(faces_mesh[face][2:-2].split())[0]]
    n[face, :] = np.cross(X1, X2)
    if np.dot(c_xyz[face, :] + n[face, :], M) > 0:
        n[face, :] = -n[face, :]
    n[face, :] = n[face, :] / np.linalg.norm(n[face, :])
    A[face] = poly_area(xyz_anders[np.int32(faces_mesh[face][2:-2].split())], n[face, :])
    po = np.empty([len(xyz_anders[np.int32(faces_mesh[face][2:-2].split())]), p_raw.shape[1]], dtype=float)
    for hut in range(len(xyz_anders[np.int32(faces_mesh[face][2:-2].split())])):
        indices = np.where(
            (xyz[:, 0] == xyz_anders[np.int32(faces_mesh[face][2:-2].split())][hut][0]) &
            (xyz[:, 1] == xyz_anders[np.int32(faces_mesh[face][2:-2].split())][hut][1]) &
            (xyz[:, 2] == xyz_anders[np.int32(faces_mesh[face][2:-2].split())][hut][2])
        )[0]
        if len(indices) > 0:
            po[hut, :] = p_raw[indices, :]
        else:
            po[hut, :] = 0.0
            print(f"Warning: No matching points found for face {face}, hut {hut}. Setting pressure to zeros.")
    p_cxyz[face, :] = np.mean(po, axis=0).T
    print(face)

os.chdir(path_base)
with open('data_geo.pckl', 'wb') as f:
    pickle.dump([A, n, c_xyz, p_cxyz], f, protocol=4)
