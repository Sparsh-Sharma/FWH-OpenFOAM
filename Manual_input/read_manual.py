#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Dec 20 12:10:42 2023

@author: shar_sp
"""
# -*- coding: utf-8 -*-
"""
Created on Sat Nov 21 14:44:46 2020

@author: WS1
"""

import numpy as np
from pylab import *
import os
#import matplotlib.pyplot as plt
import pickle
#from mpl_toolkits.mplot3d import Axes3D
import csv


#-------------------- user input
#put here these two values for your sampled patch from
#polyMesh/boundary: 
nFaces=65130
startFace=58012835
# nFaces=792868
# startFace=88896228

M=np.array([0.02,-0.005,0.00]) #random point inside of the body, to determine the direction of the surface normal vector
path_base='/scratch/ws/m0/shar_sp-LES/LEN_CNode'
file_name='p_AIRFOIL.raw' # name of the surface sampling data for every time step
patch= 'AIRFOIL' # name of the patch where surface pressure was sampled




#------------------------------------------------------
path=path_base+'/postProcessingM/surfaceSampling/'
path_points = path_base+'/constant/polyMesh/points' #path  polyMesh/points
path_faces=path_base+'/constant/polyMesh/faces' #path polyMesh/faces
path_bound=path_base+'/constant/polyMesh/boundary'#path polyMesh/vboundary
os.chdir(path)

dir = sort(os.listdir(path)) 

test = open(path + dir[0] + '/' +file_name).read()

## initialise the geometrical information from the first time step
a=int(-1)
with open(path + dir[0] + '/' + file_name, 'r') as csvfile:
    reader = csv.reader(csvfile, delimiter=' ')
    nxyz=int(next(reader, None)[-1])
    xyz= np.empty([nxyz,3], dtype=float) # skip the header 2 times because two rows for header and read number of points from first row
    next(reader, None)
    for row in reader:   #this will read each row line by line
        a=a+1
        xyz[a,0] = float(row[0])
        xyz[a,1] = float(row[1])
        xyz[a,2] = float(row[2])
print('Finished Initializing geometrical information')

# read out the pressure information for every time step for every point
p_raw = np.empty([nxyz,len(dir)], dtype=float)
b=int(-1)   
for x in range(len(dir)):  #every time step 
    a=int(-1)
    b=b+1
    with open(path + dir[x] + '/' + file_name, 'r') as csvfile:
        reader = csv.reader(csvfile, delimiter=' ')
        next(reader, None) # skip the header, 2 times because two rows for header
        next(reader, None)
        for row in reader:   #one time step at a time for all points 
            a=a+1
            p_raw[a,b] = float(row[3])
    print('Finished timestep ' +dir[x])

print('Finished Reading out the pressure time signals')

#read all the points from polyMesh/points            
points = open(path_points).readlines()
xyz_mesh = points[21:-4] #delete first 21 rows for header and the last 4
xyz_anders = np.ones([len(xyz_mesh),3], dtype=float)
for bla in range(len(xyz_mesh)):
   xyz_anders[bla,:] = np.float64(xyz_mesh[bla][1:-2].split()) #read points from strings
close(path_points)

#read all the faces from polyMesh/faces 
faces = open(path_faces).readlines()
faces_mesh = faces[startFace+20:startFace+20+nFaces] # read faces from startface until number of faces 
close(path_faces)

#-----------------------------------------
#area of polygon poly --> https://stackoverflow.com/questions/12642256/python-find-area-of-polygon-from-xyz-coordinates
def poly_area(poly,n):
    if len(poly) < 3: # not a plane - no area
        return 0
    total = [0, 0, 0]
    N = len(poly)
    for i in range(N):
        vi1 = poly[i]
        vi2 = poly[(i+1) % N]
        prod = np.cross(vi1, vi2)
        total[0] += prod[0]
        total[1] += prod[1]
        total[2] += prod[2]
#    result = np.dot(total, unit_normal(poly[0], poly[1], poly[2]))
    result = np.dot(total, n)
    return abs(result/2)
#--------------------------------------

#initialise variables
p_cxyz=np.empty([len(faces_mesh), p_raw.shape[1]], dtype=float)
c_xyz=np.empty([len(faces_mesh), 3], dtype=float)
n=np.empty([len(faces_mesh), 3], dtype=float)
A=np.empty([len(faces_mesh)], dtype=float)

xyz=np.round(xyz, decimals=6)
xyz_anders=np.round(xyz_anders, decimals=6)


#for every face: calculate centroid c_xyz, normal vector n, Area A, pressure spatial mean value p_cxyz of all points from face
for face in range(len(faces_mesh)):

    c_xyz[face,:]=mean(xyz_anders[np.int32(faces_mesh[face][2:-2].split())], axis=0) #centroid is arithmetic mean value of points
    
    #calcualte normal vector with cross product, this can be done with any three points from the list, so just take the first 3
    X1 = xyz_anders[np.int32(faces_mesh[face][2:-2].split())[1]] - xyz_anders[np.int32(faces_mesh[face][2:-2].split())[0]] #with point 1
    X2 = xyz_anders[np.int32(faces_mesh[face][2:-2].split())[2]] - xyz_anders[np.int32(faces_mesh[face][2:-2].split())[0]] #with point 2
    n[face,:] = np.cross(X1, X2) #normal vector
    
    #test if normal vector is pointing in the right direction:
    if np.dot(c_xyz[face,:]+n[face,:], M) >0:
        n[face,:]=-n[face,:] #if not, turn the vector
        
    n[face,:] = n[face,:]/np.linalg.norm(n[face,:]) #normalize to length 1
    
    #calculate area of n-point polygon with function from above
    A[face] = poly_area(xyz_anders[np.int32(faces_mesh[face][2:-2].split())], n[face,:])
    
    #look up which points belong to face and calculate the mean value for pressure
    # Look up which points belong to face and calculate the mean value for pressure
    # Look up which points belong to face and calculate the mean value for pressure
    po = np.empty([len(xyz_anders[np.int32(faces_mesh[face][2:-2].split())]), p_raw.shape[1]], dtype=float)
    for hut in range(len(xyz_anders[np.int32(faces_mesh[face][2:-2].split())])):
        # Use try-except block to handle cases where the where condition returns an empty array
        indices = np.where(
            (xyz[:, 0] == xyz_anders[np.int32(faces_mesh[face][2:-2].split())][hut][0]) &
            (xyz[:, 1] == xyz_anders[np.int32(faces_mesh[face][2:-2].split())][hut][1]) &
            (xyz[:, 2] == xyz_anders[np.int32(faces_mesh[face][2:-2].split())][hut][2])
        )[0]  # Access the first element of the array
    
        if len(indices) > 0:
            po[hut, :] = p_raw[indices, :]
        else:
            # Set the first row to zeros when no matching points are found
            po[hut, :] = 0.0
            print(f"Warning: No matching points found for face {face}, hut {hut}. Setting pressure to zeros.")

    # po = np.empty([len(xyz_anders[np.int32(faces_mesh[face][2:-2].split())]), p_raw.shape[1]], dtype=float)
    # for hut in range(len(xyz_anders[np.int32(faces_mesh[face][2:-2].split())])):
    #     po[hut,:] = p_raw[np.where( (xyz[:,0]==xyz_anders[np.int32(faces_mesh[face][2:-2].split())][hut][0]) & (xyz[:,1]==xyz_anders[np.int32(faces_mesh[face][2:-2].split())][hut][1]) & (xyz[:,2]==xyz_anders[np.int32(faces_mesh[face][2:-2].split())][hut][2]) ),:]
    p_cxyz[face,:]=mean(po, axis=0).T
    print(face)  
    

os.chdir(path_base)
f = open('data_geo.pckl', 'wb')
pickle.dump([A, n, c_xyz, p_cxyz], f, protocol=4)
f.close()

