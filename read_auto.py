#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import numpy as np
import os
import pickle
import csv

def process_openfoam_data(M):
    # Get the directory of the script
    script_dir = os.path.dirname(__file__)
    
    # Set the base path
    path_base = script_dir

    # ------------------------------------------------------
    path = os.path.join(path_base, 'postProcessing/surfaceSampling/')
    path_points = os.path.join(path_base, 'constant/polyMesh/points')
    path_faces = os.path.join(path_base, 'constant/polyMesh/faces')
    path_boundary = os.path.join(path_base, 'constant/polyMesh')
    path_system = os.path.join(path_base, 'system')

    os.chdir(path)

    dir = sorted(os.listdir(path))

    # Automatically generate file_name based on the provided patch
    # Assuming 'AIRFOIL' as the patch, you can modify this based on your requirement
    file_name = 'p_AIRFOIL.raw'
    test = open(os.path.join(path, dir[0], file_name)).read()

    # initialise the geometrical information from the first time step
    a = int(-1)
    with open(os.path.join(path, dir[0], file_name), 'r') as csvfile:
        reader = csv.reader(csvfile, delimiter=' ')
        nxyz = int(next(reader, None)[-1])
        xyz = np.empty([nxyz, 3], dtype=float)
        next(reader, None)
        for row in reader:
            a = a + 1
            xyz[a, 0] = float(row[0])
            xyz[a, 1] = float(row[1])
            xyz[a, 2] = float(row[2])
    print('Finished Initializing geometrical information')

    # read out the pressure information for every time step for every point
    p_raw = np.empty([nxyz, len(dir)], dtype=float)
    b = int(-1)
    for x in range(len(dir)):
        a = int(-1)
        b = b + 1
        with open(os.path.join(path, dir[x], file_name), 'r') as csvfile:
            reader = csv.reader(csvfile, delimiter=' ')
            next(reader, None)
            next(reader, None)
            for row in reader:
                a = a + 1
                p_raw[a, b] = float(row[3])
        print('Finished timestep ' + dir[x])

    print('Finished Reading out the pressure time signals')

    # read all the points from polyMesh/points
    points = open(path_points).readlines()
    xyz_mesh = points[21:-4]
    xyz_anders = np.ones([len(xyz_mesh), 3], dtype=float)
    for bla in range(len(xyz_mesh)):
        xyz_anders[bla, :] = np.float64(xyz_mesh[bla][1:-2].split())

    # Read the 'controlDict' file to extract the patch information
    with open(os.path.join(path_system, 'controlDict'), 'r') as control_file:
        in_surface_sampling = False
        for line in control_file:
            if 'surfaces' in line:
                in_surface_sampling = True
            elif in_surface_sampling and 'patches' in line:
                patch = line.split('(')[1].split(')')[0].strip()
                break

    # Read nFaces and startFace from polyMesh/boundary
    with open(os.path.join(path_boundary, 'boundary'), 'r') as boundary_file:
        in_patch = False
        for line in boundary_file:
            if patch in line:
                in_patch = True
            elif in_patch and 'nFaces' in line:
                nFaces = int(''.join(filter(str.isdigit, line)))
            elif in_patch and 'startFace' in line:
                startFace = int(''.join(filter(str.isdigit, line)))
                in_patch = False

    # Read all the faces from polyMesh/faces
    faces = open(path_faces).readlines()
    faces_mesh = faces[startFace + 20 : startFace + 20 + nFaces]

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

    # --------------------------------------

    p_cxyz = np.empty([len(faces_mesh), p_raw.shape[1]], dtype=float)
    c_xyz = np.empty([len(faces_mesh), 3], dtype=float)
    n = np.empty([len(faces_mesh), 3], dtype=float)
    A = np.empty([len(faces_mesh)], dtype=float)

    xyz = np.round(xyz, decimals=6)
    xyz_anders = np.round(xyz_anders, decimals=6)

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
    f = open('data_geo.pckl', 'wb')
    pickle.dump([A, n, c_xyz, p_cxyz], f, protocol=4)
    f.close()

# Example usage:
M = np.array([0.02, -0.005, 0.00])  # a random point inside the body to determine the direction of the surface normal vector.

process_openfoam_data(M)
