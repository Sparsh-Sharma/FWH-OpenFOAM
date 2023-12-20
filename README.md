# FWH-OpenFOAM

![Python Logo](https://www.python.org/static/community_logos/python-logo.png)

[![License](https://img.shields.io/badge/License-MIT-blue.svg)](LICENSE)
[![GitHub issues](https://img.shields.io/github/issues/your-username/your-repo)](https://github.com/your-username/your-repo/issues)
[![GitHub stars](https://img.shields.io/github/stars/your-username/your-repo)](https://github.com/your-username/your-repo/stargazers)

## Description

Briefly describe your project here. Highlight its purpose, features, and potential benefits.

## Table of Contents

- [Installation](#installation)
- [Usage](#usage)
- [Contributing](#contributing)
- [License](#license)

## Installation

Provide step-by-step instructions on how to install and run your project. Include any dependencies and prerequisites needed.

```bash
pip install your-package

# FWH-OpenFOAM
FWH post-processing scripts for OpenFOAM solutions

 The formulation I implemented is only for the pressure on solid surfaces though and does not include convection. 
 
To use the script, you need to sample the surface pressure during the simulation. Use this function in OpenFOAM:


```
functions
{

  surfaceSampling
    {
        type surfaces;

        // Where to load it from (if not already in solver)
        libs            ("libsampling.so");
        writeControl    writeTime;

        interpolationScheme cellPoint;
        setFormat ascii;
        surfaceFormat raw;

        // Fields to be sampled
        fields
        (
            p
        );

        surfaces
        (
                airfoil
            {
                type            patchInternalField;
                patches         ( airfoil );
                distance 0;
                interpolate     true;
                triangulate     false;
            }
        );
    }
}

```

You can simply add this to the end of your controlDict and replace "airfoil" with the name of the noSlip/solidsurface-patch from which you want to calculate the sound radiation. Just make sure that "writeTime" is a feasible sampling frequency.

The "readData_Mesh_OpenFOAM"-script reads all the geometrical information needed from the mesh and the sampled pressure data and saves them. "acousticPython6O" calculates the acoustic pressure from the saved file.
