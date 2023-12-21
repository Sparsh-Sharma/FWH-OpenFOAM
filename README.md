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




## Usage for directivity.py

1. Clone the repository:

```bash
git clone https://github.com/Sparsh-Sharma/FWH-OpenFOAM.git
cd FWH-OpenFOAM
```

2. Run the script:

```bash
python directivity.py
```

3. View the generated polar plots in the `plots/` directory.

## Parameters

- `fs`: Sampling frequency
- `c`: Speed of sound
- `rho0`: Density of air
- `dt`: Time step
- `freq_of_interest`: Frequencies for SPL calculation
- `x0_values`: Observer positions (x0 values)
- `num_angles`: Number of angles to consider
- `theta_values`: Angles in degrees

## Files

- `directivity_script.py`: Main script for calculating SPL directivity.
- `data_geo.pckl`: Pickled file containing geometry data.
- `plots/`: Directory containing the generated polar plots.

## License

This project is licensed under the MIT License - see the [LICENSE](LICENSE) file for details.

## Acknowledgments

- This script was developed by Sparsh Sharma.

Feel free to customize this template according to your project's specific details. Add more sections as needed, such as installation instructions, troubleshooting tips, or any additional acknowledgments.
