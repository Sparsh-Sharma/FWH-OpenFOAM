# FWH-OpenFOAM
FWH post-processing scripts for OpenFOAM solutions

 The formulation I implemented is only for the pressure on solid surfaces though and does not include convection. I hope this also suits your needs.

To use the script, you need to sample the surface pressure during the simulation. I use this function in OpenFOAM:

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

You can simply add this to the end of your controlDict and replace "airfoil" with the name of the noSlip/solidsurface-patch from which you want to calculate the sound radiation. Just make sure that "writeTime" is a feasible sampling frequency.
