# hfr

performs quantum chemical calculation based on Restricted Hartree-Fock-Roothaan equation

# Prerequirement

compiled with NVIDIA HPC SDK 23.7 on AlmaLinux 9.2 on amd64

# Specifications

## electron integrals

Electron integrals ( e.g. overlap integral, electron repulsion integral ) are calculated using Obara-Saika methods.

## parallel computation

Electron repulsion integrals are parallely calculated using OpenMPI.

# Example

Example input files are stored in inputs directory.

File format is specified in inputs/sample.dat

Summary of calculation result for h2o.dat is the below.

 - MolecularEnergy= -0.73214037E+02
 - ElectronEnergy= -0.90588536E+02
 - OrbitalEnergies
    - -0.2046E+02
    - -0.1688E+01
    - -0.9450E+00
    - -0.5844E+00
    - -0.5472E+00
    - 0.9919E+00
    - 0.1501E+01
