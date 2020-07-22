# Computational Analysis of Dynamic Allostery and Control in the SARS-CoV-2 Main Protease
Processed data and figure generation codes for a project on an Elastic Network Model analysis of the SARS-CoV-2 main protease (M<sup>pro</sup>, aka 3CL<sup>pro</sup>).
See [preprint](https://www.biorxiv.org/content/10.1101/2020.05.21.105965v1) in bioRxiv.

## PDB
We used [6LU7](https://www.rcsb.org/structure/6LU7) X-ray crystallographic structure of the main protease to construct apo (ligand-free), holo1 (one ligand present) and holo2 protein forms (two ligands present, native form).
Note: all water molecules (`H2O`) were removed from all protein structure PDB files.

## Software
All ENM's data was generated and processed via [Durham Dynamic Protein Toolbox](https://sourceforge.net/projects/durham-ddpt/), 'DDPT' for short, an open-source FORTRAN based protein analysis toolbox. 

## Visualization
Python scripts were used for the processed data visualization. 

## Data
This work presents four types of data:
1. Wild-type ENM Gibbs free energies 
2. Cross-correlation of motion 
3. 1-point mutational scans.
4. 2-point mutational scans.