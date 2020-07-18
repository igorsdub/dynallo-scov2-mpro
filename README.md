# Computational Analysis of Dynamic Allostery and Control in the SARS-CoV-2 Main Protease
Processed data and figures genration code for a project on an Elastic Network Model analysis of the SARS-CoV-2 Main Protease (M<sup>pro<\sup>, aka 3CLsup>pro<\sup>).
See a [prepeint](https://www.biorxiv.org/content/10.1101/2020.05.21.105965v1) at bioRxiv.

## PDB
We used wild-type [6LU7](https://www.rcsb.org/structure/6LU7) X-ray crystallographic structure to construct apo (ligand-free), holo1 (one ligand present) and holo2 protein forms (two ligands present, native form).
Note: all water molecues (`H2O`) were removed from all protein stucture PDB files.

## Software
ENMs were generated and processed in [Durham Dynamic Protein Toolbox](https://sourceforge.net/projects/durham-ddpt/), 'DDPT' for short, an open-source FORTRAN based protein analysis toolbox.

## Data
1. Cross-correlation of motion 
2. 1-point mutaional scan: relative change of fluctuation Gibbs free energy; ligand dissocaition contant ratio - in response ENM modifications.
3. 2-point mutational scan: relative change of fluctuation Gibbs free energy; ligand dissocaition contant ratio - in response ENM modifications.
 

## Code
1. Cross-correlation and distance map
2. 1-point mutational scan map (apo and allo)
3. 2-point mutational scan map (apo and allo)
