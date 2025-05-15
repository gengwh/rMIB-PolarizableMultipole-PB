# MIB-PMPB
This is the repo for the code for solving the Polarizable Multipole Poisson-Boltzmann (PMPB) model using regularized Matched Interface and Boundary (rMIB) method.


## References:
1. X. Yang, S. Zhao, W. Geng, A regularized Matched Interface and Boundary Method (rMIB) for Solving Polarizable Multipole Poisson-Boltzmann model, in preparation. 

2.  W. Geng and S. Zhao, A two-component Matched Interface and Boundary (MIB) regularization for charge singularity in implicit solvation, J. Comput. Phys., 351, 25-39 (2017).


## Usage:
1. Install the triangular surface generators MSMS (https://ccsb.scripps.edu/msms/downloads/) and/or ESES (https://github.com/WeilabMSU/ESES); Make sure the binary executive file of one or both were added into the system path (e.g. modify .bashrc file for Linux and bash_profile or .zshrc file for MacOS)

2. Clone the code and run Makefile to build executive file rMIB.

3. Modify the parameters in usrdata.in file as needed.

==Note==: dcel is the grid size; den is the density for MSMS; eps0 and eps1 are dielectric constants for solute and solvent.

==Note==: icg,isf,ibd gives the different handling of charge, interface and boundary.

icg:  0: no charge;          2: green function. 

isf:  0: exact interface;    1: MSMS interface;   4: ESES interface.  (Note: When isf=4, pqr file has to be width-controled; not space seperated)

ibd:  0: exact soln      1: kirkwood soln;        2: 1/r decay      3:exp decay. (Note: when isf=1 or 4, set ibd=2 for zero ionic strength and ibd=3 for nonzero inoic strength)

inl:  0:linear PB;           1: nonlinear PB. 

ipm:  0:point charge		1: point multipole (monopole+dipole+quadrupole)

## Point Multipole 
<!-- 1. Install Tinker: ./source/pdbxyz.x 1pgb.pdb -key 1pgb.key  
1pgb.key file : parameter amoebapro13.prm  
2. python readData.py 1pgb.xyz -->


1. Users will need to download the Tinker software package first.
This can be done using command line:
```
git clone https://github.com/TinkerTools/Tinker
```
Then, a Brookhaven PDB file downloaded from Protein Data Bank (https://www.rcsb.org/) can be converted into a Tinker **.xyz** Cartesian coordinate file using the program **pdbxyz** in Tinker. This conversion can be done under **source** subdirectory under Tinker directory using the command line:
```
pdbxyz.x PDBID.pdb -key PDBID.key
```
Here, the key file is a single line file, i.e., ``parameters amoebapro13", which can redirect to the force field file **amoebapro13.prm**. More parameter files can be found under **params** subdirectory under Tinker directory.

2. Next, a **.xyz** file needs to be converted to a dummy **.pqr** file which contains charge positions, point charges, dipole moments, and quadrupole moments information, plus a **.xyzr** file which contains the charge position and radius information. The **.xyzr** file will also be used when calling **MSMS** molecular surface software. 
The **.pqr** file can be generated using the updated python script under __src/test_proteins__ under MIB directory by the command line:
```
python readData.py PDBID.xyz
```
The **readData.py** is originally from the software PYGBE (https://github.com/cdcooper84/pygbe).
The **.xyzr** file can be generated using the python script under __src/test_proteins__ under MIB directory by the command line:
```
python tinker\_to\_xyzr.py PDBID
```
This file is also from software PYGBE (https://github.com/cdcooper84/pygbe).
