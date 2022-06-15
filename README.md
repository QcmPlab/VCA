# Variational Cluster Approximation

An equilibrium Variational Cluster Approximation code.

The code is based on:  
* SciFortran [https://github.com/aamaricci/SciFortran]  
* DMFT_Tools [https://github.com/aamaricci/DMFTtools]

--

### REQUIRES

SciFortran
DMFTtools

### INSTALLATION
```
mkdir build
cd build
cmake .. -DBUILD_TYPE=(TESTING/DEBUG) -DEXE=(FILE IN ../drivers) 
make
```
### HOW IT WORKS

The user has to provide the hopping matrices for the cluster and the lattice (in the cluster basis). Optionally, a bath can be fed to the program (type "normal" of DMFT-ED). VCA_SOLVE returns the value of the Potthoff potential for the input the user provides. Minimization procedure or plots of Ω must be implemented by the user. Periodization must be implemented by the user (examples are provided in drivers).

### NEW BATH IMPLEMENTATION

The bath levels and coupling between them and the interacting cluster are now completely generic. One must pass to the solve function two matrices: the square bath_h(Nlat_bath,Nlat_bath,Nspin,Nspin,Norb_Bath,Norb_bath) and the rectangular bath_v(Nlat,Nlat_bath,Nspin,Nspin,Norb,Norb_bath). This make every bath geometry accessible, not restricting ourselves to diagonal levels/all-to-all hybrid or replicas anymore.

--

***COPYRIGHT & LICENSING***  
Copyright 2012 -  (c) Lorenzo Crippa, Adriano Amaricci.  
All rights reserved. 

The software is provided with no license, as such it is protected by copyright.
The software is provided as it is and can be read and copied, in agreement with 
the Terms of Service of GITHUB. 
The use of the code is constrained to authors agreement.

