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

### TODOs

- [x] save weights&poles in file as gf_cluster_restart.vca
- [ ] make it possible to pass t as a (Nlso,Nlso,Nk) object (useful for
superconducting case, that is not realizable in this branch)
- [x] calculate the Self-Energy
- [x] implement G-scheme periodization
- [x] fix flags for the integration routines (actually in Scifor)
- [ ] extend to orbital-resolved case
- [x] add bath
- [x] add sample periodization routines in drivers
- [x] add finite T case

--

***COPYRIGHT & LICENSING***  
Copyright 2012 -  (c) Lorenzo Crippa, Adriano Amaricci.  
All rights reserved. 

The software is provided with no license, as such it is protected by copyright.
The software is provided as it is and can be read and copied, in agreement with 
the Terms of Service of GITHUB. 
The use of the code is constrained to authors agreement.

