# Variational Cluster Approximation

An equilibrium Variational Cluster Approximation code.

The code is based on:  
* SciFortran [https://github.com/aamaricci/SciFortran]  
* DMFT_Tools [https://github.com/aamaricci/DMFTtools]

--

### TODOs

- [x] save weights&poles in file as gf_cluster_restart.vca
- [ ] make it possible to pass t as a (Nlso,Nlso,Nk) object (useful for
superconducting case, that is not realizable in this branch)
- [ ] calculate the Self-Energy
- [ ] fix flags for the integration routines (actually in Scifor)
- [ ] extend to orbital-resolved case
- [ ] add bath
- [ ] add finite T case

--

***COPYRIGHT & LICENSING***  
Copyright 2012 - 2017 (c), Adriano Amaricci.  
All rights reserved. 

The software is provided with no license, as such it is protected by copyright.
The software is provided as it is and can be read and copied, in agreement with 
the Terms of Service of GITHUB. Use of the code is constrained to author agreement.   

