  !
  !diagonal bath hamiltonian: +energy of the bath=\sum_a=1,Norb\sum_{l=1,Nbath}\e^a_l n^a_l
  htmp=zero
  do ilat=1,size(vca_bath%e,1)    !Nlat for bath_type=normal, 1 for bath_type=hybrid
     do iorb=1,size(vca_bath%e,2) !Norb for bath_type=normal, 1 for bath_type=hybrid
        do kp=1,Nbath
           alfa=getBathStride(ilat,iorb,kp)
           htmp =htmp + vca_bath%e(ilat,iorb,1,kp)*ib(alfa)        !UP
           htmp =htmp + vca_bath%e(ilat,iorb,Nspin,kp)*ib(alfa+Ns) !DW
        enddo
     enddo
  enddo
  !
  Hmat(i,i) = Hmat(i,i) + htmp
  !
