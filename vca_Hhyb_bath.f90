  !diagonal, spin conserving:
  htmp=zero
  do ilat=1,Nlat
     do iorb=1,Norb
        do kp=1,Nbath
           ms=getBathStride(ilat,iorb,kp)
           !
           ! IMP UP <--> BATH UP
           is = imp_state_index(ilat,iorb,1)
           Jcondition = (vca_bath%v(ilat,iorb,1,kp)/=0d0).AND.(ib(is)==1).AND.(ib(ms)==0)
           if(Jcondition)then
              call c(iorb,m,k1,sg1)
              call cdg(ms,k1,k2,sg2)
              j = binary_search(H%map,k2)
              htmp = vca_bath%v(ilat,iorb,1,kp)*sg1*sg2
              Hmat(i,j) = Hmat(i,j) + htmp
           endif
           Jcondition = (vca_bath%v(ilat,iorb,1,kp)/=0d0).AND.(ib(is)==0).AND.(ib(ms)==1)
           if(Jcondition)then
              call c(ms,m,k1,sg1)
              call cdg(iorb,k1,k2,sg2)
              j=binary_search(H%map,k2)
              htmp = vca_bath%v(ilat,iorb,1,kp)*sg1*sg2
              Hmat(i,j) = Hmat(i,j) + htmp
           endif
           !
           !IMP DW <--> BATH DW
           is = imp_state_index(ilat,iorb,2)           
           Jcondition = (vca_bath%v(ilat,iorb,Nspin,kp)/=0d0).AND.(ib(is)==1).AND.(ib(ms+Ns)==0)
           if(Jcondition)then
              call c(iorb+Ns,m,k1,sg1)
              call cdg(ms+Ns,k1,k2,sg2)
              j=binary_search(H%map,k2)
              htmp=vca_bath%v(ilat,iorb,Nspin,kp)*sg1*sg2
              Hmat(i,j) = Hmat(i,j) + htmp
           endif
           Jcondition = (vca_bath%v(ilat,iorb,Nspin,kp)/=0d0).AND.(ib(is)==0).AND.(ib(ms+Ns)==1)
           if(Jcondition)then
              call c(ms+Ns,m,k1,sg1)
              call cdg(iorb+Ns,k1,k2,sg2)
              j=binary_search(H%map,k2)
              htmp=vca_bath%v(ilat,iorb,Nspin,kp)*sg1*sg2
              Hmat(i,j) = Hmat(i,j) + htmp
           endif
        enddo
     enddo
  enddo

