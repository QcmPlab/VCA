  !CLUSTER  HAMILTONIAN
  !Diagonal Elements, i.e. local part
  htmp = zero
  htmp = htmp - xmu*(sum(nup)+sum(ndw))
  !
  do ilat=1,Nlat
     do iorb=1,Norb
        htmp = htmp + impHloc(ilat,ilat,1,1,iorb,iorb)*nup(ilat,iorb)
        htmp = htmp + impHloc(ilat,ilat,Nspin,Nspin,iorb,iorb)*ndw(ilat,iorb)
     enddo
  enddo
  !
  Hmat(i,i)=Hmat(i,i) + htmp
  !
  !
  !Off-diagonal elements, i.e. non-local part
  !this loop considers only the site-orbitals off-diagonal terms
  !because ilat/iorb=jlat/jorb can not have simultaneously
  !occupation 0 and 1, as required by this if Jcondition:  
  do ilat=1,Nlat
     do jlat=1,Nlat
        do iorb=1,Norb
           do jorb=1,Norb
              !
              !UP
              is = imp_state_index(ilat,iorb,1)
              js = imp_state_index(jlat,jorb,1)
              Jcondition = (impHloc(ilat,jlat,1,1,iorb,jorb)/=0d0).AND.(ib(js)==1).AND.(ib(is)==0)
              if (Jcondition) then
                 call c(js,m,k1,sg1)
                 call cdg(is,k1,k2,sg2)
                 j = binary_search(H%map,k2)
                 htmp = impHloc(ilat,jlat,1,1,iorb,jorb)*sg1*sg2
                 !
                 Hmat(i,j) = Hmat(i,j) + htmp
                 !
              endif
              !
              !DW
              is = imp_state_index(ilat,iorb,2)
              js = imp_state_index(jlat,jorb,2)
              Jcondition = (impHloc(ilat,jlat,Nspin,Nspin,iorb,jorb)/=0d0).AND.(ib(js)==1).AND.(ib(is)==0)
              if (Jcondition) then
                 call c(js,m,k1,sg1)
                 call cdg(is,k1,k2,sg2)
                 j = binary_search(H%map,k2)
                 htmp = impHloc(ilat,jlat,Nspin,Nspin,iorb,jorb)*sg1*sg2
                 !
                 Hmat(i,j) = Hmat(i,j) + htmp
                 !
              endif
           enddo
        enddo
     enddo
  enddo
