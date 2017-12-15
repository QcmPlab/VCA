  !LOCAL INTERACTION
  !density-density interaction: same orbital, opposite spins:
  ! = \sum_\a U_\a*(n_{\a,up}*n_{\a,dw})
  htmp = zero
  do ilat=1,Nlat
     do iorb=1,Norb
        htmp = htmp + Uloc(iorb)*nup(ilat,iorb)*ndw(ilat,iorb)
     enddo
  enddo
  if(Norb>1)then
     !density-density interaction: different orbitals, opposite spins:
     ! =   U'   *     sum_{i/=j} [ n_{i,up}*n_{j,dw} + n_{j,up}*n_{i,dw} ]
     ! =  (Uloc-2*Jh)*sum_{i/=j} [ n_{i,up}*n_{j,dw} + n_{j,up}*n_{i,dw} ]
     do ilat=1,Nlat
        do iorb=1,Norb
           do jorb=iorb+1,Norb
              htmp = htmp + Ust*(nup(ilat,iorb)*ndw(ilat,jorb) + nup(ilat,jorb)*ndw(ilat,iorb))
           enddo
        enddo
     enddo
     !density-density interaction: different orbitals, parallel spins
     ! = \sum_{i<j}    U''     *[ n_{i,up}*n_{j,up} + n_{i,dw}*n_{j,dw} ]
     ! = \sum_{i<j} (Uloc-3*Jh)*[ n_{i,up}*n_{j,up} + n_{i,dw}*n_{j,dw} ]
     do ilat=1,Nlat
        do iorb=1,Norb
           do jorb=iorb+1,Norb
              htmp = htmp + (Ust-Jh)*(nup(ilat,iorb)*nup(ilat,jorb) + ndw(ilat,iorb)*ndw(ilat,jorb))
           enddo
        enddo
     enddo
  endif
  !if using the Hartree-shifted chemical potential: mu=0 for half-filling
  !sum up the contributions of hartree terms:
  if(hfmode)then
     do ilat=1,Nlat
        do iorb=1,Norb
           htmp = htmp - 0.5d0*Uloc(iorb)*(nup(ilat,iorb)+ndw(ilat,iorb)) + 0.25d0*Uloc(iorb)
        enddo
     enddo
     if(Norb>1)then
        do ilat=1,Nlat
           do iorb=1,Norb
              do jorb=iorb+1,Norb
                 htmp=htmp-0.5d0*Ust*(nup(ilat,iorb)+ndw(ilat,iorb)+nup(ilat,jorb)+ndw(ilat,jorb))+0.25d0*Ust
                 htmp=htmp-0.5d0*(Ust-Jh)*(nup(ilat,iorb)+ndw(ilat,iorb)+nup(ilat,jorb)+ndw(ilat,jorb))+0.25d0*(Ust-Jh)
              enddo
           enddo
        enddo
     endif
  endif
  !  
  ! Hmat(i,i) = Hmat(i,i) + htmp
  call sp_insert_element(spH0,htmp,i,i)
  !
  !
  !
  ! SPIN-EXCHANGE (S-E) and PAIR-HOPPING TERMS
  !    S-E: J c^+_iorb_up c^+_jorb_dw c_iorb_dw c_jorb_up  (i.ne.j) 
  !    S-E: J c^+_{iorb} c^+_{jorb+Ns} c_{iorb+Ns} c_{jorb}
  if(Norb>1.AND.Jhflag)then
     do ilat=1,Nlat
        do iorb=1,Norb
           do jorb=1,Norb
              i_up = imp_state_index(ilat,iorb,1)
              i_dw = imp_state_index(ilat,iorb,2)
              j_up = imp_state_index(ilat,jorb,1)
              j_dw = imp_state_index(ilat,jorb,2)
              Jcondition=(&
                   (iorb/=jorb).AND.&
                   (ib(j_up)==1).AND.&
                   (ib(i_dw)==1).AND.&
                   (ib(j_dw)==0).AND.&
                   (ib(i_up)==0))
              if(Jcondition)then
                 call c(j_up,m,k1,sg1)
                 call c(i_dw,k1,k2,sg2)
                 call cdg(j_dw,k2,k3,sg3)
                 call cdg(i_up,k3,k4,sg4)
                 j=binary_search(H%map,k4)
                 htmp = one*Jx*sg1*sg2*sg3*sg4
                 !
                 ! if(j/=0)Hmat(i,j) = Hmat(i,j) + htmp
                 if(j/=0)call sp_insert_element(spH0,htmp,i,j)
                 !
              endif
           enddo
        enddo
     enddo
  endif
  !
  ! PAIR-HOPPING (P-H) TERMS
  !    P-H: J c^+_iorb_up c^+_iorb_dw   c_jorb_dw   c_jorb_up  (i.ne.j) 
  !    P-H: J c^+_{iorb}  c^+_{iorb+Ns} c_{jorb+Ns} c_{jorb}
  if(Norb>1.AND.Jhflag)then
     do ilat=1,Nlat
        do iorb=1,Norb
           do jorb=1,Norb
              i_up = imp_state_index(ilat,iorb,1)
              i_dw = imp_state_index(ilat,iorb,2)
              j_up = imp_state_index(ilat,jorb,1)
              j_dw = imp_state_index(ilat,jorb,2)
              Jcondition=(&
                   (iorb/=jorb).AND.&
                   (ib(j_up)==1).AND.&
                   (ib(j_dw)==1).AND.&
                   (ib(i_dw)==0).AND.&
                   (ib(i_up)==0))
              if(Jcondition)then
                 call c(j_up,m,k1,sg1)
                 call c(j_dw,k1,k2,sg2)
                 call cdg(i_dw,k2,k3,sg3)
                 call cdg(i_up,k3,k4,sg4)
                 j=binary_search(H%map,k4)
                 htmp = one*Jp*sg1*sg2*sg3*sg4
                 !
                 ! if(j/=0)Hmat(i,j) = Hmat(i,j) + htmp
                 if(j/=0)call sp_insert_element(spH0,htmp,i,j)
                 !
              endif
           enddo
        enddo
     enddo
  endif
