  do iup=1,DimUp
     mup = Hs(1)%map(iup)
     ibup = bdecomp(mup,Ns)
     !
     !
     !> H_imp: Off-diagonal elements, i.e. non-local part. 
     !remark: iorb=jorb cant have simultaneously n=0 and n=1 (Jcondition)
     do ilat=1,Nlat
         do jlat=1,Nlat
            do iorb=1,Norb
               do jorb=1,Norb
                  is = imp_state_index(ilat,iorb,1) !CONTROLLA 1 is actually not used
                  js = imp_state_index(jlat,jorb,1)
                  Jcondition = (impHloc(ilat,jlat,1,1,iorb,jorb)/=0d0).AND.(ibup(js)==1).AND.(ibup(is)==0)
               if (Jcondition) then
                  call c(js,mup,k1,sg1)
                  call cdg(is,k1,k2,sg2)
                  jup = binary_search(Hs(1)%map,k2)
                  htmp = impHloc(ilat,jlat,1,1,iorb,jorb)*sg1*sg2
                  !
                  call sp_insert_element(spH0ups(1),htmp,iup,jup)
                  !
               endif
             enddo
          enddo
       enddo
    enddo
  enddo
     !> H_Bath: inter-orbital bath hopping contribution.
     !if(bath_type=="replica") then
     !   do kp=1,Nbath
     !      do iorb=1,Norb
     !         do jorb=1,Norb
     !            !
     !            ialfa = getBathStride(iorb,kp)
     !            ibeta = getBathStride(jorb,kp)
     !            Jcondition = &
     !                 (dmft_bath%h(1,1,iorb,jorb,kp)/=zero) &
     !                 .AND. (Nup(ibeta)==1) .AND. (Nup(ialfa)==0)
     !            !
     !            if (Jcondition)then
     !               call c(ibeta,mup,k1,sg1)
     !               call cdg(ialfa,k1,k2,sg2)
     !               jup = binary_search(Hs(1)%map,k2)
     !               htmp = dmft_bath%h(1,1,iorb,jorb,kp)*sg1*sg2
     !               !
     !               call sp_insert_element(spH0ups(1),htmp,iup,jup)
     !               !
     !            endif
     !         enddo
     !      enddo
     !   enddo
        !
     !\\\cccend if
     !
     !
     !>H_hyb: hopping terms for a given spin (imp <--> bath)
     !do iorb=1,Norb
        !do kp=1,Nbath
           !ialfa=getBathStride(iorb,kp)
           !
           !if( (diag_hybr(1,iorb,kp)/=0d0) &
            !    .AND. (nup(iorb)==1) .AND. (nup(ialfa)==0) )then              
            !  call c(iorb,mup,k1,sg1)
            !  call cdg(ialfa,k1,k2,sg2)
            !  jup = binary_search(Hs(1)%map,k2)
            !  htmp = diag_hybr(1,iorb,kp)*sg1*sg2
              !
            !  call sp_insert_element(spH0ups(1),htmp,iup,jup)
              !
           !endif
           !
           !if( (diag_hybr(1,iorb,kp)/=0d0) &
             !   .AND. (nup(iorb)==0) .AND. (nup(ialfa)==1) )then
              !call c(ialfa,mup,k1,sg1)
              !call cdg(iorb,k1,k2,sg2)
              !jup=binary_search(Hs(1)%map,k2)
              !htmp = diag_hybr(1,iorb,kp)*sg1*sg2
              !
              !call sp_insert_element(spH0ups(1),htmp,iup,jup)
              !
           !endif
        !enddo
     !enddo
     !
  !enddo


