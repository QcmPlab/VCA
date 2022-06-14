  do jdw=1,MpiQdw
     do jup=1,DimUp
        mup = Hs(1)%map(jup)
        ibup = bdecomp(mup,Ns)
        j    = jup + (jdw-1)*dimUp
        !
        !
        !> H_imp: Off-diagonal elements, i.e. non-local part. 
        !remark: iorb=jorb cant have simultaneously n=0 and n=1 (Jcondition)
        do ilat=1,Nlat
          do jlat=1,Nlat
            do iorb=1,Norb
               do jorb=1,Norb
                  is = imp_state_index(ilat,iorb)
                  js = imp_state_index(jlat,jorb)
                  Jcondition = (impHloc(ilat,jlat,1,1,iorb,jorb)/=0d0).AND.(ibup(js)==1).AND.(ibup(is)==0)
                  if (Jcondition) then
                      call c(js,mup,k1,sg1)
                      call cdg(is,k1,k2,sg2)
                      iup = binary_search(Hs(1)%map,k2)
                      htmp = impHloc(ilat,jlat,1,1,iorb,jorb)*sg1*sg2
                      i   = iup + (jdw-1)*DimUp
                      !
                      Hv(i) = Hv(i) + htmp*vin(j)
                      !
                  endif
                enddo
              enddo
            enddo
          enddo
        
  !> H_Bath: inter-orbital bath hopping contribution.
  if(Nlat_bath>0 .and. Norb_bath>0)then
    do ilat=1,Nlat_bath
      do jlat=1,Nlat_bath
        do iorb=1,Norb_bath
            do jorb=1,Norb_bath
               !
               ialfa = getBathStride(ilat,iorb)
               ibeta = getBathStride(jlat,jorb)
               Jcondition = &
                    (ialfa/=ibeta .AND. vca_bath%h(ilat,jlat,1,1,iorb,jorb)/=zero) &
                    .AND. (ibup(ibeta)==1) .AND. (ibup(ialfa)==0)
               !
               if (Jcondition)then
                  call c(ibeta,mup,k1,sg1)
                  call cdg(ialfa,k1,k2,sg2)
                  iup = binary_search(Hs(1)%map,k2)
                  htmp = vca_bath%h(ilat,jlat,1,1,iorb,jorb)*sg1*sg2
                  !
                  hv(i) = hv(i) + htmp*vin(j)
                  !
               endif
            enddo
         enddo
       enddo
     enddo
     !
     !>H_hyb: hopping terms for a given spin (imp <--> bath)
     do ilat=1,Nlat
       do iorb=1,Norb
          do jlat=1,Nlat_bath
            do jorb=1,Norb_bath
              ialfa=getBathStride(jlat,jorb) !bath site
              is = imp_state_index(ilat,iorb) !imp site
              if( (diag_hybr(ilat,jlat,1,iorb,jorb)/=0d0) &
                   .AND. (ibup(ialfa)==1) .AND. (ibup(is)==0) )then              
                 call c(ialfa,mup,k1,sg1)
                 call cdg(is,k1,k2,sg2)
                 iup = binary_search(Hs(1)%map,k2)
                 htmp = diag_hybr(ilat,jlat,1,iorb,jorb)*sg1*sg2
                 !
                 hv(i) = hv(i) + htmp*vin(j)
                 !
              endif
              if( (diag_hybr(ilat,jlat,1,iorb,jorb)/=0d0) &
                   .AND. (ibup(ialfa)==0) .AND. (ibup(is)==1) )then
                 call c(is,mup,k1,sg1)
                 call cdg(ialfa,k1,k2,sg2)
                 iup=binary_search(Hs(1)%map,k2)
                 htmp = conjg(diag_hybr(ilat,jlat,1,iorb,jorb))*sg1*sg2
                 !
                 hv(i) = hv(i) + htmp*vin(j)
                 !
             endif
            enddo
          enddo
       enddo
     enddo
  endif 
        
        
     enddo
  enddo







