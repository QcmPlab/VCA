  do idw=1,MpiQup
     do iup=1,DimDw
        mdw  = Hs(2)%map(iup)
        ibdw  = bdecomp(mdw,Ns)
        !
        i    = iup + (idw-1)*DimDw
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
                    Jcondition = &
                         (impHloc(ilat,jlat,Nspin,Nspin,iorb,jorb)/=zero) .AND. &
                         (ibdw(js)==1) .AND. (ibdw(is)==0)
                    if (Jcondition) then
                       call c(js,mdw,k1,sg1)
                       call cdg(is,k1,k2,sg2)
                       jup = binary_search(Hs(2)%map,k2)
                       jdw = idw             
                       j   = jup + (jdw-1)*DimDw
                       htmp = impHloc(ilat,jlat,Nspin,Nspin,iorb,jorb)*sg1*sg2
                       !
                       Hvt(i) = Hvt(i) + htmp*vt(j)
                       !
                    endif
                 enddo
              enddo
            enddo
          enddo
     enddo
  enddo

