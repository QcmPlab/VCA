!+------------------------------------------------------------------+
!PURPOSE  : Evaluate Green's functions
!+------------------------------------------------------------------+
subroutine lanc_build_gf_normal()
  integer :: iorb,ispin
  integer :: isite,jsite
  !
  call start_timer
  !
  call vca_allocate_time_freq_arrays()
  !
  !Spin-Orbital diagonal:
  do ispin=1,Nspin
     do iorb=1,Norb
        write(LOGfile,"(A)")"Get impG_l"//str(iorb)//"_s"//str(ispin)
        !
        do isite=1,Nlat
           !site-digonal:
           call GFmatrix_allocate(impGmatrix(isite,isite,ispin,ispin,iorb,iorb),N=2) !2= add,del exc. c^+_i|psi>
           call lanc_build_gf_normal_diag_c(isite,iorb,ispin)
           !site-off-diagonal:
           do jsite=1,Nlat
              if(isite==jsite)cycle   !this is not elegant but who cares?
              call GFmatrix_allocate(impGmatrix(isite,jsite,ispin,ispin,iorb,iorb),N=4)!4=add,del exc. (c^+_i + c^+_j)/(c^+_i +ic^+_j)|psi>
              call lanc_build_gf_normal_offdiag_c(isite,jsite,iorb,ispin)
           enddo
        enddo
        !
        !
        do isite=1,Nlat
           do jsite=1,Nlat
              if(isite==jsite)cycle
              impGmats(isite,jsite,ispin,ispin,iorb,iorb,:) = 0.5d0*(impGmats(isite,jsite,ispin,ispin,iorb,iorb,:) &
                   - (one-xi)*impGmats(isite,isite,ispin,ispin,iorb,iorb,:) - (one-xi)*impGmats(jsite,jsite,ispin,ispin,iorb,iorb,:))
              impGreal(isite,jsite,ispin,ispin,iorb,iorb,:) = 0.5d0*(impGreal(isite,jsite,ispin,ispin,iorb,iorb,:) &
                   - (one-xi)*impGreal(isite,isite,ispin,ispin,iorb,iorb,:) - (one-xi)*impGreal(jsite,jsite,ispin,ispin,iorb,iorb,:))
           enddo
        enddo
        !
     enddo
  enddo
  call vca_deallocate_time_freq_arrays()
  !
  call stop_timer
end subroutine lanc_build_gf_normal




subroutine lanc_build_gf_normal_diag_c(isite,iorb,ispin)
  integer                          :: iorb,ispin,isite
  complex(8),allocatable           :: vvinit(:)
  real(8),allocatable              :: alfa_(:),beta_(:)
  integer                          :: isector,istate
  integer                          :: idim,jsector
  integer                          :: jdim,is
  integer                          :: ib(Nlevels)
  integer                          :: m,i,j,r,numstates
  real(8)                          :: sgn,norm2,norm0
  complex(8)                       :: cnorm2
  integer                          :: Nitermax,Nlanc
  type(sector_map)                 :: HI,HJ
  !
  write(LOGfile,*)"Solving G_cluster_I"//str(isite,3)//"_J"//str(isite,3)
  !
  is = imp_state_index(isite,iorb,ispin)
  !
  do istate=1,state_list%size
     isector    =  es_return_sector(state_list,istate)
     state_e    =  es_return_energy(state_list,istate)
     state_cvec => es_return_cvector(state_list,istate)
     norm0=sqrt(dot_product(state_cvec,state_cvec))
     if(abs(norm0-1.d0)>1.d-9)stop "GS is not normalized"
     idim  = getdim(isector)
     call build_sector(isector,HI)
     !
     !ADD ONE PARTICLE:
     jsector = getCDGsector(ispin,isector)
     if(jsector/=0)then 
        jdim  = getdim(jsector)
        allocate(vvinit(jdim))
        call build_sector(jsector,HJ)
        vvinit=zero
        do m=1,idim
           i=HI%map(m)
           ib = bdecomp(i,2*Ns)
           if(ib(is)==0)then
              call cdg(is,i,r,sgn)
              j=binary_search(HJ%map,r)
              vvinit(j) = sgn*state_cvec(m)
           endif
        enddo
        deallocate(HJ%map)
        norm2=dot_product(vvinit,vvinit)
        vvinit=vvinit/sqrt(norm2)
        !
        call setup_Hv_sector(jsector)
        call buildH_c()
        !
        nlanc=min(jdim,lanc_nGFiter)
        allocate(alfa_(nlanc),beta_(nlanc))
        call sp_lanc_tridiag(spHtimesV_cc,vvinit,alfa_,beta_)
        call add_to_lanczos_gf_normal(one*norm2,state_e,alfa_,beta_,1,isite,isite,iorb,ispin,1)
        !
        call delete_Hv_sector()
        !
        deallocate(vvinit,alfa_,beta_)
        if(spH0%status)call sp_delete_matrix(spH0)
     endif
     !
     !REMOVE ONE PARTICLE:
     jsector = getCsector(ispin,isector)
     if(jsector/=0)then
        jdim  = getdim(jsector)        
        allocate(vvinit(jdim))
        call build_sector(jsector,HJ)
        vvinit=zero
        do m=1,idim
           i=HI%map(m)
           ib = bdecomp(i,2*Ns)
           if(ib(is)==1)then
              call c(is,i,r,sgn)
              j=binary_search(HJ%map,r)
              vvinit(j) = sgn*state_cvec(m)
           endif
        enddo
        deallocate(HJ%map)
        norm2=dot_product(vvinit,vvinit)
        vvinit=vvinit/sqrt(norm2)
        !
        call setup_Hv_sector(jsector)
        call buildH_c()
        !
        nlanc=min(jdim,lanc_nGFiter)
        allocate(alfa_(nlanc),beta_(nlanc))
        call sp_lanc_tridiag(spHtimesV_cc,vvinit,alfa_,beta_)
        call add_to_lanczos_gf_normal(one*norm2,state_e,alfa_,beta_,-1,isite,isite,iorb,ispin,2)
        !
        call delete_Hv_sector()
        !
        deallocate(vvinit,alfa_,beta_)
        if(spH0%status)call sp_delete_matrix(spH0)
     endif
     !
     nullify(state_cvec)
     deallocate(HI%map)
     !
  enddo
end subroutine lanc_build_gf_normal_diag_c





subroutine lanc_build_gf_normal_offdiag_c(isite,jsite,iorb,ispin)
  integer                          :: iorb,jorb,ispin,isite,jsite,isector,istate
  integer                          :: idim,jsector
  integer                          :: jdim,is,js
  integer                          :: ib(Nlevels)
  integer                          :: m,i,j,r,numstates
  real(8)                          :: sgn,norm2,norm0
  complex(8)                       :: cnorm2
  complex(8),allocatable           :: vvinit(:)
  complex(8),allocatable           :: cvinit(:)
  real(8),allocatable              :: alfa_(:),beta_(:)
  integer                          :: Nitermax,Nlanc
  type(sector_map)                 :: HI,HJ
  !
  write(LOGfile,*)"Solving G_cluster_I"//str(isite,3)//"_J"//str(jsite,3)
  !
  is = imp_state_index(isite,iorb,ispin)
  js = imp_state_index(jsite,iorb,ispin)
  !
  do istate=1,state_list%size
     isector    =  es_return_sector(state_list,istate)
     state_e    =  es_return_energy(state_list,istate)
     state_cvec => es_return_cvector(state_list,istate)
     norm0=sqrt(dot_product(state_cvec,state_cvec))
     if(abs(norm0-1.d0)>1.d-9)stop "GS is not normalized"
     !
     idim  = getdim(isector)
     call build_sector(isector,HI)
     !
     !EVALUATE (c^+_is + c^+_js)|gs>
     jsector = getCDGsector(ispin,isector)
     if(jsector/=0)then 
        jdim  = getdim(jsector)
        allocate(vvinit(jdim))
        call build_sector(jsector,HJ)
        vvinit=zero
        do m=1,idim
           i=HI%map(m)
           ib = bdecomp(i,2*Ns)
           if(ib(is)==0)then
              call cdg(is,i,r,sgn)
              j=binary_search(HJ%map,r)
              vvinit(j) = sgn*state_cvec(m)
           endif
        enddo
        do m=1,idim
           i=HI%map(m)
           ib = bdecomp(i,2*Ns)
           if(ib(js)==0)then
              call cdg(js,i,r,sgn)
              j=binary_search(HJ%map,r)
              vvinit(j) = vvinit(j) + sgn*state_cvec(m)
           endif
        enddo
        deallocate(HJ%map)
        norm2=dot_product(vvinit,vvinit)
        vvinit=vvinit/sqrt(norm2)
        !
        call setup_Hv_sector(jsector)
        call buildH_c()
        !
        nlanc=min(jdim,lanc_nGFiter)
        allocate(alfa_(nlanc),beta_(nlanc))
        call sp_lanc_tridiag(spHtimesV_cc,vvinit,alfa_,beta_)
        call add_to_lanczos_gf_normal(one*norm2,state_e,alfa_,beta_,1,isite,jsite,iorb,ispin,1)
        !
        call delete_Hv_sector()
        !
        deallocate(vvinit,alfa_,beta_)
        if(spH0%status)call sp_delete_matrix(spH0)
     endif
     !
     !EVALUATE (c_i + c_j)|gs>
     jsector = getCsector(ispin,isector)
     if(jsector/=0)then
        jdim   = getdim(jsector)
        allocate(vvinit(jdim))
        call build_sector(jsector,HJ)
        vvinit=zero
        do m=1,idim
           i=HI%map(m)
           ib = bdecomp(i,2*Ns)
           if(ib(is)==1)then
              call c(is,i,r,sgn)
              j=binary_search(HJ%map,r)
              vvinit(j) = sgn*state_cvec(m)
           endif
        enddo
        do m=1,idim
           i=HI%map(m)
           ib = bdecomp(i,2*Ns)
           if(ib(js)==1)then
              call c(js,i,r,sgn)
              j=binary_search(HJ%map,r)
              vvinit(j) = vvinit(j) + sgn*state_cvec(m)
           endif
        enddo
        deallocate(HJ%map)
        norm2=dot_product(vvinit,vvinit)
        vvinit=vvinit/sqrt(norm2)
        !
        call setup_Hv_sector(jsector)
        call buildH_c()
        !
        nlanc=min(jdim,lanc_nGFiter)
        allocate(alfa_(nlanc),beta_(nlanc))
        call sp_lanc_tridiag(spHtimesV_cc,vvinit,alfa_,beta_)
        call add_to_lanczos_gf_normal(one*norm2,state_e,alfa_,beta_,-1,isite,jsite,iorb,ispin,2)
        !
        call delete_Hv_sector()
        !
        deallocate(vvinit,alfa_,beta_)
        if(spH0%status)call sp_delete_matrix(spH0)
     endif
     !
     !EVALUATE (c^+_iorb + i*c^+_jorb)|gs>
     jsector = getCDGsector(ispin,isector)
     if(jsector/=0)then 
        jdim  = getdim(jsector)
        allocate(cvinit(jdim))
        call build_sector(jsector,HJ)
        cvinit=zero
        do m=1,idim
           i=HI%map(m)
           ib = bdecomp(i,2*Ns)
           if(ib(is)==0)then
              call cdg(is,i,r,sgn)
              j=binary_search(HJ%map,r)
              cvinit(j) = sgn*state_cvec(m)
           endif
        enddo
        do m=1,idim
           i=HI%map(m)
           ib = bdecomp(i,2*Ns)
           if(ib(js)==0)then
              call cdg(js,i,r,sgn)
              j=binary_search(HJ%map,r)
              cvinit(j) = cvinit(j) + xi*sgn*state_cvec(m)
           endif
        enddo
        deallocate(HJ%map)
        norm2=dot_product(cvinit,cvinit)
        cvinit=cvinit/sqrt(norm2)
        !
        call setup_Hv_sector(jsector)
        call buildH_c()
        !
        nlanc=min(jdim,lanc_nGFiter)
        allocate(alfa_(nlanc),beta_(nlanc))
        call sp_lanc_tridiag(spHtimesV_cc,cvinit,alfa_,beta_)
        call add_to_lanczos_gf_normal(-xi*norm2,state_e,alfa_,beta_,1,isite,jsite,iorb,ispin,3)
        !
        call delete_Hv_sector()
        !
        deallocate(cvinit,alfa_,beta_)
        if(spH0%status)call sp_delete_matrix(spH0)
     endif
     !
     !EVALUATE (c_iorb - xi*c_jorb)|gs>
     jsector = getCsector(ispin,isector)
     if(jsector/=0)then
        jdim   = getdim(jsector)
        allocate(cvinit(jdim))
        call build_sector(jsector,HJ)
        cvinit=zero
        do m=1,idim
           i=HI%map(m)
           ib = bdecomp(i,2*Ns)
           if(ib(is)==1)then
              call c(is,i,r,sgn)
              j=binary_search(HJ%map,r)
              cvinit(j) = sgn*state_cvec(m)
           endif
        enddo
        do m=1,idim
           i=HI%map(m)
           ib = bdecomp(i,2*Ns)
           if(ib(js)==1)then
              call c(js,i,r,sgn)
              j=binary_search(HJ%map,r)
              cvinit(j) = cvinit(j) - xi*sgn*state_cvec(m)
           endif
        enddo
        deallocate(HJ%map)
        norm2=dot_product(cvinit,cvinit)
        cvinit=cvinit/sqrt(norm2)
        !
        call setup_Hv_sector(jsector)
        call buildH_c()
        !
        nlanc=min(jdim,lanc_nGFiter)
        allocate(alfa_(nlanc),beta_(nlanc))
        call sp_lanc_tridiag(spHtimesV_cc,cvinit,alfa_,beta_)        
        call add_to_lanczos_gf_normal(-xi*norm2,state_e,alfa_,beta_,-1,isite,jsite,iorb,ispin,4)
        !
        call delete_Hv_sector()
        !
        deallocate(cvinit,alfa_,beta_)
        if(spH0%status)call sp_delete_matrix(spH0)
     endif
     !
     nullify(state_cvec)
     deallocate(HI%map)
     !
  enddo
  !
end subroutine lanc_build_gf_normal_offdiag_c


















!+------------------------------------------------------------------+
!PURPOSE  : 
!+------------------------------------------------------------------+
subroutine add_to_lanczos_gf_normal(vnorm2,Ei,alanc,blanc,isign,ilat,jlat,iorb,ispin,ichan)
  integer                                    :: ilat,jlat
  complex(8)                                 :: vnorm2,pesoBZ,peso
  real(8)                                    :: Ei,Egs,de
  integer                                    :: nlanc,itype
  real(8),dimension(:)                       :: alanc
  real(8),dimension(size(alanc))             :: blanc 
  integer                                    :: isign,iorb,jorb,ispin,ichan
  real(8),dimension(size(alanc),size(alanc)) :: Z
  real(8),dimension(size(alanc))             :: diag,subdiag
  integer                                    :: i,j,ierr
  complex(8)                                 :: iw
  !
  Egs = state_list%emin       !get the gs energy
  !
  Nlanc = size(alanc)
  !
  ! if((finiteT).and.(beta*(Ei-Egs).lt.200))then
  !    pesoBZ = vnorm2*exp(-beta*(Ei-Egs))/zeta_function
  ! elseif(.not.finiteT)then
  !    pesoBZ = vnorm2/zeta_function
  ! else
  !    pesoBZ=0.d0
  ! endif
  !
  pesoBZ = vnorm2/zeta_function
  if(finiteT)pesoBZ = vnorm2*exp(-beta*(Ei-Egs))/zeta_function
  !
  diag             = 0.d0
  subdiag          = 0.d0
  Z                = eye(Nlanc)
  diag(1:Nlanc)    = alanc(1:Nlanc)
  subdiag(2:Nlanc) = blanc(2:Nlanc)
  call tql2(Nlanc,diag,subdiag,Z,ierr)
  !
  call GFmatrix_allocate(impGmatrix(ilat,jlat,ispin,ispin,iorb,iorb),i=ichan,Nexc=Nlanc)
  !
  do j=1,nlanc
     de = diag(j)-Ei
     peso = pesoBZ*Z(1,j)*Z(1,j)
     !
     impGmatrix(ilat,jlat,ispin,ispin,iorb,iorb)%channel(ichan)%weight(j) = peso
     impGmatrix(ilat,jlat,ispin,ispin,iorb,iorb)%channel(ichan)%poles(j)  = isign*de
     !
     do i=1,Lmats
        iw=xi*wm(i)
        impGmats(ilat,jlat,ispin,ispin,iorb,iorb,i)=impGmats(ilat,jlat,ispin,ispin,iorb,iorb,i) + peso/(iw-isign*de)
     enddo
     do i=1,Lreal
        iw=dcmplx(wr(i),eps)
        impGreal(ilat,jlat,ispin,ispin,iorb,iorb,i)=impGreal(ilat,jlat,ispin,ispin,iorb,iorb,i) + peso/(iw-isign*de)
     enddo
  enddo
end subroutine add_to_lanczos_gf_normal











