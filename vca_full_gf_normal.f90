!+------------------------------------------------------------------+
!PURPOSE  : Evaluate Green's functions
!+------------------------------------------------------------------+
subroutine full_build_gf_normal(iorb,ispin)
  integer :: iorb,ispin
  integer :: isite,jsite
  !
  call start_timer
  !
  call vca_allocate_time_freq_arrays()
  !
  do isite=1,Nlat
     do jsite=1,Nlat
        !
        call full_build_gf_normal_diag_c(isite,jsite,iorb,ispin)
        !
     enddo
  enddo
  !
  call vca_deallocate_time_freq_arrays()
  !
  call stop_timer
end subroutine full_build_gf_normal




subroutine full_build_gf_normal_diag_c(isite,jsite,iorb,ispin)
  integer          :: isite,jsite
  integer          :: iorb,ispin
  integer          :: is,js
  integer          :: idim,isector
  integer          :: jdim,jsector
  complex(8)       :: op_mat(2)
  complex(8)       :: spectral_weight
  real(8)          :: sgn_cdg,sgn_c
  integer          :: ib(Nlevels)
  integer          :: li,rj
  integer          :: m,i,j,r,k,p
  real(8)          :: Ei,Ej
  real(8)          :: expterm,peso,de,w0
  complex(8)       :: iw
  type(sector_map) :: HI,HJ

  is = imp_state_index(isite,iorb,ispin)
  js = imp_state_index(jsite,iorb,ispin)
  do isector=1,Nsectors
     jsector=getCDGsector(ispin,isector)
     if(jsector==0)cycle
     !
     idim=getdim(isector)     !i-th sector dimension
     jdim=getdim(jsector)     !j-th sector dimension
     call build_sector(isector,HI)
     call build_sector(jsector,HJ)
     !
     do i=1,idim          !loop over the states in the i-th sect.
        do j=1,jdim       !loop over the states in the j-th sect.
           op_mat=0.d0
           expterm=exp(-beta*espace(isector)%e(i))+exp(-beta*espace(jsector)%e(j))
           if(expterm < cutoff)cycle
           !
           do li=1,idim              !loop over the component of |I> (IN state!)
              m = HI%map(li)
              ib = bdecomp(m,2*Ns)
              if(ib(isite) == 1)cycle
              call cdg(is,m,k,sgn_cdg)
              rj = binary_search(HJ%map,k)
              op_mat(1)=op_mat(1) + conjg(espace(jsector)%M(rj,j))*sgn_cdg*espace(isector)%M(li,i)
           enddo
           !
           do rj=1,jdim
              m = HJ%map(rj)
              ib = bdecomp(m,2*Ns)
              if(ib(jsite) == 0)cycle
              call c(js,m,k,sgn_c)
              li = binary_search(HI%map,k)
              op_mat(2)=op_mat(2) + conjg(espace(isector)%M(li,i))*sgn_c*espace(jsector)%M(rj,j)
           enddo
           !
           Ei=espace(isector)%e(i)
           Ej=espace(jsector)%e(j)
           de=Ej-Ei
           peso=expterm/zeta_function
           spectral_weight=peso*product(op_mat)
           !
           do m=1,Lmats
              iw=xi*wm(m)
              impGmats(isite,jsite,ispin,ispin,iorb,iorb,m)=impGmats(isite,jsite,ispin,ispin,iorb,iorb,m)+&
                   spectral_weight/(iw-de)
           enddo
           !
           do m=1,Lreal
              w0=wr(m);iw=cmplx(w0,eps)
              impGreal(isite,jsite,ispin,ispin,iorb,iorb,m)=impGreal(isite,jsite,ispin,ispin,iorb,iorb,m)+&
                   spectral_weight/(iw-de)
           enddo
           !
        enddo
     enddo
     deallocate(HI%map,HJ%map)
  enddo
end subroutine full_build_gf_normal_diag_c
