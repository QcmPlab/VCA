!+------------------------------------------------------------------+
!PURPOSE  : Evaluate Green's functions
!+------------------------------------------------------------------+
subroutine full_build_gf_normal()
  integer :: iorb,ispin
  integer :: isite,jsite
  integer :: Nexc
  !
  call start_timer
  !
  call vca_allocate_time_freq_arrays()
  !
  !Orbital diagonal:

  do ispin=1,Nspin
     Nexc = count_excitations(ispin)
     do iorb=1,Norb
        write(LOGfile,"(A)")"Get impG_l"//str(iorb)//"_s"//str(ispin)
        do isite=1,Nlat
           do jsite=1,Nlat
              !
              call GFmatrix_allocate(impGmatrix(isite,jsite,ispin,ispin,iorb,iorb),N=1)
              call GFmatrix_allocate(impGmatrix(isite,jsite,ispin,ispin,iorb,iorb),i=1,Nexc=Nexc)
              !
              call full_build_gf_normal_diag_c(isite,jsite,iorb,ispin)
              !
           enddo
        enddo
     enddo
  enddo
  !
  call vca_deallocate_time_freq_arrays()
  !
  call stop_timer
end subroutine full_build_gf_normal





!> Count number of 1p-excitations per spin-channel 
function count_excitations(ispin) result(Nexc)
  integer :: ispin
  integer :: isector,jsector
  integer :: idim,jdim
  integer :: i,j,iexc
  real(8) :: expterm
  integer :: Nexc
  iexc=0
  do isector=1,Nsectors
     jsector=getCDGsector(ispin,isector)
     if(jsector==0)cycle
     idim=getdim(isector)     !i-th sector dimension
     jdim=getdim(jsector)     !j-th sector dimension
     do i=1,idim          !loop over the states in the i-th sect.
        do j=1,jdim       !loop over the states in the j-th sect.
           expterm=exp(-beta*espace(isector)%e(i))+exp(-beta*espace(jsector)%e(j))
           if(expterm < cutoff)cycle
           iexc=iexc+1
        enddo
     enddo
  enddo
  !
  Nexc = iexc
  !
  write(LOGfile,"(A,I5,A,I2)")"Found N =",iexc," excitations with the actual cut-off"
  !
end function count_excitations



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
  integer          :: li,rj,iexc
  integer          :: m,i,j,r,k,p
  real(8)          :: Ei,Ej
  real(8)          :: expterm,peso,de,w0
  complex(8)       :: iw
  type(sector_map) :: HI,HJ
  !
  is = imp_state_index(isite,iorb,ispin)
  js = imp_state_index(jsite,iorb,ispin)
  iexc = 0
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
           iexc = iexc+1
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
           impGmatrix(isite,jsite,ispin,ispin,iorb,iorb)%channel(1)%weight(iexc) = spectral_weight
           impGmatrix(isite,jsite,ispin,ispin,iorb,iorb)%channel(1)%poles(iexc)  = de
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
