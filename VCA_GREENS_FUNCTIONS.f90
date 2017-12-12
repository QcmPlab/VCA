MODULE VCA_GREENS_FUNCTIONS
  USE VCA_VARS_GLOBAL
  USE VCA_DIAG
  USE VCA_SETUP  
  USE VCA_AUX_FUNX
  USE VCA_IO
  USE VCA_BATH_FUNCTIONS
  !
  USE SF_LINALG,    only: inv
  USE SF_CONSTANTS, only: one,xi,zero,pi
  USE SF_TIMER,     only: start_timer,stop_timer,eta
  USE SF_IOTOOLS,   only: str
  !
  implicit none
  private 


  ! !Frequency and time arrays:
  ! !=========================================================
  ! real(8),dimension(:),allocatable            :: wm,wr

  public :: build_gf_cluster

  

  
contains




  !+------------------------------------------------------------------+
  !PURPOSE  : Build the Q and Lambda matrix from Spectral sum.
  !+------------------------------------------------------------------+
  subroutine build_gf_cluster()
    impGmats=zero
    impGreal=zero
    !
    impSmats = zero
    impSreal = zero
    !
    impG0mats=zero
    impG0real=zero
    !
    !
    write(LOGfile,"(A)")"Get impurity Greens functions:"
    call build_gf_normal()
    call get_sigma_normal()
    !
    if(print_Sigma)call vca_print_impSigma()
    if(print_impG)call vca_print_impG()
    if(print_impG0)call vca_print_impG0()
    !
  end subroutine build_gf_cluster





  !+------------------------------------------------------------------+
  !PURPOSE  : Evaluate Green's functions
  !+------------------------------------------------------------------+
  subroutine build_gf_normal()
    integer :: iorb,ispin
    !
    !NORMAL: (default)    
    do ispin=1,Nspin
       do iorb=1,Norb
          write(LOGfile,"(A)")"Get impG_l"//str(iorb)//"_s"//str(ispin)
          call full_build_gf_normal(iorb,ispin)
       enddo
    enddo
    !
  end subroutine build_gf_normal





  !+------------------------------------------------------------------+
  !PURPOSE  : Evaluate Green's functions
  !+------------------------------------------------------------------+
  subroutine full_build_gf_normal(iorb,ispin)
    integer          :: iorb,ispin
    integer :: is,js
    integer          :: isite,jsite
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
    !
    call start_timer
    !
    call vca_allocate_time_freq_arrays()
    !
    do isite=1,Nlat
       do jsite=1,Nlat
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
       enddo
    enddo
    !
    call vca_deallocate_time_freq_arrays()
    !
    call stop_timer
  end subroutine full_build_gf_normal






  !+------------------------------------------------------------------+
  !                    SELF-ENERGY FUNCTIONS 
  !+------------------------------------------------------------------+
  subroutine get_sigma_normal
    integer                                                     :: i,ispin,iorb
    complex(8),dimension(Nlat,Nlat,Nspin,Nspin,Norb,Norb,Lmats) :: invG0mats,invGmats
    complex(8),dimension(Nlat,Nlat,Nspin,Nspin,Norb,Norb,Lreal) :: invG0real,invGreal
    complex(8),dimension(Nlat,Nlat)                             :: invGimp
    !
    call vca_allocate_time_freq_arrays()
    !
    invG0mats = zero
    invGmats  = zero
    invG0real = zero
    invGreal  = zero
    !
    !Get G0^-1
    invG0mats = invg0_bath_mats(dcmplx(0d0,wm(:)),vca_bath)
    invG0real = invg0_bath_real(dcmplx(wr(:),eps),vca_bath)
    !
    !
    !Get Gimp^-1
    invGmats=zero
    invGreal=zero
    do ispin=1,Nspin
       do iorb=1,Norb
          do i=1,Lmats
             invGimp = impGmats(:,:,ispin,ispin,iorb,iorb,i)
             call inv(invGimp)
             invGmats(:,:,ispin,ispin,iorb,iorb,i)=invGimp
          enddo
          !
          do i=1,Lreal
             invGimp = impGreal(:,:,ispin,ispin,iorb,iorb,i)
             call inv(invGimp)
             invGreal(:,:,ispin,ispin,iorb,iorb,i)=invGimp
          enddo
       enddo
    enddo
    !
    !Get Sigma functions: Sigma= G0^-1 - G^-1
    impSmats=zero
    impSreal=zero
    do ispin=1,Nspin
       do iorb=1,Norb
          impSmats(:,:,ispin,ispin,iorb,iorb,:) = invG0mats(:,:,ispin,ispin,iorb,iorb,:) - invGmats(:,:,ispin,ispin,iorb,iorb,:)
          impSreal(:,:,ispin,ispin,iorb,iorb,:) = invG0real(:,:,ispin,ispin,iorb,iorb,:) - invGreal(:,:,ispin,ispin,iorb,iorb,:)
       enddo
    enddo
    !
    !
    !Get G0and:
    impG0mats = g0and_bath_mats(dcmplx(0d0,wm(:)),vca_bath)
    impG0real = g0and_bath_real(dcmplx(wr(:),eps),vca_bath)
    !
    call vca_deallocate_time_freq_arrays()
    !
  end subroutine get_sigma_normal




end MODULE VCA_GREENS_FUNCTIONS













