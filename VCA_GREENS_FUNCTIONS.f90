MODULE VCA_GREENS_FUNCTIONS
  USE VCA_VARS_GLOBAL
  USE VCA_IO                     !< this contains the routine to print GF
  USE VCA_DIAG
  USE VCA_SETUP
  USE VCA_AUX_FUNX
  !
  USE SF_CONSTANTS, only:one,xi,zero,pi
  USE SF_TIMER,     only: start_timer,stop_timer,eta
  USE SF_IOTOOLS,   only: str
  USE SF_ARRAYS,    only: arange,linspace
  USE SF_MISC,      only: assert_shape
  !
  implicit none
  private 


  public :: buildGf_cluster



  !Frequency and time arrays:
  !=========================================================
  real(8),dimension(:),allocatable            :: wm,tau,wr,vm




contains



  subroutine buildgf_cluster()
    call allocate_grids()
    !
    impGmats=zero
    impGreal=zero
    !
    write(LOGfile,"(A)")"Get impurity Greens functions:"
    call build_gf_normal()
    !
    if(print_G)call vca_print_impG()
    !
    call deallocate_grids()
  end subroutine buildgf_cluster





  !+------------------------------------------------------------------+
  !PURPOSE  : Evaluate Green's functions
  !+------------------------------------------------------------------+
  subroutine build_gf_normal()
    integer :: ilat,jlat,iorb,jorb,ispin
    !
    call allocate_QLmatrix()
    !
    call start_timer
    !        
    do ispin=1,Nspin
       do ilat=1,Nlat
          do jlat=1,Nlat
             do iorb=1,Norb
                do jorb=1,Norb
                   write(LOGfile,"(A)")"Get G_i"//str(ilat,3)//"_j"//str(jlat,3)//&
                        "_l"//str(iorb)//"_m"//str(jorb)//"_s"//str(ispin)
                   call full_build_gf_QLmatrix(ilat,jlat,iorb,iorb,ispin)
                enddo
             enddo
          enddo
       enddo
    enddo
    !
    call stop_timer
    !
  end subroutine build_gf_normal






  !+------------------------------------------------------------------+
  !PURPOSE  : Spectral sum evaluation of the Green's functions
  ! Matsubara and Real-Axis
  !+------------------------------------------------------------------+
  subroutine full_build_gf_QLmatrix(ilat,jlat,iorb,jorb,ispin)
    integer                            :: ilat,jlat
    integer                            :: iorb,jorb
    integer                            :: ispin,jspin
    integer                            :: isite,jsite
    integer                            :: idim,isector
    integer                            :: jdim,jsector
    real(8)                            :: cdgOp_mat,cOp_mat
    real(8)                            :: op_mat(2),matCC
    real(8)                            :: spectral_weight
    real(8)                            :: sgn_cdg,sgn_c
    integer                            :: ib(Nlevels)
    integer                            :: n,m,p,iexc
    integer                            :: ni,mj
    integer                            :: i,j,r,k
    real(8)                            :: Ei,Ej
    real(8)                            :: expterm,de,w0
    complex(8)                         :: iw
    type(sector_map)                   :: HI,HJ
    !
    isite = state_index(ilat,ispin,iorb)
    jsite = state_index(jlat,ispin,jorb)
    !
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
             op_mat = 0d0
             matCC = 0d0
             expterm=exp(-beta*espace(isector)%e(i))+exp(-beta*espace(jsector)%e(j))
             if(expterm < cutoff)cycle
             !
             iexc=iexc+1
             !
             do ni=1,idim              !loop over the component of |I> (IN state!)
                n  = HI%map(ni)
                ib = bdecomp(n,2*Ns)
                if(ib(isite) == 1)cycle
                call cdg(isite,n,m,sgn_cdg)
                mj = binary_search(HJ%map,m)
                !
                cdgOp_mat = espace(jsector)%M(mj,j)*sgn_cdg*espace(isector)%M(ni,i)
                !
                op_mat(1)=op_mat(1) + cdgOp_mat
             enddo
             !
             do mj=1,jdim              !loop over the component of |J> (IN state!)
                m  = HJ%map(mj)
                ib = bdecomp(m,2*Ns)
                if(ib(jsite) == 0)cycle
                call c(jsite,m,n,sgn_c)
                ni = binary_search(HI%map,n)
                !
                cOp_mat   = espace(isector)%M(ni,i)*sgn_c*espace(jsector)%M(mj,j)
                !
                op_mat(2)=op_mat(2) + cOp_mat
                !
             enddo
             !
             Ei=espace(isector)%e(i)
             Ej=espace(jsector)%e(j)
             de=Ej-Ei
             spectral_weight=expterm/zeta_function*product(op_mat)
             !
             cdgQmatrix(ilat,ispin,iorb,iexc) = op_mat(1)*sqrt(expterm/zeta_function)
             cQmatrix(jlat,ispin,jorb,iexc)   = op_mat(2)*sqrt(expterm/zeta_function)
             Lmatrix(ilat,jlat,ispin,ispin,iorb,jorb,iexc) = de
             !             
             do m=1,Lmats
                iw=xi*wm(m)
                impGmats(ilat,jlat,ispin,ispin,iorb,jorb,m)=impGmats(ilat,jlat,ispin,ispin,iorb,jorb,m)+spectral_weight/(iw-de)
             enddo
             do m=1,Lreal
                iw=dcmplx(wr(m),eps)
                impGreal(ilat,jlat,ispin,ispin,iorb,jorb,m)=impGreal(ilat,jlat,ispin,ispin,iorb,jorb,m)+spectral_weight/(iw-de)
             enddo
             !
          enddo
       enddo
       call delete_sector(isector,HI)
       call delete_sector(jsector,HJ)
    enddo
    !
  end subroutine full_build_gf_QLmatrix









  !+------------------------------------------------------------------+
  !PURPOSE  : Allocate Lambda and Q matrices, storing the spectral
  ! sum terms
  !+------------------------------------------------------------------+
  subroutine allocate_QLmatrix()
    integer :: Nexc
    integer :: isector,jsector
    integer :: idim,jdim
    integer :: i,j
    real(8) :: expterm
    !Evaluate the number of one-body excitations by dry run of the spectral sum:
    Nexc=0
    do isector=1,Nsectors
       jsector=getCDGsector(1,isector)
       if(jsector==0)cycle
       idim=getdim(isector)     !i-th sector dimension
       jdim=getdim(jsector)     !j-th sector dimension
       do i=1,idim          !loop over the states in the i-th sect.
          do j=1,jdim       !loop over the states in the j-th sect.
             expterm=exp(-beta*espace(isector)%e(i))+exp(-beta*espace(jsector)%e(j))
             if(expterm < cutoff)cycle
             Nexc=Nexc+1
          enddo
       enddo
    enddo
    !
    if(allocated(cQmatrix))deallocate(cQmatrix)    
    allocate(cQmatrix(Nlat,Nspin,Norb,Nexc))
    cQmatrix=0d0
    !
    if(allocated(cdgQmatrix))deallocate(cdgQmatrix)
    allocate(cdgQmatrix(Nlat,Nspin,Norb,Nexc))
    cdgQmatrix=0d0
    !
    if(allocated(Lmatrix))deallocate(Lmatrix)
    allocate(Lmatrix(Nlat,Nlat,Nspin,Nspin,Norb,Norb,Nexc))
    Lmatrix=0d0
    !
  end subroutine allocate_QLmatrix



  !+------------------------------------------------------------------+
  !PURPOSE  : Allocate arrays and setup frequencies and times
  !+------------------------------------------------------------------+
  subroutine allocate_grids
    integer :: i
    if(.not.allocated(wm))allocate(wm(Lmats))
    if(.not.allocated(vm))allocate(vm(0:Lmats))          !bosonic frequencies
    if(.not.allocated(wr))allocate(wr(Lreal))
    if(.not.allocated(tau))allocate(tau(0:Ltau))
    wm     = pi/beta*real(2*arange(1,Lmats)-1,8)
    do i=0,Lmats
       vm(i) = pi/beta*2.d0*dble(i)
    enddo
    wr     = linspace(wini,wfin,Lreal)
    tau(0:)= linspace(0.d0,beta,Ltau+1)
  end subroutine allocate_grids


  subroutine deallocate_grids
    if(allocated(wm))deallocate(wm)
    if(allocated(vm))deallocate(vm)
    if(allocated(tau))deallocate(tau)
    if(allocated(wr))deallocate(wr)
  end subroutine deallocate_grids




end MODULE VCA_GREENS_FUNCTIONS





  ! !+------------------------------------------------------------------+
  ! !PURPOSE  : DOUBLE COMPLEX
  ! !+------------------------------------------------------------------+
  ! subroutine full_build_gf_local(ilat,iorb,ispin)
  !   integer                            :: ilat
  !   integer                            :: iorb
  !   integer                            :: ispin
  !   integer                            :: isite
  !   integer                            :: idim,isector
  !   integer                            :: jdim,jsector
  !   real(8)                            :: cdgOp_mat,cOp_mat
  !   real(8)                            :: op_mat(2),matCC
  !   real(8)                            :: spectral_weight
  !   real(8)                            :: sgn_cdg,sgn_c
  !   integer                            :: ib(Nlevels)
  !   integer                            :: n,m,p
  !   integer                            :: ni,mj
  !   integer                            :: i,j,r,k
  !   real(8)                            :: Ei,Ej
  !   real(8)                            :: expterm,de,w0
  !   complex(8)                         :: iw
  !   type(sector_map)                   :: HI,HJ
  !   !
  !   isite = state_index(ilat,ispin,iorb)
  !   !
  !   call start_timer
  !   !
  !   do isector=1,Nsectors
  !      jsector=getCDGsector(ispin,isector)
  !      if(jsector==0)cycle
  !      !
  !      ! call eta(isector,Nsectors)
  !      !
  !      idim=getdim(isector)     !i-th sector dimension
  !      jdim=getdim(jsector)     !j-th sector dimension
  !      call build_sector(isector,HI)
  !      call build_sector(jsector,HJ)
  !      !
  !      do i=1,idim          !loop over the states in the i-th sect.
  !         do j=1,jdim       !loop over the states in the j-th sect.
  !            op_mat = 0d0
  !            expterm=exp(-beta*espace(isector)%e(i))+exp(-beta*espace(jsector)%e(j))
  !            if(expterm < cutoff)cycle
  !            !
  !            do ni=1,idim              !loop over the component of |I> (IN state!)
  !               n  = HI%map(ni)
  !               !
  !               ib = bdecomp(n,2*Ns)
  !               if(ib(isite) == 1)cycle
  !               call cdg(isite,n,m,sgn_cdg)
  !               !
  !               ib = bdecomp(m,2*Ns)
  !               if(ib(isite) == 0)cycle
  !               call c(isite,m,p,sgn_c)
  !               !
  !               mj = binary_search(HJ%map,m)
  !               !
  !               cdgOp_mat = espace(jsector)%M(mj,j)*sgn_cdg*espace(isector)%M(ni,i)
  !               cOp_mat   = espace(isector)%M(ni,i)*sgn_c*espace(jsector)%M(mj,j)
  !               !
  !               op_mat(1)=op_mat(1) + cdgOp_mat
  !               op_mat(2)=op_mat(2) + cOp_mat
  !               !
  !            enddo
  !            !
  !            Ei=espace(isector)%e(i)
  !            Ej=espace(jsector)%e(j)
  !            de=Ej-Ei
  !            spectral_weight=expterm/zeta_function*product(op_mat)
  !            !
  !            do m=1,Lmats
  !               iw=xi*wm(m)
  !               impGmats(ilat,ilat,ispin,ispin,iorb,iorb,m)=impGmats(ilat,ilat,ispin,ispin,iorb,iorb,m)+spectral_weight/(iw+de)
  !            enddo
  !            do m=1,Lreal
  !               iw=dcmplx(wr(m),eps)
  !               impGreal(ilat,ilat,ispin,ispin,iorb,iorb,m)=impGreal(ilat,ilat,ispin,ispin,iorb,iorb,m)+spectral_weight/(iw+de)
  !            enddo
  !            !
  !         enddo
  !      enddo
  !      call delete_sector(isector,HI)
  !      call delete_sector(jsector,HJ)
  !   enddo
  !   !
  !   call stop_progress
  !   !
  ! end subroutine full_build_gf_local
