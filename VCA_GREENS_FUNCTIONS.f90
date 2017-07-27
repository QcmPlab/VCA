MODULE VCA_GREENS_FUNCTIONS
  USE VCA_VARS_GLOBAL
  USE VCA_IO                     !< this contains the routine to print GF
  USE VCA_DIAG
  USE VCA_SETUP
  USE VCA_AUX_FUNX
  !
  USE SF_CONSTANTS, only: one,xi,zero,pi
  USE SF_TIMER,     only: start_timer,stop_timer
  USE SF_IOTOOLS,   only: str
  !
  implicit none
  private 


  public :: build_QL_cluster


contains



  subroutine build_QL_cluster()
    integer :: ilat,jlat,iorb,jorb,ispin
    !
    call allocate_QLmatrix()
    !
    call start_timer
    !        

    do ispin=1,Nspin
       call full_build_Lmatrix(ispin)
       do ilat=1,Nlat
          do jlat=1,Nlat
             do iorb=1,Norb
                do jorb=1,Norb
                   !
                   write(LOGfile,"(A)")"Get Q_i"//str(ilat,3)//"_j"//str(jlat,3)//"_l"//str(iorb)//"_m"//str(jorb)//"_s"//str(ispin)
                   !
                   call full_build_Qmatrix(ilat,jlat,iorb,iorb,ispin)
                   !
                enddo
             enddo
          enddo
       enddo
    enddo
    !
    call stop_timer
    !
  end subroutine build_QL_cluster








  !+------------------------------------------------------------------+
  !PURPOSE  : Spectral sum evaluation of the Green's functions
  ! Matsubara and Real-Axis
  !+------------------------------------------------------------------+
  subroutine full_build_Qmatrix(ilat,jlat,iorb,jorb,ispin)
    integer                            :: ilat,jlat
    integer                            :: iorb,jorb
    integer                            :: ispin,jspin
    integer                            :: isite,jsite
    integer                            :: idim,isector
    integer                            :: jdim,jsector
    real(8)                            :: cdgOp_mat,cOp_mat
    real(8)                            :: sgn_cdg,sgn_c
    integer                            :: ib(Nlevels)
    integer                            :: n,m,p,iexc
    integer                            :: ni,mj
    integer                            :: i,j,r,k
    real(8)                            :: expterm
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
             cdgOp_mat = 0d0
             cOp_mat   = 0d0
             !
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
                cdgOp_mat = cdgOp_mat + espace(jsector)%M(mj,j)*sgn_cdg*espace(isector)%M(ni,i)
                !
             enddo
             !
             do mj=1,jdim              !loop over the component of |J> (IN state!)
                m  = HJ%map(mj)
                ib = bdecomp(m,2*Ns)
                if(ib(jsite) == 0)cycle
                call c(jsite,m,n,sgn_c)
                ni = binary_search(HI%map,n)
                !
                cOp_mat = cOp_mat + espace(isector)%M(ni,i)*sgn_c*espace(jsector)%M(mj,j)
                !
             enddo
             !
             cdgQmatrix(ilat,ispin,iorb,iexc) = cdgOp_mat*sqrt(expterm/zeta_function)
             cQmatrix(jlat,ispin,jorb,iexc)   =   cOp_mat*sqrt(expterm/zeta_function)
             !             
          enddo
       enddo
       call delete_sector(isector,HI)
       call delete_sector(jsector,HJ)
    enddo
    !
  end subroutine full_build_Qmatrix


  subroutine full_build_Lmatrix(ispin)
    integer                            :: ilat,jlat
    integer                            :: iorb,jorb
    integer                            :: ispin,jspin
    integer                            :: isite,jsite
    integer                            :: idim,isector
    integer                            :: jdim,jsector
    real(8)                            :: cdgOp_mat,cOp_mat
    real(8)                            :: sgn_cdg,sgn_c
    integer                            :: ib(Nlevels)
    integer                            :: n,m,p,iexc
    integer                            :: ni,mj
    integer                            :: i,j,r,k
    real(8)                            :: expterm
    type(sector_map)                   :: HI,HJ
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
             !
             expterm=exp(-beta*espace(isector)%e(i))+exp(-beta*espace(jsector)%e(j))
             if(expterm < cutoff)cycle
             !
             iexc=iexc+1
             !
             Lmatrix(ispin,iexc) = espace(jsector)%e(j)-espace(isector)%e(i)
             !             
          enddo
       enddo
       call delete_sector(isector,HI)
       call delete_sector(jsector,HJ)
    enddo
    !
  end subroutine full_build_Lmatrix







  !+------------------------------------------------------------------+
  !PURPOSE  : Allocate Lambda and Q matrices, storing the spectral
  ! sum terms
  !+------------------------------------------------------------------+
  subroutine allocate_QLmatrix()
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
    allocate(Lmatrix(2,Nexc))
    Lmatrix=0d0
    !
    write(LOGfile,*)"Found",Nexc," excitations with the actual cut-off"
    !
  end subroutine allocate_QLmatrix







end MODULE VCA_GREENS_FUNCTIONS














! !+------------------------------------------------------------------+
! !PURPOSE  : Evaluate Green's functions
! !+------------------------------------------------------------------+
! subroutine build_gf_normal()
!   integer :: ilat,jlat,iorb,jorb,ispin
!   !
!   call allocate_QLmatrix()
!   !
!   call start_timer
!   !        
!   do ispin=1,Nspin
!      do ilat=1,Nlat
!         do jlat=1,Nlat
!            do iorb=1,Norb
!               do jorb=1,Norb
!                  write(LOGfile,"(A)")"Get G_i"//str(ilat,3)//"_j"//str(jlat,3)//&
!                       "_l"//str(iorb)//"_m"//str(jorb)//"_s"//str(ispin)
!                  call full_build_gf_QLmatrix(ilat,jlat,iorb,iorb,ispin)
!               enddo
!            enddo
!         enddo
!      enddo
!   enddo
!   !
!   call stop_timer
!   !
! end subroutine build_gf_normal



! !+------------------------------------------------------------------+
! !PURPOSE  : Spectral sum evaluation of the Green's functions
! ! Matsubara and Real-Axis
! !+------------------------------------------------------------------+
! subroutine full_build_gf_QLmatrix(ilat,jlat,iorb,jorb,ispin)
!   integer                            :: ilat,jlat
!   integer                            :: iorb,jorb
!   integer                            :: ispin,jspin
!   integer                            :: isite,jsite
!   integer                            :: idim,isector
!   integer                            :: jdim,jsector
!   real(8)                            :: cdgOp_mat,cOp_mat
!   real(8)                            :: op_mat(2),matCC
!   real(8)                            :: spectral_weight
!   real(8)                            :: sgn_cdg,sgn_c
!   integer                            :: ib(Nlevels)
!   integer                            :: n,m,p,iexc
!   integer                            :: ni,mj
!   integer                            :: i,j,r,k
!   real(8)                            :: Ei,Ej
!   real(8)                            :: expterm,de,w0
!   complex(8)                         :: iw
!   type(sector_map)                   :: HI,HJ
!   !
!   isite = state_index(ilat,ispin,iorb)
!   jsite = state_index(jlat,ispin,jorb)
!   !
!   iexc = 0
!   do isector=1,Nsectors
!      jsector=getCDGsector(ispin,isector)
!      if(jsector==0)cycle
!      !
!      idim=getdim(isector)     !i-th sector dimension
!      jdim=getdim(jsector)     !j-th sector dimension
!      call build_sector(isector,HI)
!      call build_sector(jsector,HJ)
!      !
!      do i=1,idim          !loop over the states in the i-th sect.
!         do j=1,jdim       !loop over the states in the j-th sect.
!            op_mat = 0d0
!            matCC = 0d0
!            expterm=exp(-beta*espace(isector)%e(i))+exp(-beta*espace(jsector)%e(j))
!            if(expterm < cutoff)cycle
!            !
!            iexc=iexc+1
!            !
!            do ni=1,idim              !loop over the component of |I> (IN state!)
!               n  = HI%map(ni)
!               ib = bdecomp(n,2*Ns)
!               if(ib(isite) == 1)cycle
!               call cdg(isite,n,m,sgn_cdg)
!               mj = binary_search(HJ%map,m)
!               !
!               cdgOp_mat = espace(jsector)%M(mj,j)*sgn_cdg*espace(isector)%M(ni,i)
!               !
!               op_mat(1)=op_mat(1) + cdgOp_mat
!            enddo
!            !
!            do mj=1,jdim              !loop over the component of |J> (IN state!)
!               m  = HJ%map(mj)
!               ib = bdecomp(m,2*Ns)
!               if(ib(jsite) == 0)cycle
!               call c(jsite,m,n,sgn_c)
!               ni = binary_search(HI%map,n)
!               !
!               cOp_mat   = espace(isector)%M(ni,i)*sgn_c*espace(jsector)%M(mj,j)
!               !
!               op_mat(2)=op_mat(2) + cOp_mat
!               !
!            enddo
!            !
!            Ei=espace(isector)%e(i)
!            Ej=espace(jsector)%e(j)
!            de=Ej-Ei
!            spectral_weight=expterm/zeta_function*product(op_mat)
!            !
!            cdgQmatrix(ilat,ispin,iorb,iexc) = op_mat(1)*sqrt(expterm/zeta_function)
!            cQmatrix(jlat,ispin,jorb,iexc)   = op_mat(2)*sqrt(expterm/zeta_function)
!            Lmatrix(ilat,jlat,ispin,ispin,iorb,jorb,iexc) = de
!            !             
!            do m=1,Lmats
!               iw=xi*wm(m)
!               impGmats(ilat,jlat,ispin,ispin,iorb,jorb,m)=impGmats(ilat,jlat,ispin,ispin,iorb,jorb,m)+spectral_weight/(iw-de)
!            enddo
!            do m=1,Lreal
!               iw=dcmplx(wr(m),eps)
!               impGreal(ilat,jlat,ispin,ispin,iorb,jorb,m)=impGreal(ilat,jlat,ispin,ispin,iorb,jorb,m)+spectral_weight/(iw-de)
!            enddo
!            !
!         enddo
!      enddo
!      call delete_sector(isector,HI)
!      call delete_sector(jsector,HJ)
!   enddo
!   !
! end subroutine full_build_gf_QLmatrix
