MODULE VCA_QMATRIX
  USE VCA_VARS_GLOBAL
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


  !> Allocate and Destruct a generic Qmatrix
  public :: allocate_Qmatrix
  public :: deallocate_Qmatrix



  !> Build the Qmatrix structure for the cluster
  public :: build_Qmatrix_cluster



contains



  subroutine build_Qmatrix_cluster()
    integer :: ilat,jlat,iorb,jorb,ispin
    !
    call start_timer
    !
    call allocate_Qmatrix(Qcluster)
    !
    call full_build_Lmatrix(Qcluster)
    ! do ispin=1,Nspin
    !    !
    !    do ilat=1,Nlat
    !       do jlat=1,Nlat
    !          do iorb=1,Norb
    !             do jorb=1,Norb
    !                !
    !                write(LOGfile,"(A)")"Get Q_i"//str(ilat,3)//"_j"//str(jlat,3)//"_l"//str(iorb)//"_m"//str(jorb)//"_s"//str(ispin)
    !
    call full_build_Qmatrix(Qcluster)!,ispin,ilat,jlat,iorb,jorb)
    !
    !             enddo
    !          enddo
    !       enddo
    !    enddo
    ! enddo
    !
    call stop_timer
    !
  end subroutine build_Qmatrix_cluster









  !+------------------------------------------------------------------+
  !PURPOSE  : Build the Lambda matrix, storing 1-particle excitations
  !+------------------------------------------------------------------+
  subroutine full_build_Lmatrix(matrix)
    type(Qmatrix)    :: matrix
    integer          :: ispin
    integer          :: idim,isector
    integer          :: jdim,jsector
    integer          :: i,j,iexc
    real(8)          :: expterm
    !
    iexc = 0
    do ispin=1,Nspin
       !
       do isector=1,Nsectors
          jsector=getCDGsector(ispin,isector)
          if(jsector==0)cycle
          !
          idim=getdim(isector)     !i-th sector dimension
          jdim=getdim(jsector)     !j-th sector dimension
          !
          do i=1,idim          !loop over the states in the i-th sect.
             do j=1,jdim       !loop over the states in the j-th sect.
                !
                expterm=exp(-beta*espace(isector)%e(i))+exp(-beta*espace(jsector)%e(j))
                if(expterm < cutoff)cycle
                !
                iexc=iexc+1
                !
                matrix%poles(iexc) = espace(jsector)%e(j)-espace(isector)%e(i)
                !             
             enddo
          enddo
       enddo
       !
    enddo
    !
  end subroutine full_build_Lmatrix




  !+------------------------------------------------------------------+
  !PURPOSE  : Spectral sum evaluation of the Green's functions
  ! Matsubara and Real-Axis
  !+------------------------------------------------------------------+
  subroutine full_build_Qmatrix(matrix)!,ispin,ilat,jlat,iorb,jorb)
    type(Qmatrix)    :: matrix
    integer          :: ilat,jlat
    integer          :: iorb,jorb
    integer          :: ispin,jspin
    integer          :: isite,jsite
    integer          :: idim,isector
    integer          :: jdim,jsector
    real(8)          :: cdgOp_mat,cOp_mat
    real(8)          :: sgn_cdg,sgn_c
    integer          :: ib(Nlevels)
    integer          :: n,m,p,iexc
    integer          :: ni,mj
    integer          :: i,j,is,js
    real(8)          :: expterm
    type(sector_map) :: HI,HJ
    !
    do iorb=1,Norb
       do ilat=1,Nlat
          iexc = 0
          do ispin=1,Nspin
             isite = state_index(ilat,iorb,ispin)
             write(LOGfile,"(A,A)")&
                  "Get Q_i"//str(ilat,3)//"_l"//str(iorb)//"_s"//str(ispin),&
                  "Q^+_i"//str(ilat,3)//"_l"//str(iorb)//"_s"//str(ispin)
             do isector=1,Nsectors
                jsector=getCDGsector(ispin,isector)
                if(jsector==0)cycle
                idim=getdim(isector)     !i-th sector dimension
                jdim=getdim(jsector)     !j-th sector dimension
                call build_sector(isector,HI)
                call build_sector(jsector,HJ)
                do i=1,idim          !loop over the states in the i-th sect.
                   do j=1,jdim       !loop over the states in the j-th sect.
                      expterm=exp(-beta*espace(isector)%e(i))+exp(-beta*espace(jsector)%e(j))
                      if(expterm < cutoff)cycle
                      !
                      iexc=iexc+1
                      !
                      cdgOp_mat = 0d0
                      do ni=1,idim              !loop over the component of |I> (IN state!)
                         n  = HI%map(ni)
                         ib = bdecomp(n,2*Ns)
                         if(ib(isite) == 1)cycle
                         call cdg(isite,n,m,sgn_cdg)
                         mj = binary_search(HJ%map,m)
                         cdgOp_mat = cdgOp_mat + espace(jsector)%M(mj,j)*sgn_cdg*espace(isector)%M(ni,i)
                      enddo
                      matrix%cdg(iexc,isite) = cdgOp_mat*sqrt(expterm/zeta_function)
                      !
                      cOp_mat   = 0d0
                      do mj=1,jdim              !loop over the component of |J> (IN state!)
                         m  = HJ%map(mj)
                         ib = bdecomp(m,2*Ns)
                         if(ib(isite) == 0)cycle
                         call c(isite,m,n,sgn_c)
                         ni = binary_search(HI%map,n)
                         cOp_mat = cOp_mat + espace(isector)%M(ni,i)*sgn_c*espace(jsector)%M(mj,j)
                      enddo
                      matrix%c(isite,iexc)   =   cOp_mat*sqrt(expterm/zeta_function)
                      !             
                   enddo
                enddo
                call delete_sector(isector,HI)
                call delete_sector(jsector,HJ)
             enddo
             !
          enddo
       enddo
    enddo
    !
  end subroutine full_build_Qmatrix











  !+------------------------------------------------------------------+
  !PURPOSE  : Allocate Lambda and Q matrices, storing the spectral
  ! sum terms
  !+------------------------------------------------------------------+
  subroutine allocate_Qmatrix(matrix)
    type(Qmatrix)    :: matrix
    integer          :: Nexc
    !
    if(matrix%allocated)call deallocate_Qmatrix(matrix)
    !
    write(LOGfile,"(A)")"Allocating Qmatrix"
    !
    Nexc = Nexcitations
    !
    matrix%Nexc = Nexc
    !
    matrix%allocated = .true.
    !
    allocate(matrix%c(Nlat*Norb,Nexc))
    matrix%c=0d0
    !
    allocate(matrix%cdg(Nexc,Nlat*Norb))
    matrix%cdg=0d0
    !
    allocate(matrix%poles(Nexc))
    matrix%poles=0d0
    !
  end subroutine allocate_Qmatrix


  subroutine deallocate_Qmatrix(matrix)
    type(Qmatrix) :: matrix
    !
    if(allocated(matrix%c))deallocate(matrix%c)
    if(allocated(matrix%cdg))deallocate(matrix%cdg)
    if(allocated(matrix%poles))deallocate(matrix%poles)
    !
    matrix%Nexc = 0
    !
    matrix%allocated = .false.
    !
  end subroutine deallocate_Qmatrix




end MODULE VCA_QMATRIX













! subroutine full_build_Qmatrix(matrix,ispin,ilat,jlat,iorb,jorb)
!   type(Qmatrix)    :: matrix
!   integer          :: ilat,jlat
!   integer          :: iorb,jorb
!   integer          :: ispin,jspin
!   integer          :: isite,jsite
!   integer          :: idim,isector
!   integer          :: jdim,jsector
!   real(8)          :: cdgOp_mat,cOp_mat
!   real(8)          :: sgn_cdg,sgn_c
!   integer          :: ib(Nlevels)
!   integer          :: n,m,p,iexc
!   integer          :: ni,mj
!   integer          :: i,j,is,js
!   real(8)          :: expterm
!   type(sector_map) :: HI,HJ
!   !
!   isite = state_index(ilat,iorb,ispin)
!   jsite = state_index(jlat,jorb,ispin)
!   is = iorb + (ilat-1)*Norb
!   js = jorb + (jlat-1)*Norb
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
!            cdgOp_mat = 0d0
!            cOp_mat   = 0d0
!            !
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
!               cdgOp_mat = cdgOp_mat + espace(jsector)%M(mj,j)*sgn_cdg*espace(isector)%M(ni,i)
!               !
!            enddo
!            !
!            do mj=1,jdim              !loop over the component of |J> (IN state!)
!               m  = HJ%map(mj)
!               ib = bdecomp(m,2*Ns)
!               if(ib(jsite) == 0)cycle
!               call c(jsite,m,n,sgn_c)
!               ni = binary_search(HI%map,n)
!               !
!               cOp_mat = cOp_mat + espace(isector)%M(ni,i)*sgn_c*espace(jsector)%M(mj,j)
!               !
!            enddo
!            !
!            matrix%cdg(iexc,is) = cdgOp_mat*sqrt(expterm/zeta_function)
!            matrix%c(js,iexc)   =   cOp_mat*sqrt(expterm/zeta_function)
!            !             
!         enddo
!      enddo
!      call delete_sector(isector,HI)
!      call delete_sector(jsector,HJ)
!   enddo
!   !
! end subroutine full_build_Qmatrix
