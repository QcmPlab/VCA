module VCA_DIAG
  USE VCA_VARS_GLOBAL
  USE VCA_SETUP
  !
  USE SF_CONSTANTS
  USE SF_LINALG, only: eigh
  USE SF_TIMER,  only: start_timer,stop_timer,eta
  USE SF_IOTOOLS, only:reg,free_unit
  USE SF_MISC, only: assert_shape
  !
  implicit none
  private


  integer                      :: Hsector=0
  logical                      :: Hstatus=.false.
  type(sector_map)             :: H,Hup,Hdw


  !>build sparse hamiltonian of the sector
  public  :: buildH_c

  !>auxiliary routines
  public  :: setup_Hv_sector
  public  :: delete_Hv_sector

  !>diagonalize cluster Hamiltonian
  public  :: diagonalize_cluster



contains




  !> Set the actual Sector
  subroutine setup_Hv_sector(isector)
    integer                   :: isector
    Hsector=isector
    Hstatus=.true.
    call build_sector(isector,H)
  end subroutine setup_Hv_sector

  !> Reset the sector 
  subroutine delete_Hv_sector()
    call delete_sector(Hsector,H)
    Hsector=0
    Hstatus=.false.
  end subroutine delete_Hv_sector




  !+-------------------------------------------------------------------+
  !PURPOSE  : diagonalize the Hamiltonian in each sector 
  !+------------------------------------------------------------------+
  subroutine diagonalize_cluster
    integer                     :: nup,ndw,isector,dim
    integer                     :: sz,nt
    integer                     :: i,j,unit
    real(8),dimension(Nsectors) :: e0 
    real(8)                     :: egs
    logical                     :: Tflag
    !
    e0=1000.d0
    write(LOGfile,"(A)")"Diagonalize Cluster Hc+Hint:"
    call start_timer()
    !
    sector: do isector=1,Nsectors
       !
       Dim      = getdim(isector)
       !
       if(verbose==3)then
          nup  = getnup(isector)
          ndw  = getndw(isector)
          write(LOGfile,"(A,I4,A6,I2,A6,I2,A6,I15)")"Solving sector:",isector,", nup:",nup,", ndw:",ndw,", dim=",getdim(isector)
       elseif(verbose==1.OR.verbose==2)then
          call eta(isector,Nsectors,LOGfile)
       endif
       !
       call setup_Hv_sector(isector)
       call buildH_c(espace(isector)%M)
       call delete_Hv_sector()
       call eigh(espace(isector)%M,espace(isector)%e,'V','U')
       if(dim==1)espace(isector)%M=one
       !
       e0(isector)=minval(espace(isector)%e)
       !
    enddo sector
    !
    call stop_timer(LOGfile)
    !
    !Get the ground state energy and rescale energies
    write(LOGfile,"(A)")"DIAG summary:"
    egs=minval(e0)
    open(free_unit(unit),file='egs'//reg(file_suffix)//".vca",position='append')
    do isector=1,Nsectors
       if(abs(e0(isector)-egs)>gs_threshold)cycle
       nup  = getnup(isector)
       ndw  = getndw(isector)
       write(LOGfile,"(A,F20.12,2I4)")'Egs =',e0(isector),nup,ndw
       write(unit,"(F20.12,2I4)")e0(isector),nup,ndw
    enddo
    close(unit)
    !
    !Get the partition function Z
    zeta_function=0.d0
    forall(isector=1:Nsectors)espace(isector)%e = espace(isector)%e - egs
    do isector=1,Nsectors
       dim=getdim(isector)
       do i=1,dim
          zeta_function=zeta_function+exp(-beta*espace(isector)%e(i))
       enddo
    enddo
    omega_potential = -1d0/beta*log(zeta_function)
    write(LOGfile,"(A,F20.12)")'Z     =',zeta_function
    write(LOGfile,"(A,F20.12)")'Omega =',omega_potential
    !
    return
  end subroutine diagonalize_cluster





  !> Build the Cluster Hamiltonian (into Hmat)
  subroutine buildH_c(Hmat)
    complex(8),dimension(:,:)    :: Hmat
    integer                      :: isector
    integer,dimension(Nlevels)   :: ib
    real(8),dimension(Nlat,Norb) :: nup,ndw
    integer                      :: dim,dimUp,dimDw
    integer                      :: i
    integer                      :: m
    integer                      :: ishift,ishift_up,ishift_dw
    integer                      :: j,ms,impi
    integer                      :: ilat,jlat,iorb,jorb,ispin,jspin,is,js
    integer                      :: i_up,i_dw,j_up,j_dw
    integer                      :: kp,k1,k2,k3,k4
    integer                      :: alfa,beta
    real(8)                      :: sg1,sg2,sg3,sg4
    complex(8)                   :: htmp,htmpup,htmpdw
    logical                      :: Jcondition
    !
    if(.not.Hstatus)stop "buildH_c ERROR: Hsector NOT set"
    isector=Hsector
    !
    dim=getdim(isector)
    !
    call assert_shape(Hmat,[dim,dim],"buildH_c","Hmat")
    !
    Hmat=zero
    !
    !-----------------------------------------------!
    states: do i=1,Dim
       m = H%map(i)
       impi = i
       ib = bdecomp(m,2*Ns)
       !
       do ilat=1,Nlat
          do iorb=1,Norb
             nup(ilat,iorb)=dble(ib(imp_state_index(ilat,iorb,1)))
             ndw(ilat,iorb)=dble(ib(imp_state_index(ilat,iorb,2)))
          enddo
       enddo
       !
       include "vca_Hcluster.f90" 
       !
       include "vca_Hint.f90"
       !
       if(vca_bath%status)then
          include "vca_Hbath.f90"
          !
          include "vca_Hhyb_bath.f90"
       endif
       !
    enddo states
  end subroutine buildH_c

END MODULE VCA_DIAG




! !> Count number of 1p-excitations per spin-channel 
! subroutine enumerate_Nexcitations
!   integer          :: ispin
!   integer          :: isector,jsector
!   integer          :: idim,jdim
!   integer          :: i,j,iexc
!   real(8)          :: expterm
!   iexc=0
!   do ispin=1,Nspin
!      do isector=1,Nsectors
!         jsector=getCDGsector(ispin,isector)
!         if(jsector==0)cycle
!         idim=getdim(isector)     !i-th sector dimension
!         jdim=getdim(jsector)     !j-th sector dimension
!         do i=1,idim          !loop over the states in the i-th sect.
!            do j=1,jdim       !loop over the states in the j-th sect.
!               expterm=exp(-beta*espace(isector)%e(i))+exp(-beta*espace(jsector)%e(j))
!               if(expterm < cutoff)cycle
!               iexc=iexc+1
!            enddo
!         enddo
!      enddo
!   enddo
!   !
!   Nexcitations = iexc
!   !
!   write(LOGfile,"(A,I5,A,I2)")"Found N =",iexc," excitations with the actual cut-off"
!   !
! end subroutine enumerate_Nexcitations









