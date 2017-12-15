MODULE VCA_HAMILTONIAN_MATVEC
  USE SCIFOR, only:zero,assert_shape
  USE VCA_INPUT_VARS
  USE VCA_VARS_GLOBAL
  USE VCA_BATH_SETUP
  USE VCA_SETUP

  implicit none
  private

  !>build sparse hamiltonian of the sector
  public  :: buildH_c
  !
  !
  !>Sparse Mat-Vec product using stored sparse matrix 
  public  :: spMatVec_cc
  !
  !
  !> Related auxiliary routines:
  public  :: setup_Hv_sector
  public  :: delete_Hv_sector

  integer                      :: Hsector=0
  logical                      :: Hstatus=.false.
  type(sector_map)             :: H,Hup,Hdw

contains



  !####################################################################
  !             BUILD SPARSE HAMILTONIAN of the SECTOR
  !####################################################################
  subroutine buildH_c(Hmat)
    complex(8),dimension(:,:),optional    :: Hmat
    complex(8),dimension(:,:),allocatable :: Hredux
    integer                               :: isector
    integer,dimension(Nlevels)            :: ib
    real(8),dimension(Nlat,Norb)          :: nup,ndw
    integer                               :: dim,dimUp,dimDw
    integer                               :: i
    integer                               :: m
    integer                               :: ishift,ishift_up,ishift_dw
    integer                               :: j,ms,impi
    integer                               :: ilat,jlat,iorb,jorb,ispin,jspin,is,js
    integer                               :: i_up,i_dw,j_up,j_dw
    integer                               :: kp,k1,k2,k3,k4
    integer                               :: alfa,beta
    real(8)                               :: sg1,sg2,sg3,sg4
    complex(8)                            :: htmp,htmpup,htmpdw
    logical                               :: Jcondition
    !
    if(.not.Hstatus)stop "buildH_c ERROR: Hsector NOT set"
    isector=Hsector
    !
    if(spH0%status)call sp_delete_matrix(spH0)
    !
    dim=getdim(isector)
    !
    call sp_init_matrix(spH0,dim)
    !
    !
    !-----------------------------------------------!
    states: do i=1,Dim
       m = H%map(i)
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
    !
    if(present(Hmat))then
       call assert_shape(Hmat,[dim,dim],"buildH_c","Hmat")
       Hmat=zero
       call sp_dump_matrix(spH0,Hmat)
    endif
  end subroutine buildH_c









  !####################################################################
  !        SPARSE MAT-VEC PRODUCT USING STORED SPARSE MATRIX 
  !####################################################################
  !+------------------------------------------------------------------+
  !PURPOSE: Perform the matrix-vector product H*v used in the
  !+------------------------------------------------------------------+
  subroutine spMatVec_cc(Nloc,v,Hv)
    integer                      :: Nloc
    complex(8),dimension(Nloc)   :: v
    complex(8),dimension(Nloc)   :: Hv
    integer                      :: i
    type(sparse_element),pointer :: c
    Hv=zero
    do i=1,Nloc
       c => spH0%row(i)%root%next       
       matmul: do while(associated(c))
          Hv(i) = Hv(i) + c%cval*v(c%col)
          c => c%next
       end do matmul
    end do
    nullify(c)
  end subroutine spMatVec_cc






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



END MODULE VCA_HAMILTONIAN_MATVEC
