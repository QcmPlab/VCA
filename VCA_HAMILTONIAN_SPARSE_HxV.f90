! > BUILD SPARSE HAMILTONIAN of the SECTOR
MODULE VCA_HAMILTONIAN_SPARSE_HxV
  USE VCA_HAMILTONIAN_COMMON
  implicit none
  private

  !>build sparse hamiltonian of the sector
  !public  :: buildH_c
  public :: vca_buildh_main
  !public :: vca_buildh_orbs
  !
  !
  !>Sparse Mat-Vec product using stored sparse matrix 
  public  :: spMatVec_main
  !public :: spMatVec_orbs
  !
  !
  !> Related auxiliary routines:
  !public  :: setup_Hv_sector
  !public  :: delete_Hv_sector

  integer                              :: i,iup,idw
  integer                              :: j,jup,jdw
  integer                              :: m,mup,mdw
  integer                              :: ms,iud
  integer                              :: impi
  integer                              :: ilat,jlat,iorb,jorb,ispin,jspin,is,js
  integer                              :: kp,k1,k2,k3,k4
  integer                              :: ialfa,ibeta,indx
  real(8)                              :: sg1,sg2,sg3,sg4
  real(8)                              :: htmp,htmpup,htmpdw
  logical                              :: Jcondition
  integer                              :: Nfoo

contains



  !####################################################################
  !             BUILD SPARSE HAMILTONIAN of the SECTOR
  !####################################################################
  subroutine vca_buildh_main(isector,Hmat)
    integer                              :: isector
    real(8),dimension(:,:),optional      :: Hmat
    real(8),dimension(:,:),allocatable   :: Htmp_up,Htmp_dw,Hrdx
    integer,dimension(Ns)                :: ibup,ibdw
    integer,dimension(2*Ns_Ud)           :: Indices    ![2-2*Norb]  CONTROLLA
    integer,dimension(Ns_Ud,Ns_Orb)      :: Nups,Ndws  ![1,Ns]-[Norb,1+Nbath]  CONTROLLA
    integer,dimension(Nlat,Norb)         :: Nup,Ndw    !CONTROLLA

    !
    if(.not.Hstatus)stop "buildH_c ERROR: Hsector NOT set"
    isector=Hsector
    !
    !
    if(present(Hmat))&
         call assert_shape(Hmat,[getdim(isector), getdim(isector)],"vca_buildh_main","Hmat")
    !
    if(spH0d%status)call sp_delete_matrix(spH0d)
    !
    !Get diagonal hybridization, bath energy
    !include "VCA_HAMILTONIAN/diag_hybr_bath.f90"                 CONTROLLA
    !    
    call sp_init_matrix(spH0d,Dim)
    call sp_init_matrix(spH0dws(1),DimDw)
    call sp_init_matrix(spH0ups(1),DimUp)
    !
    !
    !-----------------------------------------------!
    !LOCAL CLUSTER HAMILTONIAN TERMS
    include "VCA_HAMILTONIAN/stored/H_cluster_local.f90"
    !
    !UP CLUSTER TERMS
    include "VCA_HAMILTONIAN/stored/H_cluster_up.f90"
    !
    !DW CLUSTER TERMS
    include "VCA_HAMILTONIAN/stored/H_cluster_dw.f90"
    !-----------------------------------------------!

    !   if(vca_bath%status)then
    !      include "vca_Hbath.f90"
          !
    !      include "vca_Hhyb_bath.f90"
    !   endif
    !

    if(present(Hmat))then
       Hmat = 0d0
       allocate(Htmp_up(DimUp,DimUp));Htmp_up=0d0
       allocate(Htmp_dw(DimDw,DimDw));Htmp_dw=0d0
       !
!#ifdef _MPI
 !      if(MpiStatus)then
 !         call sp_dump_matrix(MpiComm,spH0d,Hmat)
 !      else
          call sp_dump_matrix(spH0d,Hmat)
 !      endif
!#else
       call sp_dump_matrix(spH0d,Hmat)
!#endif
       call sp_dump_matrix(spH0ups(1),Htmp_up)
       call sp_dump_matrix(spH0dws(1),Htmp_dw)
       Hmat = Hmat + kronecker_product(eye(DimUp),Htmp_dw) ! kron(eye(DimDw),Htmp_up)
       Hmat = Hmat + kronecker_product(Htmp_up,eye(DimDw)) ! kron(Htmp_dw,eye(DimUp))
       !
       deallocate(Htmp_up,Htmp_dw)
    endif
    !
    !deallocate(diag_hybr,bath_diag)
    return

    !if(present(Hmat))then
       !call assert_shape(Hmat,[dim,dim],"buildH_c","Hmat")
       !Hmat=zero
       !call sp_dump_matrix(spH0,Hmat)
    !endif
  end subroutine vca_buildh_main


  !####################################################################
  !        TODO: ORBS
  !####################################################################






  !####################################################################
  !        SPARSE MAT-VEC PRODUCT USING STORED SPARSE MATRIX 
  !####################################################################
  !+------------------------------------------------------------------+
  !PURPOSE: Perform the matrix-vector product H*v used in the
  !+------------------------------------------------------------------+
  subroutine spMatVec_main(Nloc,v,Hv)
    integer                         :: Nloc
    real(8),dimension(Nloc)         :: v
    real(8),dimension(Nloc)         :: Hv
    real(8)                         :: val
    integer                         :: i,iup,idw,j,jup,jdw,jj
    !
    !
    Hv=0d0
    !
    do i = 1,Nloc
       do j=1,spH0d%row(i)%Size
          Hv(i) = Hv(i) + spH0d%row(i)%vals(j)*v(spH0d%row(i)%cols(j))
       enddo
    enddo
    !
    !DW:
    do iup=1,DimUp
       !
       do idw=1,DimDw
          i = iup + (idw-1)*DimUp
          do jj=1,spH0dws(1)%row(idw)%Size
             jup = iup
             jdw = spH0dws(1)%row(idw)%cols(jj)
             val = spH0dws(1)%row(idw)%vals(jj)
             j     = jup +  (jdw-1)*DimUp
             Hv(i) = Hv(i) + val*V(j)
          enddo
       enddo
       !
    enddo
    !
    !UP:
    do idw=1,DimDw
       !
       do iup=1,DimUp
          i = iup + (idw-1)*DimUp
          do jj=1,spH0ups(1)%row(iup)%Size
             jup = spH0ups(1)%row(iup)%cols(jj)
             jdw = idw
             val = spH0ups(1)%row(iup)%vals(jj)
             j =  jup + (jdw-1)*DimUp
             Hv(i) = Hv(i) + val*V(j)
          enddo
       enddo
       !
    enddo
    !
  end subroutine spMatVec_main


  !####################################################################
  !        TODO: ORBS
  !####################################################################


END MODULE VCA_HAMILTONIAN_SPARSE_HxV
