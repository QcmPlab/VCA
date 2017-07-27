module VCA_MAIN
  USE VCA_VARS_GLOBAL
  USE VCA_AUX_FUNX
  USE VCA_SETUP
  USE VCA_DIAG
  USE VCA_QMATRIX
  USE VCA_OBSERVABLES
  !
  USE SF_LINALG,  only: eigh,diag
  USE SF_IOTOOLS, only: str
  USE SF_TIMER,only: start_timer,stop_timer
  USE SF_MISC, only: assert_shape
  implicit none
  private

  !>INIT VCA SOLVER
  public :: vca_init_solver


  !> VCA DIAG CLUSTER
  public :: vca_diag_cluster


  ! !> VCA GF POLES UPDATE
  ! public :: vca_diag_system



contains



  !+-----------------------------------------------------------------------------+!
  ! PURPOSE: allocate and initialize one or multiple baths -+!
  !+-----------------------------------------------------------------------------+!
  subroutine vca_init_solver(Hloc)
    real(8),intent(in)              :: Hloc(:,:,:,:,:,:) !Hloc(Nlat,Nlat,Nspin,Nspin,Norb,Norb)
    logical,save                    :: isetup=.true.
    !
    write(LOGfile,"(A)")"INIT SOLVER FOR "//trim(file_suffix)
    !
    Nlat = size(Hloc,1)
    call assert_shape(Hloc,[Nlat,Nlat,Nspin,Nspin,Norb,Norb],"vca_init_solver","Hloc")
    !
    !Init Structure & memory
    if(isetup)call init_cluster_structure()
    !
    !Init Hcluster
    call set_Hcluster(Hloc)
    !
    if(isetup)call setup_pointers_normal
    isetup=.false.
    !
  end subroutine vca_init_solver


  !+-----------------------------------------------------------------------------+!
  !PURPOSE: Diag the cluster, reference system
  !+-----------------------------------------------------------------------------+!
  subroutine vca_diag_cluster(Hloc)
    real(8),optional,intent(in) :: Hloc(:,:,:,:,:,:) !Hloc(Nlat,Nlat,Nspin,Nspin,Norb,Norb)
    !
    Nlat = size(Hloc,1)
    call assert_shape(Hloc,[Nlat,Nlat,Nspin,Nspin,Norb,Norb],"vca_diag","Hloc")
    !
    if(present(Hloc))call set_Hcluster(Hloc)
    !
    call setup_eigenspace()
    !
    call diagonalize_cluster()         !find target states by digonalization of Hamiltonian
    call observables_cluster()         !obtain impurity observables as thermal averages.  
    call build_Qmatrix_cluster()       !build the one-particle Q and Lambda matrices, storing spectral sums
    !
    call delete_eigenspace()
  end subroutine vca_diag_cluster





  ! !+-----------------------------------------------------------------------------+!
  ! !PURPOSE: Get the spectrum of the system GF from the Cluster tiling 
  ! !+-----------------------------------------------------------------------------+!
  ! subroutine vca_diag_system(Vmat)
  !   real(8),dimension(Nlat*Nspin*Norb,Nlat*Nspin*Norb) :: Vmat
  !   integer                                            :: ilat,jlat
  !   integer                                            :: iorb,jorb
  !   integer                                            :: ispin
  !   integer                                            :: iexc,jexc
  !   integer                                            :: i,j,Nexc
  !   real(8),dimension(:,:),allocatable                 :: Mmat
  !   real(8),dimension(:),allocatable                   :: Lvec

  !   do ispin=1,Nspin
  !      !
  !      Nexc = Qcluster(ispin)%Nexc
  !      call allocate_Qmatrix(Qsystem(ispin),ispin)
  !      !
  !      allocate(Mmat(Nexc,Nexc))
  !      allocate(Lvec(Nexc))
  !      !
  !      Mmat = diag( Qcluster(ispin)%poles(:) )
  !      !
  !      ! do iexc=1,Nexc
  !      !    do jexc=1,Nexc
  !      !       !
  !      !       do i=1,Nlat*Norb
  !      !          do j=1,Nlat*Norb
  !      !             Mmat(iexc,jexc) = Mmat(iexc,jexc) + &
  !      !                  Qcluster(ispin)%cdg(iexc,i)*Vmat(i,j)*Qcluster(ispin)%c(j,jexc)
  !      !          enddo
  !      !       enddo
  !      !       !
  !      !    enddo
  !      ! enddo
  !      Mmat = matmul(Qcluster(ispin)%cdg, matmul(Vmat,Qcluster(ispin)%c))
  !      Mmat = Mmat + diag( Qcluster(ispin)%poles )
  !      !
  !      call eigh(Mmat,Lvec)
  !      !
  !      do i=1,Nlat*Norb
  !         do iexc=1,Nexc
  !            do jexc=1,Nexc
  !               Qsystem(ispin)%c(i,iexc)   = Qsystem(ispin)%c(i,iexc)   + &
  !                    Qcluster(ispin)%c(i,jexc)*Mmat(jexc,iexc)
  !               !
  !               Qsystem(ispin)%cdg(iexc,i) = Qsystem(ispin)%cdg(iexc,i) + &
  !                    Mmat(iexc,jexc)*Qcluster(ispin)%cdg(jexc,i)
  !            enddo
  !         enddo
  !      enddo
  !      Qsystem(ispin)%poles(:) = Lvec
  !      !
  !   enddo
  ! end subroutine vca_diag_system







end module VCA_MAIN
