module VCA_MAIN
  USE VCA_VARS_GLOBAL
  USE VCA_AUX_FUNX
  USE VCA_SETUP
  USE VCA_DIAG
  USE VCA_QMATRIX
  USE VCA_OBSERVABLES
  !
  USE SF_LINALG,  only: eigh,diag
  USE SF_IOTOOLS, only: str,free_unit
  USE SF_TIMER,only: start_timer,stop_timer
  USE SF_MISC, only: assert_shape
  implicit none
  private

  !>INIT VCA SOLVER
  public :: vca_init_solver


  !> VCA DIAG CLUSTER
  public :: vca_diag_cluster


  !> VCA GF POLES UPDATE
  public :: vca_diag_system


  !> VCA SFT GRAND POTENTIAL
  public :: vca_sft_potential


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





  !+-----------------------------------------------------------------------------+!
  !PURPOSE: Get the spectrum of the system GF from the Cluster tiling 
  !+-----------------------------------------------------------------------------+!
  subroutine vca_diag_system(Vmat)
    real(8),dimension(Nlat*Nspin*Norb,Nlat*Nspin*Norb) :: Vmat
    integer                                            :: ilat,jlat
    integer                                            :: iorb,jorb
    integer                                            :: ispin
    integer                                            :: iexc,jexc
    integer                                            :: i,j,Nexc
    real(8),dimension(:,:),allocatable                 :: Mmat
    real(8),dimension(:),allocatable                   :: Lvec
    !
    Nexc = Qcluster%Nexc
    call allocate_Qmatrix(Qsystem)
    !
    allocate(Mmat(Nexc,Nexc))
    allocate(Lvec(Nexc))
    !
    Mmat = matmul(Qcluster%cdg, matmul(Vmat,Qcluster%c))
    Mmat = Mmat + diag( Qcluster%poles )
    !
    call eigh(Mmat,Lvec)
    !
    Qsystem%c     = matmul( Qcluster%c, Mmat )              ![Nlos,Nexc][Nexc,Nexc]
    Qsystem%cdg   = matmul( transpose(Mmat), Qcluster%cdg)  ![Nexc,Nexc][Nexc,Nlos]
    Qsystem%poles = Lvec
    !
  end subroutine vca_diag_system




  !+-----------------------------------------------------------------------------+!
  !PURPOSE: Get the SFT grand potential
  !+-----------------------------------------------------------------------------+!
  subroutine vca_sft_potential(sft_potential)
    real(8) :: sft_potential
    real(8) :: Tr(2)
    integer :: unit,iexc,Nexc
    Nexc = Qcluster%Nexc    
    Tr=0d0
    do iexc=1,Nexc
       Tr(1)  = Tr(1) - log(1d0+exp(-beta*Qsystem%poles(iexc)))
       Tr(2)  = Tr(2) - log(1d0+exp(-beta*Qcluster%poles(iexc)))
    enddo
    sft_potential = omega_potential + Tr(1)/beta - Tr(2)/beta
    open(free_unit(unit),file="SFT_potential.vca",access='append')
    write(unit,*)sft_potential
    close(unit)
  end subroutine vca_sft_potential

end module VCA_MAIN
