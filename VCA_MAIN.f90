module VCA_MAIN
  USE VCA_VARS_GLOBAL
  USE VCA_AUX_FUNX
  USE VCA_SETUP
  USE VCA_DIAG
  USE VCA_GREENS_FUNCTIONS
  USE VCA_OBSERVABLES
  !
  USE SF_LINALG
  USE SF_IOTOOLS, only: str
  USE SF_TIMER,only: start_timer,stop_timer

  implicit none
  private

  !>INIT VCA SOLVER
  public :: vca_init_solver


  !> VCA DIAG CLUSTER
  public :: vca_diag



contains



  !+-----------------------------------------------------------------------------+!
  ! PURPOSE: allocate and initialize one or multiple baths -+!
  !+-----------------------------------------------------------------------------+!
  subroutine vca_init_solver(Hloc)
    real(8),intent(in)              :: Hloc(Nlat,Nlat,Nspin,Nspin,Norb,Norb)
    logical,save                       :: isetup=.true.
    !
    write(LOGfile,"(A)")"INIT SOLVER FOR "//trim(file_suffix)
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
  subroutine vca_diag(Hloc)
    real(8),optional,intent(in)  :: Hloc(Nlat,Nlat,Nspin,Nspin,Norb,Norb)
    !
    if(present(Hloc))call set_Hcluster(Hloc)
    !
    call setup_eigenspace()
    !
    !SOLVE THE QUANTUM IMPURITY PROBLEM:
    call diagonalize_cluster()         !find target states by digonalization of Hamiltonian
    call observables_cluster()         !obtain impurity observables as thermal averages.  
    call buildgf_cluster()             !build the one-particle impurity Green's functions  & Self-energy
    !
    call delete_eigenspace()
  end subroutine vca_diag


end module VCA_MAIN
