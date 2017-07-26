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
    logical,save                    :: isetup=.true.
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
    call build_QL_cluster()            !build the one-particle Q and Lambda matrices, storing spectral sums
    !
    call delete_eigenspace()
  end subroutine vca_diag






  subroutine vca_update_poles()
    integer                            :: Nexc
    integer                            :: ilat,jlat
    integer                            :: iorb,jorb
    integer                            :: ispin
    integer                            :: iexc,jexc
    integer                            :: i
    real(8),dimension(:,:),allocatable :: Mmat

    call vca_get_Nexc(Nexc)
    print*,Nexc

    allocate(Mmat(Nexc,Nexc))

    Mmat = 0d0
    do iexc=1,Nexc
       do jexc=1,Nexc
          !
          do ispin=1,Nspin
             do ilat=1,Nlat
                do jlat=1,Nlat
                   do iorb=1,Norb
                      do jorb=1,Norb

                      enddo
                   enddo
                enddo
             enddo
          enddo
          !
       enddo
    enddo

  end subroutine vca_update_poles

  
end module VCA_MAIN
