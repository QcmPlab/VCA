module VCA_MAIN
  USE VCA_VARS_GLOBAL
  USE VCA_AUX_FUNX
  USE VCA_SETUP
  USE VCA_DIAG
  USE VCA_OBSERVABLES
  USE VCA_GREENS_FUNCTIONS
  USE VCA_BATH_SETUP
  USE VCA_EIGENSPACE
  USE VCA_HAMILTONIAN
  USE VCA_OMEGA
  !
  USE SF_LINALG,  only: eigh,diag,kron,eye
  USE SF_IOTOOLS, only: str,free_unit
  USE SF_TIMER,only: start_timer,stop_timer
  USE SF_MISC, only: assert_shape,sort_array
  implicit none
  private


  !>INIT VCA SOLVER
  interface vca_init_solver
     module procedure :: vca_init_solver_serial
#ifdef _MPI
     module procedure :: vca_init_solver_mpi
#endif
  end interface vca_init_solver
  !>
  public :: vca_init_solver

  !>VCA SOLVER
  interface vca_solve
     module procedure :: vca_solve_serial
#ifdef _MPI
     module procedure :: vca_solve_mpi
#endif
  end interface vca_solve
  !>
  public :: vca_solve




contains



  !+-----------------------------------------------------------------------------+!
  ! PURPOSE: allocate and initialize one or multiple baths -+!
  !+-----------------------------------------------------------------------------+!
  subroutine vca_init_solver_serial(bath_h,bath_v)
    complex(8),intent(inout),optional :: bath_h(:,:,:,:,:,:)
    complex(8),intent(inout),optional :: bath_v(:,:,:,:,:,:)
    logical,save                      :: isetup=.true.
    !
    write(LOGfile,"(A)")"INIT SOLVER FOR "//trim(file_suffix)
    !
    if(present(bath_h).and.present(bath_v))then
       if(Nlat_bath<1 .or. Norb_bath<1)stop "VCA_INIT_SOLVER error: wrong bath dimensions"
       call vca_allocate_bath(vca_bath)
       call vca_init_bath(vca_bath)
    else
      write(LOGfile,"(A)") "Bath not present, setting Nbath to 0"
      Nlat_bath=0
      Norb_bath=0
    endif
    !Init Structure & memory
    if(isetup)call init_cluster_structure()
    !
    if(isetup)call setup_global
    !
    if(vca_bath%status)call vca_deallocate_bath(vca_bath)
    !
    isetup=.false.
    !
  end subroutine vca_init_solver_serial

#ifdef _MPI
  subroutine vca_init_solver_mpi(MpiComm,bath_h,bath_v)
    integer                                     :: MpiComm
    complex(8),intent(inout),optional           :: bath_h(:,:,:,:,:,:)
    complex(8),intent(inout),optional           :: bath_v(:,:,:,:,:,:)
    logical                                     :: check 
    logical,save                                :: isetup=.true.
    integer                                     :: i
    !
    call vca_set_mpicomm(MpiComm)
    !
    write(LOGfile,"(A)")"INIT SOLVER FOR "//trim(file_suffix)
    !
    if(present(bath_h).and.present(bath_v))then
       if(Nlat_bath<1 .or. Norb_bath<1)stop "VCA_INIT_SOLVER error: wrong bath dimensions"
       call vca_allocate_bath(vca_bath)
       call vca_init_bath(vca_bath)
    else
      write(LOGfile,"(A)") "Bath not present, setting Nbath to 0"
      Nlat_bath=0
      Norb_bath=0
    endif
    !Init Structure & memory
    if(isetup)call init_cluster_structure()
    !
    !Init bath:
    !call set_hloc(Hloc)
    !
    !
    if(isetup)call setup_global
    !
    if(vca_bath%status)call vca_deallocate_bath(vca_bath)
    !
    isetup=.false.
    !
    call vca_del_MpiComm()
    !
  end subroutine vca_init_solver_mpi
#endif




  !+-----------------------------------------------------------------------------+!
  !PURPOSE: Diag the cluster, reference system
  !+-----------------------------------------------------------------------------+!
  subroutine vca_solve_serial(Hloc,Hk,bath_h,bath_v)
    complex(8),intent(in),dimension(:,:,:,:,:,:)   :: Hloc ![Nlat,Nlat,Nspin,Nspin,Norb,Norb]
    complex(8),intent(in),dimension(:,:,:,:,:,:,:) :: Hk ![Nlat,Nlat,Nspin,Nspin,Norb,Norb,Nktot]
    complex(8),intent(inout),optional              :: bath_h(:,:,:,:,:,:)
    complex(8),intent(inout),optional              :: bath_v(:,:,:,:,:,:)
    !
    integer                                        :: Lk,Nsites
    integer                                        :: ilat,jlat
    integer                                        :: iorb,jorb
    integer                                        :: ispin

    integer                                        :: i,j
    real(8)                                        :: Tr,omega_integral,omegaprime
    integer                                        :: unit

    !
    write(LOGfile,"(A)")"SOLVING VCA"
    !
    if(rank(Hloc) .ne. 6) STOP "STOP: wrong cluster matrix dimensions"
    if(rank(Hk)   .ne. 7) STOP "STOP: wrong lattice matrix dimensions"
    !
    if(present(bath_h).and.present(bath_v))then
       call assert_shape(bath_h,[Nlat_bath,Nlat_bath,Nspin,Nspin,Norb_bath,Norb_bath],"vca_solve","bath_h")
       call assert_shape(bath_v,[Nlat     ,Nlat_bath,Nspin,Nspin,Norb     ,Norb_bath],"vca_solve","bath_h")
       call vca_allocate_bath(vca_bath)
       call vca_set_bath(bath_h,bath_v,vca_bath)
    endif
    !
    !
    select case(vca_sparse_H)
      case (.true.)
        spHtimesV_p => spMatVec_main
      case (.false.)
        spHtimesV_p => directMatVec_main
      case default
        stop "vca_solve_single ERROR: vca_sparse_H undefined"
    end select
    !
    !
    !GENERATE THE CLUSTER HAMILTONIAN AND THE HOPPING MATRIX FOR THE LATTICE
    !
    allocate(impHloc(Nlat,Nlat,Nspin,Nspin,Norb,Norb))
    allocate(impHk(Nlat,Nlat,Nspin,Nspin,Norb,Norb,size(Hk,7)))
    impHloc=zero
    impHk=zero
    call vca_set_Hcluster(Hloc)
    call vca_set_Hk(Hk)
    !
    !GET CLUSTER GREEN'S FUNCTION AND GROUND STATE ENERGY
    !
    call diagonalize_cluster()    !find target states by digonalization of Hamiltonian
    call build_gf_cluster()       !build the one-particle Green's functions and Self-Energies
    if(print_observables)then
      call observables_cluster()    !obtain impurity observables as thermal averages.
      call get_custom_observables()    !obtain impurity observables as thermal averages.
    endif
    !
    !CALCULATE THE VARIATIONAL GRAND POTENTIAL
    !
    omegaprime=0.d0
    omega_integral=0.d0
    !
    if(finiteT)then
      do i=1,state_list%size
        omegaprime=omegaprime+exp(-beta*(es_return_energy(state_list,i)-state_list%emin))
      enddo
      omegaprime=state_list%emin-(1.d0/beta)*log(omegaprime)
      omega_integral=frequency_integration_finite_t()
    else
      omegaprime=state_list%emin
      omega_integral=frequency_integration()
    endif
    !
    sft_potential = omegaprime-omega_integral
    !
    write(LOGfile,"(A,10f18.12,A)")"EGS PER SITE",state_list%emin/NLAT
    write(LOGfile,"(A,10f18.12,A)")"OMEGA POTENTIAL PER SITE=",(state_list%emin-omega_integral)/NLAT
    open(free_unit(unit),file="SFT_potential.vca",position='append')
    write(unit,*)sft_potential,omegaprime,-omega_integral
    close(unit)
    !
    !CLEAN UP
    !
    call es_delete_espace(state_list)
    nullify(spHtimesV_p)
    !
    !
    deallocate(impHloc)
    deallocate(impHk)
    if(vca_bath%status)call vca_deallocate_bath(vca_bath)
    !
  end subroutine vca_solve_serial


#ifdef _MPI

  subroutine vca_solve_mpi(MpiComm,Hloc,Hk,bath_h,bath_v)
    complex(8),intent(in),dimension(:,:,:,:,:,:)   :: Hloc ![Nlat,Nlat,Nspin,Nspin,Norb,Norb]
    complex(8),intent(in),dimension(:,:,:,:,:,:,:) :: Hk   ![Nlat,Nlat,Nspin,Nspin,Norb,Norb,Nktot]
    complex(8),intent(inout),optional              :: bath_h(:,:,:,:,:,:)
    complex(8),intent(inout),optional              :: bath_v(:,:,:,:,:,:)
    !
    integer                                        :: Lk,Nsites
    integer                                        :: ilat,jlat
    integer                                        :: iorb,jorb
    integer                                        :: ispin

    integer                                        :: i,j
    real(8)                                        :: Tr,omega_integral,omegaprime
    integer                                        :: unit
    integer                                        :: MpiComm
    logical                                        :: check
    !
    !SET THE LOCAL MPI COMMUNICATOR :
    !
    call vca_set_MpiComm(MpiComm)
    !
    if(rank(Hloc) .ne. 6) STOP "STOP: wrong cluster matrix dimensions"
    if(rank(Hk)   .ne. 7) STOP "STOP: wrong lattice matrix dimensions"
    !    
    !
    !
    if(present(bath_h).and.present(bath_v))then
       call assert_shape(bath_h,[Nlat_bath,Nlat_bath,Nspin,Nspin,Norb_bath,Norb_bath],"vca_solve","bath_h")
       call assert_shape(bath_v,[Nlat     ,Nlat_bath,Nspin,Nspin,Norb     ,Norb_bath],"vca_solve","bath_h")
       call vca_allocate_bath(vca_bath)
       call vca_set_bath(bath_h,bath_v,vca_bath)
    endif
    !
    select case(vca_sparse_H)
      case (.true.)
        spHtimesV_p => spMatVec_main
      case (.false.)
        spHtimesV_p => directMatVec_main
      case default
        stop "vca_solve_single ERROR: vca_sparse_H undefined"
    end select
    !
    !GENERATE THE CLUSTER HAMILTONIAN AND THE HOPPING MATRIX FOR THE LATTICE
    !
    allocate(impHloc(Nlat,Nlat,Nspin,Nspin,Norb,Norb))
    allocate(impHk(Nlat,Nlat,Nspin,Nspin,Norb,Norb,size(Hk,7)))
    impHloc=zero
    impHk=zero
    call vca_set_Hcluster(Hloc)
    call vca_set_Hk(Hk)

    !
    !GET CLUSTER GREEN'S FUNCTION AND GROUND STATE ENERGY
    !
    call diagonalize_cluster()    !find target states by digonalization of Hamiltonian
    call build_gf_cluster()       !build the one-particle Green's functions and Self-Energies
    if(print_observables)then
      call observables_cluster()       !obtain impurity observables as thermal averages.
      call get_custom_observables()    !obtain impurity observables as thermal averages.
    endif
    !call save_gfprime("gfprime",use_formatted=.true.)
    !call read_gfprime("gfprime",use_formatted=.true.)
    !call reconstruct_g()
    !
    !CALCULATE THE VARIATIONAL GRAND POTENTIAL
    !
    omegaprime=0.d0
    omega_integral=0.d0
    !
    if(finiteT)then
      do i=1,state_list%size
        omegaprime=omegaprime+exp(-beta*(es_return_energy(state_list,i)-state_list%emin))
      enddo
      omegaprime=state_list%emin-(1.d0/beta)*log(omegaprime)
      if(MPIMASTER) omega_integral=frequency_integration_finite_t()
    else
      omegaprime=state_list%emin
      if(MPIMASTER) omega_integral=frequency_integration()
    endif
    !
    if(MPIMASTER) sft_potential = omegaprime-omega_integral
    call Bcast_MPI(MpiComm,sft_potential)
    !
    if(MPIMASTER)then
      write(LOGfile,"(A,10f18.12,A)")"EGS PER SITE",omegaprime/NLAT
      write(LOGfile,"(A,10f18.12,A)")"OMEGA POTENTIAL PER SITE=",(omegaprime-omega_integral)/NLAT
      open(free_unit(unit),file="SFT_potential.vca",position='append')
      write(unit,*)sft_potential,omegaprime,-omega_integral
      close(unit)
    endif
    !
    call barrier_mpi(MpiComm)
    !
    !CLEAN UP
    !
    call es_delete_espace(state_list)
    nullify(spHtimesV_p)
    !
    deallocate(impHloc)
    deallocate(impHk)
    if(vca_bath%status)call vca_deallocate_bath(vca_bath)
    call vca_del_MpiComm()
    !
  end subroutine vca_solve_mpi
#endif


end module VCA_MAIN


! !
