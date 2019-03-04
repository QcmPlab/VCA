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
  subroutine vca_init_solver_serial(bath)
    real(8),intent(inout),optional :: bath(:)
    logical,save                   :: isetup=.true.
    !
    write(LOGfile,"(A)")"INIT SOLVER FOR "//trim(file_suffix)
    !
    if(present(bath))then
       if(.not.check_bath_dimension(bath))stop "VCA_INIT_SOLVER error: wrong bath dimensions"
       bath=0d0
       call vca_allocate_bath(vca_bath)
       call vca_init_bath(vca_bath)
       call vca_get_bath(vca_bath,bath)
    else
      write(LOGfile,"(A)") "Bath not present, setting Nbath to 0"
      Nbath=0
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
  subroutine vca_init_solver_mpi(MpiComm,bath)
    integer                                     :: MpiComm
    real(8),dimension(:),intent(inout),optional :: bath
    !complex(8),intent(in)                      :: Hloc(Nspin,Nspin,Norb,Norb)
    logical                                     :: check 
    logical,save                                :: isetup=.true.
    integer                                     :: i
    logical                                     :: MPI_MASTER=.true.
    integer                                     :: MPI_RANK
    integer                                     :: MPI_ERR
    !
    MPI_RANK   = get_Rank_MPI(MpiComm)
    MPI_MASTER = get_Master_MPI(MpiComm)
    !
    write(LOGfile,"(A)")"INIT SOLVER FOR "//trim(file_suffix)
    !
    if(present(bath))then
       if(.not.check_bath_dimension(bath))stop "VCA_INIT_SOLVER error: wrong bath dimensions"
       bath=0d0
       call vca_allocate_bath(vca_bath)
       call vca_init_bath(vca_bath)
       call vca_get_bath(vca_bath,bath)
    else
      write(LOGfile,"(A)") "Bath not present, setting Nbath to 0"
      Nbath=0
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
    call MPI_Barrier(MpiComm,MPI_ERR)
    !
  end subroutine vca_init_solver_mpi
#endif




  !+-----------------------------------------------------------------------------+!
  !PURPOSE: Diag the cluster, reference system
  !+-----------------------------------------------------------------------------+!
  subroutine vca_solve_serial(Hloc,Hk,bath)
    complex(8),intent(in),dimension(:,:,:,:,:,:)   :: Hloc ![Nlat,Nlat,Nspin,Nspin,Norb,Norb]
    complex(8),intent(in),dimension(:,:,:,:,:,:,:) :: Hk ![Nlat,Nlat,Nspin,Nspin,Norb,Norb,Nkpts**Ndim]
    real(8),intent(inout),dimension(:),optional    :: bath
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
    call assert_shape(Hloc,[Nlat,Nlat,Nspin,Nspin,Norb,Norb],"vca_solve","Hloc") 
    call assert_shape(Hk,[Nlat,Nlat,Nspin,Nspin,Norb,Norb,Nkpts**Ndim],"vca_solve","Hk") 
    !
    if(present(Bath))then
       if(.not.check_bath_dimension(bath))stop "vca_diag_solve Error: wrong bath dimensions"
       call vca_allocate_bath(vca_bath)
       call vca_set_bath(bath,vca_bath)
       call vca_write_bath(vca_bath,LOGfile)
       call vca_save_bath(vca_bath,used=.true.)
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
    !
    !GENERATE THE CLUSTER HAMILTONIAN AND THE HOPPING MATRIX FOR THE LATTICE
    !
    call vca_set_Hcluster(Hloc)
    call vca_set_Hk(Hk)
    !call embed_hcluster(Hloc)
    !call embed_hk(Hk)
    !
    !GET CLUSTER GREEN'S FUNCTION AND GROUND STATE ENERGY
    !
    call diagonalize_cluster()    !find target states by digonalization of Hamiltonian
    call build_gf_cluster()       !build the one-particle Green's functions and Self-Energies
    call observables_cluster()    !obtain impurity observables as thermal averages.
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
    write(unit,*)sft_potential
    close(unit)
    !
    !CLEAN UP
    !
    call es_delete_espace(state_list)
    nullify(spHtimesV_p)
    !
    open(free_unit(unit),file="SFT_potential.vca",position='append')
    write(unit,*)sft_potential
    close(unit)
    !
    if(vca_bath%status)call vca_deallocate_bath(vca_bath)
    !
  end subroutine vca_solve_serial


#ifdef _MPI

  subroutine vca_solve_mpi(MpiComm,Hloc,Hk,bath)
    complex(8),intent(in),dimension(:,:,:,:,:,:)   :: Hloc ![Nlat,Nlat,Nspin,Nspin,Norb,Norb]
    complex(8),intent(in),dimension(:,:,:,:,:,:,:) :: Hk   ![Nlat,Nlat,Nspin,Nspin,Norb,Norb,Nkpts**ndim]
    real(8),intent(inout),dimension(:),optional    :: bath
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
    logical                                        :: MPI_MASTER=.true.
    !
    MPI_MASTER = get_Master_MPI(MpiComm)
    !
    if(MPI_MASTER)write(LOGfile,"(A)")"SOLVING VCA"
    !
    call assert_shape(Hloc,[Nlat,Nlat,Nspin,Nspin,Norb,Norb],"vca_solve","Hloc")
    call assert_shape(Hk,[Nlat,Nlat,Nspin,Nspin,Norb,Norb,Nkpts**Ndim],"vca_solve","Hk")
    !
    if(present(Bath))then
       if(.not.check_bath_dimension(bath))stop "vca_diag_solve Error: wrong bath dimensions"
       call vca_allocate_bath(vca_bath)
       call vca_set_bath(bath,vca_bath)
       call vca_write_bath(vca_bath,LOGfile)
       call vca_save_bath(vca_bath,used=.true.)
    endif
    !SET THE LOCAL MPI COMMUNICATOR :
    call vca_set_MpiComm(MpiComm)
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
    call vca_set_Hcluster(Hloc)
    call vca_set_Hk(Hk)
    !call embed_hcluster(Hloc)
    !call embed_hk(Hk)
    !
    !GET CLUSTER GREEN'S FUNCTION AND GROUND STATE ENERGY
    !
    call diagonalize_cluster()    !find target states by digonalization of Hamiltonian
    call build_gf_cluster()       !build the one-particle Green's functions and Self-Energies
    call observables_cluster()    !obtain impurity observables as thermal averages.

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
      omega_integral=frequency_integration_finite_t()
    else
      omegaprime=state_list%emin
      omega_integral=frequency_integration()
    endif
    !
    sft_potential = omegaprime-omega_integral
    !
    if(MPI_MASTER) then
        write(LOGfile,"(A,10f18.12,A)")"EGS PER SITE",omegaprime/NLAT
        write(LOGfile,"(A,10f18.12,A)")"OMEGA POTENTIAL PER SITE=",(omegaprime-omega_integral)/NLAT
        open(free_unit(unit),file="SFT_potential.vca",position='append')
        write(unit,*)sft_potential,omegaprime,-omega_integral
        close(unit)
    endif
    !
    !CLEAN UP
    !
    call es_delete_espace(state_list)
    nullify(spHtimesV_p)
    !
    if(vca_bath%status)call vca_deallocate_bath(vca_bath)
    call vca_del_MpiComm()
    !
  end subroutine vca_solve_mpi
#endif


end module VCA_MAIN








! Nexc     = Qcluster%Nexc
! Nexc_sys = Ncopies*Nexc
! call allocate_Qmatrix(Qsystem)
! !
! print*,"Nexc_sys=",Nexc_sys
! allocate(Mmat(Nexc_sys,Nexc_sys))
! allocate(Lvec(Nexc_sys));Lvec=zero
! !
! Nc = vca_get_cluster_dimension(present(bath))
! !
! Mmat=zero
! !
! do icopy=1,Ncopies
!    do jcopy=1,Ncopies
!       i1 = 1 + (icopy-1)*Nc
!       i2 = icopy*Nc
!       j1 = 1 + (jcopy-1)*Nc
!       j2 = jcopy*Nc
!       !
!       ii1 = 1 + (icopy-1)*Nexc
!       ii2 = icopy*Nexc
!       jj1 = 1 + (jcopy-1)*Nexc
!       jj2 = jcopy*Nexc
!       !
!       Mmat(ii1:ii2,jj1:jj2)  = matmul(Qcluster%cdg,  matmul(Vmat(i1:i2,j1:j2),Qcluster%c))
!    enddo
! enddo
! !
! Mmat = Mmat + kron(eye(Ncopies),diag( Qcluster%poles ))
! !
! call eigh(Mmat,Lvec)!>from here on M-->S
! !
! Qsystem%poles = Lvec
! call sort_array(Qcluster%poles)
! Lvec=zero
! do iexc=1,Nexc
!    do icopy=1,Ncopies
!       i = icopy + (iexc-1)*Ncopies
!       Lvec(i) = Qcluster%poles(iexc)
!    enddo
! enddo
! !
! Tr=0d0
! do i=1,Nexc_sys
!    !
!    arg1 = beta*Lvec(i)
!    arg2 = beta*Qsystem%poles(i)
!    if(arg1 < -20d0 .OR. arg2 < -20d0)then
!       Tr  = Tr - beta*( Lvec(i) - Qsystem%poles(i) )
!    else
!       Tr  = Tr + log( (1d0 + exp(-beta*Lvec(i)))/(1d0 + exp(-beta*Qsystem%poles(i))) )
!    endif
! enddo
! !
