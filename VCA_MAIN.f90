module VCA_MAIN
  USE VCA_VARS_GLOBAL
  USE VCA_AUX_FUNX
  USE VCA_SETUP
  USE VCA_DIAG
  USE VCA_OBSERVABLES
  USE VCA_GREENS_FUNCTIONS
  USE VCA_BATH_SETUP
  USE VCA_EIGENSPACE
  USE VCA_HAMILTONIAN_MATVEC
  !
  USE SF_LINALG,  only: eigh,diag,kron,eye
  USE SF_IOTOOLS, only: str,free_unit
  USE SF_TIMER,only: start_timer,stop_timer
  USE SF_MISC, only: assert_shape,sort_array
  implicit none
  private


  !>INIT VCA SOLVER
  public :: vca_init_solver


  !>VCA SOLVER
  public :: vca_solve




contains



  !+-----------------------------------------------------------------------------+!
  ! PURPOSE: allocate and initialize one or multiple baths -+!
  !+-----------------------------------------------------------------------------+!
  subroutine vca_init_solver(bath)
    real(8),intent(inout),optional :: bath(:)
    logical,save                   :: isetup=.true.
    !
    write(LOGfile,"(A)")"INIT SOLVER FOR "//trim(file_suffix)
    !
    !Init Structure & memory
    if(isetup)call init_cluster_structure()
    !
    if(present(bath))then
       if(.not.check_bath_dimension(bath))stop "VCA_INIT_SOLVER error: wrong bath dimensions"
       bath=0d0
       call vca_allocate_bath(vca_bath)
       call vca_init_bath(vca_bath)
       call vca_get_bath(vca_bath,bath)
    endif
    !
    if(isetup)call setup_pointers_normal
    !
    if(vca_bath%status)call vca_deallocate_bath(vca_bath)
    !
    isetup=.false.
    !
  end subroutine vca_init_solver




  !+-----------------------------------------------------------------------------+!
  !PURPOSE: Diag the cluster, reference system
  !+-----------------------------------------------------------------------------+!
  subroutine vca_solve(Hloc,Vq,bath)
    complex(8),intent(in),dimension(:,:,:,:,:,:) :: Hloc ![Nlat,Nlat,Nspin,Nspin,Norb,Norb]
    real(8),intent(in),dimension(:,:,:),optional :: Vq   ![Nlat*Nspin*Norb,Nlat*Nspin*Norb,Lk]
    real(8),intent(inout),dimension(:),optional  :: bath
    !
    integer                                      :: Lk,Nsites
    integer                                      :: ilat,jlat
    integer                                      :: iorb,jorb
    integer                                      :: ispin

    integer                                      :: i,j
    real(8)                                      :: Tr
    integer                                      :: unit

    !
    write(LOGfile,"(A)")"SOLVING VCA"
    !
    call assert_shape(Hloc,[Nlat,Nlat,Nspin,Nspin,Norb,Norb],"vca_solve","Hloc")
    ! if(present(Vq))Lk     = size(Vq,3)
    ! call assert_shape(Vq,[Nlat*Nspin*Norb,Nlat*Nspin*Norb,Lk],"vca_solve","Vq")
    !
    if(present(Bath))then
       if(.not.check_bath_dimension(bath))stop "vca_diag_solve Error: wrong bath dimensions"
       call vca_allocate_bath(vca_bath)
       call vca_set_bath(bath,vca_bath)
       call vca_write_bath(vca_bath,LOGfile)
       call vca_save_bath(vca_bath,used=.true.)
    endif
    !
    ! select case(vca_sparse_H)
    ! case (.true.)
    spHtimesV_cc => spMatVec_cc
    ! case (.false.)
    !    spHtimesV_cc => directMatVec_cc
    ! case default
    !    stop "vca_solve_single ERROR: vca_sparse_H undefined"
    ! end select
    !
    call vca_set_Hcluster(Hloc)
    !
    call setup_eigenspace()
    !
    call diagonalize_cluster()    !find target states by digonalization of Hamiltonian
    call build_gf_cluster()       !build the one-particle Green's functions and Self-Energies
    call observables_cluster()    !obtain impurity observables as thermal averages.
    !
    call delete_eigenspace()
    call es_delete_espace(state_list)
    nullify(spHtimesV_cc)


    sft_potential = omega_potential 
    open(free_unit(unit),file="SFT_potential.vca",position='append')
    write(unit,*)sft_potential
    close(unit)
    !
    if(vca_bath%status)call vca_deallocate_bath(vca_bath)
    !
  end subroutine vca_solve

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
