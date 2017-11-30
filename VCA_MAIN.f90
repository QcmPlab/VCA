module VCA_MAIN
  USE VCA_VARS_GLOBAL
  USE VCA_AUX_FUNX
  USE VCA_SETUP
  USE VCA_DIAG
  USE VCA_QMATRIX
  USE VCA_OBSERVABLES
  USE VCA_BATH_SETUP
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
  subroutine vca_solve(Hloc,Vmat,bath)
    real(8),intent(in)                    :: Hloc(:,:,:,:,:,:)
    real(8),dimension(:,:)                :: Vmat
    real(8),intent(inout),optional        :: bath(:)
    !
    integer                               :: ilat,jlat
    integer                               :: iorb,jorb
    integer                               :: ispin
    integer                               :: iexc,jexc
    integer                               :: i,j,Nexc,Nvmat,Nexc_sys,Nc
    integer :: icopy,jcopy,i1,i2,j1,j2,ii1,ii2,jj1,jj2
    complex(8),dimension(:,:),allocatable :: Mmat
    real(8),dimension(:),allocatable      :: Lvec
    real(8)                               :: earg
    !
    real(8)                               :: Tr,arg1,arg2
    integer                               :: unit
    !
    write(LOGfile,"(A)")"SOLVING VCA"
    !
    Nlat = size(Hloc,1)
    call assert_shape(Hloc,[Nlat,Nlat,Norb,Norb,Nspin,Nspin],"vca_solve","Hloc")
    Nvmat = vca_get_system_dimension(present(bath))
    call assert_shape(Vmat,[Nvmat,Nvmat],"vca_solve","Vmat")
    !
    if(present(Bath))then
       if(.not.check_bath_dimension(bath))stop "vca_diag_solve Error: wrong bath dimensions"
       call vca_allocate_bath(vca_bath)
       call vca_set_bath(bath,vca_bath)
       call vca_write_bath(vca_bath,LOGfile)
       call vca_save_bath(vca_bath,used=.true.)
    endif
    !
    call vca_set_Hcluster(Hloc)
    !
    call setup_eigenspace()
    !
    call diagonalize_cluster()    !find target states by digonalization of Hamiltonian
    call observables_cluster()    !obtain impurity observables as thermal averages.
    call build_Qmatrix_cluster()  !build the one-particle Q and Lambda matrices
    !
    call delete_eigenspace()
    !
    Nexc     = Qcluster%Nexc
    Nexc_sys = Ncopies*Nexc
    call allocate_Qmatrix(Qsystem)
    !    
    allocate(Mmat(Nexc_sys,Nexc_sys))
    allocate(Lvec(Nexc_sys));Lvec=zero
    !
    Nc = vca_get_cluster_dimension(present(bath))
    !
    Mmat=zero
    !
    do icopy=1,Ncopies
       do jcopy=1,Ncopies
          i1 = 1 + (icopy-1)*Nc
          i2 = icopy*Nc
          j1 = 1 + (jcopy-1)*Nc
          j2 = jcopy*Nc
          !
          ii1 = 1 + (icopy-1)*Nexc
          ii2 = icopy*Nexc
          jj1 = 1 + (jcopy-1)*Nexc
          jj2 = jcopy*Nexc
          !
          Mmat(ii1:ii2,jj1:jj2)  = matmul(Qcluster%cdg,  matmul(Vmat(i1:i2,j1:j2),Qcluster%c))
       enddo
    enddo
    !
    Mmat = Mmat + kron(eye(Ncopies),diag( Qcluster%poles ))
    !
    call eigh(Mmat,Lvec)!>from here on M-->S
    !
    !>
    ! Qsystem%c     = matmul( Qcluster%c, Mmat )                     ![Nlos,Nexc][Nexc,Nexc]
    ! Qsystem%cdg   = matmul( conjg(transpose(Mmat)), Qcluster%cdg)  ![Nexc,Nexc][Nexc,Nlos]
    Qsystem%poles = Lvec
    call sort_array(Qcluster%poles)
    !
    Tr=0d0
    do iexc=1,Nexc
       do icopy=1,Ncopies
          i = icopy + (iexc-1)*Ncopies
          !
          arg1 = beta*Qcluster%poles(iexc)
          arg2 = beta*Qsystem%poles(i)
          if(arg1 < -20d0 .OR. arg2 < -20d0)then
             Tr  = Tr - beta*( Qcluster%poles(iexc) - Qsystem%poles(i) )
          else
             Tr  = Tr + log( (1d0 + exp(-beta*Qcluster%poles(iexc)))/(1d0 + exp(-beta*Qsystem%poles(i))) )
          endif
       enddo
    enddo
    sft_potential = omega_potential + Tr/beta - 1d0/beta*log(dble(Ncopies))
    open(free_unit(unit),file="SFT_potential.vca",position='append')
    write(unit,*)sft_potential
    close(unit)
    !
    if(vca_bath%status)call vca_deallocate_bath(vca_bath)
    !
  end subroutine vca_solve

end module VCA_MAIN



