module VCA_MAIN
  USE VCA_VARS_GLOBAL
  USE VCA_AUX_FUNX
  USE VCA_SETUP
  USE VCA_DIAG
  USE VCA_QMATRIX
  USE VCA_OBSERVABLES
  USE VCA_BATH_SETUP
  !
  USE SF_LINALG,  only: eigh,diag
  USE SF_IOTOOLS, only: str,free_unit
  USE SF_TIMER,only: start_timer,stop_timer
  USE SF_MISC, only: assert_shape,sort_array
  implicit none
  private

  !>INIT VCA SOLVER
  public :: vca_init_solver

  !>VCA SOLVER
  public :: vca_solve


  ! !> VCA DIAG CLUSTER
  ! public :: vca_diag_cluster


  ! !> VCA GF POLES UPDATE
  ! public :: vca_diag_system


  ! !> VCA SFT GRAND POTENTIAL
  ! public :: vca_sft_potential


contains



  !+-----------------------------------------------------------------------------+!
  ! PURPOSE: allocate and initialize one or multiple baths -+!
  !+-----------------------------------------------------------------------------+!
  subroutine vca_init_solver(Hloc,bath)    
    real(8),intent(in)             :: Hloc(:,:,:,:,:,:)
    real(8),intent(inout),optional :: bath(:)
    logical,save                   :: isetup=.true.
    !
    write(LOGfile,"(A)")"INIT SOLVER FOR "//trim(file_suffix)
    !
    Nlat = size(Hloc,1)
    call assert_shape(Hloc,[Nlat,Nlat,Nspin,Nspin,Norb,Norb],"vca_init_solver","Hloc")
    !
    if(present(bath))then
       if(.not.check_bath_dimension(bath))stop "vca_init_solver error: wrong bath dimensions"
       bath=0d0
       call allocate_vca_bath(vca_bath)
       call init_vca_bath(vca_bath)
       call get_vca_bath(vca_bath,bath)
    endif
    !
    !Init Structure & memory
    if(isetup)call init_cluster_structure()
    !
    !Init Hcluster
    call set_Hcluster(Hloc)
    !
    if(isetup)call setup_pointers_normal
    !
    if(vca_bath%status)call deallocate_vca_bath(vca_bath)
    !
    isetup=.false.
    !
  end subroutine vca_init_solver




  !+-----------------------------------------------------------------------------+!
  !PURPOSE: Diag the cluster, reference system
  !+-----------------------------------------------------------------------------+!
  subroutine vca_solve(Hloc,Vmat,bath)
    real(8),optional,intent(in)                        :: Hloc(:,:,:,:,:,:)
    real(8),dimension(Nlat*Nspin*Norb,Nlat*Nspin*Norb) :: Vmat
    real(8),intent(inout),optional                     :: bath(:)
    !
    integer                                            :: ilat,jlat
    integer                                            :: iorb,jorb
    integer                                            :: ispin
    integer                                            :: iexc,jexc
    integer                                            :: i,j,Nexc
    complex(8),dimension(:,:),allocatable              :: Mmat
    real(8),dimension(:),allocatable                   :: Lvec
    real(8)                                            :: earg
    !
    real(8)                                            :: Tr,arg1,arg2
    integer                                            :: unit
    !
    Nlat = size(Hloc,1)
    call assert_shape(Hloc,[Nlat,Nlat,Nspin,Nspin,Norb,Norb],"vca_diag","Hloc")
    !
    if(present(Bath))then
       if(.not.check_bath_dimension(bath))stop "vca_diag_solve Error: wrong bath dimensions"
       call allocate_vca_bath(vca_bath)
       call set_vca_bath(bath,vca_bath)
       call write_vca_bath(vca_bath,LOGfile)
       call save_vca_bath(vca_bath,used=.true.)
    endif
    !
    if(present(Hloc))call set_Hcluster(Hloc)
    !
    call setup_eigenspace()
    !
    call diagonalize_cluster()    !find target states by digonalization of Hamiltonian
    call observables_cluster()    !obtain impurity observables as thermal averages.  
    call build_Qmatrix_cluster()  !build the one-particle Q and Lambda matrices
    !
    call delete_eigenspace()
    !
    if(vca_bath%status)call deallocate_vca_bath(vca_bath)
    !
    Nexc = Qcluster%Nexc
    call allocate_Qmatrix(Qsystem)
    !
    allocate(Mmat(Nexc,Nexc))
    allocate(Lvec(Nexc))
    !   
    Mmat = diag( Qcluster%poles ) + matmul(Qcluster%cdg, matmul(Vmat,Qcluster%c))
    !
    call eigh(Mmat,Lvec)!>from here on M-->S
    !
    Qsystem%c     = matmul( Qcluster%c, Mmat )                     ![Nlos,Nexc][Nexc,Nexc]
    Qsystem%cdg   = matmul( conjg(transpose(Mmat)), Qcluster%cdg)  ![Nexc,Nexc][Nexc,Nlos]
    Qsystem%poles = Lvec
    call sort_array(Qcluster%poles)
    !
    Nexc = Qcluster%Nexc    
    Tr=0d0
    do iexc=1,Nexc
       arg1 = beta*Qcluster%poles(iexc)
       arg2 = beta*Qsystem%poles(iexc)
       if(arg1 < -20d0 .OR. arg2 < -20d0)then
          Tr  = Tr - beta*( Qcluster%poles(iexc) - Qsystem%poles(iexc) )
       else
          Tr  = Tr + log( (1d0 + exp(-beta*Qcluster%poles(iexc)))/(1d0 + exp(-beta*Qsystem%poles(iexc))) )
       endif
    enddo
    sft_potential = omega_potential + Tr(3)/beta
    open(free_unit(unit),file="SFT_potential.vca",access='append')
    write(unit,*)sft_potential
    close(unit)
  end subroutine vca_solve

end module VCA_MAIN



