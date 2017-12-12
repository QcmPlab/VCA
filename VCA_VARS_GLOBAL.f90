MODULE VCA_VARS_GLOBAL
  USE VCA_INPUT_VARS
  !
  USE SF_CONSTANTS
  USE SF_ARRAYS,  only: arange,linspace
  USE SF_IOTOOLS, only:str
  implicit none


  !-------------------- EFFECTIVE BATH STRUCTURE ----------------------!
  type effective_bath
     real(8),dimension(:,:,:,:),allocatable        :: e     !local energies [1//Nlat][1//Norb][Nspin][Nbath]
     real(8),dimension(:,:,:,:),allocatable        :: v     !spin-keep hyb. [Nlat][Norb][Nspin][Nbath]
     logical                                       :: status=.false.
  end type effective_bath


  !LOG UNITS
  !=========================================================
  integer,save                                    :: LOGfile=6


  !Ns       =              Number of levels per spin
  !Nlevels  = 2*Ns       = Total Number  of levels
  !Nsectors =              Number of sectors
  !=========================================================
  integer                                         :: Ns
  integer                                         :: Nlevels
  integer                                         :: Nsectors


  !non-interacting cluster Hamiltonian
  !=========================================================
  complex(8),dimension(:,:,:,:,:,:),allocatable   :: impHloc ![Nlat][Nlat][Nspin][Nspin][Norb][Norb]



  !Some maps between sectors and full Hilbert space (pointers)
  !=========================================================
  integer,allocatable,dimension(:,:)              :: getsector
  integer,allocatable,dimension(:,:)              :: getCsector
  integer,allocatable,dimension(:,:)              :: getCDGsector
  integer,allocatable,dimension(:)                :: getDim,getDimUp,getDimDw
  integer,allocatable,dimension(:)                :: getNup,getNdw
  integer,allocatable,dimension(:,:,:)            :: getBathStride


  !Effective Bath used in the VCA code (this is opaque to user)
  !=========================================================
  type(effective_bath)                            :: vca_bath


  !Eigenvalues,Eigenvectors FULL DIAGONALIZATION
  !Hamiltonian eig-space structure
  !=========================================================
  type full_espace
     real(8),dimension(:),pointer                 :: e
     complex(8),dimension(:,:),pointer            :: M
  end type full_espace
  type(full_espace),dimension(:),allocatable      :: espace


  !Partition function, Omega potential, SFT potential
  !=========================================================
  real(8)                                         :: zeta_function
  real(8)                                         :: omega_potential
  real(8)                                         :: sft_potential


  !Cluster Green's functions
  !(Nlat,Nlat,Nspin,Nspin,Norb,Norb,:)
  !=========================================================
  complex(8),allocatable,dimension(:,:,:,:,:,:,:) :: impGmats ![Nlat][Nlat][Nspin][Nspin][Norb][Norb][L]
  complex(8),allocatable,dimension(:,:,:,:,:,:,:) :: impGreal ![Nlat][Nlat][Nspin][Nspin][Norb][Norb][L]
  !
  complex(8),allocatable,dimension(:,:,:,:,:,:,:) :: impG0mats ![Nlat][Nlat][Nspin][Nspin][Norb][Norb][L]
  complex(8),allocatable,dimension(:,:,:,:,:,:,:) :: impG0real ![Nlat][Nlat][Nspin][Nspin][Norb][Norb][L]
  !
  complex(8),allocatable,dimension(:,:,:,:,:,:,:) :: impSmats ![Nlat][Nlat][Nspin][Nspin][Norb][Norb][L]
  complex(8),allocatable,dimension(:,:,:,:,:,:,:) :: impSreal ![Nlat][Nlat][Nspin][Nspin][Norb][Norb][L]



  !Cluster local observables:
  !=========================================================
  real(8),dimension(:,:),allocatable              ::  imp_dens    ![Nlat][Norb]
  real(8),dimension(:,:),allocatable              ::  imp_dens_up ![Nlat][Norb]
  real(8),dimension(:,:),allocatable              ::  imp_dens_dw ![Nlat][Norb]
  real(8),dimension(:,:),allocatable              ::  imp_docc    ![Nlat][Norb]



  !Suffix string attached to the output files.
  !=========================================================
  character(len=64)                               :: file_suffix=""


  !SECTOR-TO-FOCK SPACE STRUCTURE
  !=========================================================
  type sector_map
     integer,dimension(:),allocatable             :: map
  end type sector_map


  !Frequency and time arrays:
  !=========================================================
  real(8),dimension(:),allocatable                :: wm,tau,wr,vm



contains



  !> Get stride position in the one-particle many-body space 
  function index_stride_los(ilat,iorb,ispin) result(indx)
    integer :: ilat
    integer :: iorb
    integer :: ispin
    integer :: indx
    indx = iorb + (ilat-1)*Norb + (ispin-1)*Norb*Nlat
  end function index_stride_los



  !+------------------------------------------------------------------+
  !PURPOSE  : Allocate arrays and setup frequencies and times
  !+------------------------------------------------------------------+
  subroutine vca_allocate_time_freq_arrays
    integer :: i
    call vca_deallocate_time_freq_arrays()
    allocate(wm(Lmats))
    allocate(vm(0:Lmats))          !bosonic frequencies
    allocate(wr(Lreal))
    allocate(tau(0:Ltau))
    wm     = pi/beta*(2*arange(1,Lmats)-1)
    do i=0,Lmats
       vm(i) = pi/beta*2*dble(i)
    enddo
    wr     = linspace(wini,wfin,Lreal)
    tau(0:)= linspace(0d0,beta,Ltau+1)
  end subroutine vca_allocate_time_freq_arrays
  !
  subroutine vca_deallocate_time_freq_arrays
    if(allocated(wm))deallocate(wm)
    if(allocated(vm))deallocate(vm)
    if(allocated(tau))deallocate(tau)
    if(allocated(wr))deallocate(wr)
  end subroutine vca_deallocate_time_freq_arrays


END MODULE VCA_VARS_GLOBAL
