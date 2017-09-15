MODULE VCA_VARS_GLOBAL
  USE VCA_INPUT_VARS
  !
  USE SF_CONSTANTS
  USE SF_IOTOOLS, only:str
  implicit none



  !-------------------- EFFECTIVE BATH STRUCTURE ----------------------!
  type effective_bath
     real(8),dimension(:,:,:,:),allocatable        :: e     !local energies [Nlat][Nspin][1//Norb][Nbath]
     real(8),dimension(:,:,:,:),allocatable        :: v     !spin-keep hyb. [Nlat][Nspin][Norb][Nbath]
     logical                                     :: status=.false.
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

  !Other System dimension
  !=========================================================  
  integer                                         :: Nlat
  integer                                         :: Nexcitations

  !non-interacting cluster Hamiltonian
  !=========================================================
  real(8),dimension(:,:,:,:,:,:),allocatable      :: impHloc ![Nlat][Nlat][Norb][Norb][Nspin][Nspin]



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
     real(8),dimension(:,:),pointer               :: M
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
  complex(8),allocatable,dimension(:,:,:,:,:,:,:) :: impGmats ![Nlat][Nlat][Norb][Norb][Nspin][Nspin][L]
  complex(8),allocatable,dimension(:,:,:,:,:,:,:) :: impGreal ![Nlat][Nlat][Norb][Norb][Nspin][Nspin][L]


  !Q and Lambda VCA matrices:
  !=============================================== ==========
  type Qmatrix
     logical                            :: allocated
     integer                            :: Nexc
     real(8),dimension(:,:),allocatable :: c     ![Nlat*Norb*Nspin][Nexc]
     real(8),dimension(:,:),allocatable :: cdg   ![Nlat*Norb*Nspin][Nexc]
     real(8),dimension(:),allocatable   :: poles ![Nexc] diagonal matrix
  end type Qmatrix
  type(Qmatrix)                                   :: Qcluster
  type(Qmatrix)                                   :: Qsystem


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


END MODULE VCA_VARS_GLOBAL
