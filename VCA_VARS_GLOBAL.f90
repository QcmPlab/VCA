MODULE VCA_VARS_GLOBAL
  USE VCA_INPUT_VARS
  USE VCA_SPARSE_MATRIX
  !
  USE SF_CONSTANTS
  USE SF_ARRAYS,  only: arange,linspace
  USE SF_IOTOOLS, only: str
#ifdef _MPI
  USE MPI
  USE SF_MPI
#endif
  implicit none


  !-------------------- EFFECTIVE BATH STRUCTURE ----------------------!
  type effective_bath
     real(8),dimension(:,:,:,:),allocatable        :: e     !local energies [1//Nlat][1//Norb][Nspin][Nbath]
     real(8),dimension(:,:,:,:),allocatable        :: v     !spin-keep hyb. [Nlat][Norb][Nspin][Nbath]
     logical                                       :: status=.false.
  end type effective_bath

  !------------------ FULL HAMILTONIAN STRUCTURE ---------------------!
  type full_espace
     real(8),dimension(:),pointer      :: e
     complex(8),dimension(:,:),pointer :: M
  end type full_espace


  !---------------- SECTOR-TO-FOCK SPACE STRUCTURE -------------------!
  type sector_map
     integer,dimension(:),allocatable :: map
     logical                          :: status=.false.
  end type sector_map

  interface map_allocate
     module procedure :: map_allocate_scalar
     module procedure :: map_allocate_vector
  end interface map_allocate

  interface map_deallocate
     module procedure :: map_deallocate_scalar
     module procedure :: map_deallocate_vector
  end interface map_deallocate

  !-------------- GMATRIX FOR FAST EVALUATION OF GF ------------------!
  !note that we use a single Qmatrix here which must be intended as
  !the product Q^+.Q, corresponding to the weights of the GF, and a
  !component corresponding to the poles. 
  type GFspectrum
     real(8),dimension(:),allocatable :: weight
     real(8),dimension(:),allocatable :: poles
  end type GFspectrum
  type GFmatrix
     type(GFspectrum),dimension(:),allocatable :: channel
  end type GFmatrix
  interface GFmatrix_allocate
     module procedure :: allocate_GFmatrix_Nchan
     module procedure :: allocate_GFmatrix_Nexc
  end interface GFmatrix_allocate


  !------------------ ABTRACT INTERFACES PROCEDURES ------------------!
  !SPARSE MATRIX-VECTOR PRODUCTS USED IN ED_MATVEC
  !dbleMat*dbleVec
  abstract interface
     subroutine dd_sparse_HxV(Nloc,v,Hv)
       integer                 :: Nloc
       real(8),dimension(Nloc) :: v
       real(8),dimension(Nloc) :: Hv
     end subroutine dd_sparse_HxV
  end interface




 !-------------------------- ED  VARIABLES --------------------------!

  !SIZE OF THE PROBLEM
  !=========================================================
  integer,save                                       :: Ns       !Number of levels per spin
  integer,save                                       :: Nsectors !Number of sectors
  integer,save                                       :: Ns_orb
  integer,save                                       :: Ns_ud


  !non-interacting cluster Hamiltonian and full system hopping matrix
  !=========================================================
  complex(8),dimension(:,:,:,:,:,:),allocatable   :: impHloc ![Nlat][Nlat][Nspin][Nspin][Norb][Norb]
  complex(8),dimension(:,:,:,:,:,:,:),allocatable   :: impHk   ![Nlat][Nlat][Nspin][Nspin][Norb][Norb]



  !Some maps between sectors and full Hilbert space (pointers)  CHECK!!!!!!!!!!!!!!!!!!!!!
  !=========================================================
  integer,allocatable,dimension(:,:)              :: getsector
  integer,allocatable,dimension(:,:,:)            :: getCsector
  integer,allocatable,dimension(:,:,:)              :: getCDGsector
  integer,allocatable,dimension(:)                :: getDim,getDimUp,getDimDw
  integer,allocatable,dimension(:)                :: getNup,getNdw
  integer,allocatable,dimension(:,:,:)            :: getBathStride
  integer,allocatable,dimension(:,:)              :: impIndex
  logical,allocatable,dimension(:)                :: twin_mask
  logical,allocatable,dimension(:)                :: sectors_mask

  !Effective Bath used in the VCA code (this is opaque to user)
  !=========================================================
  type(effective_bath)                            :: vca_bath


  !Eigenvalues,Eigenvectors FULL DIAGONALIZATION
  !=========================================================
  type(full_espace),dimension(:),allocatable        :: espace


  !Sparse matrix for Lanczos diagonalization.
  !=========================================================  
  type(sparse_matrix_csr)                            :: spH0d !diagonal part
  type(sparse_matrix_csr)                            :: spH0nd !non-diagonal part
  type(sparse_matrix_csr),dimension(:),allocatable   :: spH0ups,spH0dws !reduced UP and DW parts
  procedure(dd_sparse_HxV),pointer                   :: spHtimesV_p=>null()



  !Variables for DIAGONALIZATION
  !=========================================================  
  integer,allocatable,dimension(:)                  :: neigen_sector
  !--------------- LATTICE WRAP VARIABLES -----------------!  
  integer,allocatable,dimension(:,:)                 :: neigen_sectorii
  integer,allocatable,dimension(:)                   :: neigen_totalii
  logical                                            :: trim_state_list=.false.


  !Partition function, Omega potential, SFT potential
  !=========================================================
  real(8)                                           :: zeta_function
  real(8)                                           :: omega_potential
  real(8)                                           :: sft_potential


  !Cluster Green's functions
  !(Nlat,Nlat,Nspin,Nspin,Norb,Norb,:)
  !=========================================================
  complex(8),allocatable,dimension(:,:,:,:,:,:,:)   :: impGmats ![Nlat][Nlat][Nspin][Nspin][Norb][Norb][L]
  complex(8),allocatable,dimension(:,:,:,:,:,:,:)   :: impGreal ![Nlat][Nlat][Nspin][Nspin][Norb][Norb][L]
  !
  complex(8),allocatable,dimension(:,:,:,:,:,:,:)   :: impG0mats ![Nlat][Nlat][Nspin][Nspin][Norb][Norb][L]
  complex(8),allocatable,dimension(:,:,:,:,:,:,:)   :: impG0real ![Nlat][Nlat][Nspin][Nspin][Norb][Norb][L]
  !
  complex(8),allocatable,dimension(:,:,:,:,:,:,:)   :: impSmats ![Nlat][Nlat][Nspin][Nspin][Norb][Norb][L]
  complex(8),allocatable,dimension(:,:,:,:,:,:,:)   :: impSreal ![Nlat][Nlat][Nspin][Nspin][Norb][Norb][L]
  !
  type(GFmatrix),allocatable,dimension(:,:,:,:,:,:) :: impGmatrix

  !Spin Susceptibilities
  !=========================================================
  real(8),allocatable,dimension(:,:)                 :: spinChi_tau
  complex(8),allocatable,dimension(:,:)              :: spinChi_w
  complex(8),allocatable,dimension(:,:)              :: spinChi_iv

  !Cluster local observables:
  !=========================================================
  real(8),dimension(:,:),allocatable                ::  imp_dens    ![Nlat][Norb]
  real(8),dimension(:,:),allocatable                ::  imp_dens_up ![Nlat][Norb]
  real(8),dimension(:,:),allocatable                ::  imp_dens_dw ![Nlat][Norb]
  real(8),dimension(:,:),allocatable                ::  imp_docc    ![Nlat][Norb]



  !Suffix string attached to the output files.
  !=========================================================
  character(len=64)                                 :: file_suffix=""
  logical                                           :: offdiag_gf_flag=.false.


  !This is the internal Mpi Communicator and variables.
  !=========================================================
#ifdef _MPI
  integer                                            :: MpiComm_Global=MPI_UNDEFINED
  integer                                            :: MpiComm=MPI_UNDEFINED
#endif
  integer                                            :: MpiGroup_Global=MPI_GROUP_NULL
  integer                                            :: MpiGroup=MPI_GROUP_NULL
  logical                                            :: MpiStatus=.false.
  logical                                            :: MpiMaster=.true.
  integer                                            :: MpiRank=0
  integer                                            :: MpiSize=1
  integer,allocatable,dimension(:)                   :: MpiMembers
  integer                                            :: mpiQup=0
  integer                                            :: mpiRup=0
  integer                                            :: mpiQdw=0
  integer                                            :: mpiRdw=0
  integer                                            :: mpiQ=0
  integer                                            :: mpiR=0
  integer                                            :: mpiIstart
  integer                                            :: mpiIend
  integer                                            :: mpiIshift
  logical                                            :: mpiAllThreads=.true.


  !Frequency and time arrays:
  !=========================================================
  !real(8),dimension(:),allocatable                  :: wm,tau,wr,vm



  !flag for finite temperature calculation
  !=========================================================
  logical                                           :: finiteT 


contains



  !> Get stride position in the one-particle many-body space 
  function index_stride_lso(ilat,iorb,ispin) result(indx)
    integer :: ilat
    integer :: iorb
    integer :: ispin
    integer :: indx
    indx = iorb + (ilat-1)*Norb + (ispin-1)*Norb*Nlat
  end function index_stride_lso



  !>Allocate and Deallocate Hilbert space maps (sector<-->Fock)
  subroutine map_allocate_scalar(H,N)
    type(sector_map) :: H
    integer          :: N
    if(H%status) call map_deallocate_scalar(H)
    allocate(H%map(N))
    H%status=.true.
  end subroutine map_allocate_scalar
  !
  subroutine map_allocate_vector(H,N)
    type(sector_map),dimension(:)       :: H
    integer,dimension(size(H))          :: N
    integer                             :: i
    do i=1,size(H)
       call map_allocate_scalar(H(i),N(i))
    enddo
  end subroutine map_allocate_vector



  subroutine map_deallocate_scalar(H)
    type(sector_map) :: H
    if(.not.H%status)then
       write(*,*) "WARNING map_deallocate_scalar: H is not allocated"
       return
    endif
    if(allocated(H%map))deallocate(H%map)
    H%status=.false.
  end subroutine map_deallocate_scalar
  !
  subroutine map_deallocate_vector(H)
    type(sector_map),dimension(:) :: H
    integer                       :: i
    do i=1,size(H)
       call map_deallocate_scalar(H(i))
    enddo
  end subroutine map_deallocate_vector
 
  !=========================================================

 subroutine vca_set_MpiComm(comm)
#ifdef _MPI
    integer :: comm,ierr
    MpiComm_Global = comm
    MpiComm        = MpiComm_Global
    MpiStatus      = .true.
    MpiSize        = get_Size_MPI(MpiComm_Global)
    MpiRank        = get_Rank_MPI(MpiComm_Global)
    MpiMaster      = get_Master_MPI(MpiComm_Global)
    call Mpi_Comm_group(MpiComm_Global,MpiGroup_Global,ierr)
#else
    integer,optional :: comm
#endif
  end subroutine vca_set_MpiComm

  subroutine vca_del_MpiComm()
#ifdef _MPI
    MpiComm_Global = MPI_UNDEFINED
    MpiComm        = MPI_UNDEFINED
    MpiStatus      = .false.
    MpiSize        = 1
    MpiRank        = 0
    MpiMaster      = .true.
#endif
  end subroutine vca_del_MpiComm

  !Allocate the channels in GFmatrix structure
  subroutine allocate_gfmatrix_Nchan(self,N)
    type(GFmatrix) :: self
    integer        :: N
    if(allocated(self%channel))deallocate(self%channel)
    allocate(self%channel(N))
  end subroutine allocate_gfmatrix_Nchan

  !Allocate the Excitations spectrum at a given channel
  subroutine allocate_gfmatrix_Nexc(self,i,Nexc)
    type(GFmatrix) :: self
    integer        :: i
    integer        :: Nexc
    if(allocated(self%channel(i)%weight))deallocate(self%channel(i)%weight)
    if(allocated(self%channel(i)%poles))deallocate(self%channel(i)%poles)
    allocate(self%channel(i)%weight(Nexc))
    allocate(self%channel(i)%poles(Nexc))
  end subroutine allocate_gfmatrix_Nexc
 



END MODULE VCA_VARS_GLOBAL
