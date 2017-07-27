MODULE VCA_VARS_GLOBAL
  USE VCA_INPUT_VARS
  !
  USE SF_CONSTANTS
  USE SF_IOTOOLS, only:str
  implicit none



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
  integer                                         :: Nexc


  !non-interacting cluster Hamiltonian
  !=========================================================
  real(8),dimension(:,:,:,:,:,:),allocatable      :: impHloc ![Nlat][Nlat][Nspin][Nspin][Norb][Norb]


  !Some maps between sectors and full Hilbert space (pointers)
  !=========================================================
  integer,allocatable,dimension(:,:)              :: getsector
  integer,allocatable,dimension(:,:)              :: getCsector
  integer,allocatable,dimension(:,:)              :: getCDGsector
  integer,allocatable,dimension(:)                :: getDim,getDimUp,getDimDw
  integer,allocatable,dimension(:)                :: getNup,getNdw


  !Eigenvalues,Eigenvectors FULL DIAGONALIZATION
  !Hamiltonian eig-space structure
  !=========================================================
  type full_espace
     real(8),dimension(:),pointer                 :: e
     real(8),dimension(:,:),pointer               :: M
  end type full_espace
  type(full_espace),dimension(:),allocatable      :: espace


  !Partition function
  !=========================================================
  real(8)                                         :: zeta_function
  real(8)                                         :: omega_potential


  !Cluster Green's functions
  !(Nlat,Nlat,Nspin,Nspin,Norb,Norb,:)
  !=========================================================
  complex(8),allocatable,dimension(:,:,:,:,:,:,:) :: impGmats ![Nlat][Nlat][Nspin][Nspin][Norb][Norb][L]
  complex(8),allocatable,dimension(:,:,:,:,:,:,:) :: impGreal ![Nlat][Nlat][Nspin][Nspin][Norb][Norb][L]


  !Q and Lambda VCA matrices:
  !=========================================================
  real(8),allocatable,dimension(:,:,:,:)          :: cQmatrix   ![Nlat][Nspin][Norb][Nexcitations]
  real(8),allocatable,dimension(:,:,:,:)          :: cdgQmatrix ![Nlat][Nspin][Norb][Nexcitations]
  real(8),allocatable,dimension(:,:)                :: Lmatrix    ![Nspin][Nexcitations]
  ! real(8),allocatable,dimension(:,:,:,:,:,:,:)    :: Lmatrix    ![Nlat][Nlat][Nspin][Nspin][Norb][Norb][Nexcitations]


  !Cluster local observables:
  !=========================================================
  real(8),dimension(:,:),allocatable              ::  imp_dens ![Nlat][Norb]
  real(8),dimension(:,:),allocatable              ::  imp_dens_up ![Nlat][Norb]
  real(8),dimension(:,:),allocatable              ::  imp_dens_dw ![Nlat][Norb]
  real(8),dimension(:,:),allocatable              ::  imp_docc ![Nlat][Norb]


  !Suffix string attached to the output files.
  !=========================================================
  character(len=64)                               :: file_suffix=""


  !SECTOR-TO-FOCK SPACE STRUCTURE
  !=========================================================
  type sector_map
     integer,dimension(:),allocatable             :: map
  end type sector_map

  interface map_allocate
     module procedure                             :: map_allocate_scalar
     module procedure                             :: map_allocate_vector
  end interface map_allocate

  interface map_deallocate
     module procedure                             :: map_deallocate_scalar
     module procedure                             :: map_deallocate_vector
  end interface map_deallocate




contains





  subroutine map_allocate_scalar(H,N)
    type(sector_map)                              :: H
    integer                                       :: N
    allocate(H%map(N))
  end subroutine map_allocate_scalar
  !
  subroutine map_allocate_vector(H,N)
    type(sector_map),dimension(:)                 :: H
    integer,dimension(size(H))                    :: N
    integer                                       :: i
    do i=1,size(H)
       allocate(H(i)%map(N(i)))
    enddo
  end subroutine map_allocate_vector


  subroutine map_deallocate_scalar(H)
    type(sector_map)                              :: H
    deallocate(H%map)
  end subroutine map_deallocate_scalar
  !
  subroutine map_deallocate_vector(H)
    type(sector_map),dimension(:)                 :: H
    integer                                       :: i
    do i=1,size(H)
       deallocate(H(i)%map)
    enddo
  end subroutine map_deallocate_vector











END MODULE VCA_VARS_GLOBAL
