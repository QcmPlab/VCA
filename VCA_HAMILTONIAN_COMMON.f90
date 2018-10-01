MODULE VCA_HAMILTONIAN_COMMON
  USE SF_MISC,    only: assert_shape
  USE SF_LINALG,  only: kronecker_product,eye
  USE VCA_INPUT_VARS
  USE VCA_VARS_GLOBAL
!  USE VCA_BATH
  USE VCA_SETUP
  implicit none


  !> MPI local variables (shared)
  ! #ifdef _MPI
  !   integer                                   :: MpiComm=MPI_UNDEFINED
  ! #else
  !   integer                                   :: MpiComm=0
  ! #endif
  !   logical                                   :: MpiStatus=.false.
  !   logical                                   :: MpiMaster=.true.
  !   integer                                   :: MpiIerr
  !   integer                                   :: MpiRank=0
  !   integer                                   :: MpiSize=1
  !   integer                                   :: MpiQ=1
  !   integer                                   :: MpiQup=1
  !   integer                                   :: MpiQdw=1
  !   integer                                   :: MpiR=0
  !   integer                                   :: MpiRup=0
  !   integer                                   :: MpiRdw=0
  !   integer                                   :: MpiIstart
  !   integer                                   :: MpiIend
  !   integer                                   :: MpiIshift
  !
  integer                                   :: Dim
  integer                                   :: DimUp
  integer                                   :: DimDw
  integer,allocatable,dimension(:)          :: DimUps
  integer,allocatable,dimension(:)          :: DimDws
  !
  integer                                   :: Hsector=0
  logical                                   :: Hstatus=.false.
  type(sector_map),dimension(:),allocatable :: Hs




end MODULE VCA_HAMILTONIAN_COMMON




