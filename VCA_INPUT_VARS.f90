MODULE VCA_INPUT_VARS
  USE SF_PARSE_INPUT
  USE SF_IOTOOLS, only:str
  USE SF_VERSION
  implicit none

  !GIT VERSION
  include "revision.inc"  !this file is generated at compilation time in the Makefile


  !INPUT VARIABLES (to be exported in MAIN module)
  !input variables
  !=========================================================
  integer              :: Norb         !# of lattice orbitals per site
  integer              :: Nspin        !# spin degeneracy (max 2)
  integer              :: Nlat         !# size of cluster
  integer              :: Nbath        !# of bath sites (per orbital or not depending on bath_type)
  real(8),dimension(2) :: Uloc         !local interactions
  real(8)              :: Ust          !intra-orbitals interactions
  real(8)              :: Jh           !J_Hund: Hunds' coupling constant 
  real(8)              :: Jx           !J_X: coupling constant for the spin-eXchange interaction term
  real(8)              :: Jp           !J_P: coupling constant for the Pair-hopping interaction term 
  real(8)              :: xmu          !chemical potential
  real(8)              :: beta         !inverse temperature
  real(8)              :: eps          !broadening
  real(8)              :: wini,wfin    !
  logical              :: Jhflag       !spin-exchange and pair-hopping flag.
  logical              :: HFmode       !flag for HF interaction form U(n-1/2)(n-1/2) VS Unn
  real(8)              :: cutoff       !cutoff for spectral summation
  real(8)              :: gs_threshold !Energy threshold for ground state degeneracy loop up
  real(8)              :: sb_field     !symmetry breaking field
  logical              :: print_Sigma  !flag to print impurity Green`s functions
  logical              :: print_impG   !flag to print impurity Green`s functions
  logical              :: print_impG0  !flag to print impurity Green`s functions
  logical              :: diag_twin    !flag to reduce (T) or not (F,default) the number of visited sector using twin symmetry.
  ! character(len=7)     :: bath_type           !flag to set bath type: normal (1bath/imp), hybrid(1bath)
  !
  character(len=6)     :: vca_method          !flag to set ED method: full (full Diagonalization) OR lanc (Lanczos based T=0/T>0 diagonalization)
  !
  real(8)              :: lanc_tolerance      !Tolerance for the Lanczos iterations as used in Arpack and plain lanczos. 
  integer              :: lanc_niter          !Max number of Lanczos iterations
  integer              :: lanc_niter_spectrum !Max number of iterations in the spectrum annealing
  integer              :: lanc_ngfiter        !Max number of iteration in resolvant tri-diagonalization
  integer              :: lanc_ncv_factor     !Set the size of the block used in Lanczos-Arpack by multiplying the required Neigen (Ncv=lanc_ncv_factor*Neigen+lanc_ncv_add)
  integer              :: lanc_ncv_add        !Adds up to the size of the block to prevent it to become too small (Ncv=lanc_ncv_factor*Neigen+lanc_ncv_add)
  integer              :: lanc_nstates_sector !Max number of required eigenvalues per sector
  integer              :: lanc_nstates_total  !Max number of states hold in the finite T calculation
  integer              :: lanc_nstates_step   !Number of states added at each step to determine the optimal spectrum size at finite T
  integer              :: lanc_dim_threshold  !Min dimension threshold to use Lanczos determination of the spectrum rather than Lapack based exact diagonalization.
  real(8)              :: lanc_spectrum_threshold  !Threshold for the spectrum annealing error.
  logical              :: vca_sparse_H         !flag to select  storage of sparse matrix H (mem--, cpu++) if TRUE, or direct on-the-fly H*v product (mem++, cpu--) if FALSE
  !

  real(8)              :: nread        !fixed density. if 0.d0 fixed chemical potential calculation.
  real(8)              :: nerr         !fix density threshold. a loop over from 1.d-1 to required nerr is performed
  real(8)              :: ndelta       !initial chemical potential step
  integer              :: niter        !
  integer              :: verbose      !


  !Some parameters for function dimension:
  !=========================================================
  integer              :: Lmats
  integer              :: Lreal
  integer              :: Ltau





contains




  !+-------------------------------------------------------------------+
  !PURPOSE  : READ THE INPUT FILE AND SETUP GLOBAL VARIABLES
  !+-------------------------------------------------------------------+
  subroutine vca_read_input(INPUTunit)
    character(len=*) :: INPUTunit
    !
    !DEFAULT VALUES OF THE PARAMETERS:
    call parse_input_variable(Norb,"NORB",INPUTunit,default=1,comment="Number of orbitals per cluster site.")
    call parse_input_variable(Nspin,"NSPIN",INPUTunit,default=1,comment="Number of spin degeneracy")
    call parse_input_variable(Nlat,"NLAT",INPUTunit,default=1,comment="Number of cluster copies tiling the system")
    call parse_input_variable(Nbath,"NBATH",INPUTunit,default=0,comment="Number of bath sites:(normal=>Nbath per orb)(hybrid=>Nbath total)")
    !
    call parse_input_variable(uloc,"ULOC",INPUTunit,default=[2.d0,0.d0,0.d0],comment="Values of the local interaction per orbital")
    call parse_input_variable(ust,"UST",INPUTunit,default=0.d0,comment="Value of the inter-orbital interaction term")
    call parse_input_variable(Jh,"JH",INPUTunit,default=0.d0,comment="Hunds coupling")
    call parse_input_variable(jhflag,"JHFLAG",INPUTunit,default=.false.,comment="Flag to include full rotational invariant terms: SE, PH.")
    call parse_input_variable(Jx,"JX",INPUTunit,default=0.d0,comment="Spin-Exchance coupling")
    call parse_input_variable(Jp,"JP",INPUTunit,default=0.d0,comment="Pair-Hopping coupling")
    !
    call parse_input_variable(beta,"BETA",INPUTunit,default=1000.d0,comment="Inverse temperature, at T=0 it is used as a IR cut-off.")
    call parse_input_variable(xmu,"XMU",INPUTunit,default=0.d0,comment="Chemical potential. If HFMODE=T, xmu=0 indicates half-filling condition.")
    call parse_input_variable(hfmode,"HFMODE",INPUTunit,default=.true.,comment="Flag to set the Hartree form of the interaction (n-1/2). see xmu.")
    !
    call parse_input_variable(Lmats,"LMATS",INPUTunit,default=5000,comment="Number of Matsubara frequencies.")
    call parse_input_variable(Lreal,"LREAL",INPUTunit,default=5000,comment="Number of Real-axis frequencies.")
    call parse_input_variable(Ltau,"LTAU",INPUTunit,default=1000,comment="Number of Imaginary time points.")
    call parse_input_variable(wini,"WINI",INPUTunit,default=-5.d0,comment="Smallest real-axis frequency")
    call parse_input_variable(wfin,"WFIN",INPUTunit,default=5.d0,comment="Largest real-axis frequency")
    call parse_input_variable(eps,"EPS",INPUTunit,default=0.01d0,comment="Broadening on the real-axis.")
    call parse_input_variable(print_Sigma,"PRINT_SIGMA",INPUTunit,default=.true.,comment="Flag to print impurity Self-energy")
    call parse_input_variable(print_impG,"PRINT_IMPG",INPUTunit,default=.true.,comment="Flag to print impurity interacting Greens function")
    call parse_input_variable(print_impG0,"PRINT_IMPG0",INPUTunit,default=.true.,comment="Flag to print impurity non-interacting Greens function")
    ! call parse_input_variable(bath_type,"BATH_TYPE",INPUTunit,default='normal',comment="flag to set bath type: normal (1bath/imp), hybrid(1bath)")
    !
    call parse_input_variable(vca_method,"VCA_METHOD",INPUTunit,default="lanc",comment="flag to set ED method: full (full Diagonalization) OR lanc (Lanczos based T=0/T>0 diagonalization)")
    !
    call parse_input_variable(diag_twin,"DIAG_TWIN",INPUTunit,default=.false.,comment="flag to reduce (T) or not (F,default) the number of visited sector using twin symmetry.")
    !
    call parse_input_variable(vca_sparse_H,"VCA_SPARSE_H",INPUTunit,default=.true.,comment="flag to select  storage of sparse matrix H (mem--, cpu++) if TRUE, or direct on-the-fly H*v product (mem++, cpu--) if FALSE ")
    call parse_input_variable(lanc_nstates_sector,"LANC_NSTATES_SECTOR",INPUTunit,default=6,comment="Initial number of states per sector to be determined.")
    call parse_input_variable(lanc_nstates_total,"LANC_NSTATES_TOTAL",INPUTunit,default=1,comment="Initial number of total states to be determined.")
    call parse_input_variable(lanc_nstates_step,"LANC_NSTATES_STEP",INPUTunit,default=2,comment="Number of states added to the spectrum at each step.")
    call parse_input_variable(lanc_ncv_factor,"LANC_NCV_FACTOR",INPUTunit,default=3,comment="Set the size of the block used in Lanczos-Arpack by multiplying the required Neigen (Ncv=lanc_ncv_factor*Neigen+lanc_ncv_add)")
    call parse_input_variable(lanc_ncv_add,"LANC_NCV_ADD",INPUTunit,default=5,comment="Adds up to the size of the block to prevent it to become too small (Ncv=lanc_ncv_factor*Neigen+lanc_ncv_add)")
    call parse_input_variable(lanc_niter,"LANC_NITER",INPUTunit,default=512,comment="Number of Lanczos iteration in spectrum determination.")
    call parse_input_variable(lanc_niter_spectrum,"LANC_NITER_SPECTRUM",INPUTunit,default=10,comment="Number of iterations in spectrum annealing.")
    call parse_input_variable(lanc_ngfiter,"LANC_NGFITER",INPUTunit,default=200,comment="Number of Lanczos iteration in GF determination. Number of momenta.")
    call parse_input_variable(lanc_tolerance,"LANC_TOLERANCE",INPUTunit,default=0.000000000001d0,comment="Tolerance for the Lanczos iterations as used in Arpack and plain lanczos.")
    call parse_input_variable(lanc_dim_threshold,"LANC_DIM_THRESHOLD",INPUTunit,default=256,comment="Min dimension threshold to use Lanczos determination of the spectrum rather than Lapack based exact diagonalization.")
    call parse_input_variable(lanc_spectrum_threshold,"LANC_SPECTRUM_THRESHOLD",INPUTunit,default=1d-5,comment="Threshold for the spectrum annealing error.")
    !
    call parse_input_variable(nread,"NREAD",INPUTunit,default=0.d0,comment="Objective density for fixed density calculations.")
    call parse_input_variable(nerr,"NERR",INPUTunit,default=1.d-4,comment="Error threshold for fixed density calculations.")
    call parse_input_variable(ndelta,"NDELTA",INPUTunit,default=0.1d0,comment="Initial step for fixed density calculations.")
    !
    call parse_input_variable(cutoff,"CUTOFF",INPUTunit,default=1.d-9,comment="Spectrum cut-off, used to determine the number states to be retained.")
    call parse_input_variable(gs_threshold,"GS_THRESHOLD",INPUTunit,default=1.d-9,comment="Energy threshold for ground state degeneracy loop up")
    call parse_input_variable(verbose,"VERBOSE",INPUTunit,default=3,comment="Verbosity level: 0=almost nothing --> 3:all.")
    !
    !
    Ltau=max(int(beta),Ltau)
    call save_input_file(INPUTunit)
    call scifor_version()
    call code_version(revision)
  end subroutine vca_read_input


END MODULE VCA_INPUT_VARS
