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
  integer                                       :: Nlat         !# of cluster sites
  integer                                       :: Norb         !# of lattice orbitals per site
  integer                                       :: Nspin        !# spin degeneracy (max 2)
  real(8),dimension(2)                          :: Uloc         !local interactions
  real(8)                                       :: Ust          !intra-orbitals interactions
  real(8)                                       :: Jh           !J_Hund: Hunds' coupling constant 
  real(8)                                       :: Jx           !J_X: coupling constant for the spin-eXchange interaction term
  real(8)                                       :: Jp           !J_P: coupling constant for the Pair-hopping interaction term 
  real(8)                                       :: xmu          !chemical potential
  real(8)                                       :: beta         !inverse temperature
  real(8)                                       :: eps          !broadening
  real(8)                                       :: wini,wfin    !
  logical                                       :: Jhflag       !spin-exchange and pair-hopping flag.
  logical                                       :: HFmode       !flag for HF interaction form U(n-1/2)(n-1/2) VS Unn
  real(8)                                       :: cutoff       !cutoff for spectral summation
  real(8)                                       :: gs_threshold !Energy threshold for ground state degeneracy loop up
  real(8)                                       :: sb_field     !symmetry breaking field
  logical                                       :: print_G      !flag to print impurity Green`s functions
  logical                                       :: diag_twin    !flag to reduce (T) or not (F,default) the number of visited sector using twin symmetry.
  real(8)                                       :: nread        !fixed density. if 0.d0 fixed chemical potential calculation.
  real(8)                                       :: nerr         !fix density threshold. a loop over from 1.d-1 to required nerr is performed
  real(8)                                       :: ndelta       !initial chemical potential step
  integer                                       :: niter        !
  integer                                       :: verbose      !


  !Some parameters for function dimension:
  !=========================================================
  integer                                       :: Lmats
  integer                                       :: Lreal
  integer                                       :: Ltau





contains



  !+-------------------------------------------------------------------+
  !PURPOSE  : READ THE INPUT FILE AND SETUP GLOBAL VARIABLES
  !+-------------------------------------------------------------------+
  subroutine vca_read_input(INPUTunit)
    character(len=*) :: INPUTunit
    integer          :: i


    !DEFAULT VALUES OF THE PARAMETERS:
    call parse_input_variable(Nlat,"NLAT",INPUTunit,default=1,comment="Number of cluster sites.")
    call parse_input_variable(Norb,"NORB",INPUTunit,default=1,comment="Number of orbitals per cluster site.")
    call parse_input_variable(Nspin,"NSPIN",INPUTunit,default=1,comment="Number of spin degeneracy")
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
    call parse_input_variable(print_G,"PRINT_G",INPUTunit,default=.true.,comment="Flag to print impurity Greens function")
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