MODULE VCA_INPUT_VARS
  USE SF_VERSION
  USE SF_PARSE_INPUT
  USE SF_IOTOOLS, only:str
  implicit none

  !GIT VERSION
  include "revision.inc"  !this file is generated at compilation time in the Makefile

  !input variables
  !=========================================================
  integer              :: Nlat         !# of cluster sites
  integer              :: Norb         !# of lattice orbitals per site
  integer              :: Nspin        !# spin degeneracy (max 2)
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
  logical              :: print_G      !flag to print impurity Green`s functions
  logical              :: diag_twin    !flag to reduce (T) or not (F,default) the number of visited sector using twin symmetry.
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


  !LOG AND HLOC UNITS
  !=========================================================
  character(len=100)   :: HLOCfile
  integer,save         :: LOGfile



contains


  !+-------------------------------------------------------------------+
  !PURPOSE  : READ THE INPUT FILE AND SETUP GLOBAL VARIABLES
  !+-------------------------------------------------------------------+
  subroutine vca_read_input(INPUTunit,comm)
#ifdef _MPI
    USE MPI
    USE SF_MPI
#endif
    character(len=*) :: INPUTunit
    integer,optional :: comm
    logical          :: master=.true.
    integer          :: i,rank=0
#ifdef _MPI
    if(present(comm))then
       master=get_Master_MPI(comm)
       rank  =get_Rank_MPI(comm)
    endif
#endif


    !DEFAULT VALUES OF THE PARAMETERS:
    call parse_input_variable(Norb,"NORB",INPUTunit,default=1,comment="Number of impurity orbitals.")
    call parse_input_variable(Nspin,"NSPIN",INPUTunit,default=1,comment="Number of spin degeneracy (max 2)")
    call parse_input_variable(uloc,"ULOC",INPUTunit,default=[2.d0,0.d0,0.d0],comment="Values of the local interaction per orbital (max 3)")
    call parse_input_variable(ust,"UST",INPUTunit,default=0.d0,comment="Value of the inter-orbital interaction term")
    call parse_input_variable(Jh,"JH",INPUTunit,default=0.d0,comment="Hunds coupling")
    call parse_input_variable(Jx,"JX",INPUTunit,default=0.d0,comment="S-E coupling")
    call parse_input_variable(Jp,"JP",INPUTunit,default=0.d0,comment="P-H coupling")
    call parse_input_variable(beta,"BETA",INPUTunit,default=1000.d0,comment="Inverse temperature, at T=0 is used as a IR cut-off.")
    call parse_input_variable(xmu,"XMU",INPUTunit,default=0.d0,comment="Chemical potential. If HFMODE=T, xmu=0 indicates half-filling condition.")
    call parse_input_variable(print_Sigma,"PRINT_SIGMA",INPUTunit,default=.true.,comment="flag to print impurity Self-energies")
    call parse_input_variable(print_G,"PRINT_G",INPUTunit,default=.true.,comment="flag to print impurity Greens function")
    call parse_input_variable(print_G0,"PRINT_G0",INPUTunit,default=.true.,comment="flag to print non-interacting impurity Greens function")
    call parse_input_variable(Lmats,"LMATS",INPUTunit,default=5000,comment="Number of Matsubara frequencies.")
    call parse_input_variable(Lreal,"LREAL",INPUTunit,default=5000,comment="Number of real-axis frequencies.")
    call parse_input_variable(Ltau,"LTAU",INPUTunit,default=1000,comment="Number of imaginary time points.")
    call parse_input_variable(nread,"NREAD",INPUTunit,default=0.d0,comment="Objective density for fixed density calculations.")
    call parse_input_variable(nerr,"NERR",INPUTunit,default=1.d-4,comment="Error threshold for fixed density calculations.")
    call parse_input_variable(ndelta,"NDELTA",INPUTunit,default=0.1d0,comment="Initial step for fixed density calculations.")
    call parse_input_variable(wini,"WINI",INPUTunit,default=-5.d0,comment="Smallest real-axis frequency")
    call parse_input_variable(wfin,"WFIN",INPUTunit,default=5.d0,comment="Largest real-axis frequency")
    call parse_input_variable(jhflag,"JHFLAG",INPUTunit,default=.false.,comment="Flag to include full SU(2) invariant term: spin-flip, pair-hopping.")
    call parse_input_variable(hfmode,"HFMODE",INPUTunit,default=.true.,comment="Flag to set the Hartree form of the interaction (n-1/2). see xmu.")
    call parse_input_variable(eps,"EPS",INPUTunit,default=0.01d0,comment="Broadening on the real-axis.")
    call parse_input_variable(cutoff,"CUTOFF",INPUTunit,default=1.d-9,comment="Spectrum cut-off, used to determine the number states to be retained.")
    call parse_input_variable(gs_threshold,"GS_THRESHOLD",INPUTunit,default=1.d-9,comment="Energy threshold for ground state degeneracy loop up")
    call parse_input_variable(HLOCfile,"HLOCFILE",INPUTunit,default="inputHLOC.in",comment="File read the input local H.")
    call parse_input_variable(LOGfile,"LOGFILE",INPUTunit,default=6,comment="LOG unit.")
    call parse_input_variable(file_suffix,"FILE_SUFFIX",INPUTunit,default=".vca",comment="Suffix in the output files.")        
    call parse_input_variable(verbose,"VERBOSE",INPUTunit,default=3,comment="Verbosity level: 0=almost nothing --> 3:all.")

#ifdef _MPI
    if(present(comm))then
       if(.not.master)then
          LOGfile=1000-rank
          open(LOGfile,file="stdOUT.rank"//str(rank)//".vca")
          do i=1,get_Size_MPI(comm)
             if(i==rank)write(*,"(A,I0,A,I0)")"Rank ",rank," writing to unit: ",LOGfile
          enddo
       endif
    endif
#endif
    !
    !
    Ltau=max(int(beta),Ltau)
    if(master)then
       call save_input_file(INPUTunit)
       call scifor_version()
       call code_version(revision)
    endif
    call substring_delete(file_suffix,".vca")
    call substring_delete(Hfile,".restart")
    call substring_delete(Hfile,".vca")
  end subroutine vca_read_input




  subroutine substring_delete (s,sub)
    !! S_S_DELETE2 recursively removes a substring from a string.
    !    The remainder is left justified and padded with blanks.
    !    The substitution is recursive, so
    !    that, for example, removing all occurrences of "ab" from
    !    "aaaaabbbbbQ" results in "Q".
    !  Parameters:
    !    Input/output, character ( len = * ) S, the string to be transformed.
    !    Input, character ( len = * ) SUB, the substring to be removed.
    !    Output, integer ( kind = 4 ) IREP, the number of occurrences of
    !    the substring.
    integer          :: ihi
    integer          :: irep
    integer          :: loc
    integer          :: nsub
    character(len=*) ::  s
    integer          :: s_length
    character(len=*) :: sub
    s_length = len ( s )
    nsub = len ( sub )
    irep = 0
    ihi = s_length
    do while ( 0 < ihi )
       loc = index ( s(1:ihi), sub )
       if ( loc == 0 ) then
          return
       end if
       irep = irep + 1
       call s_chop ( s, loc, loc+nsub-1 )
       ihi = ihi - nsub
    end do
    return
  end subroutine substring_delete

  subroutine s_chop ( s, ilo, ihi )
    !! S_CHOP "chops out" a portion of a string, and closes up the hole.
    !  Example:
    !    S = 'Fred is not a jerk!'
    !    call s_chop ( S, 9, 12 )
    !    S = 'Fred is a jerk!    '
    !  Parameters:
    !    Input/output, character ( len = * ) S, the string to be transformed.
    !    Input, integer ( kind = 4 ) ILO, IHI, the locations of the first and last
    !    characters to be removed.
    integer               ::ihi
    integer               ::ihi2
    integer               ::ilo
    integer               ::ilo2
    character ( len = * ) :: s
    integer               ::s_length
    s_length = len ( s )
    ilo2 = max ( ilo, 1 )
    ihi2 = min ( ihi, s_length )
    if ( ihi2 < ilo2 ) then
       return
    end if
    s(ilo2:s_length+ilo2-ihi2-1) = s(ihi2+1:s_length)
    s(s_length+ilo2-ihi2:s_length) = ' '
    return
  end subroutine s_chop


END MODULE VCA_INPUT_VARS
