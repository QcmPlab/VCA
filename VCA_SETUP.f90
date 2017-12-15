MODULE VCA_SETUP
  USE VCA_VARS_GLOBAL
  USE VCA_AUX_FUNX
  !
  USE SF_TIMER
  USE SF_IOTOOLS, only:free_unit,reg,file_length
  implicit none
  private


  public :: init_cluster_structure
  !
  public :: setup_eigenspace
  public :: delete_eigenspace
  !
  public :: setup_pointers_normal
  !
  public :: build_sector
  public :: delete_sector
  !
  public :: imp_state_index
  !
  public :: bdecomp
  public :: bjoin
  public :: flip_state
  !
  public :: c,cdg
  !
  public :: binary_search
  !
  public :: twin_sector_order
  public :: get_twin_sector


contains





  !+------------------------------------------------------------------+
  !PURPOSE  : Init ED structure and calculation
  !+------------------------------------------------------------------+
  subroutine init_cluster_structure()
    logical                                           :: control
    real(8),dimension(Nspin,Nspin,Norb,Norb)          :: reHloc         !local hamiltonian, real part 
    real(8),dimension(Nspin,Nspin,Norb,Norb)          :: imHloc         !local hamiltonian, imag part
    integer                                           :: i,dim_sector_max(2),iorb,jorb,ispin,jspin
    integer                                           :: isector,in,shift
    logical                                           :: MPI_MASTER=.true.
    !
    ! select case(bath_type)
    ! case default
    Ns = (Nbath+1)*Nlat*Norb !Norb per site plus Nbath per orb per site
    ! case ('hybrid')
    !    Ns = Nbath+Nlat*Norb     !Norb per site plus shared Nbath sites
    ! end select
    !
    Nlevels  = 2*Ns
    !
    Nsectors = (Ns+1)*(Ns+1) !nup=0:Ns;ndw=0:Ns
    !
    dim_sector_max=0
    dim_sector_max(1)=get_normal_sector_dimension(Ns/2)
    dim_sector_max(2)=get_normal_sector_dimension(Ns-Ns/2)
    !
    write(LOGfile,"(A)")"Summary:"
    write(LOGfile,"(A)")"--------------------------------------------"
    write(LOGfile,"(A,I15)") '# of levels/spin      = ',Ns
    write(LOGfile,"(A,I15)") 'Total size            = ',Nlevels
    write(LOGfile,"(A,I15)") '# of sites            = ',Nlat
    write(LOGfile,"(A,I15)") '# of orbitals         = ',Norb
    write(LOGfile,"(A,I15)")'# of bath              = ',Nbath
    write(LOGfile,"(A,2I15)")'Fock space size       = ',2**Ns*2**Ns
    write(LOGfile,"(A,2I15)")'Largest Sector        = ',dim_sector_max
    write(LOGfile,"(A,I15)") 'Number of sectors     = ',Nsectors
    write(LOGfile,"(A)")"--------------------------------------------"
    !
    !>CHECKS:
    ! if(bath_type/='normal')stop "VCA ERROR: bath_type != normal is not yet supported. ask developers"
    if(Nspin>1)stop "VCA ERROR: Nspin > 1 is not yet supported. Uncomment this line in VCA_SETUP to use it anyway"
    if(Norb>1)stop "VCA ERROR: Norb > 1 is not yet supported. Uncomment this line in VCA_SETUP to use it anyway"
    !    
    allocate(impHloc(Nlat,Nlat,Nspin,Nspin,Norb,Norb))
    impHloc=zero
    !
    !
    !Allocate indexing arrays
    allocate(getDim(Nsectors));getDim=0
    !
    allocate(getDimUp(Nsectors),getDimDw(Nsectors));getDimUp=0;getDimDw=0
    allocate(getNup(Nsectors),getNdw(Nsectors));getNup=0;getNdw=0
    !
    !
    allocate(getSector(0:Ns,0:Ns))
    getSector=0
    !
    allocate(getCsector(2,Nsectors));getCsector=0
    allocate(getCDGsector(2,Nsectors));getCDGsector=0
    !
    allocate(getBathStride(Nlat,Norb,Nbath));getBathStride=0
    allocate(twin_mask(Nsectors));
    allocate(neigen_sector(Nsectors))
    !
    !check finiteT
    finiteT=.true.              !assume doing finite T per default
    if(lanc_nstates_total==1)then     !is you only want to keep 1 state
       finiteT=.false.          !set to do zero temperature calculations
       write(LOGfile,"(A)")"Required Lanc_nstates_total=1 => set T=0 calculation"
    endif
    !
    !
    !check whether lanc_nstates_sector and lanc_states are even (we do want to keep doublet among states)
    if(finiteT)then
       if(mod(lanc_nstates_sector,2)/=0)then
          lanc_nstates_sector=lanc_nstates_sector+1
          write(LOGfile,"(A,I10)")"Increased Lanc_nstates_sector:",lanc_nstates_sector
       endif
       if(mod(lanc_nstates_total,2)/=0)then
          lanc_nstates_total=lanc_nstates_total+1
          write(LOGfile,"(A,I10)")"Increased Lanc_nstates_total:",lanc_nstates_total
       endif
    endif
    !
    !
    if(finiteT)then
       write(LOGfile,"(A)")"Lanczos FINITE temperature calculation:"
       write(LOGfile,"(A,I3)")"Nstates x Sector = ", lanc_nstates_sector
       write(LOGfile,"(A,I3)")"Nstates   Total  = ", lanc_nstates_total
       call sleep(1)
    else
       write(LOGfile,"(A)")"Lanczos ZERO temperature calculation:"
       call sleep(1)
    endif

    !
    !CHECKS:
    if(Nspin>2)stop "ED ERROR: Nspin > 2 is currently not supported"
    if(Norb>2)stop "ED ERROR: Norb > 2 is currently not supported"
    !
    if(nread/=0.d0)then
       i=abs(floor(log10(abs(nerr)))) !modulus of the order of magnitude of nerror
    endif
    !
    !
    !
    !allocate functions
    allocate(impSmats(Nlat,Nlat,Nspin,Nspin,Norb,Norb,Lmats))
    allocate(impSreal(Nlat,Nlat,Nspin,Nspin,Norb,Norb,Lreal))
    impSmats=zero
    impSreal=zero
    !
    allocate(impGmats(Nlat,Nlat,Nspin,Nspin,Norb,Norb,Lmats))
    allocate(impGreal(Nlat,Nlat,Nspin,Nspin,Norb,Norb,Lreal))
    impGmats=zero
    impGreal=zero
    !
    !allocate functions
    allocate(impG0mats(Nlat,Nlat,Nspin,Nspin,Norb,Norb,Lmats))
    allocate(impG0real(Nlat,Nlat,Nspin,Nspin,Norb,Norb,Lreal))
    impG0mats=zero
    impG0real=zero
    !
    !allocate observables
    allocate(imp_dens(Nlat,Norb),imp_docc(Nlat,Norb),imp_dens_up(Nlat,Norb),imp_dens_dw(Nlat,Norb))
    imp_dens=0d0
    imp_docc=0d0
    imp_dens_up=0d0
    imp_dens_dw=0d0
    !
  end subroutine init_cluster_structure







  !+------------------------------------------------------------------+
  !PURPOSE: SETUP THE GLOBAL POINTERS FOR THE ED CALCULAIONS.
  !+------------------------------------------------------------------+
  !
  !NORMAL CASE
  !
  subroutine setup_pointers_normal
    integer                          :: i,in,dim,isector,jsector
    integer                          :: nup,ndw,jup,jdw,iorb,ilat
    integer                          :: unit,status,istate,stride
    logical                          :: IOfile
    integer                          :: anint
    real(8)                          :: adouble
    integer                          :: list_len
    integer,dimension(:),allocatable :: list_sector
    isector=0
    do nup=0,Ns
       do ndw=0,Ns
          isector=isector+1
          getSector(nup,ndw)=isector
          getNup(isector)=nup
          getNdw(isector)=ndw
          dim = get_normal_sector_dimension(nup,ndw)
          getDim(isector)=dim
          getDimUp(isector)=get_normal_sector_dimension(nup)
          getDimDw(isector)=get_normal_sector_dimension(ndw)
       enddo
    enddo


    inquire(file="state_list"//reg(file_suffix)//".restart",exist=IOfile)
    if(IOfile)then
       write(LOGfile,"(A)")"Restarting from a state_list file:"
       list_len=file_length("state_list"//reg(file_suffix)//".restart")
       allocate(list_sector(list_len))
       open(free_unit(unit),file="state_list"//reg(file_suffix)//".restart",status="old")
       read(unit,*)!read comment line
       status=0
       do while(status>=0)
          read(unit,"(i3,f18.12,2x,ES19.12,1x,2i3,3x,i3,i10)",iostat=status)istate,adouble,adouble,nup,ndw,isector,anint
          list_sector(istate)=isector
          if(nup/=getnup(isector).OR.ndw/=getndw(isector))&
               stop "setup_pointers_normal error: nup!=getnup(isector).OR.ndw!=getndw(isector) "
       enddo
       close(unit)
       lanc_nstates_total = list_len
       do isector=1,Nsectors
          neigen_sector(isector) = max(1,count(list_sector==isector))
       enddo
    else
       do isector=1,Nsectors
          neigen_sector(isector) = min(getdim(isector),lanc_nstates_sector)   !init every sector to required eigenstates
       enddo
    endif
    !
    !
    twin_mask=.true.
    if(diag_twin)then
       ! stop "WARNING: In this updated version with Nup-Ndw factorization the twin-sectors have not been tested!!"
       do isector=1,Nsectors
          nup=getnup(isector)
          ndw=getndw(isector)
          if(nup<ndw)twin_mask(isector)=.false.
       enddo
       write(LOGfile,"(A,I4,A,I4)")"Looking into ",count(twin_mask)," sectors out of ",Nsectors
    endif



    !normal:
    !|imp_up>|bath_up> * |imp_dw>|bath_dw>
    !
    !|imp_sigma>=
    !|(1..Na)_1
    !  ...
    ! (1..Na)_Nl; <-- Norb*Nlat
    ! ([1..Nb]_1...[1..Nb]_Na)_1
    !  ...
    ! ([1..Nb]_1...[1..Nb]_Na)_Nl> <-- Nbath*Norb*Nlat
    stride=Nlat*Norb
    ! select case(bath_type)
    ! case default
    do ilat=1,Nlat
       do iorb=1,Norb
          do i=1,Nbath
             getBathStride(ilat,iorb,i) = i + &
                  (iorb-1)*Nbath + (ilat-1)*Norb*Nbath + stride
          enddo
       enddo
    enddo
    ! case ('hybrid')
    !hybrid:
    !|(1..Na)_1
    !  ...
    ! (1..Na)_Nl; <-- Norb*Nlat
    ! [1..Nb]>    <-- Nbath
    !    do i=1,Nbath
    !       getBathStride(1:Nlat,1:Norb,i) = i + stride
    !    enddo
    ! end select
    !
    getCsector=0
    do isector=1,Nsectors
       nup=getnup(isector);ndw=getndw(isector)
       jup=nup-1;jdw=ndw;if(jup < 0)cycle
       jsector=getsector(jup,jdw)
       getCsector(1,isector)=jsector
    enddo
    !
    do isector=1,Nsectors
       nup=getnup(isector);ndw=getndw(isector)
       jup=nup;jdw=ndw-1;if(jdw < 0)cycle
       jsector=getsector(jup,jdw)
       getCsector(2,isector)=jsector
    enddo
    !
    getCDGsector=0
    do isector=1,Nsectors
       nup=getnup(isector);ndw=getndw(isector)
       jup=nup+1;jdw=ndw;if(jup > Ns)cycle
       jsector=getsector(jup,jdw)
       getCDGsector(1,isector)=jsector
    enddo
    !
    do isector=1,Nsectors
       nup=getnup(isector);ndw=getndw(isector)
       jup=nup;jdw=ndw+1;if(jdw > Ns)cycle
       jsector=getsector(jup,jdw)
       getCDGsector(2,isector)=jsector
    enddo
  end subroutine setup_pointers_normal







  !+------------------------------------------------------------------+
  !PURPOSE  : 
  !+------------------------------------------------------------------+
  subroutine setup_eigenspace
    integer :: isector,dim,jsector
    if(allocated(espace)) deallocate(espace)
    allocate(espace(1:Nsectors))
    do isector=1,Nsectors
       dim=getdim(isector);if(dim==0)stop "setup_eigenspace: dim==0!"
       allocate(espace(isector)%e(dim))
       allocate(espace(isector)%M(dim,dim))
    enddo
  end subroutine setup_eigenspace


  !+------------------------------------------------------------------+
  !PURPOSE  : 
  !+------------------------------------------------------------------+
  subroutine delete_eigenspace
    integer :: isector
    if(allocated(espace))then
       do isector=1,size(espace)
          deallocate(espace(isector)%e)
          deallocate(espace(isector)%M)
       end do
       deallocate(espace)
    endif
  end subroutine delete_eigenspace




  !+------------------------------------------------------------------+
  !PURPOSE  : return the dimension of a sector
  !+------------------------------------------------------------------+
  !NORMAL
  function get_normal_sector_dimension(nup,ndw) result(dim)
    integer :: nup
    integer,optional :: ndw
    integer :: dim,dimup,dimdw
    if(present(ndw))then
       dimup = binomial(Ns,nup)    !this ensures better evaluation of the dimension
       dimdw = binomial(Ns,ndw)    !as it avoids large numbers
    else
       dimup = binomial(Ns,nup)
       dimdw = 1
    endif
    dim=dimup*dimdw
  end function get_normal_sector_dimension





  !+------------------------------------------------------------------+
  !PURPOSE  : constructs the sectors by storing the map to the 
  !states i\in Hilbert_space from the states count in H_sector.
  !+------------------------------------------------------------------+
  subroutine build_sector(isector,Hup)
    integer                                      :: isector
    type(sector_map)                             :: Hup
    integer                                      :: nup,ndw,sz,nt,twoJz
    integer                                      :: nup_,ndw_,sz_,nt_
    integer                                      :: twoSz_,twoLz_
    integer                                      :: i,ibath,iorb
    integer                                      :: iup,idw
    integer                                      :: dim
    integer                                      :: ivec(Ns),jvec(Ns)
    nup = getNup(isector)
    ndw = getNdw(isector)
    dim = getDim(isector)
    call map_allocate(Hup,dim)
    dim=0
    do idw=0,2**Ns-1
       jvec  = bdecomp(idw,Ns)
       ndw_  = sum(jvec)
       if(ndw_ /= ndw)cycle
       do iup=0,2**Ns-1
          ivec  = bdecomp(iup,Ns)
          nup_  = sum(ivec)
          if(nup_ /= nup)cycle
          dim      = dim+1
          Hup%map(dim) = iup + idw*2**Ns
       enddo
    enddo
  end subroutine build_sector



  subroutine delete_sector(isector,Hup)!,Hdw)
    integer                   :: isector
    type(sector_map)          :: Hup
    ! type(sector_map),optional :: Hdw
    call map_deallocate(Hup)
  end subroutine delete_sector






  !> Find position in the state vector for a given lattice-spin-orbital position for the cluster (no bath considered)
  !normal:
  !|imp_up>|bath_up> * |imp_dw>|bath_dw>
  !
  !|imp_sigma>=
  !|(1..Na)_1
  !  ...
  ! (1..Na)_Nl; <-- Norb*Nlat
  ! ([1..Nb]_1...[1..Nb]_Na)_1
  !  ...
  ! ([1..Nb]_1...[1..Nb]_Na)_Nl> <-- Nbath*Norb*Nlat
  function imp_state_index(ilat,iorb,ispin) result(indx)
    integer :: ilat
    integer :: ispin
    integer :: iorb
    integer :: indx
    indx = iorb + (ilat-1)*Norb + (ispin-1)*(Nbath+1)*Norb*Nlat
  end function imp_state_index




  !+-------------------------------------------------------------------+
  !PURPOSE: input state |in> of the basis and calculates 
  !   |out>=C_pos|in>  OR  |out>=C^+_pos|in> ; 
  !   the sign of |out> has the phase convention, pos labels the sites
  !+-------------------------------------------------------------------+
  subroutine c(pos,in,out,fsgn)
    integer,intent(in)    :: pos
    integer,intent(in)    :: in
    integer,intent(inout) :: out
    real(8),intent(inout) :: fsgn    
    integer               :: l
    if(.not.btest(in,pos-1))stop "C error: C_i|...0_i...>"
    fsgn=1d0
    do l=1,pos-1
       if(btest(in,l-1))fsgn=-fsgn
    enddo
    out = ibclr(in,pos-1)
  end subroutine c


  subroutine cdg(pos,in,out,fsgn)
    integer,intent(in)    :: pos
    integer,intent(in)    :: in
    integer,intent(inout) :: out
    real(8),intent(inout) :: fsgn    
    integer               :: l
    if(btest(in,pos-1))stop "C^+ error: C^+_i|...1_i...>"
    fsgn=1d0
    do l=1,pos-1
       if(btest(in,l-1))fsgn=-fsgn
    enddo
    out = ibset(in,pos-1)
  end subroutine cdg






  !##################################################################
  !##################################################################
  !TWIN SECTORS ROUTINES:
  !##################################################################
  !##################################################################

  !+------------------------------------------------------------------+
  !PURPOSE  : Build the re-ordering map to go from sector A(nup,ndw)
  ! to its twin sector B(ndw,nup), with nup!=ndw. 
  !+------------------------------------------------------------------+
  subroutine twin_sector_order(isector,order)
    integer                          :: isector
    integer,dimension(:)             :: order
    type(sector_map)                 :: H,Hup,Hdw
    integer                          :: i,dim
    dim = getdim(isector)
    if(size(Order)/=dim)stop "twin_sector_order error: wrong dimensions of *order* array"
    !- build the map from the A-sector to \HHH
    !- get the list of states in \HHH corresponding to sector B twin of A
    !- return the ordering of B-states in \HHH with respect to those of A
    call build_sector(isector,H)
    do i=1,dim
       Order(i)=flip_state(H%map(i))
    enddo
    call sort_array(Order)
    deallocate(H%map)
    !
  end subroutine twin_sector_order



  !+------------------------------------------------------------------+
  !PURPOSE  : Flip an Hilbert space state m=|{up}>|{dw}> into:
  !
  ! normal: j=|{dw}>|{up}>  , nup --> ndw
  ! superc: j=|{dw}>|{up}>  , sz  --> -sz
  ! nonsu2: j=|{!up}>|{!dw}>, n   --> 2*Ns-n
  !+------------------------------------------------------------------+
  function flip_state(m,n) result(j)
    integer          :: m
    integer,optional :: n
    integer          :: j
    integer          :: ivec(2*Ns),foo(2*Ns),ivup(Ns),ivdw(Ns)
    ! Ivup = bdecomp(m,Ns)
    ! Ivdw = bdecomp(n,Ns)
    Ivec = bdecomp(m,2*Ns)
    foo(1:Ns)     =Ivec(Ns+1:2*Ns)!Ivdw
    foo(Ns+1:2*Ns)=Ivec(1:Ns)     !Ivup
    !
    j = bjoin(foo,2*Ns)
    !
  end function flip_state



  !+------------------------------------------------------------------+
  !PURPOSE  : get the twin of a given sector (the one with opposite 
  ! quantum numbers): 
  ! nup,ndw ==> ndw,nup (spin-exchange)
  ! sz      ==> -sz     (total spin flip)
  ! n       ==> 2*Ns-n  (particle hole)
  !+------------------------------------------------------------------+
  function get_twin_sector(isector) result(jsector)
    integer,intent(in) :: isector
    integer :: jsector
    integer :: iup,idw,in,isz
    iup=getnup(isector)
    idw=getndw(isector)
    jsector=getsector(idw,iup)
  end function get_twin_sector





  !##################################################################
  !##################################################################
  !AUXILIARY COMPUTATIONAL ROUTINES ARE HERE BELOW:
  !##################################################################
  !##################################################################

  !+------------------------------------------------------------------+
  !PURPOSE  : input a state |i> and output a vector ivec(Nlevels)
  !with its binary decomposition
  !(corresponds to the decomposition of the number i-1)
  !+------------------------------------------------------------------+
  function bdecomp(i,Ntot) result(ivec)
    integer :: Ntot,ivec(Ntot),l,i
    logical :: busy
    !this is the configuration vector |1,..,Ns,Ns+1,...,Ntot>
    !obtained from binary decomposition of the state/number i\in 2^Ntot
    do l=0,Ntot-1
       busy=btest(i,l)
       ivec(l+1)=0
       if(busy)ivec(l+1)=1
    enddo
  end function bdecomp



  !+------------------------------------------------------------------+
  !PURPOSE  : input a vector ib(Nlevels) with the binary sequence 
  ! and output the corresponding state |i>
  !(corresponds to the recomposition of the number i-1)
  !+------------------------------------------------------------------+
  function bjoin(ib,Ntot) result(i)
    integer                 :: Ntot
    integer,dimension(Ntot) :: ib
    integer                 :: i,j
    i=0
    do j=0,Ntot-1
       i=i+ib(j+1)*2**j
    enddo
  end function bjoin




  !+------------------------------------------------------------------+
  !PURPOSE  : calculate the factorial of an integer N!=1.2.3...(N-1).N
  !+------------------------------------------------------------------+
  recursive function factorial(n) result(f)
    integer            :: f
    integer,intent(in) :: n
    if(n<=0)then
       f=1
    else
       f=n*factorial(n-1)
    end if
  end function factorial



  !+------------------------------------------------------------------+
  !PURPOSE  : calculate the binomial factor n1 over n2
  !+------------------------------------------------------------------+
  function binomial(n1,n2) result(nchoos)
    real(8) :: xh
    integer :: n1,n2,i
    integer nchoos
    xh = 1.d0
    if(n2<0) then
       nchoos = 0
       return
    endif
    if(n2==0) then
       nchoos = 1
       return
    endif
    do i = 1,n2
       xh = xh*dble(n1+1-i)/dble(i)
    enddo
    nchoos = int(xh + 0.5d0)
  end function binomial



  !+------------------------------------------------------------------+
  !PURPOSE : binary search of a value in an array
  !+------------------------------------------------------------------+
  recursive function binary_search(a,value) result(bsresult)
    integer,intent(in) :: a(:), value
    integer            :: bsresult, mid
    mid = size(a)/2 + 1
    if (size(a) == 0) then
       bsresult = 0        ! not found
       !stop "binary_search error: value not found"
    else if (a(mid) > value) then
       bsresult= binary_search(a(:mid-1), value)
    else if (a(mid) < value) then
       bsresult = binary_search(a(mid+1:), value)
       if (bsresult /= 0) then
          bsresult = mid + bsresult
       end if
    else
       bsresult = mid      ! SUCCESS!!
    end if
  end function binary_search






  !+------------------------------------------------------------------+
  !PURPOSE : sort array of integer using random algorithm
  !+------------------------------------------------------------------+
  subroutine sort_array(array)
    integer,dimension(:),intent(inout)      :: array
    integer,dimension(size(array))          :: order
    integer                                 :: i
    forall(i=1:size(array))order(i)=i
    call qsort_sort( array, order, 1, size(array) )
    array=order
  contains
    recursive subroutine qsort_sort( array, order, left, right )
      integer, dimension(:)                 :: array
      integer, dimension(:)                 :: order
      integer                               :: left
      integer                               :: right
      integer                               :: i
      integer                               :: last
      if ( left .ge. right ) return
      call qsort_swap( order, left, qsort_rand(left,right) )
      last = left
      do i = left+1, right
         if ( compare(array(order(i)), array(order(left)) ) .lt. 0 ) then
            last = last + 1
            call qsort_swap( order, last, i )
         endif
      enddo
      call qsort_swap( order, left, last )
      call qsort_sort( array, order, left, last-1 )
      call qsort_sort( array, order, last+1, right )
    end subroutine qsort_sort
    !---------------------------------------------!
    subroutine qsort_swap( order, first, second )
      integer, dimension(:)                 :: order
      integer                               :: first, second
      integer                               :: tmp
      tmp           = order(first)
      order(first)  = order(second)
      order(second) = tmp
    end subroutine qsort_swap
    !---------------------------------------------!
    function qsort_rand( lower, upper )
      implicit none
      integer                               :: lower, upper
      real(8)                               :: r
      integer                               :: qsort_rand
      call random_number(r)
      qsort_rand =  lower + nint(r * (upper-lower))
    end function qsort_rand
    function compare(f,g)
      integer                               :: f,g
      integer                               :: compare
      compare=1
      if(f<g)compare=-1
    end function compare
  end subroutine sort_array



end MODULE VCA_SETUP
