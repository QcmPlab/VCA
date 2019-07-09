program vca_bhz_2d
  USE SCIFOR
  USE DMFT_TOOLS
  USE MPI
  USE VCA
  !
  !System parameters
  implicit none
  integer                                         :: Nlso
  integer                                         :: Nx,Ny
  integer                                         :: ilat,jlat
  real(8)                                         :: ts,ts_var,Mh,Mh_var,lambdauser,lambdauser_var
  real(8)                                         :: M,M_var,t,t_var,lambda,lambda_var,mu,mu_var
  !Bath
  integer                                         :: Nb
  real(8),allocatable                             :: Bath(:)
  !Matrices:
  real(8),allocatable,dimension(:)                :: wm,wr
  complex(8),allocatable,dimension(:,:,:,:,:,:)   :: t_prime
  complex(8),allocatable,dimension(:,:,:,:,:,:)   :: observable_matrix
  complex(8),allocatable,dimension(:,:,:,:,:,:,:) :: h_k
  complex(8),allocatable,dimension(:,:,:,:,:)     :: gfmats_local,gfmats_periodized     ![Nspin][Nspin][Norb][Norb][L]
  complex(8),allocatable,dimension(:,:,:,:,:)     :: gfreal_local,gfreal_periodized     ![Nspin][Nspin][Norb][Norb][L]
  complex(8),allocatable,dimension(:,:,:,:,:)     :: Smats_local,Smats_periodized       ![Nspin][Nspin][Norb][Norb][L]
  complex(8),allocatable,dimension(:,:,:,:,:)     :: Sreal_local,Sreal_periodized       ![Nspin][Nspin][Norb][Norb][L]
  complex(8),allocatable,dimension(:,:,:)         :: Smats_periodized_lso,Sreal_periodized_lso    
  !Utility variables:
  integer                                         :: unit
  integer                                         :: comm,rank
  integer                                         :: iloop,jloop,nloop
  integer                                         :: iii,jjj,kkk
  logical                                         :: master,wloop,wmin,MULTIMAX
  logical                                         :: usez
  logical                                         :: print_mats,print_real
  character(len=6)                                :: scheme
  character(len=16)                               :: finput
  real(8)                                         :: omegadummy,observable_dummy
  real(8),dimension(:),allocatable                :: ts_array_x,ts_array_y,params
  real(8),dimension(:,:),allocatable              :: omega_grid
  real(8),allocatable,dimension(:,:)              :: kgrid_test,kpath_test
  real(8),dimension(3)                            :: df
  !
  !MPI INIT
  !
  call init_MPI()
  comm = MPI_COMM_WORLD
  call StartMsg_MPI(comm)
  rank = get_Rank_MPI(comm)
  master = get_Master_MPI(comm)
  !
  !PARSE INPUT VARIABLES
  !
  call parse_cmd_variable(finput,"FINPUT",default='inputVCA.conf')
  call parse_input_variable(ts,"ts",finput,default=0.5d0,comment="Hopping parameter (units of epsilon)")
  call parse_input_variable(Mh,"Mh",finput,default=3d0,comment="Field splitting (units of epsilon)")
  call parse_input_variable(lambdauser,"lambda",finput,default=0.3d0,comment="Spin/orbit coupling (units of epsilon)")
  call parse_input_variable(ts_var,"ts_Var",finput,default=0.5d0,comment="variational hopping parameter (units of epsilon)")
  call parse_input_variable(Mh_var,"Mh_Var",finput,default=3d0,comment="variational field splitting (units of epsilon)")
  call parse_input_variable(lambdauser_var,"lambda_var",finput,default=0.3d0,comment="variational spin/orbit coupling (units of epsilon)")
  call parse_input_variable(Nx,"Nx",finput,default=2,comment="Number of sites along X")
  call parse_input_variable(Ny,"Ny",finput,default=2,comment="Number of sites along Y")
  call parse_input_variable(nloop,"NLOOP",finput,default=100)
  call parse_input_variable(wloop,"WLOOP",finput,default=.false.)
  call parse_input_variable(wmin,"WMIN",finput,default=.false.,comment="T: includes global minimization")
  call parse_input_variable(scheme,"SCHEME",finput,default="g")
  call parse_input_variable(print_mats,"PRINT_MATS",finput,default=.true.)
  call parse_input_variable(print_real,"PRINT_REAL",finput,default=.true.)
  call parse_input_variable(scheme,"SCHEME",finput,default="g")
  call parse_input_variable(usez,"USEZ",finput,default=.false.)
  !
  call vca_read_input(trim(finput),comm)
  !
  !
  !Add DMFT CTRL Variables:
  call add_ctrl_var(Nlat,"NLAT")
  call add_ctrl_var(Norb,"norb")
  call add_ctrl_var(Nspin,"nspin")
  call add_ctrl_var(beta,"beta")
  call add_ctrl_var(xmu,"xmu")
  call add_ctrl_var(wini,'wini')
  call add_ctrl_var(wfin,'wfin')
  call add_ctrl_var(eps,"eps")
  !
  !SET BATH
  !
  Nb=vca_get_bath_dimension()
  allocate(Bath(Nb))
  !
  !SET CLUSTER DIMENSIONS (ASSUME SQUARE CLUSTER):
  !
  Nlat=Nx**Ndim
  Ny=Nx
  Nlso = Nlat*Norb*Nspin
  !
  !SET LATTICE PARAMETERS (GLOBAL VARIABLES FOR THE DRIVER):
  !
  t=ts
  M=(2.d0*t)*Mh
  lambda=(2.d0*t)*lambdauser
  mu=0.d0*t
  !
  !ALLOCATE VECTORS:
  !
  if(.not.allocated(wm))allocate(wm(Lmats))
  if(.not.allocated(wr))allocate(wr(Lreal))
  if(.not.allocated(params))allocate(params(3))
  wm     = pi/beta*real(2*arange(1,Lmats)-1,8)
  wr     = linspace(wini,wfin,Lreal)
  !
  !INITIALIZE SOLVER:
  !
  call vca_init_solver(comm,bath)
  print_impG=.false.
  print_impG0=.false.
  print_Sigma=.false.
  print_observables=.false.
  MULTIMAX=.false.
  !
  !SOLVE INTERACTING PROBLEM:
  ! 
  if(wmin)then
    !
    !
    !INITIALIZE VARIABLES TO THE LATTICE VALUES
    !
    params=[t,M,lambda]
    !
    call minimize_parameters(params,0.5d0)
    !call fmin_brent(params,0.2d0)
    !
    print_Sigma=.true.
    print_observables=.true.
    omegadummy=solve_vca_multi(params)
    !
    write(*,"(A,F15.9,A,3F15.9)")bold_green("FOUND STATIONARY POINT "),omegadummy,bold_green(" AT "),t_var,m_var,lambda_var
    write(*,"(A)")""
    !
    call solve_Htop_new()
    !
  elseif(wloop)then
    !
    allocate(ts_array_x(Nloop))
    allocate(omega_grid(Nloop,Nloop))
    !
    ts_array_x = linspace(0.01d0,1.0d0,Nloop)
    do iloop=1,Nloop
        omega_grid(iloop,1)=solve_vca_multi([ts_var,Mh_var,lambda,ts_array_x(iloop)])
    enddo
    !
    call splot("sft_Omega_loopVSts.dat",ts_array_x,omega_grid(:,1))
    !!
  else
    print_observables=.true.
    omegadummy=solve_vca_multi([ts_var,Mh_Var,lambdauser_var])
    !print*,"calculate gradient"
    !call fdjac_1n_func(solve_vca_multi,[ts_var,Mh_Var,lambdauser_var],df)
    !print*,"gradient is", df
    !
    write(*,"(A,F15.9,A,3F15.9)")bold_green("OMEGA IS "),omegadummy,bold_green(" AT "),ts_var,Mh_Var,lambdauser_var
    !
    allocate(observable_matrix(Nlat,Nlat,Nspin,Nspin,Norb,Norb))
    !
    !SET OBSERVABLE MATRIX: AS AN EXAMPLE, THE TWO ORBITAL OCCUPATIONS
    !
    observable_matrix=zero
    do iii=1,Nlat
          observable_matrix(iii,iii,1,1,1,1)=1.d0
          observable_matrix(iii,iii,2,2,1,1)=1.d0
    enddo
    call observables_lattice(observable_matrix,observable_dummy)
    print*,"User-requested observable with value store is ",observable_dummy
    !
    observable_matrix=zero
    do iii=1,Nlat
          observable_matrix(iii,iii,1,1,2,2)=1.d0
          observable_matrix(iii,iii,2,2,2,2)=1.d0
    enddo
    call observables_lattice(observable_matrix)
  endif
  !
  !PRINT LOCAL GF AND SIGMA
  !
  call solve_Htop_new()
  !
  if(allocated(wm))deallocate(wm)
  if(allocated(wr))deallocate(wr)
  if(allocated(params))deallocate(params)
  !
  call finalize_MPI()
  !
contains

  !+------------------------------------------------------------------+
  !PURPOSE  : solve the model
  !+------------------------------------------------------------------+

  function solve_vca_multi(pars) result(Omega)
    integer                      :: iy,ik
    real(8)                      :: Vij,Eij,deltae
    real(8),dimension(:)         :: pars
    real(8),dimension(Nbath)     :: evector,vvector,tmp
    logical                      :: invert
    real(8)                      :: Omega
    !
    !SET VARIATIONAL PARAMETERS (GLOBAL VARIABLES FOR THE DRIVER):
    !
    t_var=pars(1)  
    M_var=pars(2)
    lambda_var=pars(3)
    mu_var=0.d0*t_var
    Eij=0.d0
    deltae=0.1d0
    Vij=pars(4)
    if(NBATH>1)then
      tmp=linspace(0.d0,deltae,Nbath)
    else
      tmp=0.5d0*deltae
    endif
    !
    do iy=1,Nbath
      evector(iy)=Eij+tmp(iy)-0.5d0*deltae
      vvector(iy)=Vij
    enddo
    do iy=1,Nspin
      do ik=1,Norb
        call set_bath_component(bath,1,iy,ik,e_component=evector)
        call set_bath_component(bath,1,iy,ik,v_component=vvector)
      enddo
    enddo

    !
    print*,""
    print*,"Variational parameters:"
    print*,"t      = ",t_var
    print*,"M      = ",m_var
    print*,"lambda = ",lambda_var
    print*,"Bath E = ",Eij
    print*,"Bath V = ",Vij
    print*,"Lattice parameters:"
    print*,"t      = ",t
    print*,"M      = ",m
    print*,"lambda = ",lambda
    call generate_tcluster()
    call generate_hk()
    call vca_solve(comm,t_prime,h_k,bath)
    call vca_get_sft_potential(omega)
    !
    if(MULTIMAX)omega=-omega
    !    
    print*,""
    !
  end function solve_vca_multi


  function solve_vca_t(variable) result(Omega)
    logical                      :: invert
    real(8)                      :: Omega,variable
    !
    !SET VARIATIONAL PARAMETERS (GLOBAL VARIABLES FOR THE DRIVER):
    !
    t_var=variable 
    !
    print*,""
    print*,"Variational parameters:"
    print*,"t      = ",t_var
    print*,"M      = ",m_var
    print*,"lambda = ",lambda_var
    print*,"Lattice parameters:"
    print*,"t      = ",t
    print*,"M      = ",m
    print*,"lambda = ",lambda
    call generate_tcluster()
    call generate_hk()
    call vca_solve(comm,t_prime,h_k)
    call vca_get_sft_potential(omega)
    !
    if(MULTIMAX)omega=-omega
    !    
    print*,""
    !
  end function solve_vca_t

  function solve_vca_m(variable) result(Omega)
    logical                      :: invert
    real(8)                      :: Omega,variable
    !
    !SET VARIATIONAL PARAMETERS (GLOBAL VARIABLES FOR THE DRIVER):
    !
    M_var=variable
    !
    print*,""
    print*,"Variational parameters:"
    print*,"t      = ",t_var
    print*,"M      = ",m_var
    print*,"lambda = ",lambda_var
    print*,"Lattice parameters:"
    print*,"t      = ",t
    print*,"M      = ",m
    print*,"lambda = ",lambda
    call generate_tcluster()
    call generate_hk()
    call vca_solve(comm,t_prime,h_k)
    call vca_get_sft_potential(omega)
    !
    if(MULTIMAX)omega=-omega
    !    
    print*,""
    !
  end function solve_vca_m

  function solve_vca_l(variable) result(Omega)
    logical                      :: invert
    real(8)                      :: Omega,variable
    !
    !SET VARIATIONAL PARAMETERS (GLOBAL VARIABLES FOR THE DRIVER):
    !
    lambda_var=variable
    !
    print*,""
    print*,"Variational parameters:"
    print*,"t      = ",t_var
    print*,"M      = ",m_var
    print*,"lambda = ",lambda_var
    print*,"Lattice parameters:"
    print*,"t      = ",t
    print*,"M      = ",m
    print*,"lambda = ",lambda
    call generate_tcluster()
    call generate_hk()
    call vca_solve(comm,t_prime,h_k)
    call vca_get_sft_potential(omega)
    !
    if(MULTIMAX)omega=-omega
    !    
    print*,""
    !
  end function solve_vca_l

 
  !+------------------------------------------------------------------+
  !PURPOSE:  multidimensional finder of stationary points
  !+------------------------------------------------------------------+
  subroutine minimize_parameters(v,radius)
    real(8),dimension(:),allocatable          :: v,l,lold,u,uold,parvec
    integer,dimension(:),allocatable          :: nbd
    real(8)                                   :: radius     
    integer                                   :: i,iprint_         
    !
    allocate ( nbd(size(v)), parvec(size(v)), l(size(v)), u(size(v)), lold(size(v)), uold(size(v)) )
    !
    !INITIALIZE FLAGS
    !
    iprint_=1
    !if(verbose .ge. 1)iprint_=1
    MULTIMAX=.false.
    !
    !INITIALIZE PARAMETERS VECTOR AND BOUNDARIES
    !
    parvec=v
    !
    do i=1,size(v)
      nbd(i) = 2
      l(i)   = parvec(i)-radius
      u(i)   = parvec(i)+radius
    enddo
    lold=l
    uold=u
    !
    write(*,"(A)")""
    write(*,"(A)")bold_red("LOOKING FOR MINIMUMS")
    !
    !FIND LOCAL MINIMA
    !
    call fmin_bfgs(solve_vca_multi,parvec,l,u,nbd,factr=1.d8,iprint=iprint_,nloop=Nloop)
    !
    !RESET MINIMIZER FOR VARIABLES ON THE BORDER
    !
    do i=1,size(v)
      if((abs(parvec(i)-l(i)) .lt. 1.d-6) .or. (abs(parvec(i)-u(i)) .lt. 1.d-6))then
        parvec(i)=v(i)
      MULTIMAX=.true.
      else
        l(i)=parvec(i)
        u(i)=parvec(i)
      endif
    enddo
    !
    !IF NEEDED FIND MAXIMUMS/SADDLE POINTS
    !
    if (MULTIMAX) then
      write(*,"(A)")""
      write(*,"(A)")bold_red("LOOKING FOR MAXIMUMS")
      call fmin_bfgs(solve_vca_multi,parvec,l,u,nbd,factr=1.d8,iprint=iprint_,nloop=Nloop)
      do i=1,size(v)
        if((abs(parvec(i)-lold(i)) .lt. 1.d-6) .or. (abs(parvec(i)-uold(i)) .lt. 1.d-6))stop "STATIONARY POINT NOT FOUND!"
      enddo
      !
      MULTIMAX=.false.
      v=parvec
    endif
    !
    deallocate(nbd,parvec,l,u,lold,uold)
    !
  end subroutine minimize_parameters



  subroutine fmin_brent(v,radius)
    real(8),dimension(3)                      :: v
    real(8)                                   :: radius 
    integer                                   :: i        
    !
    !
    !INITIALIZE FLAGS
    !
    t_var=v(1)
    M_var=v(2)
    lambda_var=v(3)
    mu_var=0.d0*t_var
    MULTIMAX=.false.
    !
    print*,"MINIMIZE T"
    !
    call  brent(solve_vca_t,v(1),[t_var-radius,t_var+radius])   
    if((abs(t_var-v(1)-radius) .lt. 1.d-5) .or. (abs(t_var-v(1)+radius) .lt. 1.d-5) )then
      MULTIMAX=.true.
      call  brent(solve_vca_t,v(1),[t_var-radius,t_var+radius])
      if((abs(t_var-v(1)-radius) .lt. 1.d-5) .or. (abs(t_var-v(1)+radius) .lt. 1.d-5) )STOP "error on minimizing t"
    endif
    t_var=v(1)
    !
    print*,"MINIMIZE M"
    !
    MULTIMAX=.true.
    call  brent(solve_vca_m,v(2),[M_var-radius,M_var+radius])   
    if((abs(M_var-v(2)-radius) .lt. 1.d-5) .or. (abs(M_var-v(2)+radius) .lt. 1.d-5) )then
      MULTIMAX=.false.
      call  brent(solve_vca_m,v(2),[M_var-radius,M_var+radius])
      if((abs(M_var-v(2)-radius) .lt. 1.d-5) .or. (abs(M_var-v(2)+radius) .lt. 1.d-5) )STOP "error on minimizing t"
    endif
    M_var=v(2)
    !
    print*,"MINIMIZE LAMBDA"
    !
    MULTIMAX=.true.
    call  brent(solve_vca_l,v(3),[lambda_var-radius,lambda_var+radius])   
    if((abs(lambda_var-v(3)-radius) .lt. 1.d-5) .or. (abs(lambda_var-v(3)+radius) .lt. 1.d-5) )then
      MULTIMAX=.false.
      call  brent(solve_vca_l,v(3),[lambda_var-radius,lambda_var+radius])
      if((abs(lambda_var-v(3)-radius) .lt. 1.d-5) .or. (abs(lambda_var-v(3)+radius) .lt. 1.d-5) )STOP "error on minimizing t"
    endif
    lambda_var=v(3)
    MULTIMAX=.false.
    call fdjac_1n_func(solve_vca_multi,v,df)
    print*,"norm of gradient is", dot_product(df,df)
    !
  end subroutine fmin_brent


  subroutine fdjac_1n_func(funcv,x,fjac,epsfcn)
    implicit none
    interface 
       function funcv(x)
         implicit none
         real(8),dimension(:) :: x
         real(8)              :: funcv
       end function funcv
    end interface
    integer          ::  n
    real(8)          ::  x(:)
    real(8)          ::  fvec
    real(8)          ::  fjac(size(x))
    real(8),optional ::  epsfcn
    real(8)          ::  eps,eps_
    real(8)          ::  epsmch
    real(8)          ::  h,temp
    real(8)          ::  wa1
    real(8)          ::  wa2
    integer          :: i,j,k
    real(8)          :: df_eps=1.d-3
    n=size(x)
    fjac=0.d0
    eps_= df_eps; if(present(epsfcn))eps_=epsfcn
    epsmch = epsilon(epsmch)
    eps  = sqrt(max(eps_,epsmch))
    !  Evaluate the function
    fvec = funcv(x)
    do j=1,n
       temp = x(j)
       h    = eps*abs(temp)
       if(h==0.d0) h = eps
       x(j) = temp + h
       wa1  = funcv(x)
       x(j) = temp
       fjac(j) = (wa1 - fvec)/h
    enddo
  end subroutine fdjac_1n_func


  !+------------------------------------------------------------------+
  !PURPOSE  : generate hopping matrices
  !+------------------------------------------------------------------+


  subroutine generate_tcluster()
    integer                                      :: ilat,jlat,ispin,iorb,jorb
    complex(8),dimension(Nlso,Nlso)              :: H0
    character(len=64)                            :: file_
    integer                                      :: unit
    file_ = "tcluster_matrix.dat"
    !
    if(allocated(t_prime))deallocate(t_prime)
    allocate(t_prime(Nlat,Nlat,Nspin,Nspin,Norb,Norb))
    t_prime=zero
    !
    do ispin=1,Nspin
      do ilat=1,Nlat
        t_prime(ilat,ilat,ispin,ispin,:,:)= t_m(m_var)
        if(ilat<Nx)then
          t_prime(ilat,ilat,ispin,ispin,:,:)= t_x(t_var,lambda_var,ispin)
        endif
        if(ilat>1)then
          t_prime(ilat,ilat,ispin,ispin,:,:)= dconjg(transpose(t_x(t_var,lambda_var,ispin)))
        endif
      enddo
    enddo
  end subroutine generate_tcluster


 function tk(kpoint) result(hopping_matrix)
    integer                                                                 :: ilat,jlat,ispin,iorb,jorb,i,j
    real(8),dimension(Ndim),intent(in)                                      :: kpoint
    complex(8),dimension(Nlat,Nlat,Nspin,Nspin,Norb,Norb)                   :: hopping_matrix
    !
    hopping_matrix=zero
    !
    do ispin=1,Nspin
      do ilat=1,Nx
        hopping_matrix(ilat,ilat,ispin,ispin,:,:)= t_m(m)
        if(ilat<Nx)then
          hopping_matrix(ilat,ilat,ispin,ispin,:,:)= t_x(t,lambda,ispin)
        endif
        if(ilat>1)then
          hopping_matrix(ilat,ilat,ispin,ispin,:,:)= dconjg(transpose(t_x(t,lambda,ispin)))
        endif
      enddo
    enddo
    !
    !
    do ispin=1,Nspin
      do ilat=1,Nx
        hopping_matrix(ilat,ilat,ispin,ispin,:,:)=hopping_matrix(ilat,ilat,ispin,ispin,:,:) + dconjg(transpose(t_x(t,lambda,ispin)))*exp(xi*kpoint(1)*Nx)
        hopping_matrix(ilat,ilat,ispin,ispin,:,:)=hopping_matrix(ilat,ilat,ispin,ispin,:,:) + t_x(t,lambda,ispin)*exp(-xi*kpoint(1)*Nx)
        !
        hopping_matrix(ilat,ilat,ispin,ispin,:,:)=hopping_matrix(ilat,ilat,ispin,ispin,:,:) + transpose(t_y(t,lambda))*exp(xi*kpoint(2)*Ny)
        hopping_matrix(ilat,ilat,ispin,ispin,:,:)=hopping_matrix(ilat,ilat,ispin,ispin,:,:) + t_y(t,lambda)*exp(-xi*kpoint(2)*Ny)
      enddo
    enddo
    !
    ! 
  end function tk

  subroutine generate_hk()
    integer                                      :: ik,ii,ispin,iorb,unit,jj
    real(8),dimension(Nkpts**Ndim,Ndim)          :: kgrid
    real(8),dimension(Nlso,Nlso)                 :: H0
    character(len=64)                            :: file_
    file_ = "tlattice_matrix.dat"
    !
    call TB_build_kgrid([Nkpts,Nkpts],kgrid)
    !Reduced Brillouin Zone
    kgrid=kgrid/Nx 
    !
    if(allocated(h_k))deallocate(h_k)
    allocate(h_k(Nlat,Nlat,Nspin,Nspin,Norb,Norb,Nkpts**ndim)) 
    h_k=zero
    !
    do ik=1,Nkpts**ndim
        !
        h_k(:,:,:,:,:,:,ik)=tk(kgrid(ik,:))
        !
    enddo
  end subroutine generate_hk


!AUXILLIARY HOPPING MATRIX CONSTRUCTORS

  function t_m(mass) result(tmpmat)
    complex(8),dimension(Norb,Norb) :: tmpmat
    real(8)                         :: mass
    !
    tmpmat=zero
    tmpmat=mass*pauli_sigma_z
    !
  end function t_m

  function t_x(hop1,hop2,spinsign) result(tmpmat)
    complex(8),dimension(Norb,Norb) :: tmpmat
    real(8)                         :: hop1,hop2,sz
    integer                         :: spinsign
    !
    tmpmat=zero
    sz=(-1.d0)**(spinsign+1)
    tmpmat=-hop1*pauli_sigma_z+0.5d0*sz*xi*hop2*pauli_sigma_x
    !
  end function t_x

  function t_y(hop1,hop2) result(tmpmat)
    complex(8),dimension(Norb,Norb) :: tmpmat
    real(8)                         :: hop1,hop2
    !
    tmpmat=zero
    tmpmat=-hop1*pauli_sigma_z
    tmpmat(1,2)=-hop2*0.5d0
    tmpmat(2,1)=hop2*0.5d0
    !
  end function t_y


  !+------------------------------------------------------------------+
  !PURPOSE  : PRINT HAMILTONIAN ALONG PATH
  !+------------------------------------------------------------------+

  function hk_bhz_clusterbase(kpoint,N) result(hopping_matrix_lso)
    integer                                                       :: N,ilat,jlat
    real(8),dimension(:)                                          :: kpoint
    real(8),dimension(Ndim)                                       :: kpoint_
    complex(8),dimension(N,N)                                     :: hopping_matrix_lso
    complex(8),dimension(Nlat,Nlat,Nspin,Nspin,Norb,Norb)         :: hopping_matrix_big
    complex(8),dimension(Nlat,Nlat,Nspin,Nspin,Norb,Norb,Lmats)   :: tmpSigmaMat
    real(8)                                                       :: energy_scale
    if(N.ne.Nlat*Nspin*Norb)stop "error N: wrong dimension"
    !
    hopping_matrix_lso=zero
    hopping_matrix_big=zero
    !
    !
    call vca_get_sigma_matsubara(tmpSigmaMat)
    hopping_matrix_big=tk(kpoint)+DREAL(TmpSigmaMat(:,:,:,:,:,:,1))
    !
    !
    hopping_matrix_lso=vca_nnn2lso_reshape(Hopping_matrix_big,Nlat,Nspin,Norb)
    !
  end function hk_bhz_clusterbase



  subroutine solve_Htop_new(kpath_)
    integer                                  :: i,j
    integer                                  :: Npts,Nkpath
    type(rgb_color),dimension(:),allocatable :: colors
    real(8),dimension(:,:),optional          :: kpath_
    real(8),dimension(:,:),allocatable       :: kpath
    character(len=64)                        :: file
    !
    Nkpath=100
    !
    if(present(kpath_))then
       if(master)write(LOGfile,*)"Build H(k) BHZ along a given path:"
       Npts = size(kpath_,1)
       allocate(kpath(Npts,size(kpath_,2)))
       kpath=kpath_
    else
       if(master)write(LOGfile,*)"Build H(k) BHZ along the path GXMG:"
       Npts = 4
       allocate(kpath(Npts,2))
       kpath(1,:)=[0.d0,0.d0]
       kpath(2,:)=[pi/Nx,pi/Nx]
       kpath(3,:)=[pi/Nx,0.d0]
       kpath(4,:)=[0.d0,0.d0]
    endif
    allocate(colors(Nlat*Nspin*Norb))
    colors = gray99
    colors(1) = red1
    colors(2) = blue1
    colors(3) = red1
    colors(4) = blue1
   !
   file="Eig_Htop_clusterbase.nint"
   if(master) call TB_Solve_model(hk_bhz_clusterbase,Nlat*Nspin*Norb,kpath,Nkpath,&   
         colors_name=colors,&
         points_name=[character(len=20) :: 'G', 'M', 'X', 'G'],&
         file=reg(file))
  end subroutine solve_Htop_new



  !+------------------------------------------------------------------+
  !Auxilliary functions
  !+------------------------------------------------------------------+


 
  !SET THE BATH DELTA FUNCTION

  function set_delta(freq,vps,eps) result(DELTA)
    complex(8),allocatable,dimension(:,:,:,:,:,:)               :: DELTA ![Nlat][Nlat][Nspin][Nspin][Norb][Norb]
    complex(8)                                                  :: freq
    real(8),dimension(:)                                        :: vps,eps
    integer                                                     :: ispin,iorb,ilat
    !
    allocate(DELTA(Nlat,Nlat,Nspin,Nspin,Norb,Norb))
    DELTA=zero
    !
    if (Nbath .ne. 0)then
      do ilat=1,Nlat
        do ispin=1,Nspin
           do iorb=1,Norb
             DELTA(ilat,ilat,ispin,ispin,iorb,iorb)=sum( vps(:)*vps(:)/(freq - eps(:)+XMU) )
           enddo
        enddo
      enddo
    endif
  end function set_delta


end program vca_bhz_2d









