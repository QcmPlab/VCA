program vca_bhz_2d
  USE SCIFOR
  USE DMFT_TOOLS
  USE MPI
  USE VCA
  !
  !System parameters
  implicit none
  integer                                         :: Nlso
  integer                                         :: Nx,Ny,Ndim
  integer,dimension(2)                            :: Nkpts
  integer                                         :: ilat,jlat
  real(8)                                         :: ts,ts_var,Mh,Mh_var,lambdauser,lambdauser_var,random_1,random_2,random_3
  real(8)                                         :: M,M_var,t,t_var,lambda,lambda_var,mu,mu_var
  !Bath
  complex(8),dimension(:,:,:,:,:,:),allocatable   :: bath_h,bath_v
  !Matrices:
  real(8),allocatable,dimension(:)                :: wm,wr
  complex(8),allocatable,dimension(:,:,:,:,:,:)   :: t_prime
  complex(8),allocatable,dimension(:,:,:,:,:,:)   :: observable_matrix
  complex(8),allocatable,dimension(:,:,:,:,:,:,:) :: h_k
  complex(8),allocatable,dimension(:,:,:)         :: hk
  complex(8),allocatable,dimension(:,:,:,:,:,:,:) :: Gmats,Greal,Smats,Sreal
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
  real(8),dimension(:),allocatable                :: ts_array_x,ts_array_y,params,omega_array
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
  call parse_input_variable(Nkpts,"Nkpts",finput,default=[10,10],comment="Number of k-points along each direction")
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
  !
  call vca_read_input(trim(finput),comm)
  !
  call naming_convention()
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
  !SET CLUSTER DIMENSIONS (ASSUME SQUARE CLUSTER):
  !
  Ndim=size(Nkpts)
  Nlat=Nx*Ny
  Nlso = Nlat*Norb*Nspin
  !
  if(allocated(bath_h))deallocate(bath_h)
  if(allocated(bath_v))deallocate(bath_v)
  if(allocated(t_prime))deallocate(t_prime)
  allocate(t_prime(Nlat,Nlat,Nspin,Nspin,Norb,Norb))
  allocate(bath_h(Nlat_bath,Nlat_bath,Nspin,Nspin,Norb_bath,Norb_bath))
  allocate(bath_v(Nlat     ,Nlat_bath,Nspin,Nspin,Norb     ,Norb_bath))
  allocate(Smats(Nlat,Nlat,Nspin,Nspin,Norb,Norb,Lmats),Sreal(Nlat,Nlat,Nspin,Nspin,Norb,Norb,Lreal))
  allocate(Greal(Nlat,Nlat,Nspin,Nspin,Norb,Norb,Lreal))
  !
  !ALLOCATE VECTORS:
  !
  if(.not.allocated(wm))allocate(wm(Lmats))
  if(.not.allocated(wr))allocate(wr(Lreal))
  wm     = pi/beta*real(2*arange(1,Lmats)-1,8)
  wr     = linspace(wini,wfin,Lreal)
  !
  !INITIALIZE SOLVER:
  !
  call vca_init_solver(comm,bath_h,bath_v)
  print_impG=.true.
  print_impG0=.true.
  print_Sigma=.true.
  !
  !LATTICE PARAMETERS
  !
  !fixed to cpt configuration
  Mh_var=Mh
  Ts_var=ts
  lambdauser_var=lambdauser
  !
  t=ts
  M=mh
  lambda=lambdauser
  !
  !SOLVE INTERACTING PROBLEM:
  ! 
  if(wloop)then
    !
    allocate(ts_array_x(Nloop))
    allocate(omega_array(Nloop))
    !
    !
    ts_array_x = linspace(0.05d0,0.7d0,Nloop)

    do iloop=1,Nloop
        omega_array(iloop)=solve_vca([ts_var,Mh_var,lambdauser_var,0d0,ts_array_x(iloop)])
    enddo
    !
    call splot("sft_Omega_loopVSts.dat",ts_array_x,omega_array)
  endif
  !
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

  function solve_vca(pars) result(Omega)
    real(8),dimension(:)             :: pars
    logical                          :: invert
    real(8)                          :: Omega,E,V
    !
    !SET PARAMETERS (GLOBAL VARIABLES FOR THE DRIVER):
    !
    t_var=pars(1)
    M_var=pars(2)
    lambda_var=pars(3)
    !
    E=pars(4)
    V=pars(5)
    !
    mu_var=XMU
    mu=xmu
    !
    print*,""
    print*,"Cluster parameters:"
    print*,"t      = ",t_var
    print*,"M      = ",m_var
    print*,"lambda = ",lambda_var
    print*,"Lattice parameters:"
    print*,"t      = ",t
    print*,"M      = ",m
    print*,"lambda = ",lambda
    !
    call construct_bath(E,V)
    call generate_tcluster()
    call generate_hk()
    call vca_solve(comm,t_prime,h_k,bath_h,bath_v)
    call vca_get_sft_potential(omega)
    !
    !
  end function solve_vca
  
  !+------------------------------------------------------------------+
  !PURPOSE  : generate hopping matrices
  !+------------------------------------------------------------------+

  subroutine construct_bath(eps,v)
    real(8)                 :: eps,v
    integer                 :: i,ib,o,ob,ispin
    !
    bath_h=zero
    bath_v=zero
    !
    do ispin=1,Nspin
      do ob=1,Norb_bath
        do ib=1,Nlat_bath
          bath_h(ib,ib,ispin,ispin,ob,ob)=(-1d0)**(ob+1)*M_var+eps
          do i=1,Nlat
            bath_v(ib,ib,ispin,ispin,ob,ob)=v
          enddo
        enddo
      enddo
    enddo
  end subroutine construct_bath


  subroutine generate_tcluster()
    integer                                                       :: ilat,jlat,ispin,iorb,jorb,ind1,ind2
    complex(8),dimension(Nlat,Nlat,Nspin,Nspin,2,2)               :: t_tmp
    !
    t_prime=zero
    t_tmp=zero
    !
    do ispin=1,Nspin
      do ilat=1,Nx
        do jlat=1,Ny
          ind1=indices2N([ilat,jlat])
          t_tmp(ind1,ind1,ispin,ispin,:,:)= t_m(m_var)
          if(ilat<Nx)then
            ind2=indices2N([ilat+1,jlat])
            t_tmp(ind2,ind1,ispin,ispin,:,:)= t_x(t_var,lambda_var,ispin)
          endif
          if(ilat>1)then
            ind2=indices2N([ilat-1,jlat])
            t_tmp(ind2,ind1,ispin,ispin,:,:)= conjg(transpose(t_x(t_var,lambda_var,ispin)))
          endif
          if(jlat<Ny)then
            ind2=indices2N([ilat,jlat+1])
            t_tmp(ind2,ind1,ispin,ispin,:,:)= t_y(t_var,lambda_var)
          endif
          if(jlat>1)then
            ind2=indices2N([ilat,jlat-1])
            t_tmp(ind2,ind1,ispin,ispin,:,:)= conjg(transpose(t_y(t_var,lambda_var)))
          endif
        enddo
      enddo
    enddo
    t_prime=t_tmp(:,:,:,:,1:Norb,1:Norb)
    !
  end subroutine generate_tcluster


 function tk(kpoint) result(hopping_matrix)
    integer                                                                 :: ilat,jlat,ispin,iorb,jorb,i,j,ind1,ind2
    real(8),dimension(Ndim),intent(in)                                      :: kpoint
    complex(8),dimension(Nlat,Nlat,Nspin,Nspin,Norb,Norb)                   :: hopping_matrix
    complex(8),dimension(Nlat,Nlat,Nspin,Nspin,2,2)                         :: t_tmp
    !
    hopping_matrix=zero
    t_tmp=zero
    !
    do ilat=1,Nx
      do jlat=1,Ny
        do ispin=1,Nspin
          ind1=indices2N([ilat,jlat])
          t_tmp(ind1,ind1,ispin,ispin,:,:)= t_m(m)
          if(ilat<Nx)then
            ind2=indices2N([ilat+1,jlat])
            t_tmp(ind2,ind1,ispin,ispin,:,:)= t_x(t,lambda,ispin)
          endif
          if(ilat>1)then
            ind2=indices2N([ilat-1,jlat])
            t_tmp(ind2,ind1,ispin,ispin,:,:)= conjg(transpose(t_x(t,lambda,ispin)))
          endif
          if(jlat<Ny)then
            ind2=indices2N([ilat,jlat+1])
            t_tmp(ind2,ind1,ispin,ispin,:,:)= t_y(t,lambda)
          endif
          if(jlat>1)then
            ind2=indices2N([ilat,jlat-1])
            t_tmp(ind2,ind1,ispin,ispin,:,:)= conjg(transpose(t_y(t,lambda)))
          endif
        enddo
      enddo
    enddo
    !
    do ispin=1,Nspin
      do ilat=1,Ny
        ind1=indices2N([1,ilat])
        ind2=indices2N([Nx,ilat])
        t_tmp(ind2,ind1,ispin,ispin,:,:)=t_tmp(ind2,ind1,ispin,ispin,:,:) +conjg(transpose(t_x(t,lambda,ispin)))*exp(xi*kpoint(1)*Nx)
        t_tmp(ind1,ind2,ispin,ispin,:,:)=t_tmp(ind1,ind2,ispin,ispin,:,:) +t_x(t,lambda,ispin)*exp(-xi*kpoint(1)*Nx)
      enddo
      do ilat=1,Nx
        ind1=indices2N([ilat,1])
        ind2=indices2N([ilat,Ny])
        t_tmp(ind2,ind1,ispin,ispin,:,:)=t_tmp(ind2,ind1,ispin,ispin,:,:) +conjg(transpose(t_y(t,lambda)))*exp(xi*kpoint(2)*Ny)
        t_tmp(ind1,ind2,ispin,ispin,:,:)=t_tmp(ind1,ind2,ispin,ispin,:,:) +t_y(t,lambda)*exp(-xi*kpoint(2)*Ny)
      enddo
    enddo
    !
    hopping_matrix=t_tmp(:,:,:,:,1:Norb,1:Norb)
    ! 
  end function tk

  subroutine generate_hk()
    integer                                      :: ik,ii,ispin,iorb,unit,jj
    real(8),dimension(product(Nkpts),Ndim)       :: kgrid
    real(8),dimension(Nlso,Nlso)                 :: H0
    real(8),dimension(2)                         :: e1,e2,bk1,bk2
    real(8)                                      :: bklen
    !
    e1 = [1d0, 0d0]
    e2 = [0d0, 1d0]
    call TB_set_ei(eix=e1,eiy=e2)
    bklen=2d0*pi
    bk1=bklen*[1d0, 0d0]
    bk2=bklen*[0d0, 1d0]
    call TB_set_bk(bkx=bk1,bky=bk2)
!
    call TB_build_kgrid(Nkpts,kgrid)
    kgrid(:,1)=kgrid(:,1)/Nx
    kgrid(:,2)=kgrid(:,2)/Ny
    !
    if(allocated(h_k))deallocate(h_k)
    allocate(h_k(Nlat,Nlat,Nspin,Nspin,Norb,Norb,product(Nkpts))) 
    h_k=zero
    !
    do ik=1,product(Nkpts)
        !
        h_k(:,:,:,:,:,:,ik)=tk(kgrid(ik,:))
        !
    enddo
    !
  end subroutine generate_hk


  !+------------------------------------------------------------------+
  !H block functions
  !+------------------------------------------------------------------+

  function t_m(mass) result(tmpmat)
    complex(8),dimension(2,2) :: tmpmat
    real(8)                   :: mass
    !
    tmpmat=mass*pauli_sigma_z
    !
  end function t_m

  function t_x(hop1,hop2,spinsign) result(tmpmat)
    complex(8),dimension(2,2) :: tmpmat
    real(8)                   :: hop1,hop2,sz
    integer                   :: spinsign
    !
    tmpmat=zero
    sz=(-1.d0)**(spinsign+1)
    tmpmat=-hop1*pauli_sigma_z+0.5d0*sz*xi*hop2*pauli_sigma_x
    !
  end function t_x

  function t_y(hop1,hop2) result(tmpmat)
    complex(8),dimension(2,2) :: tmpmat
    real(8)                   :: hop1,hop2
    !
    tmpmat=zero
    tmpmat=-hop1*pauli_sigma_z
    tmpmat(1,2)=-hop2*0.5d0
    tmpmat(2,1)=hop2*0.5d0
    !
  end function t_y


  !+------------------------------------------------------------------+
  !Auxilliary functions
  !+------------------------------------------------------------------+

   function indices2N(indices) result(N)
      integer,dimension(2)         :: indices
      integer                      :: N,i
      !
      !
      N=Nx*(indices(2)-1)+indices(1)
   end function indices2N

   function N2indices(N) result(indices) 
      integer,dimension(2)         :: indices
      integer                      :: N,i
      !
      indices(1)=mod(N,Nx)
      if(indices(1)==0)then
         indices(1)=Nx
         indices(2)=(N-Nx)/Nx+1
      else
         indices(2)=N/Nx+1
      endif
   end function N2indices

   subroutine naming_convention()
      integer                       :: i,j
      integer,dimension(Nx,Ny)      :: matrix
      !
      do j=1,Ny
         do i=1,Nx
            matrix(i,j)=indices2N([i,j])
         enddo
      enddo
      !
      write(LOGfile,"(A)")"The unique index of each site (on the cartesian plane) is as follows:"
      write(LOGfile,"(A)")" "
      do j=1,Ny
         write(LOGfile,"(20(I2,2x))")(matrix(i,Ny+1-j),i =1,Nx)
      enddo
      write(LOGfile,"(A)")" "
   end subroutine naming_convention


   function lso2nnn(Hlso) result(Hnnn)
      complex(8),dimension(Nlat*Nspin*Norb,Nlat*Nspin*Norb) :: Hlso
      complex(8),dimension(Nlat,Nlat,Nspin,Nspin,Norb,Norb) :: Hnnn
      integer                                               :: ilat,jlat
      integer                                               :: iorb,jorb
      integer                                               :: ispin,jspin
      integer                                               :: is,js
      Hnnn=zero
      do ilat=1,Nlat
         do jlat=1,Nlat
            do ispin=1,Nspin
               do jspin=1,Nspin
                  do iorb=1,Norb
                     do jorb=1,Norb
                        is = iorb + (ilat-1)*Norb + (ispin-1)*Norb*Nlat
                        js = jorb + (jlat-1)*Norb + (jspin-1)*Norb*Nlat
                        Hnnn(ilat,jlat,ispin,jspin,iorb,jorb) = Hlso(is,js)
                     enddo
                  enddo
               enddo
            enddo
         enddo
      enddo
   end function lso2nnn


   function nnn2lso(Hnnn) result(Hlso)
      complex(8),dimension(Nlat,Nlat,Nspin,Nspin,Norb,Norb) :: Hnnn
      complex(8),dimension(Nlat*Nspin*Norb,Nlat*Nspin*Norb) :: Hlso
      integer                                               :: ilat,jlat
      integer                                               :: iorb,jorb
      integer                                               :: ispin,jspin
      integer                                               :: is,js
      Hlso=zero
      do ilat=1,Nlat
         do jlat=1,Nlat
            do ispin=1,Nspin
               do jspin=1,Nspin
                  do iorb=1,Norb
                     do jorb=1,Norb
                        is = iorb + (ilat-1)*Norb + (ispin-1)*Norb*Nlat
                        js = jorb + (jlat-1)*Norb + (jspin-1)*Norb*Nlat
                        Hlso(is,js) = Hnnn(ilat,jlat,ispin,jspin,iorb,jorb)
                     enddo
                  enddo
               enddo
            enddo
         enddo
      enddo
   end function nnn2lso
   
   function so2nn(Hso) result(Hnn)
     complex(8),dimension(Nspin*Norb,Nspin*Norb) :: Hso
     complex(8),dimension(Nspin,Nspin,Norb,Norb) :: Hnn
     integer                                     :: iorb,ispin,is
     integer                                     :: jorb,jspin,js
     Hnn=zero
     do ispin=1,Nspin
        do jspin=1,Nspin
           do iorb=1,Norb
              do jorb=1,Norb
                  is = iorb + (ispin-1)*Norb  !spin-orbit stride
                 js = jorb + (jspin-1)*Norb  !spin-orbit stride
                 Hnn(ispin,jspin,iorb,jorb) = Hso(is,js)
              enddo
           enddo
        enddo
     enddo
   end function so2nn
   !
   function nn2so(Hnn) result(Hso)
     complex(8),dimension(Nspin,Nspin,Norb,Norb) :: Hnn
     complex(8),dimension(Nspin*Norb,Nspin*Norb) :: Hso
     integer                                     :: iorb,ispin,is
     integer                                     :: jorb,jspin,js
     Hso=zero
     do ispin=1,Nspin
        do jspin=1,Nspin
           do iorb=1,Norb
              do jorb=1,Norb
                 is = iorb + (ispin-1)*Norb  !spin-orbit stride
                 js = jorb + (jspin-1)*Norb  !spin-orbit stride
                 Hso(is,js) = Hnn(ispin,jspin,iorb,jorb)
              enddo
           enddo
        enddo
     enddo
   end function nn2so


end program vca_bhz_2d









!OLD FUNCTIONS

!  function solve_vca_square(tij) result(Omega)
!    real(8)                      :: tij
!    real(8)                      :: Vij,Eij
!    real(8)                      :: Omega
!    !
!   t=0.5d0
!    t_var=tij
!    !
!    mu=0.d0*t
!    M=(2.d0*t)*Mh
!    lambda=(2.d0*t)*0.3d0
!    !  
!    mu_var=0.d0*t
!    M_var=(2.d0*t)*Mh
!    lambda_var=(2.d0*t)*0.3d0
!    !
!    print*,""
!    print*,"------ Doing for ",tij," ------"
!    call generate_tcluster()
!    call generate_hk()
!    call vca_solve(comm,t_prime,h_k)
!    call vca_get_sft_potential(omega)
!    print*,""
!    print*,"------ DONE ------"
!    print*,""
!    !
!  end function solve_vca_square
!
!  subroutine solve_Htop(kpath_)
!    integer                                  :: i,j
!    integer                                  :: Npts,Nkpath
!    type(rgb_color),dimension(:),allocatable :: colors
!    real(8),dimension(:,:),optional          :: kpath_
!    real(8),dimension(:,:),allocatable       :: kpath
!    character(len=64)                        :: file
!    
!    !This routine build the H(k) along the GXMG path in BZ,
!    !Hk(k) is constructed along this path.
!    Nkpath=100
!    !
!    if(present(kpath_))then
!       if(master)write(LOGfile,*)"Build H(k) BHZ along a given path:"
!       Npts = size(kpath_,1)
!       allocate(kpath(Npts,size(kpath_,2)))
!       kpath=kpath_
!    else
!       if(master)write(LOGfile,*)"Build H(k) BHZ along the path GXMG:"
!       Npts = 4
!       allocate(kpath(Npts,2))
!       kpath(1,:)=[0.d0,0.d0]
!       kpath(2,:)=[pi,pi]
!       kpath(3,:)=[pi,0.d0]
!       kpath(4,:)=[0.d0,0.d0]
!    endif
!   allocate(colors(Nspin*Norb))
!   colors=[red1,blue1,red1,blue1]
!   !
!   file="Eig_Htop_kSigma.nint"
!   if(master) call TB_Solve_model(hk_bhz,Nspin*Norb,kpath,Nkpath,&   
!         colors_name=colors,&
!         points_name=[character(len=20) :: 'G', 'M', 'X', 'G'],&
!         file=reg(file))
!   file="Eig_Htop_localSigma.nint"
!   if(master) call TB_Solve_model(hk_bhz_local,Nspin*Norb,kpath,Nkpath,&   
!         colors_name=colors,&
!         points_name=[character(len=20) :: 'G', 'M', 'X', 'G'],&
!         file=reg(file))
!  end subroutine solve_Htop
!
!
!  function tk_wrapper(kpoint,N) result(hopping_matrix_lso)
!    integer                                         :: N
!    real(8),dimension(:)                            :: kpoint
!    complex(8),dimension(N,N)                       :: hopping_matrix_lso
!    if(N.ne.Nlso)stop "error N: wrong dimension"
!    !
!    hopping_matrix_lso=vca_nnn2lso_reshape(tk(kpoint),Nlat,Nspin,Norb)
!    !
!  end function tk_wrapper
!
!
!  subroutine solve_hk_GXMG(kpath_)
!    integer                                  :: i,j
!    integer                                  :: Npts,Nkpath
!    type(rgb_color),dimension(:),allocatable :: colors
!    real(8),dimension(:,:),optional          :: kpath_
!    real(8),dimension(:,:),allocatable       :: kpath
!    character(len=64)                        :: file
!    
!    !This routine build the H(k) along the GXMG path in BZ,
!    !Hk(k) is constructed along this path.
!    Nkpath=500 
!    !
!    if(present(kpath_))then
!       if(master)write(LOGfile,*)"Build H(k) BHZ along a given path:"
!       Npts = size(kpath_,1)
!       allocate(kpath(Npts,size(kpath_,2)))
!       kpath=kpath_
!       file="Eig_path.nint"
!    else
!       if(master)write(LOGfile,*)"Build H(k) BHZ along the path GXMG:"
!       Npts = 4
!       allocate(kpath(Npts,3))
!       kpath(1,:)=[0.d0,0.d0,0.0d0]
!       kpath(2,:)=[pi/Nx,pi/ny,0.0d0]
!       kpath(3,:)=[pi/Nx,0.d0,0.0d0]
!       kpath(4,:)=[0.d0,0.d0,0.0d0]
!       file="Eigenbands.nint"
!    endif
!   allocate(colors(Nlso))
!   colors=[red1,blue1,red1,blue1]
!   
!   if(master) call TB_Solve_model(tk_wrapper,Nlso,kpath,Nkpath,&
!         colors_name=colors,&
!         points_name=[character(len=20) :: 'G', 'M', 'X', 'G'],&
!         file=reg(file))
!  end subroutine solve_hk_GXMG
!
!
!
!  function hk_bhz(kpoint,N) result(hopping_matrix_lso)
!    integer                                         :: N
!    real(8),dimension(:)                            :: kpoint
!    real(8),dimension(Ndim)                         :: kpoint_
!    complex(8),dimension(N,N)                       :: hopping_matrix_lso
!    real(8)                                         :: energy_scale
!    if(N.ne.Nspin*Norb)stop "error N: wrong dimension"
!    !
!    energy_scale=2.d0*t_var
!    hopping_matrix_lso=zero
!    hopping_matrix_lso= (M_var-energy_scale*(cos(kpoint(1))+cos(kpoint(2))))*kron_pauli( pauli_sigma_0, pauli_tau_z)+&
!                                                                 lambda_var*sin(kpoint(1))*kron_pauli( pauli_sigma_z, pauli_tau_x)+&
!                                                                 lambda_var*sin(kpoint(2))*kron_pauli( pauli_sigma_0, pauli_tau_y)
!    if(scheme=="g")then
!      call build_sigma_g_scheme(kpoint)
!    else
!      call periodize_sigma_scheme(kpoint)
!    endif
!    !hopping_matrix_lso=vca_nn2so_reshape(gfmats_periodized(:,:,:,:,1),Nspin,Norb)
!    !call inv(hopping_matrix_lso)
!    hopping_matrix_lso=hopping_matrix_lso+DREAL(smats_periodized_lso(:,:,1))
!    hopping_matrix_lso=matmul(Zrenorm(smats_periodized_lso(:,:,1)),DREAL(hopping_matrix_lso))
!    !
!  end function hk_bhz
!
!  function hk_bhz_local(kpoint,N) result(hopping_matrix_lso)
!    integer                                         :: N
!    real(8),dimension(:)                            :: kpoint
!    complex(8),dimension(N,N)                       :: hopping_matrix_lso
!    real(8)                                         :: energy_scale
!    if(N.ne.Nspin*Norb)stop "error N: wrong dimension"
!    !
!    energy_scale=2.d0*t_var
!    hopping_matrix_lso=zero
!    hopping_matrix_lso= (M_var-energy_scale*(cos(kpoint(1))+cos(kpoint(2))))*kron_pauli( pauli_sigma_0, pauli_tau_z)+&
!                                                                 lambda_var*sin(kpoint(1))*kron_pauli( pauli_sigma_z, pauli_tau_x)+&
!                                                                 lambda_var*sin(kpoint(2))*kron_pauli( pauli_sigma_0, pauli_tau_y)
!    hopping_matrix_lso=hopping_matrix_lso+DREAL(vca_nn2so_reshape(smats_local(:,:,:,:,1),Nspin,Norb))
!    hopping_matrix_lso=matmul(Zrenorm(vca_nn2so_reshape(smats_local(:,:,:,:,1),Nspin,Norb)),hopping_matrix_lso)
!    !
!  end function hk_bhz_local
