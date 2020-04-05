!--------------------------------------------------------------
!Solve the 2d BHZ model with N bath sites, the only variational
!parameter is currently V.
!______________________________________________________________

program vca_bhz_2d_bath
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
  real(8)                                         :: ts,ts_var,Mh,Mh_var,lambdauser,lambdauser_var
  real(8)                                         :: M,M_var,t,t_var,lambda,lambda_var,mu,mu_var
  real(8)                                         :: bath_e,bath_v
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
  !
  !SET CLUSTER DIMENSIONS (ASSUME SQUARE CLUSTER):
  !
  Ndim=size(Nkpts)
  Nlat=Nx*Ny
  Nlso = Nlat*Norb*Nspin
  !
  !SET BATH
  !
  Nb=vca_get_bath_dimension()
  allocate(Bath(Nb))
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
  
  !CUSTOM OBSERVABLE: KINETIC ENERGY
  allocate(observable_matrix(Nlat,Nlat,Nspin,Nspin,Norb,Norb))
  observable_matrix=zero
  observable_matrix(1,1,1,1,1,1)=one
  call init_custom_observables(1,product(Nkpts))
  call add_custom_observable("test",observable_matrix)
  
  MULTIMAX=.false.
  !
  !SOLVE INTERACTING PROBLEM:
  !
  if(wmin)then
    !
    bath_v=0.4
    if(master)print*,"Guess:",bath_v
    call  brent_(solve_vca_single,bath_v,[0.05d0,0.7d0])
    if(master)print*,"Result ts : ",bath_v
    omegadummy=solve_vca_single(bath_v)
    if(master)write(*,"(A,F15.9,A,3F15.9)")bold_green("FOUND STATIONARY POINT "),omegadummy,bold_green(" AT V = "),bath_v
    !
  elseif(wloop)then
    !
    allocate(ts_array_x(Nloop))
    allocate(omega_grid(Nloop,Nloop))
    !
    ts_array_x = linspace(0.05d0,0.7d0,Nloop)
    do iloop=1,Nloop
        omega_grid(iloop,1)=solve_vca_single(ts_array_x(iloop))
    enddo
    !
    call splot("sft_Omega_loopVSts.dat",ts_array_x,omega_grid(:,1))
    !
  endif
  !
  !PRINT LOCAL GF AND SIGMA
  !
  !call solve_Htop_new()
  !call get_local_gf()
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


  function solve_vca_single(x) result(Omega_)
    real(8)                   :: x,Omega_
    !
    Omega_=solve_vca_multi([t,M,lambda,2*M,x])
    !
  end function solve_vca_single

  function solve_vca_multi(pars) result(Omega)
    integer                      :: ix,iy,ik,iq,iz
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
    deltae=pars(4)
    Vij=pars(5)
    !
    mu_var=0.d0*t_var
    Eij=0.d0
    !
    !if(NBATH>1)then
    !  tmp=linspace(0.d0,deltae,Nbath)
    !else
    !  tmp=0.5d0*deltae
    !endif
    !
    do iy=1,Nbath
      !evector(iy)=Eij+tmp(iy)-0.5d0*deltae
      evector(iy)=Eij
      vvector(iy)=Vij
    enddo
    do ix=1,Nx
      do iq=1,Ny
        iz=indices2N([ix,iq])
        do iy=1,Nspin
          do ik=1,Norb
            call set_bath_component(bath,iz,iy,ik,e_component=evector*((-1)**(ik+2)))
            call set_bath_component(bath,iz,iy,ik,v_component=vvector)
          enddo
        enddo
      enddo
    enddo

    !
    if(master)then
       print*,""
       print*,"Variational parameters:"
       print*,"t      = ",t_var
       print*,"M      = ",m_var
       print*,"lambda = ",lambda_var
       print*,"Lattice parameters:"
       print*,"t      = ",t
       print*,"M      = ",m
       print*,"lambda = ",lambda
    endif
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


  !+------------------------------------------------------------------+
  !PURPOSE  : generate hopping matrices (assume Ny=1)
  !+------------------------------------------------------------------+


  subroutine generate_tcluster()
    integer                                      :: ilat,jlat,ispin,iorb,jorb,ind1,ind2
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
      do ilat=1,Nx
        do jlat=1,Ny
          ind1=indices2N([ilat,jlat])
          t_prime(ind1,ind1,ispin,ispin,:,:)= t_m(m_var)
          if(ilat<Nx)then
            ind2=indices2N([ilat+1,jlat])
            t_prime(ind1,ind2,ispin,ispin,:,:)= t_x(t_var,lambda_var,ispin)
          endif
          if(ilat>1)then
            ind2=indices2N([ilat-1,jlat])
            t_prime(ind1,ind2,ispin,ispin,:,:)= dconjg(transpose(t_x(t_var,lambda_var,ispin)))
          endif
          if(jlat<Ny)then
            ind2=indices2N([ilat,jlat+1])
            t_prime(ind1,ind2,ispin,ispin,:,:)= t_y(t_var,lambda_var)
          endif
          if(jlat>1)then
            ind2=indices2N([ilat,jlat-1])
            t_prime(ind1,ind2,ispin,ispin,:,:)= transpose(t_y(t_var,lambda_var))
          endif
        enddo
      enddo
    enddo
    !
  end subroutine generate_tcluster


 function tk(kpoint) result(hopping_matrix)
    integer                                                                 :: ilat,jlat,ispin,iorb,jorb,i,j,ind1,ind2
    real(8),dimension(Ndim),intent(in)                                      :: kpoint
    complex(8),dimension(Nlat,Nlat,Nspin,Nspin,Norb,Norb)                   :: hopping_matrix
    !
    hopping_matrix=zero
    !
    do ispin=1,Nspin
      do ilat=1,Nx
        do jlat=1,Ny
          ind1=indices2N([ilat,jlat])
          hopping_matrix(ind1,ind1,ispin,ispin,:,:)= t_m(m)
          if(ilat<Nx)then
            ind2=indices2N([ilat+1,jlat])
            hopping_matrix(ind1,ind2,ispin,ispin,:,:)= t_x(t,lambda,ispin)
          endif
          if(ilat>1)then
            ind2=indices2N([ilat-1,jlat])
            hopping_matrix(ind1,ind2,ispin,ispin,:,:)= dconjg(transpose(t_x(t,lambda,ispin)))
          endif
          if(jlat<Ny)then
            ind2=indices2N([ilat,jlat+1])
            hopping_matrix(ind1,ind2,ispin,ispin,:,:)= t_y(t,lambda)
          endif
          if(jlat>1)then
            ind2=indices2N([ilat,jlat-1])
            hopping_matrix(ind1,ind2,ispin,ispin,:,:)= transpose(t_y(t,lambda))
          endif
        enddo
      enddo
    enddo
    !
    !
    do ispin=1,Nspin
      do ilat=1,Ny
        ind1=indices2N([1,ilat])
        ind2=indices2N([Nx,ilat])
        hopping_matrix(ind1,ind2,ispin,ispin,:,:)=hopping_matrix(ind1,ind2,ispin,ispin,:,:) + dconjg(transpose(t_x(t,lambda,ispin)))*exp(xi*kpoint(1)*Nx)
        hopping_matrix(ind2,ind1,ispin,ispin,:,:)=hopping_matrix(ind2,ind1,ispin,ispin,:,:) + t_x(t,lambda,ispin)*exp(-xi*kpoint(1)*Nx)
      enddo
      do ilat =1,Nx
        ind1=indices2N([ilat,1])
        ind2=indices2N([ilat,Ny])
        hopping_matrix(ind1,ind2,ispin,ispin,:,:)=hopping_matrix(ind1,ind2,ispin,ispin,:,:) + transpose(t_y(t,lambda))*exp(xi*kpoint(2)*Ny)
        hopping_matrix(ind2,ind1,ispin,ispin,:,:)=hopping_matrix(ind2,ind1,ispin,ispin,:,:) + t_y(t,lambda)*exp(-xi*kpoint(2)*Ny)
      enddo
    enddo
    !
    ! 
  end function tk

  subroutine generate_hk()
    integer                                      :: ik,ii,ispin,iorb,unit,jj
    real(8),dimension(product(Nkpts),Ndim)          :: kgrid
    !
    call TB_build_kgrid(Nkpts,kgrid)
    !Reduced Brillouin Zone
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

 !BRENT MINIMIZING FUNCTION

  subroutine brent_(func,xmin,brack,tol,niter)
    interface
       function func(x)
         real(8) :: x
         real(8) :: func
       end function func
    end interface
    real(8),intent(inout)         :: xmin
    real(8),dimension(:),optional :: brack
    real(8),optional              :: tol
    integer,optional              :: niter
    real(8)                       :: tol_
    integer                       :: niter_
    integer                       :: iter
    real(8)                       :: ax,xx,bx,fa,fx,fb,fret
    !
    tol_=1d-6;if(present(tol))tol_=tol
    Niter_=200;if(present(Niter))Niter_=Niter
    !
    if(present(brack))then
       select case(size(brack))
       case(1)
          stop "Brent error: calling brent with size(brack)==1. None or two points are necessary."
       case(2)
          ax = brack(1)
          xx = brack(2)
       case (3)
          ax = brack(1)
          xx = brack(2)
          bx = brack(3)
       end select
    else
       ax=0d0
       xx=1d0
    endif
    fret=brent_optimize(ax,xx,bx,func,tol_,niter_,xmin)
  end subroutine brent_
  !



  function brent_optimize(ax,bx,cx,func,tol,itmax,xmin)
    real(8), intent(in)  :: ax,bx,cx,tol
    real(8), intent(out) :: xmin
    real(8)              :: brent_optimize
    integer              :: itmax
    real(8), parameter   :: cgold=0.3819660d0,zeps=1.0d-3*epsilon(ax)
    integer              :: iter
    real(8)              :: a,b,d,e,etemp,fu,fv,fw,fx,p,q,r,tol1,tol2,u,v,w,x,xm
    interface
       function func(x)
         real(8) :: x
         real(8) :: func
       end function func
    end interface
    a=min(ax,cx)
    b=max(ax,cx)
    v=bx
    w=v
    x=v
    e=0.d0
    fx=func(x)
    fv=fx
    fw=fx
    do iter=1,itmax
       xm=0.5d0*(a+b)
       tol1=tol*abs(x)+zeps
       tol2=2.0*tol1
       if (abs(x-xm) <= (tol2-0.5d0*(b-a))) then
          xmin=x
          brent_optimize=fx
          return
       end if
       if (abs(e) > tol1) then
          r=(x-w)*(fx-fv)
          q=(x-v)*(fx-fw)
          p=(x-v)*q-(x-w)*r
          q=2.d0*(q-r)
          if (q > 0.d0) p=-p
          q=abs(q)
          etemp=e
          e=d
          if (abs(p) >= abs(0.5d0*q*etemp) .or. &
               p <= q*(a-x) .or. p >= q*(b-x)) then
             e=merge(a-x,b-x, x >= xm )
             d=cgold*e
          else
             d=p/q
             u=x+d
             if (u-a < tol2 .or. b-u < tol2) d=sign(tol1,xm-x)
          end if
       else
          e=merge(a-x,b-x, x >= xm )
          d=cgold*e
       end if
       u=merge(x+d,x+sign(tol1,d), abs(d) >= tol1 )
       fu=func(u)
       if (fu <= fx) then
          if (u >= x) then
             a=x
          else
             b=x
          end if
          call shft(v,w,x,u)
          call shft(fv,fw,fx,fu)
       else
          if (u < x) then
             a=u
          else
             b=u
          end if
          if (fu <= fw .or. w == x) then
             v=w
             fv=fw
             w=u
             fw=fu
          else if (fu <= fv .or. v == x .or. v == w) then
             v=u
             fv=fu
          end if
       end if
    end do
    !pause 'brent: exceed maximum iterations'

  end function brent_optimize

    subroutine shft(a,b,c,d)
      real(8), intent(out) :: a
      real(8), intent(inout) :: b,c
      real(8), intent(in) :: d
      a=b
      b=c
      c=d
    end subroutine shft



end program vca_bhz_2d_bath









