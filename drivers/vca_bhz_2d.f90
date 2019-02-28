program vca_square_bath
  USE SCIFOR
  USE DMFT_TOOLS
  USE MPI
  USE VCA
  !
  implicit none
  integer                                         :: Nlso,Nsys
  integer                                         :: ilat,jlat
  integer                                         :: iloop
  integer                                         :: ix,iy,ik
  logical                                         :: converged
  real(8)                                         :: wband
  !Bath
  real(8),allocatable                             :: Bath(:)
  integer                                         :: Nb
  !The local hybridization function:
  real(8),dimension(:,:),allocatable              :: Tsys,Tref,Vmat,Htb,Mmat,dens
  real(8),allocatable,dimension(:)                :: wm,wr
  complex(8)                                      :: iw
  complex(8),allocatable,dimension(:,:,:,:,:,:,:) :: Gmats,Greal
  complex(8),allocatable,dimension(:,:,:,:,:,:)   :: t_prime
  complex(8),allocatable,dimension(:,:,:,:,:,:,:) :: h_k
  complex(8),allocatable,dimension(:,:,:,:,:,:,:) :: t_k
  character(len=16)                               :: finput
  real(8)                                         :: ts,hopping
  integer                                         :: Nx,Ny,Lx,Ly,Rx,Ry,SAMPLING
  integer                                         :: unit
  integer                                         :: comm,rank
  logical                                         :: master,wloop,wmin
  integer                                         :: nloop
  character(len=6)                                :: scheme
  !MINIMIZATION ROUTINES
  real(8)                                         :: ts_dummy1,ts_dummy2,DUMMY1,DUMMY2,dummy3,ts_dummy3,omegadummy
  !VARIATIONAL PARAMETERS
  real(8)                                         :: M,M_var,t,t_var,lambda,lambda_var,mu,mu_var
  !
  real(8),dimension(:),allocatable                :: ts_array,omega_array
  integer,dimension(1)                            :: min_loc
  complex(8),allocatable,dimension(:,:,:,:,:)     :: gfmats_periodized ![Nspin][Nspin][Norb][Norb][Lmats]
  complex(8),allocatable,dimension(:,:,:,:,:)     :: gfreal_periodized ![Nspin][Nspin][Norb][Norb][Lreal]
  complex(8),allocatable,dimension(:,:,:,:,:)     :: Smats_periodized         ![Nspin][Nspin][Norb][Norb][Lmats]
  complex(8),allocatable,dimension(:,:,:,:,:)     :: Sreal_periodized         ![Nspin][Nspin][Norb][Norb][Lreal]
  complex(8),allocatable,dimension(:,:,:,:,:)     :: gtest_mats,gtest_Real,sigmatest_mats,sigmatest_real
  real(8),allocatable,dimension(:,:)              :: kgrid_test,kpath_test
  integer                                         :: iii,jjj,kkk

  call init_MPI()
  comm = MPI_COMM_WORLD
  call StartMsg_MPI(comm)
  rank = get_Rank_MPI(comm)
  master = get_Master_MPI(comm)

  call parse_cmd_variable(finput,"FINPUT",default='inputVCA.conf')
  call parse_cmd_variable(SAMPLING,"SAMPLING",default=100)
  call parse_input_variable(ts,"ts",finput,default=1d0)
  call parse_input_variable(Nx,"Nx",finput,default=2,comment="Number of sites along X")
  call parse_input_variable(Ny,"Ny",finput,default=2,comment="Number of sites along Y")
  call parse_input_variable(nloop,"NLOOP",finput,default=100)
  call parse_input_variable(hopping,"HOPPING",finput,default=1.d0)
  call parse_input_variable(wloop,"WLOOP",finput,default=.false.)
  call parse_input_variable(scheme,"SCHEME",finput,default="g")
  call parse_input_variable(wmin,"wmin",finput,default=.false.,comment="T: includes global minimization")
  !
  call vca_read_input(trim(finput),comm)
  !
  !
  Nlat=Nx**Ndim
  Ny=Nx
  Nlso = Nlat*Norb*Nspin
  !
  M=1.5d0
  M_var=1.5d0
  t=1.d0
  t_var=1.d0
  lambda=0.3d0
  lambda_var=1.d0
  mu=1.d0
  mu_var=1.d0
  mu=0.d0
  mu_var=0.d0

  !Add DMFT CTRL Variables:
  call add_ctrl_var(Nlat,"NLAT")
  call add_ctrl_var(Norb,"norb")
  call add_ctrl_var(Nspin,"nspin")
  call add_ctrl_var(beta,"beta")
  call add_ctrl_var(xmu,"xmu")
  call add_ctrl_var(wini,'wini')
  call add_ctrl_var(wfin,'wfin')
  call add_ctrl_var(eps,"eps")


  if(.not.allocated(wm))allocate(wm(Lmats))
  if(.not.allocated(wr))allocate(wr(Lreal))
  wm     = pi/beta*real(2*arange(1,Lmats)-1,8)
  wr     = linspace(wini,wfin,Lreal)

  
  !Nb=vca_get_bath_dimension()
  !allocate(Bath(Nb))
  !call vca_init_solver(comm,bath)
  call vca_init_solver(comm)

  !call generate_tcluster()
  !call generate_hk()
  !call solve_hk_GXMG()

  !STOP

  if(wmin)then
    allocate(ts_array(Nloop))
    allocate(omega_array(Nloop))
    ts_array(1)=0.7/Nloop*1
    omega_array(1)=solve_vca_square(ts_array(1))
    ts_array(2)=0.7/Nloop*2
    omega_array(2)=solve_vca_square(ts_array(2))
    do iii =3,Nloop
      ts_array(iii)=0.7/Nloop*iii
      omega_array(iii)=solve_vca_square(ts_array(iii))
      if(omega_array(iii-1)<omega_array(iii) .and. omega_array(iii-1)<omega_array(iii-2))then
        print*,"FOUND MININUM, REFINE"
        call  brent_(solve_vca_square,ts,[ts_array(iii-2),ts_array(iii-1),ts_array(iii)],tol=1.d-9) 
        omegadummy=solve_vca_square(ts)
        open(free_unit(unit),file="minmax.vca",position='append')
        write(unit,*)ts,omegadummy
        close(unit)
        print*,"GOING ON"
      elseif(omega_array(iii-1)>omega_array(iii) .and. omega_array(iii-1)>omega_array(iii-2))then
        print*,"FOUND MAXIMUM, REFINE"
        call  brent_(solve_vca_square_max,ts,[ts_array(iii-2),ts_array(iii-1),ts_array(iii)],tol=1.d-9) 
        omegadummy=solve_vca_square(ts)
        open(free_unit(unit),file="minmax.vca",position='append')
        write(unit,*)ts,omegadummy
        close(unit)
        print*,"GOING ON"
      endif
    enddo
     call splot("sft_Omega_loopVSts.dat",ts_array,omega_array)
     !call get_local_gf()
     !call get_Akw()
  else if(wloop)then
    allocate(ts_array(Nloop))
    allocate(omega_array(Nloop))
    ts_array = linspace(0.05d0,0.7d0,Nloop)
    !
    do iloop=1,Nloop
       omega_array(iloop)=solve_vca_square(ts_array(iloop))
    enddo
    !
    call splot("sft_Omega_loopVSts.dat",ts_array,omega_array)
    min_loc = minloc(omega_array)
    write(800,*)min_loc,ts_array(min_loc(1)),omega_array(min_loc(1))
  else
    do ix=1,Nlat
      do iy=1,Nspin
        do ik=1,Norb
          call set_bath_component(bath,ix,iy,ik,e_component=[Uloc(1)/2])
          call set_bath_component(bath,ix,iy,ik,v_component=[ts])
        enddo
      enddo
    enddo
    call generate_tcluster()
    call generate_hk()
    call vca_solve(comm,t_prime,h_k,bath)
    !call get_local_gf()
    !call get_Akw()
  endif

  !MINIMIZATION:

  if(allocated(wm))deallocate(wm)
  if(allocated(wr))deallocate(wr)

  call finalize_MPI()





contains

  !+------------------------------------------------------------------+
  !PURPOSE  : solve the model
  !+------------------------------------------------------------------+

  function solve_vca_square_max(tij_) result(Omega_)
    real(8)                      :: tij_,Omega_

    !
      Omega_=-solve_vca_square(tij_)
    !
  end function solve_vca_square_max



  function solve_vca_square(tij) result(Omega)
    real(8)                      :: tij
    real(8)                      :: Vij,Eij
    real(8)                      :: Omega
    !
    !
    !Vij=tij
    !Eij=Uloc(1)/2d0
    t_var=tij
    print*,""
    print*,"------ Doing for ",tij," ------"
    call generate_tcluster()
    call generate_hk()
    !BATH VARIATIONAL SETUP
    !do ix=1,Nlat
    !  do iy=1,Nspin
    !    do ik=1,Norb       
    !      call set_bath_component(bath,ix,iy,ik,e_component=[Eij])
    !      call set_bath_component(bath,ix,iy,ik,v_component=[Vij])
    !    enddo
    !  enddo
    !enddo
    !call vca_solve(comm,t_prime,h_k,bath)
    call vca_solve(comm,t_prime,h_k)
    call vca_get_sft_potential(omega)
    print*,""
    print*,"------ DONE ------"
    print*,""
    !
  end function solve_vca_square


  !+------------------------------------------------------------------+
  !PURPOSE  : generate hopping matrices
  !+------------------------------------------------------------------+


  subroutine generate_tcluster()
    integer                                      :: ilat,jlat,ispin,iorb,jorb,ind1,ind2,asdlol
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
        do jlat=1,Nx
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
          if(jlat<Nx)then
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
    H0=vca_nnn2lso_reshape(t_prime,Nlat,Nspin,Norb)
    !
    open(free_unit(unit),file=trim(file_))
    do ilat=1,Nlat*Nspin*Norb
       write(unit,"(5000(F5.2,1x))")(REAL(H0(ilat,jlat)),jlat=1,Nlat*Nspin*Norb)
    enddo
    write(unit,*)"                  "
    do ilat=1,Nlat*Nspin*Norb
       write(unit,"(5000(F5.2,1x))")(IMAG(H0(ilat,jlat)),jlat=1,Nlat*Nspin*Norb)
    enddo
    close(unit)
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
        do jlat=1,Nx
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
          if(jlat<Nx)then
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
    do ilat=1,Nx
      do ispin=1,Nspin
        do iorb=1,Norb
          do jorb=1,Norb
            ind1=indices2N([1,ilat])
            ind2=indices2N([Nx,ilat])
            hopping_matrix(ind1,ind2,ispin,ispin,:,:)=hopping_matrix(ind1,ind2,ispin,ispin,:,:) + dconjg(transpose(t_x(t,lambda,ispin)))*exp(xi*kpoint(1)*Nx)
            hopping_matrix(ind2,ind1,ispin,ispin,:,:)=hopping_matrix(ind2,ind1,ispin,ispin,:,:) + t_x(t,lambda,ispin)*exp(-xi*kpoint(1)*Nx)
            !
            ind1=indices2N([ilat,1])
            ind2=indices2N([ilat,Ny])
            hopping_matrix(ind1,ind2,ispin,ispin,:,:)=hopping_matrix(ind1,ind2,ispin,ispin,:,:) + transpose(t_y(t,lambda))*exp(xi*kpoint(2)*Ny)
            hopping_matrix(ind2,ind1,ispin,ispin,:,:)=hopping_matrix(ind2,ind1,ispin,ispin,:,:) + t_y(t,lambda)*exp(-xi*kpoint(2)*Ny)
          enddo
        enddo
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
    kgrid=kgrid/Nx !!!!!DIVIDI OGNI K PER NUMERO SITI in quella direzione, RBZ
    if(allocated(h_k))deallocate(h_k)
    allocate(h_k(Nlat,Nlat,Nspin,Nspin,Norb,Norb,Nkpts**ndim)) 
    h_k=zero
    !
    do ik=1,Nkpts**ndim
        !
        h_k(:,:,:,:,:,:,ik)=tk(kgrid(ik,:))
        !
    enddo
    H0=vca_nnn2lso_reshape(tk([0.3d0,0.6d0]),Nlat,Nspin,Norb)
    !
    open(free_unit(unit),file=trim(file_))
    do ilat=1,Nlat*Nspin*Norb
       write(unit,"(5000(F5.2,1x))")(H0(ilat,jlat),jlat=1,Nlat*Nspin*Norb)
    enddo
    close(unit)    
  end subroutine generate_hk


!hopping matrices

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


  function tk_wrapper(kpoint,N) result(hopping_matrix_lso)
    integer                                         :: N
    real(8),dimension(:)                            :: kpoint
    complex(8),dimension(N,N)                       :: hopping_matrix_lso
    if(N.ne.Nlso)stop "error N: wrong dimension"
    !
    hopping_matrix_lso=vca_nnn2lso_reshape(tk(kpoint),Nlat,Nspin,Norb)
    !
  end function tk_wrapper


  subroutine solve_hk_GXMG(kpath_)
    integer                                  :: i,j
    integer                                  :: Npts,Nkpath
    type(rgb_color),dimension(:),allocatable :: colors
    real(8),dimension(:,:),optional          :: kpath_
    real(8),dimension(:,:),allocatable       :: kpath
    character(len=64)                        :: file
    
    !This routine build the H(k) along the GXMG path in BZ,
    !Hk(k) is constructed along this path.
    Nkpath=500 
    !
    if(present(kpath_))then
       if(master)write(LOGfile,*)"Build H(k) BHZ along a given path:"
       Npts = size(kpath_,1)
       allocate(kpath(Npts,size(kpath_,2)))
       kpath=kpath_
       file="Eig_path.nint"
    else
       if(master)write(LOGfile,*)"Build H(k) BHZ along the path GXMG:"
       Npts = 4
       allocate(kpath(Npts,3))
       kpath(1,:)=kpoint_Gamma/2
       kpath(2,:)=kpoint_M1/2
       kpath(3,:)=kpoint_X1/2
       kpath(4,:)=kpoint_Gamma/2
       file="Eigenbands.nint"
    endif
   allocate(colors(Nlso))
   colors = black
   
   if(master) call TB_Solve_model(tk_wrapper,Nlso,kpath,Nkpath,&
         colors_name=colors,&
         points_name=[character(len=20) :: 'G', 'M', 'X', 'G'],&
         file=reg(file))
  end subroutine solve_hk_GXMG


 
  !+------------------------------------------------------------------+
  !Auxilliary functions
  !+------------------------------------------------------------------+


  function indices2N(indices) result(N)
    integer,dimension(Ndim)      :: indices
    integer                      :: N,i
    !
    !
    N=1
    N=N+(indices(1)-1)*Nx+(indices(2)-1)
  end function indices2N

  function N2indices(N_) result(indices) ! FIXME: only for 2d
    integer,dimension(Ndim)      :: indices
    integer                      :: N,i,N_
    !
    N=N_-1
    indices(2)=mod(N,Nx)+1
    indices(1)=N/Nx+1
  end function N2indices


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
    tol_=1d-9;if(present(tol))tol_=tol
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



end program vca_square_bath

