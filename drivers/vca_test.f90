program vca_test
  USE SCIFOR
  USE DMFT_TOOLS
  USE MPI
  USE VCA
  !
  implicit none
  integer                                         :: Nlso,Nsys
  integer                                         :: ilat,jlat
  integer                                         :: iloop
  integer                                         :: ix,iy
  logical                                         :: converged
  real(8)                                         :: wband

  !The local hybridization function:
  real(8),dimension(:,:),allocatable              :: Tsys,Tref,Vmat,Htb,Mmat,dens
  real(8),allocatable,dimension(:)                :: wm,wr
  complex(8)                                      :: iw
  complex(8),allocatable,dimension(:,:,:,:,:,:,:) :: Gmats,Greal
  complex(8),allocatable,dimension(:,:,:,:,:,:)   :: t_prime
  complex(8),allocatable,dimension(:,:,:,:,:,:,:) :: t_k
  character(len=16)                               :: finput
  real(8)                                         :: ts
  integer                                         :: Nx,Ny,Lx,Ly,Rx,Ry
  integer                                         :: unit
  integer                                         :: comm,rank
  logical                                         :: master,wloop
  integer                                         :: nloop
   !FIX!!!!!!!!!!!!!!!
  real(8)                                         :: mu,t,t_var,mu_var
  real(8),dimension(:),allocatable   :: ts_array,omega_array
  integer,dimension(1)               :: min_loc

  call init_MPI()
  comm = MPI_COMM_WORLD
  call StartMsg_MPI(comm)
  rank = get_Rank_MPI(comm)
  master = get_Master_MPI(comm)

  call parse_cmd_variable(finput,"FINPUT",default='inputVCA.conf')
  call parse_input_variable(ts,"ts",finput,default=1d0)
  call parse_input_variable(Rx,"Rx",finput,default=1,comment="Ratio L/Lc=Rx along X-directions, aka # of copies along X")
  call parse_input_variable(Ry,"Ry",finput,default=1,comment="Ratio L/Lc=Ry along Y-directions, aka # of copies along Y")
  call parse_input_variable(nloop,"NLOOP",finput,default=100)
  call parse_input_variable(wloop,"WLOOP",finput,default=.true.)
  !
  call vca_read_input(trim(finput),comm)


  !if(Norb/=1)stop "Norb != 1"
  !if(Nspin/=1)stop "Nspin != 1"
  !
  Nx=2
  Ny=2
  Ndim=2
  Lx   = Rx*Nx
  Ly   = Ry*Ny
  !
  Nlat = Nx*Ny
  Nsys = Lx*Ly
  !
  Nlso = Nlat*Norb*Nspin

  t_var=1.0d0
  t=1.0d0
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


  allocate(Tsys(Nsys,Nsys))
  allocate(Tref(Nsys,Nsys))
  allocate(Vmat(Nsys,Nsys))
  allocate(Htb(Nlat,Nlat))


  call vca_init_solver(comm)

  !htb=Htb_square_lattice(2,2,1.d0)
  !t_prime=vca_lso2nnn_reshape(htb,Nlat,Nspin,Norb)
  call generate_hcluster()
  call generate_t_k()
  call vca_solve(comm,t_prime,t_k)
  !call print_2DLattice_Structure(dcmplx(htb),[Nx,Ny],1,1,file="Htb_vecchio")
  !call print_2DLattice_Structure(vca_nnn2lso_reshape(t_prime,Nlat,Nspin,Norb),[Nx,Ny],1,1,file="Htb_nuovo")

  if(wloop)then
    allocate(ts_array(Nloop))
    allocate(omega_array(Nloop))
    ts_array = linspace(1.d-2,5.d0,Nloop)
    !
    do iloop=1,Nloop
       omega_array(iloop)=solve_vca_square(ts_array(iloop))
    enddo
    !
    call splot("sft_Omega_loopVSts.dat",ts_array,omega_array)
    min_loc = minloc(omega_array)
    write(800,*)min_loc,ts_array(min_loc(1)),omega_array(min_loc(1))
  endif


  call finalize_MPI()





contains

  !+------------------------------------------------------------------+
  !PURPOSE  : solve the model
  !+------------------------------------------------------------------+



  function solve_vca_square(tij) result(Omega)
    real(8)                      :: tij
    real(8)                      :: Omega
    !
    !
    t_var=tij
    print*,""
    print*,"------ t = ",t_var,"------"
    call generate_hcluster()
    call generate_t_k()
    call vca_solve(comm,t_prime,t_k)
    call vca_get_sft_potential(omega)
    print*,""
    print*,"------ DONE ------"
    print*,""
    !
  end function solve_vca_square


  !+------------------------------------------------------------------+
  !PURPOSE  : generate test hopping matrices
  !+------------------------------------------------------------------+

  subroutine generate_hcluster()   ! FIXME: CAREFUL, I WILL NEED BOTH t' AND THE FULL H (W, interaction terms)
    integer                          :: ii,ispin,iorb,i,j
    character(len=64)                :: file_
    real(8),dimension(4,4)           :: H0
    integer                          :: unit
    file_ = "Hcluster_matrix.dat"
    !
    if(allocated(t_prime))deallocate(t_prime)
    allocate(t_prime(Nlat,Nlat,Nspin,Nspin,Norb,Norb))
    t_prime=zero
    !
    t_prime(1,1,1,1,1,1)= -mu_var
    t_prime(2,2,1,1,1,1)= -mu_var
    t_prime(3,3,1,1,1,1)= -mu_var
    t_prime(4,4,1,1,1,1)= -mu_var
    !
    t_prime(1,2,1,1,1,1)= -t_var
    t_prime(2,1,1,1,1,1)= -t_var
    !
    t_prime(2,3,1,1,1,1)= -t_var
    t_prime(3,2,1,1,1,1)= -t_var
    !
    t_prime(3,4,1,1,1,1)= -t_var
    t_prime(4,3,1,1,1,1)= -t_var
    !
    t_prime(4,1,1,1,1,1)= -t_var
    t_prime(1,4,1,1,1,1)= -t_var
    !
    H0=vca_nnn2lso_reshape(t_prime,Nlat,Nspin,Norb)
    open(free_unit(unit),file=trim(file_))
    do i=1,4
       write(unit,"(5000(F5.2,1x))")(H0(i,j),j=1,4)
    enddo
    close(unit)
  end subroutine generate_hcluster


  subroutine generate_t_k()
    integer                         :: ik,ii,ispin,iorb,unit
    real(8),dimension(2)            :: kpoint
    character(len=64)                :: file_
    real(8),dimension(Nkpts**2,2)   :: kgrid
    !
    file_ = "Hk_matrix.dat"
    call TB_build_kgrid([Nkpts,Nkpts],kgrid)
    ! 
    kgrid=kgrid/2.d0 !!!!!DIVIDI PER NUMERO SITI, RBZ
    if(allocated(t_k))deallocate(t_k)
    allocate(t_k(Nlat,Nlat,Nspin,Nspin,Norb,Norb,size(kgrid,1)))
    t_k=zero
    !
    do ik=1,size(kgrid,1)
        !
        kpoint=kgrid(ik,:)
        !
        t_k(1,1,1,1,1,1,ik)=-mu
        t_k(2,2,1,1,1,1,ik)=-mu
        t_k(3,3,1,1,1,1,ik)=-mu
        t_k(4,4,1,1,1,1,ik)=-mu
        !
        t_k(1,2,1,1,1,1,ik)= -t
        t_k(2,1,1,1,1,1,ik)= -t
        !
        t_k(2,3,1,1,1,1,ik)= -t
        t_k(3,2,1,1,1,1,ik)= -t
        !
        t_k(3,4,1,1,1,1,ik)= -t
        t_k(4,3,1,1,1,1,ik)= -t
        !
        t_k(4,1,1,1,1,1,ik)= -t
        t_k(1,4,1,1,1,1,ik)= -t
        !
        t_k(1,4,1,1,1,1,ik)=t_k(1,4,1,1,1,1,ik)-t*exp(xi*kpoint(2))
        t_k(4,1,1,1,1,1,ik)=t_k(4,1,1,1,1,1,ik)-t*exp(-xi*kpoint(2))
        t_k(1,2,1,1,1,1,ik)=t_k(1,2,1,1,1,1,ik)-t*exp(xi*kpoint(1))
        t_k(2,1,1,1,1,1,ik)=t_k(2,1,1,1,1,1,ik)-t*exp(-xi*kpoint(1))
        !
        t_k(2,3,1,1,1,1,ik)=t_k(2,3,1,1,1,1,ik)-t*exp(-xi*kpoint(2))
        t_k(3,2,1,1,1,1,ik)=t_k(3,2,1,1,1,1,ik)-t*exp(xi*kpoint(2))
        t_k(3,4,1,1,1,1,ik)=t_k(3,4,1,1,1,1,ik)-t*exp(xi*kpoint(1))
        t_k(4,3,1,1,1,1,ik)=t_k(4,3,1,1,1,1,ik)-t*exp(-xi*kpoint(1))
    enddo       
  end subroutine generate_t_k


end program vca_test



