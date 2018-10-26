program vca_chain1d
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
  logical                                         :: master,wloop,wmin
  integer                                         :: nloop
  real(8)                                         :: mu,t,t_var,mu_var
  real(8),dimension(:),allocatable                :: ts_array,omega_array
  integer,dimension(1)                            :: min_loc

  call init_MPI()
  comm = MPI_COMM_WORLD
  call StartMsg_MPI(comm)
  rank = get_Rank_MPI(comm)
  master = get_Master_MPI(comm)

  call parse_cmd_variable(finput,"FINPUT",default='inputVCA.conf')
  call parse_input_variable(ts,"ts",finput,default=1d0)
  call parse_input_variable(wloop,"wloop",finput,default=.false.,comment="T: includes loop over ts")
  call parse_input_variable(wmin,"wmin",finput,default=.false.,comment="T: includes global minimization")
  call parse_input_variable(nloop,"NLOOP",finput,default=5)
  !
  call vca_read_input(trim(finput),comm)
  
  

  !Add DMFTtools CTRL Variables:

  call add_ctrl_var(Nlat,"NLAT")
  call add_ctrl_var(Norb,"norb")
  call add_ctrl_var(Nspin,"nspin")
  call add_ctrl_var(beta,"beta")
  call add_ctrl_var(xmu,"xmu")
  call add_ctrl_var(wini,'wini')
  call add_ctrl_var(wfin,'wfin')
  call add_ctrl_var(eps,"eps")

 
  if(Norb/=1)stop "Norb != 1"
  if(Nspin/=1)stop "Nspin != 1"
  !
  ! FIXME: 
  t_var=1.d0
  t=1.0d0
  mu=0.d0
  mu_var=0.d0


  !INIT SOLVER:

  call vca_init_solver(comm)
  
  !LOOP:

  if(wloop)then
     allocate(ts_array(Nloop))
     allocate(omega_array(Nloop))
     ts_array = linspace(0.5d0,2.d0,Nloop)
     do iloop=1,Nloop
        omega_array(iloop)=solve_vca1d(ts_array(iloop))
     enddo
     call splot("sft_Omega_loopVSts.dat",ts_array,omega_array)
  else
    call generate_tcluster()
    call generate_tk()
    call vca_solve(comm,t_prime,t_k)
  endif

  !MINIMIZATION:

  if(wmin)then
     print*,"Guess:",ts
     call  brent(solve_vca1d,ts,[0.5d0,2d0])
     print*,"Result ts : ",ts
     stop
  endif



!+------------------------------------------------------------------+
!FUNCTIONS
!+------------------------------------------------------------------+

contains



  function solve_vca1d(tij) result(Omega)
    real(8)                      :: tij
    real(8)                      :: Omega
    !
    t_var=tij
    print*,""
    print*,"----- INIT -----"
    !
    print*,"T_VAR = ",t_var
    call generate_tcluster()
    call generate_tk()
    call vca_solve(comm,t_prime,t_k)
    call vca_get_sft_potential(omega)
    !
    print*,"------ DONE ------"
    print*,""
    !
  end function solve_vca1d




  subroutine generate_tcluster()
    integer                          :: ii,ispin,iorb,i,j
    character(len=64)                :: file_
    real(8),dimension(Nlat,Nlat)     :: H0
    integer                          :: unit
    file_ = "tcluster_matrix.dat"
    !
    if(allocated(t_prime))deallocate(t_prime)
    allocate(t_prime(Nlat,Nlat,Nspin,Nspin,Norb,Norb))
    t_prime=zero
    !
    do ii=1,Nlat
        t_prime(ii,ii,1,1,1,1)= -mu_var
        !
        if(ii>1)t_prime(ii,ii-1,1,1,1,1)= -t_var
        if(ii<Nlat)t_prime(ii,ii+1,1,1,1,1)= -t_var
    enddo
    !
    H0=vca_nnn2lso_reshape(t_prime,Nlat,Nspin,Norb)
    !
    open(free_unit(unit),file=trim(file_))
    do i=1,Nlat
       write(unit,"(5000(F5.2,1x))")(H0(i,j),j=1,Nlat)
    enddo
    close(unit)
  end subroutine generate_tcluster




  subroutine generate_tk()
    integer                             :: ik,ii,ispin,iorb,unit,jj
    real(8),dimension(Ndim)             :: kpoint
    character(len=64)                   :: file_
    real(8),dimension(Nkpts,1)          :: kgrid
    complex(8),dimension(Nlat,Nlat)     :: H0
    !
    file_ = "tk_matrix_at_some_k.dat"
    !
    call TB_build_kgrid([Nkpts],kgrid)
    kgrid=kgrid/Nlat !!!!!DIVIDI OGNI K PER NUMERO SITI, RBZ
    !
    if(allocated(t_k))deallocate(t_k)
    allocate(t_k(Nlat,Nlat,Nspin,Nspin,Norb,Norb,size(kgrid,1)))  
    t_k=zero
    !
    do ik=1,size(kgrid,1)
        !
        kpoint=kgrid(ik,1)
        !
        do ii=1,Nlat
            t_k(ii,ii,1,1,1,1,ik)= -mu
            !
            if(ii>1)t_k(ii,ii-1,1,1,1,1,ik)= -t
            if(ii<Nlat)t_k(ii,ii+1,1,1,1,1,ik)= -t
        enddo
        !
        t_k(1,Nlat,1,1,1,1,ik)=t_k(1,Nlat,1,1,1,1,ik)-t*exp(xi*kpoint(1)*Nlat)
        t_k(Nlat,1,1,1,1,1,ik)=t_k(Nlat,1,1,1,1,1,ik)-t*exp(-xi*kpoint(1)*Nlat)
        !
    enddo
    H0=vca_nnn2lso_reshape(t_k(:,:,:,:,:,:,10),Nlat,Nspin,Norb)
    !
    open(free_unit(unit),file=trim(file_))
    do ii=1,Nlat
       write(unit,"(5000(F5.2,1x))")(H0(ii,jj),jj=1,Nlat)
    enddo
    close(unit)       
  end subroutine generate_tk


 

end program vca_chain1d



