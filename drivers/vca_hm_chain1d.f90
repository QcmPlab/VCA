program vca_chain1d
  USE SCIFOR
  USE DMFT_TOOLS
  USE MPI
  USE VCA
  !
  implicit none
  integer                                         :: Nlso,Nsys,Ndim
  integer,dimension(1)                            :: Nkpts
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
  complex(8),allocatable,dimension(:,:,:,:,:,:,:) :: h_k
  character(len=16)                               :: finput
  real(8)                                         :: ts
  integer                                         :: Nx,Ny,Lx,Ly,Rx,Ry,SAMPLING
  integer                                         :: unit
  integer                                         :: comm,rank
  logical                                         :: master,wloop,wmin
  integer                                         :: nloop
  real(8)                                         :: mu,t,t_var,mu_var,omegadummy
  real(8),dimension(:),allocatable                :: ts_array,omega_array
  integer,dimension(1)                            :: min_loc
  complex(8),allocatable,dimension(:,:,:,:,:)     :: gfmats_periodized        ![Nspin][Nspin][Norb][Norb][Lmats]
  complex(8),allocatable,dimension(:,:,:,:,:)     :: gfreal_periodized        ![Nspin][Nspin][Norb][Norb][Lreal]
  complex(8),allocatable,dimension(:,:,:,:,:)     :: Smats_periodized         ![Nspin][Nspin][Norb][Norb][Lmats]
  complex(8),allocatable,dimension(:,:,:,:,:)     :: Sreal_periodized         ![Nspin][Nspin][Norb][Norb][Lreal]


  call init_MPI()
  comm = MPI_COMM_WORLD
  call StartMsg_MPI(comm)
  rank = get_Rank_MPI(comm)
  master = get_Master_MPI(comm)

  call parse_cmd_variable(finput,"FINPUT",default='inputVCA.conf')
  call parse_input_variable(Nkpts,"Nkpts",finput,default=[10],comment="Number of k-points along each direction")
  call parse_cmd_variable(SAMPLING,"SAMPLING",default=100)
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
  !
  !
  ! FIXME:
  Ndim=size(Nkpts) 
  Nlso=Nlat*Nspin*Norb
  t_var=ts
  t=1.0d0
  mu=0.d0
  mu_var=0.d0

  if(.not.allocated(wm))allocate(wm(Lmats))
  if(.not.allocated(wr))allocate(wr(Lreal))
  wm     = pi/beta*real(2*arange(1,Lmats)-1,8)
  wr     = linspace(wini,wfin,Lreal)

  !INIT SOLVER:

  call vca_init_solver(comm)
    print_observables=.false.
  
  !LOOP:

  if(wloop)then
     allocate(ts_array(Nloop))
     allocate(omega_array(Nloop))
     ts_array = linspace(0.2d0,2d0,Nloop)
     do iloop=1,Nloop
        omega_array(iloop)=solve_vca1d(ts_array(iloop))
     enddo
     call splot("sft_Omega_loopVSts.dat",ts_array,omega_array)
  else
    call generate_tcluster()
    call generate_hk()
    call vca_solve(comm,t_prime,h_k)
  endif

  !MINIMIZATION:

  if(wmin)then
     print*,"Guess:",ts
     call  brent(solve_vca1d,ts,[0.5d0,2d0])
     print*,"Result ts : ",ts
     stop
  else
    omegadummy=solve_vca1d(ts)
  endif
  !
  if(allocated(wm))deallocate(wm)
  if(allocated(wr))deallocate(wr)
  
  call finalize_MPI()

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
    call generate_hk()
    call vca_solve(comm,t_prime,h_k)
    call vca_get_sft_potential(omega)
    !
    print*,"------ DONE ------"
    print*,""
    !
  end function solve_vca1d




  subroutine generate_tcluster()
    integer                          :: ii,ispin,iorb,i,j
    character(len=64)                :: file_
    complex(8),dimension(Nlso,Nlso)  :: H0
    integer                          :: unit
    file_ = "tcluster_matrix.dat"
    !
    if(allocated(t_prime))deallocate(t_prime)
    allocate(t_prime(Nlat,Nlat,Nspin,Nspin,Norb,Norb))
    t_prime=zero
    !
    do ii=1,Nlat
      do ispin=1,Nspin
        do iorb=1,Norb
          t_prime(ii,ii,1,1,iorb,iorb)= -mu_var
          !
          if(ii>1)t_prime(ii,ii-1,ispin,ispin,iorb,iorb)= -t_var
          if(ii<Nlat)t_prime(ii,ii+1,ispin,ispin,iorb,iorb)= -t_var
        enddo
      enddo
    enddo
    !
    !
  end subroutine generate_tcluster




  subroutine generate_hk()
    integer                                      :: ik,ii,ispin,iorb,unit,jj
    real(8),dimension(product(Nkpts),1)          :: kgrid
    real(8),dimension(1)                         :: e1,bk1
    real(8)                                      :: bklen
    !
    e1 = [1d0]
    call TB_set_ei(eix=e1)
    bklen=2d0*pi
    bk1=bklen*[1d0]
    call TB_set_bk(bkx=bk1)
    call TB_build_kgrid(Nkpts,kgrid)
    kgrid=kgrid/Nlat !!!!!DIVIDI OGNI K PER NUMERO SITI, RBZ
    !
    !
    if(allocated(h_k))deallocate(h_k)
    allocate(h_k(Nlat,Nlat,Nspin,Nspin,Norb,Norb,size(kgrid,1)))  
    h_k=zero
    !
    do ik=1,size(kgrid,1)
        !
        h_k(:,:,:,:,:,:,ik)=tk(kgrid(ik,1))
        !
    enddo
    !    
  end subroutine generate_hk


  function tk(kpoint) result(hopping_matrix)
    integer                                                               :: ilat,ispin,iorb,unit,jj
    real(8),dimension(Ndim),intent(in)                                    :: kpoint
    complex(8),dimension(Nlat,Nlat,Nspin,Nspin,Norb,Norb)                 :: hopping_matrix
    !
    do ispin=1,Nspin
      do iorb=1,Norb
        do ilat=1,Nlat
          hopping_matrix(ilat,ilat,1,1,iorb,iorb)= -mu
          !
          if(ilat>1)hopping_matrix(ilat,ilat-1,ispin,ispin,iorb,iorb)= -t
          if(ilat<Nlat)hopping_matrix(ilat,ilat+1,ispin,ispin,iorb,iorb)= -t
        enddo
      enddo
    enddo
      !
    do ispin=1,Nspin
      do iorb=1,Norb  
        hopping_matrix(1,Nlat,ispin,ispin,iorb,iorb)=hopping_matrix(1,Nlat,ispin,ispin,iorb,iorb)-t*exp(xi*kpoint(1)*Nlat)
        hopping_matrix(Nlat,1,ispin,ispin,iorb,iorb)=hopping_matrix(Nlat,1,ispin,ispin,iorb,iorb)-t*exp(-xi*kpoint(1)*Nlat)
      enddo
    enddo
    ! 
  end function tk


end program vca_chain1d



