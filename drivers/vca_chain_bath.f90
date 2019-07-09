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
  character(len=16)                               :: finput
  real(8)                                         :: ts
  integer                                         :: Nx,Ny,Lx,Ly,Rx,Ry,SAMPLING
  integer                                         :: unit
  integer                                         :: comm,rank
  logical                                         :: master,wloop,wmin
  integer                                         :: nloop
  real(8)                                         :: mu,t,t_var,mu_var
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

  Nb=vca_get_bath_dimension()
  allocate(Bath(Nb))
  call vca_init_solver(comm,bath)
 
  if(Norb/=1)stop "Norb != 1"
  Nlso=Nlat*Nspin*Norb
  !
  ! FIXME: 
  t_var=1.d0
  t=1.0d0
  mu=0.d0
  mu_var=0.d0

  if(.not.allocated(wm))allocate(wm(Lmats))
  if(.not.allocated(wr))allocate(wr(Lreal))
  wm     = pi/beta*real(2*arange(1,Lmats)-1,8)
  wr     = linspace(wini,wfin,Lreal)

  !INIT SOLVER:

  call vca_init_solver(comm,bath)
  
  !LOOP:

  if(wloop)then
     allocate(ts_array(Nloop))
     allocate(omega_array(Nloop))
     ts_array = linspace(0.05d0,0.5d0,Nloop)
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
  endif
  !
  !
  !
 
  call finalize_MPI()

!+------------------------------------------------------------------+
!FUNCTIONS
!+------------------------------------------------------------------+

contains



  function solve_vca1d(Vij) result(Omega)
    real(8)                      :: Vij,Eij
    real(8)                      :: Omega
    integer                      :: ix,iy,ik
    !
    print*,""
    print*,"----- INIT -----"
    !
    print*,"V_VAR = ",Vij
    Eij=0.d0
    call generate_tcluster()
    call generate_hk()
    !BATH VARIATIONAL SETUP
    do iy=1,Nspin
      do ik=1,Norb
        call set_bath_component(bath,1,iy,ik,e_component=[Eij])
        call set_bath_component(bath,1,iy,ik,v_component=[Vij])
        call set_bath_component(bath,Nlat,iy,ik,e_component=[Eij])
        call set_bath_component(bath,Nlat,iy,ik,v_component=[Vij])     
        do ix=2,Nlat-1   
          call set_bath_component(bath,ix,iy,ik,e_component=[Eij])
          call set_bath_component(bath,ix,iy,ik,v_component=[0.d0])
        enddo
      enddo
    enddo
    call vca_solve(comm,t_prime,h_k,bath)
    call vca_get_sft_potential(omega)
    !
    print*,"------ DONE ------"
    print*,""
    !
  end function solve_vca1d




  subroutine generate_tcluster()
    integer                          :: ii,ispin,iorb,i,j,jj
    character(len=64)                :: file_
    real(8),dimension(Nlso,Nlso)     :: H0
    integer                          :: unit
    file_ = "tcluster_matrix.dat"
    !
    if(allocated(t_prime))deallocate(t_prime)
    allocate(t_prime(Nlat,Nlat,Nspin,Nspin,Norb,Norb))
    t_prime=zero
    !
    do ispin=1,Nspin
    do ii=1,Nlat
          t_prime(ii,ii,ispin,ispin,1,1)= -mu_var
          !
          if(ii>1)t_prime(ii,ii-1,ispin,ispin,1,1)= -t_var
          if(ii<Nlat)t_prime(ii,ii+1,ispin,ispin,1,1)= -t_var
      enddo
    enddo
    !
    H0=vca_nnn2lso_reshape(t_prime,Nlat,Nspin,Norb)
    !
    open(free_unit(unit),file=trim(file_))
    do ii=1,Nlso
       write(unit,"(5000(F5.2,1x))")(H0(ii,jj),jj=1,Nlso)
    enddo
    close(unit)
  end subroutine generate_tcluster




  subroutine generate_hk()
    integer                             :: ik,ii,ispin,iorb,unit,jj
    real(8),dimension(Nkpts,1)          :: kgrid
    real(8),dimension(Nlso,Nlso)        :: H0
    character(len=64)                   :: file_
    file_ = "tlattice_matrix.dat"
    !
    call TB_build_kgrid([Nkpts],kgrid)
    kgrid=kgrid/Nlat !!!!!DIVIDI OGNI K PER NUMERO SITI, RBZ
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
    H0=vca_nnn2lso_reshape(tk([0.3d0,0.6d0]),Nlat,Nspin,Norb)
    !
    open(free_unit(unit),file=trim(file_))
    do ii=1,Nlat*Nspin*Norb
       write(unit,"(5000(F5.2,1x))")(H0(ii,jj),jj=1,Nlat*Nspin*Norb)
    enddo
    close(unit)    
  end subroutine generate_hk


  function tk(kpoint) result(hopping_matrix)
    integer                                                               :: ilat,ispin,iorb,unit,jj
    real(8),dimension(Ndim),intent(in)                                    :: kpoint
    complex(8),dimension(Nlat,Nlat,Nspin,Nspin,Norb,Norb)                 :: hopping_matrix
    !
    hopping_matrix=zero
    do ispin=1,Nspin
      do ilat=1,Nlat
          hopping_matrix(ilat,ilat,ispin,ispin,1,1)= -mu
          !
          if(ilat>1)hopping_matrix(ilat,ilat-1,ispin,ispin,1,1)= -t
          if(ilat<Nlat)hopping_matrix(ilat,ilat+1,ispin,ispin,1,1)= -t
      enddo
      !
      hopping_matrix(1,Nlat,ispin,ispin,1,1)=hopping_matrix(1,Nlat,ispin,ispin,1,1)-t*exp(xi*kpoint(1)*Nlat)
      hopping_matrix(Nlat,1,ispin,ispin,1,1)=hopping_matrix(Nlat,1,ispin,ispin,1,1)-t*exp(-xi*kpoint(1)*Nlat)
    enddo
    ! 
  end function tk


end program vca_chain1d



