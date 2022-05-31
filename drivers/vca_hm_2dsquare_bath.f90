program vca_square_bath
  USE SCIFOR
  USE DMFT_TOOLS
  USE MPI
  USE VCA
  !
  implicit none
  integer                                         :: Nlso,Nsys,Ndim
  integer,dimension(2)                            :: Nkpts
  integer                                         :: ilat,jlat
  integer                                         :: iloop,jloop
  integer                                         :: ix,iy,ik
  logical                                         :: converged
  real(8)                                         :: wband
  !Bath
  complex(8),dimension(:,:,:,:,:,:),allocatable   :: bath_h,bath_v
  !The local hybridization function:
  real(8),dimension(:,:),allocatable              :: Tsys,Tref,Vmat,Htb,Mmat,dens
  real(8),allocatable,dimension(:)                :: wm,wr
  complex(8)                                      :: iw
  complex(8),allocatable,dimension(:,:,:,:,:,:,:) :: Gmats,Greal
  complex(8),allocatable,dimension(:,:,:,:,:,:)   :: t_prime
  complex(8),allocatable,dimension(:,:,:,:,:,:,:) :: h_k
  complex(8),allocatable,dimension(:,:,:,:,:,:,:) :: t_k
  character(len=16)                               :: finput
  real(8)                                         :: ts
  integer                                         :: Nx,Ny,Lx,Ly,Rx,Ry
  integer                                         :: unit
  integer                                         :: comm,rank
  logical                                         :: master,wloop,wmin
  integer                                         :: nloop
  character(len=6)                                :: scheme
   !FIX!!!!!!!!!!!!!!!
  real(8)                                         :: mu,t,t_var,mu_var,omegadummy
  real(8),dimension(:),allocatable                :: ts_array,omega_array
  real(8),dimension(:),allocatable                :: ts_array_x,ts_array_y,params
  real(8),dimension(:,:),allocatable              :: omega_grid
  integer,dimension(1)                            :: min_loc
  complex(8),allocatable,dimension(:,:,:,:,:)     :: gfmats_periodized ![Nspin][Nspin][Norb][Norb][Lmats]
  complex(8),allocatable,dimension(:,:,:,:,:)     :: gfreal_periodized ![Nspin][Nspin][Norb][Norb][Lreal]
  complex(8),allocatable,dimension(:,:,:,:,:)     :: Smats_periodized         ![Nspin][Nspin][Norb][Norb][Lmats]
  complex(8),allocatable,dimension(:,:,:,:,:)     :: Sreal_periodized         ![Nspin][Nspin][Norb][Norb][Lreal]
  complex(8),allocatable,dimension(:,:,:,:,:)     :: gtest_mats,gtest_Real,sigmatest_mats,sigmatest_real
  real(8),allocatable,dimension(:,:)              :: kgrid_test,kpath_test

  call init_MPI()
  comm = MPI_COMM_WORLD
  call StartMsg_MPI(comm)
  rank = get_Rank_MPI(comm)
  master = get_Master_MPI(comm)

  call parse_cmd_variable(finput,"FINPUT",default='inputVCA.conf')
  call parse_input_variable(ts,"ts",finput,default=1d0)
  call parse_input_variable(Nx,"Nx",finput,default=2,comment="Number of sites along X")
  call parse_input_variable(Ny,"Ny",finput,default=2,comment="Number of sites along Y")
  call parse_input_variable(Nkpts,"Nkpts",finput,default=[10,10],comment="Number of k-points along each direction")
  call parse_input_variable(nloop,"NLOOP",finput,default=100)
  call parse_input_variable(wloop,"WLOOP",finput,default=.false.)
  call parse_input_variable(scheme,"SCHEME",finput,default="g")
  call parse_input_variable(wmin,"wmin",finput,default=.false.,comment="T: includes global minimization")
  !
  call vca_read_input(trim(finput),comm)
  !
  !
  Ndim=size(Nkpts)
  Ny=Nx
  Nlat=Nx*Ny
  Nlso = Nlat*Norb*Nspin
  Nlat_bath=Nlat
  Norb_bath=Norb
  call naming_convention()
  !
  t_var=1.0d0
  t=1.0d0
  mu=0.d0
  mu_var=0.d0
  bandwidth=2.d0*Ndim*(2*t) !(2 times sum over dim of 2*t*max (cosine))

  !Add DMFT CTRL Variables:
  call add_ctrl_var(Nlat,"NLAT")
  call add_ctrl_var(Norb,"norb")
  call add_ctrl_var(Nlat_bath,"NLAT_bath")
  call add_ctrl_var(Norb_bath,"norb_bath")
  call add_ctrl_var(Nspin,"nspin")
  call add_ctrl_var(beta,"beta")
  call add_ctrl_var(xmu,"xmu")
  call add_ctrl_var(wini,'wini')
  call add_ctrl_var(wfin,'wfin')
  call add_ctrl_var(eps,"eps")

  if(allocated(bath_h))deallocate(bath_h)
  if(allocated(bath_v))deallocate(bath_v)
  allocate(bath_h(Nlat_bath,Nlat_bath,Nspin,Nspin,Norb_bath,Norb_bath))
  allocate(bath_v(Nlat     ,Nlat_bath,Nspin,Nspin,Norb     ,Norb_bath))

  
  call vca_init_solver(comm,bath_h,bath_v)
    

  if(wloop)then
    allocate(ts_array(Nloop))
    allocate(omega_array(Nloop))
    !
    ts_array = linspace(0.05d0,0.7d0,Nloop)
    do iloop=1,Nloop
      omega_array(iloop)=solve_vca_square(ts_array(iloop))
      open(free_unit(unit),file="TEST.vca",position='append')
      write(unit,*)ts_array(iloop),omega_array(iloop)
      close(unit)
    enddo
    !
    call splot("sft_Omega_loopVSts.dat",ts_array,omega_array)
    !
  endif
  

  omegadummy=solve_vca_square(0.7d0)
  call finalize_MPI()

contains

  !+------------------------------------------------------------------+
  !PURPOSE  : solve the model
  !+------------------------------------------------------------------+



  function solve_vca_square(V) result(Omega)
    real(8)                      :: V,E
    real(8)                      :: Omega
    !
    !
    t_var=1.0d0
    E=0.d0
    mu_var=0.d0
    print*,""
    print*,"------ Doing for ", V," ------"
    call construct_bath(E,V)
    call generate_tcluster()
    call generate_hk()
    call vca_solve(comm,t_prime,h_k,bath_h,bath_v)
    call vca_get_sft_potential(omega)
    print*,""
    print*,"------ DONE ------"
    print*,""
    !
  end function solve_vca_square


  !+------------------------------------------------------------------+
  !PURPOSE  : generate test hopping matrices
  !+------------------------------------------------------------------+

  subroutine generate_tcluster()
    integer                                      :: ilat,jlat,ispin,iorb,ind1,ind2
    !
    if(allocated(t_prime))deallocate(t_prime)
    allocate(t_prime(Nlat,Nlat,Nspin,Nspin,Norb,Norb))
    t_prime=zero
    !
    do ispin=1,Nspin
      do ilat=1,Nx
        do jlat=1,Ny
          ind1=indices2N([ilat,jlat])
          t_prime(ind1,ind1,ispin,ispin,1,1)= -mu_var
          if(ilat<Nx)then
            ind2=indices2N([ilat+1,jlat])
            t_prime(ind2,ind1,ispin,ispin,1,1)= -t_var
          endif
          if(ilat>1)then
            ind2=indices2N([ilat-1,jlat])
            t_prime(ind2,ind1,ispin,ispin,1,1)= -t_var
          endif
          if(jlat<Ny)then
            ind2=indices2N([ilat,jlat+1])
            t_prime(ind2,ind1,ispin,ispin,1,1)= -t_var
          endif
          if(jlat>1)then
            ind2=indices2N([ilat,jlat-1])
            t_prime(ind2,ind1,ispin,ispin,1,1)= -t_var
          endif
        enddo
      enddo
    enddo
    !
  end subroutine generate_tcluster




  subroutine generate_hk()
    integer                                      :: ik,ii,ispin,iorb,unit,jj
    real(8),dimension(product(Nkpts),Ndim)       :: kgrid
    real(8),dimension(Nlso,Nlso)                 :: H0
    real(8),dimension(2)                        :: e1,e2,bk1,bk2
    real(8)                                     :: bklen
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
  end subroutine generate_hk


  function tk(kpoint) result(hopping_matrix)
    integer                                                                 :: ilat,jlat,ispin,iorb,i,j,ind1,ind2
    real(8),dimension(Ndim),intent(in)                                      :: kpoint
    complex(8),dimension(Nlat,Nlat,Nspin,Nspin,Norb,Norb)                   :: hopping_matrix
    !
    hopping_matrix=zero
    !
    do ilat=1,Nx
      do jlat=1,Ny
        do ispin=1,Nspin
          ind1=indices2N([ilat,jlat])
          hopping_matrix(ind1,ind1,ispin,ispin,1,1)= -mu
          if(ilat<Nx)then
            ind2=indices2N([ilat+1,jlat])
            hopping_matrix(ind2,ind1,ispin,ispin,1,1)= -t
          endif
          if(ilat>1)then
            ind2=indices2N([ilat-1,jlat])
            hopping_matrix(ind2,ind1,ispin,ispin,1,1)= -t
          endif
          if(jlat<Ny)then
            ind2=indices2N([ilat,jlat+1])
            hopping_matrix(ind2,ind1,ispin,ispin,1,1)= -t
          endif
          if(jlat>1)then
            ind2=indices2N([ilat,jlat-1])
            hopping_matrix(ind2,ind1,ispin,ispin,1,1)= -t
          endif
        enddo
      enddo
    enddo
    !
    do ispin=1,Nspin
      do ilat=1,Ny
        ind1=indices2N([1,ilat])
        ind2=indices2N([Nx,ilat])
        hopping_matrix(ind2,ind1,ispin,ispin,1,1)=hopping_matrix(ind2,ind1,ispin,ispin,1,1) -t*exp(xi*kpoint(1)*Nx)
        hopping_matrix(ind1,ind2,ispin,ispin,1,1)=hopping_matrix(ind1,ind2,ispin,ispin,1,1) -t*exp(-xi*kpoint(1)*Nx)
      enddo
      do ilat=1,Nx
        ind1=indices2N([ilat,1])
        ind2=indices2N([ilat,Ny])
        hopping_matrix(ind2,ind1,ispin,ispin,1,1)=hopping_matrix(ind2,ind1,ispin,ispin,1,1) -t*exp(xi*kpoint(2)*Ny)
        hopping_matrix(ind1,ind2,ispin,ispin,1,1)=hopping_matrix(ind1,ind2,ispin,ispin,1,1) -t*exp(-xi*kpoint(2)*Ny)
      enddo
    enddo
    ! 
  end function tk
  
  
  subroutine construct_bath(eps,v)
    real(8)                 :: eps,v
    integer                 :: i,ib,o,ob
    !
    bath_h=zero
    bath_v=zero
    do ib=1,Nlat_bath
      do ob=1,Norb_bath
        bath_h(ib,ib,1,1,ob,ob)         = xmu-eps
        bath_h(ib,ib,Nspin,Nspin,ob,ob) = xmu-eps
        bath_v(ib,ib,1,1,ob,ob)         = v
        bath_v(ib,ib,Nspin,Nspin,ob,ob) = v
      enddo
    enddo
  end subroutine construct_bath



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

end program vca_square_bath



