program vca_bhz_2d
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
  logical                                         :: converged,usez
  real(8)                                         :: wband
  !Bath
  real(8),allocatable                             :: Bath(:)
  integer                                         :: Nb,COUNTER
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
  real(8)                                         :: M,M_var,t,t_var,lambda,lambda_var,mu,mu_var,Mh
  !
  real(8),dimension(:),allocatable                :: ts_array,omega_array
  integer,dimension(1)                            :: min_loc
  complex(8),allocatable,dimension(:,:,:,:,:)     :: gfmats_local,gfmats_periodized ![Nspin][Nspin][Norb][Norb][Lmats]
  complex(8),allocatable,dimension(:,:,:,:,:)     :: Smats_local,Smats_periodized         ![Nspin][Nspin][Norb][Norb][Lmats]
  complex(8),allocatable,dimension(:,:,:)         :: Smats_periodized_lso     
  real(8),allocatable,dimension(:,:)              :: kgrid_test,kpath_test
  integer                                         :: iii,jjj,kkk
  !
  call init_MPI()
  comm = MPI_COMM_WORLD
  call StartMsg_MPI(comm)
  rank = get_Rank_MPI(comm)
  master = get_Master_MPI(comm)
  !
  call parse_cmd_variable(finput,"FINPUT",default='inputVCA.conf')
  call parse_cmd_variable(SAMPLING,"SAMPLING",default=100)
  call parse_input_variable(ts,"ts",finput,default=1d0)
  call parse_input_variable(Mh,"Mh",finput,default=3d0)
  call parse_input_variable(Nx,"Nx",finput,default=2,comment="Number of sites along X")
  call parse_input_variable(Ny,"Ny",finput,default=2,comment="Number of sites along Y")
  call parse_input_variable(nloop,"NLOOP",finput,default=100)
  call parse_input_variable(hopping,"HOPPING",finput,default=1.d0)
  call parse_input_variable(wloop,"WLOOP",finput,default=.false.)
  call parse_input_variable(scheme,"SCHEME",finput,default="g")
  call parse_input_variable(wmin,"wmin",finput,default=.false.,comment="T: includes global minimization")
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
  Nlat=Nx**Ndim
  Ny=Nx
  Nlso = Nlat*Norb*Nspin
  COUNTER=0
  !
  !
  if(.not.allocated(wm))allocate(wm(Lmats))
  if(.not.allocated(wr))allocate(wr(Lreal))
  wm     = pi/beta*real(2*arange(1,Lmats)-1,8)
  wr     = linspace(wini,wfin,Lreal)
  !
  call vca_init_solver(comm)
  !
  if(wloop)then
    allocate(ts_array(Nloop))
    allocate(omega_array(Nloop))
    ts_array = linspace(0.4d0,1.5d0,Nloop)
    !
    do iloop=1,Nloop
       omega_array(iloop)=solve_vca_square(ts_array(iloop))
    enddo
    !
    call splot("sft_Omega_loopVSts.dat",ts_array,omega_array)
    min_loc = minloc(omega_array)
    write(800,*)min_loc,ts_array(min_loc(1)),omega_array(min_loc(1))
  else
    omegadummy=solve_vca_square(ts)
    if(scheme=="g")then
      call get_local_gf()
      !call solve_Htop()
      call solve_Htop_new()
      call build_z2_indices() 
    else
      call get_local_sigma()
      !call solve_Htop()
      call solve_Htop_new()
      call build_z2_indices() 
    endif
  endif
  !
  !MINIMIZATION:
  !
  if(allocated(wm))deallocate(wm)
  if(allocated(wr))deallocate(wr)
  !
  call finalize_MPI()
  !
contains

  !+------------------------------------------------------------------+
  !PURPOSE  : solve the model
  !+------------------------------------------------------------------+


  function solve_vca_square(tij) result(Omega)
    real(8)                      :: tij
    real(8)                      :: Vij,Eij
    real(8)                      :: Omega
    !
    t=0.5d0
    t_var=tij
    !
    mu=0.d0*t
    M=(2.d0*t)*Mh
    lambda=(2.d0*t)*0.3d0
    !  
    mu_var=0.d0*t
    M_var=(2.d0*t)*Mh
    lambda_var=(2.d0*t)*0.3d0
    !
    print*,""
    print*,"------ Doing for ",tij," ------"
    call generate_tcluster()
    call generate_hk()
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
    !
    do ispin=1,Nspin
      do ilat=1,Nx
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
    H0=vca_nnn2lso_reshape(tk([0.d0,0.d0]),Nlat,Nspin,Norb)
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
       kpath(1,:)=[0.d0,0.d0,0.0d0]
       kpath(2,:)=[pi/Nx,pi/ny,0.0d0]
       kpath(3,:)=[pi/Nx,0.d0,0.0d0]
       kpath(4,:)=[0.d0,0.d0,0.0d0]
       file="Eigenbands.nint"
    endif
   allocate(colors(Nlso))
   colors=[red1,blue1,red1,blue1]
   
   if(master) call TB_Solve_model(tk_wrapper,Nlso,kpath,Nkpath,&
         colors_name=colors,&
         points_name=[character(len=20) :: 'G', 'M', 'X', 'G'],&
         file=reg(file))
  end subroutine solve_hk_GXMG




  function hk_bhz(kpoint,N) result(hopping_matrix_lso)
    integer                                         :: N
    real(8),dimension(:)                            :: kpoint
    real(8),dimension(Ndim)                         :: kpoint_
    complex(8),dimension(N,N)                       :: hopping_matrix_lso
    real(8)                                         :: energy_scale
    if(N.ne.Nspin*Norb)stop "error N: wrong dimension"
    !
    energy_scale=2.d0*t_var
    hopping_matrix_lso=zero
    hopping_matrix_lso= (M_var-energy_scale*(cos(kpoint(1))+cos(kpoint(2))))*kron_pauli( pauli_sigma_0, pauli_tau_z)+&
                                                                 lambda_var*sin(kpoint(1))*kron_pauli( pauli_sigma_z, pauli_tau_x)+&
                                                                 lambda_var*sin(kpoint(2))*kron_pauli( pauli_sigma_0, pauli_tau_y)
    print*,COUNTER
    COUNTER=COUNTER+1
    if(scheme=="g")then
      call build_sigma_g_scheme(kpoint)
    else
      call periodize_sigma_scheme(kpoint)
    endif
    !hopping_matrix_lso=vca_nn2so_reshape(gfmats_periodized(:,:,:,:,1),Nspin,Norb)
    !call inv(hopping_matrix_lso)
    hopping_matrix_lso=hopping_matrix_lso+DREAL(smats_periodized_lso(:,:,1))
    hopping_matrix_lso=matmul(Zrenorm(smats_periodized_lso(:,:,1)),DREAL(hopping_matrix_lso))
    !
  end function hk_bhz

  function hk_bhz_local(kpoint,N) result(hopping_matrix_lso)
    integer                                         :: N
    real(8),dimension(:)                            :: kpoint
    complex(8),dimension(N,N)                       :: hopping_matrix_lso
    real(8)                                         :: energy_scale
    if(N.ne.Nspin*Norb)stop "error N: wrong dimension"
    !
    energy_scale=2.d0*t_var
    hopping_matrix_lso=zero
    hopping_matrix_lso= (M_var-energy_scale*(cos(kpoint(1))+cos(kpoint(2))))*kron_pauli( pauli_sigma_0, pauli_tau_z)+&
                                                                 lambda_var*sin(kpoint(1))*kron_pauli( pauli_sigma_z, pauli_tau_x)+&
                                                                 lambda_var*sin(kpoint(2))*kron_pauli( pauli_sigma_0, pauli_tau_y)
    print*,COUNTER
    COUNTER=COUNTER+1
    hopping_matrix_lso=hopping_matrix_lso+DREAL(vca_nn2so_reshape(smats_local(:,:,:,:,1),Nspin,Norb))
    hopping_matrix_lso=matmul(Zrenorm(vca_nn2so_reshape(smats_local(:,:,:,:,1),Nspin,Norb)),hopping_matrix_lso)
    !
  end function hk_bhz_local


  function hk_bhz_clusterbase(kpoint,N) result(hopping_matrix_lso)
    integer                                                       :: N,ilat,jlat
    real(8),dimension(:)                                          :: kpoint
    real(8),dimension(Ndim)                                       :: kpoint_,ind1,ind2
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


  subroutine solve_Htop(kpath_)
    integer                                  :: i,j
    integer                                  :: Npts,Nkpath
    type(rgb_color),dimension(:),allocatable :: colors
    real(8),dimension(:,:),optional          :: kpath_
    real(8),dimension(:,:),allocatable       :: kpath
    character(len=64)                        :: file
    
    !This routine build the H(k) along the GXMG path in BZ,
    !Hk(k) is constructed along this path.
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
       kpath(2,:)=[pi,pi]
       kpath(3,:)=[pi,0.d0]
       kpath(4,:)=[0.d0,0.d0]
    endif
   allocate(colors(Nspin*Norb))
   colors=[red1,blue1,red1,blue1]
   !
   file="Eig_Htop_kSigma.nint"
   if(master) call TB_Solve_model(hk_bhz,Nspin*Norb,kpath,Nkpath,&   
         colors_name=colors,&
         points_name=[character(len=20) :: 'G', 'M', 'X', 'G'],&
         file=reg(file))
   file="Eig_Htop_localSigma.nint"
   if(master) call TB_Solve_model(hk_bhz_local,Nspin*Norb,kpath,Nkpath,&   
         colors_name=colors,&
         points_name=[character(len=20) :: 'G', 'M', 'X', 'G'],&
         file=reg(file))
  end subroutine solve_Htop


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
    colors(5) = red1
    colors(9) = red1
    colors(13) = red1
    colors(2) = blue1
    colors(6) = blue1
    colors(10) = blue1
    colors(14) = blue1
    colors(3) = red1
    colors(7) = red1
    colors(11) = red1
    colors(15) = red1
    colors(4) = blue1
    colors(8) = blue1
    colors(12) = blue1
    colors(16) = blue1
   !
   file="Eig_Htop_clusterbase.nint"
   if(master) call TB_Solve_model(hk_bhz_clusterbase,Nlat*Nspin*Norb,kpath,Nkpath,&   
         colors_name=colors,&
         points_name=[character(len=20) :: 'G', 'M', 'X', 'G'],&
         file=reg(file))
  end subroutine solve_Htop_new


  !+------------------------------------------------------------------+
  !Periodization functions G-SCHEME
  !+------------------------------------------------------------------+
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


  subroutine periodize_g_scheme(kpoint)
    integer                                                     :: ilat,jlat,ispin,iorb,jorb,ii
    real(8),dimension(Ndim)                                     :: kpoint,ind1,ind2
    complex(8),allocatable,dimension(:,:,:,:,:,:,:)             :: gfprime ![Nlat][Nlat][Nspin][Nspin][Norb][Norb][Lmats]
    complex(8),allocatable,dimension(:,:)                       :: gfTempMat ![Nlso][Nlso]
    complex(8),allocatable,dimension(:,:)                       :: Vk_lso ![Nlso][Nlso]
    complex(8),allocatable,dimension(:,:,:,:,:,:,:)             :: gfmats_unperiodized ![Nlat][Nlat][Nspin][Nspin][Norb][Norb][Lmats]
    !
    !
    allocate(gfmats_unperiodized(Nlat,Nlat,Nspin,Nspin,Norb,Norb,Lmats))
    allocate(gfprime(Nlat,Nlat,Nspin,Nspin,Norb,Norb,Lmats))
    allocate(gfTempMat(Nlat*Nspin*Norb,Nlat*Nspin*Norb))
    allocate(Vk_lso(Nlat*Nspin*Norb,Nlat*Nspin*Norb))
    if(.not.allocated(gfmats_periodized))allocate(gfmats_periodized(Nspin,Nspin,Norb,Norb,Lmats))
    !
    gfmats_unperiodized=zero
    gfprime=zero
    gfTempMat=zero
    Vk_lso=zero
    gfmats_periodized=zero
    call vca_get_gimp_matsubara(gfprime)
    !
    do ii=1,Lmats         
      Vk_lso=vca_nnn2lso_reshape(tk(kpoint)-t_prime-set_delta(xi*wm(ii),[ts],[Uloc(1)/2]),Nlat,Nspin,Norb)    !!!CORREGGI I PARAMETRI DEL BAGNO
      gfTempMat=vca_nnn2lso_reshape(gfprime(:,:,:,:,:,:,ii),Nlat,Nspin,Norb)
      call inv(gfTempMat)
      gfTempMat=gfTempMat-Vk_lso
      call inv(gfTempMat)
      gfmats_unperiodized(:,:,:,:,:,:,ii)=vca_lso2nnn_reshape(gfTempMat,Nlat,Nspin,Norb)
    enddo
    !
    do ilat=1,Nlat
      ind1=N2indices(ilat)   
      do jlat=1,Nlat
        ind2=N2indices(jlat)
          gfmats_periodized=gfmats_periodized+exp(-xi*dot_product(kpoint,ind1-ind2))*gfmats_unperiodized(ilat,jlat,:,:,:,:,:)/Nlat
      enddo
    enddo
    !   
  end subroutine periodize_g_scheme



  subroutine build_sigma_g_scheme(kpoint)
    integer                                                       :: i,ispin,iorb,ii
    real(8)                                                       :: energy_Scale
    real(8),dimension(Ndim)                                       :: kpoint
    complex(8),dimension(Nspin*Norb,Nspin*Norb,Lmats)             :: invG0mats_lso,invGmats_lso ![Nlso][Nlso]
    !
    !
    if(.not.allocated(Smats_periodized_lso))allocate(Smats_periodized_lso(Nspin*Norb,Nspin*Norb,Lmats))
    invG0mats_lso=zero
    invGmats_lso=zero
    Smats_periodized_lso  = zero
    energy_scale=2.d0*t_var
    !
    !Get G0^-1
    do ii=1,Lmats
        invG0mats_lso(:,:,ii) = (xi*wm(ii)+xmu)*eye(Nspin*Norb)  + (M_var-energy_scale*(cos(kpoint(1))+cos(kpoint(2))))*kron_pauli( pauli_sigma_0, pauli_tau_z)+&
                                                                 lambda_var*sin(kpoint(1))*kron_pauli( pauli_sigma_z, pauli_tau_x)+&
                                                                 lambda_var*sin(kpoint(2))*kron_pauli( pauli_sigma_0, pauli_tau_y)        ! FIXME: ad-hoc solution
    enddo
    !  
    !Get Gimp^-1
    call periodize_g_scheme(kpoint)
    do ii=1,Lmats  !FIXMEEEEEEEEEEEEEEEEEEEEEEEEEEEEE
      invGmats_lso(:,:,ii)=vca_nn2so_reshape(gfmats_periodized(:,:,:,:,ii),Nspin,Norb)
      call inv(invGmats_lso(:,:,ii))
    enddo
    !Get Sigma functions: Sigma= G0^-1 - G^-1
    Smats_periodized_lso = invG0mats_lso - invGmats_lso
    !
  end subroutine build_sigma_g_scheme


  subroutine get_local_gf()
    integer                                         :: ik,ispin,iorb,jorb
    character(len=30)                               :: suffix
    !
    if(master)then
    !
      if(.not.allocated(Gfmats_local))allocate(Gfmats_local(Nspin,Nspin,Norb,Norb,Lmats))
      if(.not.allocated(Smats_local))allocate(Smats_local(Nspin,Nspin,Norb,Norb,Lmats))
      gfmats_local=zero
      smats_local=zero
      !
      allocate(kgrid_test(Nkpts**ndim,Ndim)) 
      call TB_build_kgrid([Nkpts,Nkpts],kgrid_test)
      !
      print*,"Calculating Gloc and Sigma ",scheme," scheme"  
      !
      call start_timer
      !
      do ik=1,Nkpts**ndim
        call build_sigma_g_scheme(kgrid_test(ik,:))  !also periodizes g
        do ix=1,Lmats
          gfmats_local(:,:,:,:,ix)=gfmats_local(:,:,:,:,ix)+gfmats_periodized(:,:,:,:,ix)/(Nkpts**Ndim)
          smats_local(:,:,:,:,ix)=smats_local(:,:,:,:,ix)+vca_so2nn_reshape(Smats_periodized_lso(:,:,ix),Nspin,Norb)/(Nkpts**Ndim)
        enddo
        call eta(ik,Nkpts**ndim)
      enddo
      !
      call stop_timer
      !
      do iorb=1,Norb
       do jorb=1,Norb
         do ispin=1,Nspin
           suffix="_l"//str(iorb)//str(jorb)//"_s"//str(ispin)//"_"//str(scheme)//"_scheme"
           call splot("perG"//reg(suffix)//"_iw.vca"   ,wm,gfmats_local(ispin,ispin,iorb,jorb,:))
           call splot("perSigma"//reg(suffix)//"_iw.vca"   ,wm,Smats_local(ispin,ispin,iorb,jorb,:))
          enddo
        enddo
      enddo
    endif
  !
  end subroutine get_local_gf


  !+------------------------------------------------------------------+
  !Periodization functions SIGMA-SCHEME
  !+------------------------------------------------------------------+


 subroutine periodize_sigma_scheme(kpoint)
    integer                                                     :: ilat,jlat,ispin,iorb,jorb,ii
    real(8),dimension(Ndim)                                     :: kpoint,ind1,ind2
    complex(8),dimension(Nlat,Nlat,Nspin,Nspin,Norb,Norb,Lmats) :: tmpSigmaMat
    !
    if(.not.allocated(Smats_periodized))allocate(Smats_periodized(Nspin,Nspin,Norb,Norb,Lmats))
    if(.not.allocated(Smats_periodized_lso))allocate(Smats_periodized_lso(Nspin*Norb,Nspin*Norb,Lmats))
    !
    Smats_periodized=zero
    Smats_periodized_lso=zero
    !
    call vca_get_sigma_matsubara(tmpSigmaMat)
    do ii=1,Lmats
      do ilat=1,Nlat
        ind1=N2indices(ilat)   
        do jlat=1,Nlat
          ind2=N2indices(jlat)
          Smats_periodized(:,:,:,:,ii)=Smats_periodized(:,:,:,:,ii)+exp(-xi*dot_product(kpoint,ind1-ind2))*tmpSigmaMat(ilat,jlat,:,:,:,:,ii)/Nlat
        enddo
      enddo
    enddo
    !
    do ii=1,Lmats
      Smats_periodized_lso(:,:,ii)=vca_nn2so_reshape(Smats_periodized(:,:,:,:,ii),Nspin,Norb)
    enddo
    !
  end subroutine periodize_sigma_scheme


  subroutine get_local_sigma()
    integer                                         :: ik,ispin,iorb,jorb
    character(len=30)                               :: suffix
    !
    if(master)then
    !
    if(.not.allocated(Smats_local))allocate(Smats_local(Nspin,Nspin,Norb,Norb,Lmats))
    smats_local=zero
    !
    allocate(kgrid_test(Nkpts**ndim,Ndim)) 
    call TB_build_kgrid([Nkpts,Nkpts],kgrid_test)
    !
    print*,"Calculating SigmaLoc sigma scheme"  
    !
    call start_timer
    !
    do ik=1,Nkpts**ndim
      call periodize_sigma_scheme(kgrid_test(ik,:))
      do ix=1,Lmats
        smats_local(:,:,:,:,ix)=smats_local(:,:,:,:,ix)+Smats_periodized(:,:,:,:,ix)/(Nkpts**Ndim)
      enddo
      call eta(ik,Nkpts**ndim)
    enddo
    !
    call stop_timer
    !
    do iorb=1,Norb
     do jorb=1,Norb
       do ispin=1,Nspin
         suffix="_l"//str(iorb)//str(jorb)//"_s"//str(ispin)//"_"//str(scheme)//"_scheme"
         call splot("perSigma"//reg(suffix)//"_iw.vca"   ,wm,Smats_local(ispin,ispin,iorb,jorb,:))
        enddo
      enddo
    enddo
  endif
  !
end subroutine get_local_sigma

 
  !+------------------------------------------------------------------+
  !Auxilliary functions
  !+------------------------------------------------------------------+


  subroutine build_z2_indices()
    integer                                                         :: unit
    integer                                                         :: z2
    !
    !Evaluate the Z2 index:
    z2 = z2_number(reshape( [ [0,0] , [0,1] , [1,0] , [1,1] ] , [2,4])*pi/Nx)
    unit=free_unit()
    open(unit,file="z2_invariant.ed")
    write(unit,*)z2
    close(unit)
  end subroutine build_z2_indices

  function z2_number(ktrims) result(z2)
    real(8),dimension(:,:),intent(in)                                    :: ktrims
    complex(8),dimension(Nlat*Nspin*Norb,Nlat*Nspin*Norb,size(Ktrims,2)) :: Htrims
    complex(8),dimension(size(Ktrims,2))                                 :: Delta
    integer                                                              :: z2
    integer                                                              :: i,j,Ntrim,itrim,Nocc
    real(8),dimension(Nlat*Nspin*Norb)                                   :: Eigval
    !
    Ntrim=size(Ktrims,2)
    !
    do itrim=1,Ntrim
       Htrims(:,:,itrim) = hk_bhz_clusterbase(Ktrims(:,itrim),Nlat*Nspin*Norb)
       call eigh(Htrims(:,:,itrim),Eigval)
       Delta(itrim)=one
       do Nocc=1,Nlat*Nspin*Norb/2
         if(abs(Htrims(2,Nocc,itrim))>0.125 .or. abs(Htrims(4,Nocc,itrim))>0.125 .or. abs(Htrims(6,Nocc,itrim))>0.125 .or. abs(Htrims(8,Nocc,itrim))>0.125&
          .or. abs(Htrims(10,Nocc,itrim))>0.125 .or. abs(Htrims(12,Nocc,itrim))>0.125 .or. abs(Htrims(14,Nocc,itrim))>0.125 .or. abs(Htrims(16,Nocc,itrim))>0.125) then
          Delta(itrim)=Delta(itrim)*xi
          print*,"lol",Delta(itrim)
          else
          print*,"asd", Delta(itrim)
          print*,Htrims(:,Nocc,itrim)
          endif
       enddo
      PRINT*,"NEWTRIM"
    enddo
    !
    z2=product(Delta(:))
    if(z2>0)then
       z2=0
    else
       z2=1
    end if
    print*,"Z2",z2
    !
  end function z2_number

  function Zrenorm(sigma) result(Zmats)
    !
    integer                                         :: i
    complex(8),dimension(Nspin*Norb,Nspin*Norb)     :: Zmats,Sigma
    !
    Zmats=zero
    Sigma=zero
    !
    if (usez) then
      do i=1,Nspin*Norb
      Zmats(i,i)  = 1.d0/abs( 1.d0 +  abs(dimag(sigma(i,i))/(pi/beta)) )
      end do
    else
      Zmats=eye(Nspin*Norb)
    endif 
    !
  end function Zrenorm

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


end program vca_bhz_2d

