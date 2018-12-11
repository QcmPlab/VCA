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

 
  if(Norb/=1)stop "Norb != 1"
  if(Nspin/=1)stop "Nspin != 1"
  !
  ! FIXME: 
  t_var=1.531105953d0
  t=1.0d0
  mu=0.d0
  mu_var=0.d0

  if(.not.allocated(wm))allocate(wm(Lmats))
  if(.not.allocated(wr))allocate(wr(Lreal))
  wm     = pi/beta*real(2*arange(1,Lmats)-1,8)
  wr     = linspace(wini,wfin,Lreal)

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
  !call get_local_gf()
  call get_Akw()
  !call get_dos()
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




  subroutine generate_hk()
    integer                             :: ik,ii,ispin,iorb,unit,jj
    real(8),dimension(Nkpts,1)          :: kgrid
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
    !    
  end subroutine generate_hk


  function tk(kpoint) result(hopping_matrix)
    integer                                                               :: ilat,ispin,iorb,unit,jj
    real(8),dimension(Ndim),intent(in)                                    :: kpoint
    complex(8),dimension(Nlat,Nlat,Nspin,Nspin,Norb,Norb)                 :: hopping_matrix
    !
    do ilat=1,Nlat
        hopping_matrix(ilat,ilat,1,1,1,1)= -mu
        !
        if(ilat>1)hopping_matrix(ilat,ilat-1,1,1,1,1)= -t
        if(ilat<Nlat)hopping_matrix(ilat,ilat+1,1,1,1,1)= -t
    enddo
    !
    hopping_matrix(1,Nlat,1,1,1,1)=hopping_matrix(1,Nlat,1,1,1,1)-t*exp(xi*kpoint(1)*Nlat)
    hopping_matrix(Nlat,1,1,1,1,1)=hopping_matrix(Nlat,1,1,1,1,1)-t*exp(-xi*kpoint(1)*Nlat)
    ! 
  end function tk


  !+------------------------------------------------------------------+
  !Periodization functions G-SCHEME
  !+------------------------------------------------------------------+


  subroutine periodize_g_scheme(kpoint)
    integer                                                     :: ilat,jlat,ispin,iorb,ii
    real(8),dimension(Ndim)                                     :: kpoint
    complex(8),allocatable,dimension(:,:,:,:,:,:)               :: Xgfprime ![Nlat][Nlat][Nspin][Nspin][Norb][Norb]
    complex(8),allocatable,dimension(:,:)                       :: Xgfprime_lso ![Nlso][Nlso]
    complex(8),allocatable,dimension(:,:)                       :: XVk_lso ![Nlso][Nlso]
    complex(8),allocatable,dimension(:,:,:,:,:,:,:)             :: Xgfreal_unperiodized![Nlat][Nlat][Nspin][Nspin][Norb][Norb][Lmats]
    complex(8),allocatable,dimension(:,:,:,:,:,:,:)             :: Xgfmats_unperiodized ![Nlat][Nlat][Nspin][Nspin][Norb][Norb][Lreal]
    if(master)then
    !
    allocate(Xgfmats_unperiodized(Nlat,Nlat,Nspin,Nspin,Norb,Norb,Lmats))
    allocate(Xgfreal_unperiodized(Nlat,Nlat,Nspin,Nspin,Norb,Norb,Lreal))
    allocate(Xgfprime(Nlat,Nlat,Nspin,Nspin,Norb,Norb))
    allocate(Xgfprime_lso(Nlat*Nspin*Norb,Nlat*Nspin*Norb))
    allocate(XVk_lso(Nlat*Nspin*Norb,Nlat*Nspin*Norb))
    if(.not.allocated(gfmats_periodized))allocate(gfmats_periodized(Nspin,Nspin,Norb,Norb,Lmats))
    if(.not.allocated(gfreal_periodized))allocate(gfreal_periodized(Nspin,Nspin,Norb,Norb,Lreal))
    !
    Xgfmats_unperiodized=zero
    Xgfreal_unperiodized=zero
    Xgfprime=zero
    Xgfprime_lso=zero
    gfmats_periodized=zero
    gfreal_periodized=zero
    XVk_lso=zero
    !
    !
    XVk_lso=vca_nnn2lso_reshape(tk(kpoint)-t_prime,Nlat,Nspin,Norb)
    !
    do ii=1,Lmats  
      call vca_gf_cluster(xi*wm(ii),Xgfprime)
      Xgfprime_lso=vca_nnn2lso_reshape(Xgfprime,Nlat,Nspin,Norb)
      call inv(Xgfprime_lso)
      Xgfprime_lso=Xgfprime_lso-XVk_lso
      call inv(Xgfprime_lso)
      Xgfmats_unperiodized(:,:,:,:,:,:,ii)=vca_lso2nnn_reshape(Xgfprime_lso,Nlat,Nspin,Norb)
    enddo
    !
    do ii=1,Lreal   
      call vca_gf_cluster(dcmplx(wr(ii),eps),Xgfprime)
      Xgfprime_lso=vca_nnn2lso_reshape(Xgfprime,Nlat,Nspin,Norb)
      call inv(Xgfprime_lso)
      Xgfprime_lso=Xgfprime_lso-XVk_lso
      call inv(Xgfprime_lso)
      Xgfreal_unperiodized(:,:,:,:,:,:,ii)=vca_lso2nnn_reshape(Xgfprime_lso,Nlat,Nspin,Norb)
    enddo
    !
    !
    do ii=1,Lmats
      do ilat=1,Nlat
        do jlat=1,Nlat
          gfmats_periodized(:,:,:,:,ii)=gfmats_periodized(:,:,:,:,ii)+exp(-xi*kpoint(1)*(ilat-jlat))*Xgfmats_unperiodized(ilat,jlat,:,:,:,:,ii)/Nlat
        enddo
      enddo
    enddo
    do ii=1,Lreal   
      do ilat=1,Nlat
        do jlat=1,Nlat
          gfreal_periodized(:,:,:,:,ii)=gfreal_periodized(:,:,:,:,ii)+exp(-xi*kpoint(1)*(ilat-jlat))*Xgfreal_unperiodized(ilat,jlat,:,:,:,:,ii)/Nlat
        enddo
      enddo
    enddo
    !
    deallocate(Xgfmats_unperiodized)
    deallocate(Xgfreal_unperiodized)
    deallocate(XVk_lso)
    endif   
  end subroutine periodize_g_scheme



  subroutine build_sigma_g_scheme(kpoint)
    integer                                                     :: i,ispin,iorb,ii
    real(8),dimension(Ndim)                                     :: kpoint
    complex(8),dimension(Nspin,Nspin,Norb,Norb,Lmats)           :: invertedG0mats,invertedGmats
    complex(8),dimension(Nspin,Nspin,Norb,Norb,Lreal)           :: invertedG0real,invertedGreal
    if(master)then
    !
    if(.not.allocated(Smats_periodized))allocate(Smats_periodized(Nspin,Nspin,Norb,Norb,Lmats))
    if(.not.allocated(Sreal_periodized))allocate(Sreal_periodized(Nspin,Nspin,Norb,Norb,Lreal))
    invertedG0mats = zero
    invertedGmats  = zero
    invertedG0real = zero
    invertedGreal  = zero
    Smats_periodized  = zero
    Sreal_periodized  = zero
    !
    !Get G0^-1
    !invertedG0mats = invG0_bath_mats(dcmplx(0d0,wm(:)),vca_bath)
    !invertedG0real = invG0_bath_real(dcmplx(wr(:),eps),vca_bath)
    !
    !Get G0^-1
    do ispin=1,Nspin
      do iorb=1,Norb
        do ii=1,Lmats
            invertedG0mats(ispin,ispin,iorb,iorb,ii) = (xi*wm(ii)+xmu)  -(-2.d0*t*cos(kpoint(1)))                ! FIXME: ad-hoc solution
        enddo
        do ii=1,Lreal
            invertedG0real(ispin,ispin,iorb,iorb,ii) = (wr(ii)+xmu)  -(-2.d0*t*cos(kpoint(1)))                ! FIXME: ad-hoc solution
        enddo
      enddo
    enddo
    !
    !Get Gimp^-1
    call periodize_g_scheme(kpoint)
    do ispin=1,Nspin
      do iorb=1,Norb
         invertedGmats(ispin,ispin,iorb,iorb,:) = one/gfmats_periodized(ispin,ispin,iorb,iorb,:)
         invertedGreal(ispin,ispin,iorb,iorb,:) = one/gfreal_periodized(ispin,ispin,iorb,iorb,:)
      enddo
    enddo
    !Get Sigma functions: Sigma= G0^-1 - G^-1
    Smats_periodized=zero
    Sreal_periodized=zero
    do ispin=1,Nspin
      do iorb=1,Norb
         Smats_periodized(ispin,ispin,iorb,iorb,:) = invertedG0mats(ispin,ispin,iorb,iorb,:) - invertedGmats(ispin,ispin,iorb,iorb,:)
         Sreal_periodized(ispin,ispin,iorb,iorb,:) = invertedG0real(ispin,ispin,iorb,iorb,:) - invertedGreal(ispin,ispin,iorb,iorb,:)
      enddo
    enddo
    !
    endif
  end subroutine build_sigma_g_scheme

  !+------------------------------------------------------------------+
  !Periodization functions SIGMA-SCHEME
  !+------------------------------------------------------------------+


 subroutine periodize_sigma_scheme(kpoint)
    integer                                                     :: ilat,jlat,ispin,iorb,ii
    real(8),dimension(Ndim)                                     :: kpoint,ind1,ind2
    complex(8),allocatable,dimension(:,:,:,:,:,:)               :: gfprime ![Nlat][Nlat][Nspin][Nspin][Norb][Norb]
    complex(8),allocatable,dimension(:,:)                       :: gfprime_lso,invertedG0cluster ![Nlso][Nlso]
    complex(8),allocatable,dimension(:,:,:,:,:,:,:)             :: Sreal_unperiodized![Nlat][Nlat][Nspin][Nspin][Norb][Norb][Lmats]
    complex(8),allocatable,dimension(:,:,:,:,:,:,:)             :: Smats_unperiodized ![Nlat][Nlat][Nspin][Nspin][Norb][Norb][Lreal]
    !
    if(master)then
    allocate(Smats_unperiodized(Nlat,Nlat,Nspin,Nspin,Norb,Norb,Lmats))
    allocate(Sreal_unperiodized(Nlat,Nlat,Nspin,Nspin,Norb,Norb,Lreal))
    allocate(gfprime(Nlat,Nlat,Nspin,Nspin,Norb,Norb))
    allocate(gfprime_lso(Nlat*Nspin*Norb,Nlat*Nspin*Norb))
    allocate(invertedG0cluster(Nlat*Nspin*Norb,Nlat*Nspin*Norb))
    if(.not.allocated(Smats_periodized))allocate(Smats_periodized(Nspin,Nspin,Norb,Norb,Lmats))
    if(.not.allocated(Sreal_periodized))allocate(Sreal_periodized(Nspin,Nspin,Norb,Norb,Lreal))
    !
    Smats_unperiodized=zero
    Sreal_unperiodized=zero
    gfprime=zero
    gfprime_lso=zero
    invertedG0cluster=zero
    Smats_periodized=zero
    Sreal_periodized=zero
    !
    !
    !
    do ii=1,Lmats    
        invertedG0cluster=(xi*wm(ii)+xmu)*eye(Nlat*Nspin*Norb)-vca_nnn2lso_reshape(t_prime,Nlat,Nspin,Norb)
        call vca_gf_cluster(xi*wm(ii),gfprime)
        gfprime_lso=vca_nnn2lso_reshape(gfprime,Nlat,Nspin,Norb)
        call inv(gfprime_lso)
        gfprime_lso=invertedG0cluster-gfprime_lso
        Smats_unperiodized(:,:,:,:,:,:,ii)=vca_lso2nnn_reshape(gfprime_lso,Nlat,Nspin,Norb)
    enddo
    !
    do ii=1,Lreal
        invertedG0cluster=(wr(ii)+xmu)*eye(Nlat*Nspin*Norb)-vca_nnn2lso_reshape(t_prime,Nlat,Nspin,Norb)   
        call vca_gf_cluster(dcmplx(wr(ii),eps),gfprime)
        gfprime_lso=vca_nnn2lso_reshape(gfprime,Nlat,Nspin,Norb)
        call inv(gfprime_lso)
        gfprime_lso=invertedG0cluster-gfprime_lso
        Sreal_unperiodized(:,:,:,:,:,:,ii)=vca_lso2nnn_reshape(gfprime_lso,Nlat,Nspin,Norb)
    enddo
    !
    do ii=1,Lmats
        do ilat=1,Nlat
          do jlat=1,Nlat
            Smats_periodized(:,:,:,:,ii)=Smats_periodized(:,:,:,:,ii)+exp(-xi*kpoint(1)*(ilat-jlat))*Smats_unperiodized(ilat,jlat,:,:,:,:,ii)/Nlat
          enddo
        enddo
    enddo
    !
    do ii=1,Lreal   
     do ilat=1,Nlat
        do jlat=1,Nlat
          Sreal_periodized(:,:,:,:,ii)=Sreal_periodized(:,:,:,:,ii)+exp(-xi*kpoint(1)*(ilat-jlat))*Sreal_unperiodized(ilat,jlat,:,:,:,:,ii)/Nlat
        enddo
      enddo
    enddo
    !
    deallocate(Smats_unperiodized)
    deallocate(Sreal_unperiodized)
    deallocate(gfprime)
    deallocate(gfprime_lso)
    deallocate(invertedG0cluster)
    endif   
  end subroutine periodize_sigma_scheme


  subroutine build_g_sigma_scheme(kpoint)
    integer                                                     :: i,ispin,iorb,ii
    real(8),dimension(Ndim)                                     :: kpoint
    complex(8),dimension(Nspin,Nspin,Norb,Norb,Lmats)           :: invertedG0mats
    complex(8),dimension(Nspin,Nspin,Norb,Norb,Lreal)           :: invertedG0real
    complex(8),allocatable,dimension(:,:)                       :: tmpmat ![Nso][Nso]
    !
    if(master)then
    if(.not.allocated(gfmats_periodized))allocate(gfmats_periodized(Nspin,Nspin,Norb,Norb,Lmats))
    if(.not.allocated(gfreal_periodized))allocate(gfreal_periodized(Nspin,Nspin,Norb,Norb,Lreal))
    allocate(tmpmat(Nspin*Norb,Nspin*Norb))
    invertedG0mats = zero
    invertedG0real = zero
    gfmats_periodized  = zero
    gfreal_periodized  = zero
    tmpmat = zero
    !
    !Get G0^-1
    !invertedG0mats = invg0_bath_mats(dcmplx(0d0,wm(:)),vca_bath)
    !invertedG0real = invg0_bath_real(dcmplx(wr(:),eps),vca_bath)
    !
    !Get G0^-1
    do ispin=1,Nspin
      do iorb=1,Norb
        do ii=1,Lmats
            invertedG0mats(ispin,ispin,iorb,iorb,ii) = (xi*wm(ii)+xmu)  -(-2.d0*t*cos(kpoint(1)))             ! FIXME: ad-hoc solution
        enddo
        do ii=1,Lreal
            invertedG0real(ispin,ispin,iorb,iorb,ii) = (wr(ii)+xmu)   -(-2.d0*t*cos(kpoint(1)))               ! FIXME: ad-hoc solution
        enddo
      enddo
    enddo
    call periodize_sigma_scheme(kpoint)
    !Get G: G^-1= G0-Sigma
    do ispin=1,Nspin
      do iorb=1,Norb
         gfmats_periodized(ispin,ispin,iorb,iorb,:) = invertedG0mats(ispin,ispin,iorb,iorb,:) - Smats_periodized(ispin,ispin,iorb,iorb,:)
         gfreal_periodized(ispin,ispin,iorb,iorb,:) = invertedG0real(ispin,ispin,iorb,iorb,:) - Sreal_periodized(ispin,ispin,iorb,iorb,:)
      enddo
    enddo
    !
    do ii=1,Lmats
        tmpmat=vca_nn2so_reshape(gfmats_periodized(:,:,:,:,ii),Nspin,Norb)
        call inv(tmpmat)
        gfmats_periodized(:,:,:,:,ii)=vca_so2nn_reshape(tmpmat,Nspin,Norb)
    enddo
    !
    do ii=1,Lmats
        tmpmat=vca_nn2so_reshape(gfreal_periodized(:,:,:,:,ii),Nspin,Norb)
        call inv(tmpmat)
        gfreal_periodized(:,:,:,:,ii)=vca_so2nn_reshape(tmpmat,Nspin,Norb)
    enddo
    !
    deallocate(tmpmat)
    endif
  end subroutine build_g_sigma_scheme

  !---------------------------------------------------------------------
  !PURPOSE: GET local GF
  !---------------------------------------------------------------------
  subroutine get_local_gf()
    integer                                         :: ik,ispin,iorb
    character(len=30)                               :: suffix
    complex(8),allocatable,dimension(:,:,:,:,:)     :: gtest_mats,gtest_real,sigmatest_mats,sigmatest_real
    !
    if(master)then
    !
      allocate(gtest_real(Nspin,Nspin,Norb,Norb,Lmats))
      allocate(sigmatest_real(Nspin,Nspin,Norb,Norb,Lmats))
      allocate(gtest_mats(Nspin,Nspin,Norb,Norb,Lmats))
      allocate(sigmatest_mats(Nspin,Nspin,Norb,Norb,Lmats))
      gtest_mats=zero
      sigmatest_mats=zero
      gtest_real=zero
      sigmatest_real=zero
    !
      print*,"Printing local GF"
      call start_timer
      do ik=1,SAMPLING
          !call periodize_g_scheme([(2*pi/SAMPLING)*ik])
          call build_sigma_g_scheme([(2*pi/SAMPLING)*ik])
          do ix=1,Lmats
              gtest_mats(:,:,:,:,ix)=gtest_mats(:,:,:,:,ix)+gfmats_periodized(:,:,:,:,ix)/SAMPLING
              sigmatest_mats(:,:,:,:,ix)=sigmatest_mats(:,:,:,:,ix)+Smats_periodized(:,:,:,:,ix)/SAMPLING
          enddo
          do ix=1,Lreal
              gtest_real(:,:,:,:,ix)=gtest_real(:,:,:,:,ix)+gfreal_periodized(:,:,:,:,ix)/SAMPLING
              sigmatest_real(:,:,:,:,ix)=sigmatest_real(:,:,:,:,ix)+Sreal_periodized(:,:,:,:,ix)/SAMPLING
          enddo
          call eta(ik,SAMPLING)
      enddo
      call stop_timer
      !
      do iorb=1,Norb
        do ispin=1,Nspin
          suffix="_l"//str(iorb)//str(iorb)//"_s"//str(ispin)
          call splot("Gloc"//reg(suffix)//"_iw.vca"   ,wm,gtest_mats(ispin,ispin,iorb,iorb,:))
          call splot("Gloc"//reg(suffix)//"_realw.vca",wr,gtest_real(ispin,ispin,iorb,iorb,:))
          call splot("Sigmaloc"//reg(suffix)//"_iw.vca"   ,wm,sigmatest_mats(ispin,ispin,iorb,iorb,:))
          call splot("Sigmaloc"//reg(suffix)//"_realw.vca",wr,sigmatest_real(ispin,ispin,iorb,iorb,:))
        enddo
      enddo
    endif
end subroutine get_local_gf

  !---------------------------------------------------------------------
  !PURPOSE: GET A(k,w)
  !---------------------------------------------------------------------
  subroutine get_Akw()
    integer                                                 :: ik=0,iw
    integer,parameter                                       :: Lw=10000
    integer                                                 :: iorb,jorb
    integer                                                 :: ispin,jspin
    real(8)                                                 :: akrange
    complex(8),dimension(:,:),allocatable                   :: Akreal
    complex(8),dimension(:,:,:,:,:,:),allocatable           :: Gkreal,Gkreal_interpolated
    real(8),dimension(:)                                    :: Kgrid(SAMPLING),w_int(Lw)
    character(len=30)                                       :: suffix
    !
    if(master)then
    !
    print*,"Build A(k,w)"
    !
    allocate(Gkreal(SAMPLING,Nspin,Nspin,Norb,Norb,Lreal));Gkreal=zero
    allocate(Gkreal_interpolated(SAMPLING,Nspin,Nspin,Norb,Norb,Lw));Gkreal_interpolated=zero
    allocate(Akreal(SAMPLING,Lw));Akreal=zero
    akrange=4d0
    w_int  = linspace(-akrange,akrange,Lw)
    !
    do ik=1,SAMPLING
        kgrid(ik)=pi/SAMPLING*(ik-1)
    enddo
    !
    print*,"Retrieving G(k,w)"
    call start_timer
    do ik=1,SAMPLING
        !call build_g_sigma_scheme(kgrid(ik))
        call periodize_g_scheme(kgrid(ik))
        Gkreal(ik,:,:,:,:,:)=gfreal_periodized(:,:,:,:,:)
        call eta(ik,SAMPLING)
    enddo
    call stop_timer
    !
    print*,"Interpolating G(k,w)"
    do ik=1,SAMPLING
      do ispin=1,Nspin
        do iorb=1,Norb
          call cubic_spline(wr,Gkreal(ik,ispin,ispin,iorb,iorb,:),w_int,Gkreal_interpolated(ik,ispin,ispin,iorb,iorb,:))
        enddo
      enddo
    enddo
    !
    do ispin=1,Nspin
      do iorb=1,Norb
        Akreal = Akreal - dimag(Gkreal_interpolated(:,ispin,ispin,iorb,iorb,:))/pi/Nspin/Norb
      enddo
    enddo
    print*,"Printing"
    call splot3d("Akw_real.dat",kgrid,w_int,Akreal) 
    endif
end subroutine get_Akw



  !---------------------------------------------------------------------
  !PURPOSE: GET DOS(w)
  !---------------------------------------------------------------------
  subroutine get_dos()
    integer                                                 :: ik=0,iw
    integer                                                 :: iorb,jorb
    integer                                                 :: ispin,jspin
    complex(8),dimension(:),allocatable                     :: dosreal
    complex(8),dimension(:,:,:,:,:),allocatable             :: Greal
    real(8),dimension(:),allocatable                        :: Kgrid
    character(len=30)                                       :: suffix
    !
    if(master)then
    !
    print*,"Build DOS"
    !
    allocate(Greal(Nspin,Nspin,Norb,Norb,Lreal));Greal=zero
    allocate(dosreal(Lreal));dosreal=zero
    allocate(kgrid(SAMPLING))
    !
    do ik=1,SAMPLING
        kgrid(ik)=pi/SAMPLING*(ik-1)
    enddo
    !
    call start_timer
    do ik=1,SAMPLING
        call periodize_g_scheme(kgrid(ik))
        Greal=Greal+gfreal_periodized/SAMPLING
        call eta(ik,SAMPLING)
    enddo
    call stop_timer
    !
    do ispin=1,Nspin
      do iorb=1,Norb
        dosreal = dosreal - 2*dimag(Greal(ispin,ispin,iorb,iorb,:))/(Nspin*Norb)
      enddo
    enddo
    !wr is already allocated
    call splot("dos_real.dat",wr,dosreal) 
    endif
end subroutine get_dos


end program vca_chain1d



