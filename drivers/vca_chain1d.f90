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
  integer                                         :: ix,iy,ik
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
  integer                                         :: Nx,Ny,Lx,Ly,Rx,Ry
  integer                                         :: unit
  integer                                         :: comm,rank
  logical                                         :: master,wloop,wmin
  integer                                         :: nloop
  real(8)                                         :: mu,t,t_var,mu_var
  real(8),dimension(:),allocatable                :: ts_array,omega_array
  integer,dimension(1)                            :: min_loc
  complex(8),allocatable,dimension(:,:,:,:,:)     :: gfmats_periodized ![Nspin][Nspin][Norb][Norb][Lmats]
  complex(8),allocatable,dimension(:,:,:,:,:)     :: gfreal_periodized ![Nspin][Nspin][Norb][Norb][Lreal]
  complex(8),allocatable,dimension(:,:,:,:,:)     :: impSmats_periodized         ![Nspin][Nspin][Norb][Norb][Lmats]
  complex(8),allocatable,dimension(:,:,:,:,:)     :: impSreal_periodized         ![Nspin][Nspin][Norb][Norb][Lreal]
  complex(8),allocatable,dimension(:,:,:,:,:)     :: gtest_mats,gtest_Real,sigmatest_mats,sigmatest_real


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
  t_var=1.531105953d0
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
  allocate(gtest_real(Nspin,Nspin,Norb,Norb,Lmats))
  allocate(sigmatest_real(Nspin,Nspin,Norb,Norb,Lmats))
  allocate(gtest_mats(Nspin,Nspin,Norb,Norb,Lmats))
  allocate(sigmatest_mats(Nspin,Nspin,Norb,Norb,Lmats))
  gtest_mats=zero
  sigmatest_mats=zero
  gtest_real=zero
  sigmatest_real=zero
  !
  do ik=1,Nkpts
      !print*,ik
      !call periodize_g_scheme([(2*pi/Nkpts)*ik])
      call build_sigma_g_scheme([(2*pi/Nkpts)*ik])
      do ix=1,Lmats
          gtest_mats(:,:,:,:,ix)=gtest_mats(:,:,:,:,ix)+gfmats_periodized(:,:,:,:,ix)/Nkpts
          sigmatest_mats(:,:,:,:,ix)=sigmatest_mats(:,:,:,:,ix)+impSmats_periodized(:,:,:,:,ix)/Nkpts
      enddo
      do ix=1,Lreal
          gtest_real(:,:,:,:,ix)=gtest_real(:,:,:,:,ix)+gfreal_periodized(:,:,:,:,ix)/Nkpts
          sigmatest_real(:,:,:,:,ix)=sigmatest_real(:,:,:,:,ix)+impSreal_periodized(:,:,:,:,ix)/Nkpts
      enddo
  enddo
  gfmats_periodized=gtest_mats
  impSmats_periodized=sigmatest_mats
  gfreal_periodized=gtest_real
  impSreal_periodized=sigmatest_real
  call print_periodized()

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


  subroutine periodize_g_scheme(kpoint)
    integer                                                     :: ilat,jlat,ispin,iorb,ii
    real(8),dimension(Ndim)                                     :: kpoint
    complex(8),allocatable,dimension(:,:,:,:,:,:)               :: gfprime ![Nlat][Nlat][Nspin][Nspin][Norb][Norb]
    complex(8),allocatable,dimension(:,:)                       :: gfprime_lso ![Nlso][Nlso]
    complex(8),allocatable,dimension(:,:)                       :: Vk_lso ![Nlso][Nlso]
    complex(8),allocatable,dimension(:,:,:,:,:,:,:)             :: gfreal_unperiodized![Nlat][Nlat][Nspin][Nspin][Norb][Norb][Lmats]
    complex(8),allocatable,dimension(:,:,:,:,:,:,:)             :: gfmats_unperiodized ![Nlat][Nlat][Nspin][Nspin][Norb][Norb][Lreal]
    !
    if(.not.allocated(wm))allocate(wm(Lmats))
    if(.not.allocated(wr))allocate(wr(Lreal))
    wm     = pi/beta*real(2*arange(1,Lmats)-1,8)
    wr     = linspace(wini,wfin,Lreal)
    !
    allocate(gfmats_unperiodized(Nlat,Nlat,Nspin,Nspin,Norb,Norb,Lmats))
    allocate(gfreal_unperiodized(Nlat,Nlat,Nspin,Nspin,Norb,Norb,Lreal))
    allocate(gfprime(Nlat,Nlat,Nspin,Nspin,Norb,Norb))
    allocate(gfprime_lso(Nlat*Nspin*Norb,Nlat*Nspin*Norb))
    allocate(Vk_lso(Nlat*Nspin*Norb,Nlat*Nspin*Norb))
    if(.not.allocated(gfmats_periodized))allocate(gfmats_periodized(Nspin,Nspin,Norb,Norb,Lmats))
    if(.not.allocated(gfreal_periodized))allocate(gfreal_periodized(Nspin,Nspin,Norb,Norb,Lreal))
    !
    gfmats_unperiodized=zero
    gfreal_unperiodized=zero
    gfprime=zero
    gfmats_periodized=zero
    gfreal_periodized=zero
    Vk_lso=zero
    gfprime_lso=zero
    !
    !
    Vk_lso=vca_nnn2lso_reshape(tk(kpoint)-t_prime,Nlat,Nspin,Norb)
    !
    do ii=1,Lmats  
    gfprime_lso=zero  
        call vca_gf_cluster(xi*wm(ii),gfprime)
        gfprime_lso=vca_nnn2lso_reshape(gfprime,Nlat,Nspin,Norb)
        call inv(gfprime_lso)
        gfprime_lso=gfprime_lso-Vk_lso
        call inv(gfprime_lso)
        gfmats_unperiodized(:,:,:,:,:,:,ii)=vca_lso2nnn_reshape(gfprime_lso,Nlat,Nspin,Norb)
    enddo
    !
    do ii=1,Lreal   
    gfprime_lso=zero 
        call vca_gf_cluster(dcmplx(wr(ii),eps),gfprime)
        gfprime_lso=vca_nnn2lso_reshape(gfprime,Nlat,Nspin,Norb)
        call inv(gfprime_lso)
        gfprime_lso=gfprime_lso-Vk_lso
        call inv(gfprime_lso)
        gfreal_unperiodized(:,:,:,:,:,:,ii)=vca_lso2nnn_reshape(gfprime_lso,Nlat,Nspin,Norb)
    enddo
    !
    do ii=1,Lmats
      do ilat=1,Nlat
          do jlat=1,Nlat
           gfmats_periodized(:,:,:,:,ii)=gfmats_periodized(:,:,:,:,ii)+exp(-xi*kpoint(1)*(ilat-jlat))*gfmats_unperiodized(ilat,jlat,:,:,:,:,ii)/Nlat
          enddo
      enddo
    enddo
    !
    do ii=1,Lreal   
     do ilat=1,Nlat
        do jlat=1,Nlat
          gfreal_periodized(:,:,:,:,ii)=gfreal_periodized(:,:,:,:,ii)+exp(-xi*kpoint(1)*(ilat-jlat))*gfreal_unperiodized(ilat,jlat,:,:,:,:,ii)/Nlat
        enddo
      enddo
    enddo
    !
    if(allocated(wm))deallocate(wm)
    if(allocated(wr))deallocate(wr) 
    !   
  end subroutine periodize_g_scheme



  subroutine build_sigma_g_scheme(kpoint)
    integer                                                     :: i,ispin,iorb,ii
    real(8),dimension(Ndim)                                     :: kpoint
    complex(8),dimension(Nspin,Nspin,Norb,Norb,Lmats)           :: invG0mats,invGmats
    complex(8),dimension(Nspin,Nspin,Norb,Norb,Lreal)           :: invG0real,invGreal
    !
    if(.not.allocated(wm))allocate(wm(Lmats))
    if(.not.allocated(wr))allocate(wr(Lreal))
    wm     = pi/beta*real(2*arange(1,Lmats)-1,8)
    wr     = linspace(wini,wfin,Lreal)
    !
    if(.not.allocated(impSmats_periodized))allocate(impSmats_periodized(Nspin,Nspin,Norb,Norb,Lmats))
    if(.not.allocated(impSreal_periodized))allocate(impSreal_periodized(Nspin,Nspin,Norb,Norb,Lreal))
    invG0mats = zero
    invGmats  = zero
    invG0real = zero
    invGreal  = zero
    impSmats_periodized  = zero
    impSreal_periodized  = zero
    !
    !Get G0^-1
    !invG0mats = invg0_bath_mats(dcmplx(0d0,wm(:)),vca_bath)
    !invG0real = invg0_bath_real(dcmplx(wr(:),eps),vca_bath)
    !
    !Get G0^-1
    do ispin=1,Nspin
      do iorb=1,Norb
        do ii=1,Lmats
            invG0mats(ispin,ispin,iorb,iorb,ii) = (xi*wm(ii)+xmu)  -(-2.d0*t*cos(kpoint(1)))                ! FIXME: ad-hoc solution
        enddo
        do ii=1,Lreal
            invG0real(ispin,ispin,iorb,iorb,ii) = (wr(ii)+xmu)  -(-2.d0*t*cos(kpoint(1)))                ! FIXME: ad-hoc solution
        enddo
      enddo
    enddo
    !
    !Get Gimp^-1
    call periodize_g_scheme(kpoint)
    do ispin=1,Nspin
      do iorb=1,Norb
         invGmats(ispin,ispin,iorb,iorb,:) = one/gfmats_periodized(ispin,ispin,iorb,iorb,:)
         invGreal(ispin,ispin,iorb,iorb,:) = one/gfreal_periodized(ispin,ispin,iorb,iorb,:)
      enddo
    enddo
    !Get Sigma functions: Sigma= G0^-1 - G^-1
    impSmats_periodized=zero
    impSreal_periodized=zero
    do ispin=1,Nspin
      do iorb=1,Norb
         impSmats_periodized(ispin,ispin,iorb,iorb,:) = invG0mats(ispin,ispin,iorb,iorb,:) - invGmats(ispin,ispin,iorb,iorb,:)
         impSreal_periodized(ispin,ispin,iorb,iorb,:) = invG0real(ispin,ispin,iorb,iorb,:) - invGreal(ispin,ispin,iorb,iorb,:)
      enddo
    enddo
    !
    !
    if(allocated(wm))deallocate(wm)
    if(allocated(wr))deallocate(wr)
    !
  end subroutine build_sigma_g_scheme

  !+------------------------------------------------------------------+
  !                         PRINT 
  !+------------------------------------------------------------------+  

  subroutine print_periodized()
    character(len=64) :: suffix
    integer           :: iorb,ispin
    !
    if(.not.allocated(wm))allocate(wm(Lmats))
    if(.not.allocated(wr))allocate(wr(Lreal))
    wm     = pi/beta*real(2*arange(1,Lmats)-1,8)
    wr     = linspace(wini,wfin,Lreal)
    !
    do iorb=1,Norb
     do ispin=1,Nspin
        suffix="_l"//str(iorb)//str(iorb)//"_s"//str(ispin)
        call splot("perG"//reg(suffix)//"_iw.vca"   ,wm,gfmats_periodized(ispin,ispin,iorb,iorb,:))
        call splot("perG"//reg(suffix)//"_realw.vca",wr,gfreal_periodized(ispin,ispin,iorb,iorb,:))
        call splot("perSigma"//reg(suffix)//"_iw.vca"   ,wm,impSmats_periodized(ispin,ispin,iorb,iorb,:))
        call splot("perSigma"//reg(suffix)//"_realw.vca",wr,impSreal_periodized(ispin,ispin,iorb,iorb,:))
      enddo
    enddo
    !
    if(allocated(wm))deallocate(wm)
    if(allocated(wr))deallocate(wr)
    !
  end subroutine print_periodized

 

end program vca_chain1d



