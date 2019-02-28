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
   !FIX!!!!!!!!!!!!!!!
  real(8)                                         :: mu,t,t_var,mu_var,ts_dummy1,ts_dummy2,DUMMY1,DUMMY2,dummy3,ts_dummy3,omegadummy
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
  Nlso = (Nx**Ndim)*Norb*Nspin
  Nlat=Nx**Ndim
  Ny=Nx
  !
  t_var=hopping
  t=hopping
  mu=0.d0
  mu_var=0.d0
  bandwidth=2.d0*Ndim*(2*t) !(2 times sum over dim of 2*t*max (cosine))

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

  
  Nb=vca_get_bath_dimension()
  allocate(Bath(Nb))
  call vca_init_solver(comm,bath)
 

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
     call get_local_gf()
     call get_Akw()
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
    call get_local_gf()
    call get_Akw()
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
    Vij=tij
    Eij=Uloc(1)/2d0
    print*,""
    print*,"------ Doing for ",tij," ------"
    call generate_tcluster()
    call generate_hk()
    !BATH VARIATIONAL SETUP
    do ix=1,Nlat
      do iy=1,Nspin
        do ik=1,Norb       
          call set_bath_component(bath,ix,iy,ik,e_component=[Eij])
          call set_bath_component(bath,ix,iy,ik,v_component=[Vij])
        enddo
      enddo
    enddo
    call vca_solve(comm,t_prime,h_k,bath)
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
    real(8),dimension(Nlso,Nlso)                 :: H0
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
          t_prime(ind1,ind1,ispin,ispin,1,1)= -mu_var
          if(ilat<Nx)then
            ind2=indices2N([ilat+1,jlat])
            t_prime(ind1,ind2,ispin,ispin,1,1)= -t_var
          endif
          if(ilat>1)then
            ind2=indices2N([ilat-1,jlat])
            t_prime(ind1,ind2,ispin,ispin,1,1)= -t_var
          endif
          if(jlat<Nx)then
            ind2=indices2N([ilat,jlat+1])
            t_prime(ind1,ind2,ispin,ispin,1,1)= -t_var

          endif
          if(jlat>1)then
            ind2=indices2N([ilat,jlat-1])
            t_prime(ind1,ind2,ispin,ispin,1,1)= -t_var
          endif
        enddo
      enddo
    enddo
    !
    H0=vca_nnn2lso_reshape(t_prime,Nlat,Nspin,Norb)
    !
    open(free_unit(unit),file=trim(file_))
    do ilat=1,Nlat*Nspin*Norb
       write(unit,"(5000(F5.2,1x))")(H0(ilat,jlat),jlat=1,Nlat*Nspin*Norb)
    enddo
    close(unit)
  end subroutine generate_tcluster




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
            hopping_matrix(ind1,ind2,ispin,ispin,1,1)= -t
          endif
          if(ilat>1)then
            ind2=indices2N([ilat-1,jlat])
            hopping_matrix(ind1,ind2,ispin,ispin,1,1)= -t
          endif
          if(jlat<Ny)then
            ind2=indices2N([ilat,jlat+1])
            hopping_matrix(ind1,ind2,ispin,ispin,1,1)= -t
          endif
          if(jlat>1)then
            ind2=indices2N([ilat,jlat-1])
            hopping_matrix(ind1,ind2,ispin,ispin,1,1)= -t
          endif
        enddo
      enddo
    enddo
    !
    do ilat=1,Nx
      do ispin=1,Nspin
        ind1=indices2N([1,ilat])
        ind2=indices2N([Nx,ilat])
        hopping_matrix(ind1,ind2,ispin,ispin,1,1)=hopping_matrix(ind1,ind2,ispin,ispin,1,1) -t*exp(xi*kpoint(2)*Nx)
        hopping_matrix(ind2,ind1,ispin,ispin,1,1)=hopping_matrix(ind2,ind1,ispin,ispin,1,1) -t*exp(-xi*kpoint(2)*Nx)
        !
        ind1=indices2N([ilat,1])
        ind2=indices2N([ilat,Ny])
        hopping_matrix(ind1,ind2,ispin,ispin,1,1)=hopping_matrix(ind1,ind2,ispin,ispin,1,1) -t*exp(xi*kpoint(1)*Ny)
        hopping_matrix(ind2,ind1,ispin,ispin,1,1)=hopping_matrix(ind2,ind1,ispin,ispin,1,1) -t*exp(-xi*kpoint(1)*Ny)
      enddo
    enddo
    ! 
  end function tk


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
    do ilat=1,Nlat
      do ispin=1,Nspin
         do iorb=1,Norb
           DELTA(ilat,ilat,ispin,ispin,iorb,iorb)=sum( vps(:)*vps(:)/(freq - eps(:)+XMU) )
         enddo
      enddo
    enddo
  end function set_delta


 subroutine periodize_g_scheme(kpoint)
    integer                                                     :: ilat,jlat,ispin,iorb,ii
    real(8),dimension(Ndim)                                     :: kpoint,ind1,ind2
    complex(8),allocatable,dimension(:,:,:,:,:,:)               :: gfprime ![Nlat][Nlat][Nspin][Nspin][Norb][Norb]

    complex(8),allocatable,dimension(:,:)                       :: gfprime_lso ![Nlso][Nlso]
    complex(8),allocatable,dimension(:,:)                       :: Vk_lso ![Nlso][Nlso]
    complex(8),allocatable,dimension(:,:,:,:,:,:,:)             :: gfreal_unperiodized![Nlat][Nlat][Nspin][Nspin][Norb][Norb][Lmats]
    complex(8),allocatable,dimension(:,:,:,:,:,:,:)             :: gfmats_unperiodized ![Nlat][Nlat][Nspin][Nspin][Norb][Norb][Lreal]
    !
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
    gfprime_lso=zero
    Vk_lso=zero
    gfmats_periodized=zero
    gfreal_periodized=zero
    !
    !
    !
    do ii=1,Lmats
      Vk_lso=vca_nnn2lso_reshape(tk(kpoint)-t_prime-set_delta(xi*wm(ii),[ts],[Uloc(1)/2]),Nlat,Nspin,Norb)    
      call vca_gf_cluster(xi*wm(ii),gfprime)
      gfprime_lso=vca_nnn2lso_reshape(gfprime,Nlat,Nspin,Norb)
      call inv(gfprime_lso)
      gfprime_lso=gfprime_lso-Vk_lso
      call inv(gfprime_lso)
      gfmats_unperiodized(:,:,:,:,:,:,ii)=vca_lso2nnn_reshape(gfprime_lso,Nlat,Nspin,Norb)
    enddo
    !
    do ii=1,Lreal  
      Vk_lso=vca_nnn2lso_reshape(tk(kpoint)-t_prime-set_delta(dcmplx(wr(ii),eps),[ts],[Uloc(1)/2]),Nlat,Nspin,Norb)      
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
        ind1=N2indices(ilat)        
        do jlat=1,Nlat
          ind2=N2indices(jlat)
          gfmats_periodized(:,:,:,:,ii)=gfmats_periodized(:,:,:,:,ii)+exp(-xi*dot_product(kpoint,ind1-ind2))*gfmats_unperiodized(ilat,jlat,:,:,:,:,ii)/Nlat
        enddo
      enddo
    enddo
    !
    do ii=1,Lreal   
      do ilat=1,Nlat
        ind1=N2indices(ilat)        
        do jlat=1,Nlat
          ind2=N2indices(jlat)
          gfreal_periodized(:,:,:,:,ii)=gfreal_periodized(:,:,:,:,ii)+exp(-xi*dot_product(kpoint,ind1-ind2))*gfreal_unperiodized(ilat,jlat,:,:,:,:,ii)/Nlat
        enddo
      enddo
    enddo
    !
    !   
  end subroutine periodize_g_scheme



 subroutine build_sigma_g_scheme(kpoint)
    integer                                                     :: i,ispin,iorb,ii
    real(8),dimension(Ndim)                                     :: kpoint
    complex(8),dimension(Nspin,Nspin,Norb,Norb,Lmats)           :: invG0mats,invGmats
    complex(8),dimension(Nspin,Nspin,Norb,Norb,Lreal)           :: invG0real,invGreal
    !
    !
    if(.not.allocated(Smats_periodized))allocate(Smats_periodized(Nspin,Nspin,Norb,Norb,Lmats))
    if(.not.allocated(Sreal_periodized))allocate(Sreal_periodized(Nspin,Nspin,Norb,Norb,Lreal))
    invG0mats = zero
    invGmats  = zero
    invG0real = zero
    invGreal  = zero
    Smats_periodized  = zero
    Sreal_periodized  = zero
    !
    !Get G0^-1
    do ispin=1,Nspin
      do iorb=1,Norb
        do ii=1,Lmats
            invG0mats(ispin,ispin,iorb,iorb,ii) = (xi*wm(ii)+xmu)  + 2.d0*t_var*(cos(kpoint(1)) + cos(kpoint(2)))             ! FIXME: ad-hoc solution
        enddo
        do ii=1,Lreal
            invG0real(ispin,ispin,iorb,iorb,ii) = (wr(ii)+xmu)   + 2.d0*t_var*(cos(kpoint(1)) + cos(kpoint(2)))               ! FIXME: ad-hoc solution
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
    Smats_periodized=zero
    Sreal_periodized=zero
    do ispin=1,Nspin
      do iorb=1,Norb
         Smats_periodized(ispin,ispin,iorb,iorb,:) = invG0mats(ispin,ispin,iorb,iorb,:) - invGmats(ispin,ispin,iorb,iorb,:)
         Sreal_periodized(ispin,ispin,iorb,iorb,:) = invG0real(ispin,ispin,iorb,iorb,:) - invGreal(ispin,ispin,iorb,iorb,:)
      enddo
    enddo
    !
  end subroutine build_sigma_g_scheme


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
      allocate(kgrid_test(Nkpts**ndim,Ndim)) 
      call TB_build_kgrid([Nkpts,Nkpts],kgrid_test)
      print*,"Calculating Gloc and Sigma ",scheme," scheme"  
      call start_timer
      do ik=1,Nkpts**ndim
          !if (scheme == "g" ) then
            call build_sigma_g_scheme(kgrid_test(ik,:))  !also periodizes g
          !else
            !call build_g_sigma_scheme(kgrid_test(ik,:))  !also periodizes Sigma
          !endif
          do ix=1,Lmats
              gtest_mats(:,:,:,:,ix)=gtest_mats(:,:,:,:,ix)+gfmats_periodized(:,:,:,:,ix)/(Nkpts**Ndim)
              sigmatest_mats(:,:,:,:,ix)=sigmatest_mats(:,:,:,:,ix)+Smats_periodized(:,:,:,:,ix)/(Nkpts**Ndim)
          enddo
          do ix=1,Lreal
              gtest_real(:,:,:,:,ix)=gtest_real(:,:,:,:,ix)+gfreal_periodized(:,:,:,:,ix)/(Nkpts**Ndim)
              sigmatest_real(:,:,:,:,ix)=sigmatest_real(:,:,:,:,ix)+Sreal_periodized(:,:,:,:,ix)/(Nkpts**Ndim)
          enddo
          call eta(ik,Nkpts**ndim)
      enddo
      call stop_timer
      gfmats_periodized=gtest_mats
      Smats_periodized=sigmatest_mats
      gfreal_periodized=gtest_real
      Sreal_periodized=sigmatest_real
      do iorb=1,Norb
       do ispin=1,Nspin
          suffix="_l"//str(iorb)//str(iorb)//"_s"//str(ispin)//"_"//str(scheme)//"_scheme"
          call splot("perG"//reg(suffix)//"_iw.vca"   ,wm,gfmats_periodized(ispin,ispin,iorb,iorb,:))
          call splot("perG"//reg(suffix)//"_realw.vca",wr,gfreal_periodized(ispin,ispin,iorb,iorb,:))
          call splot("perSigma"//reg(suffix)//"_iw.vca"   ,wm,Smats_periodized(ispin,ispin,iorb,iorb,:))
          call splot("perSigma"//reg(suffix)//"_realw.vca",wr,Sreal_periodized(ispin,ispin,iorb,iorb,:))
        enddo
      enddo
    endif
  !
end subroutine get_local_gf


  !---------------------------------------------------------------------
  !PURPOSE: GET A(k,w)
  !---------------------------------------------------------------------
  subroutine get_Akw()
    integer                                                 :: ik=0,iw
    integer                                                 :: Lk,Nkpath=100
    integer                                                 :: iorb,jorb
    integer                                                 :: ispin,jspin
    complex(8),dimension(:,:),allocatable                   :: Akreal
    complex(8),dimension(:,:,:,:,:,:),allocatable           :: Gkreal
    real(8),dimension(:,:),allocatable                      :: Kpath,Kgrid
    character(len=30)                                       :: suffix
  !
    if(master)then
    wr     = linspace(wini,wfin,Lreal)
  !
    print*,"Build A(k,w)"
    !
    allocate(kpath(3,3))
    kpath(1,:)=[0,0,-1]*pi
    kpath(2,:)=[0,0,0]*pi
    kpath(3,:)=[0,0,1]*pi
    !
    Lk=(size(kpath,1)-1)*Nkpath
    Nlso=Nlat*Nspin*Norb
    !
    allocate(Gkreal(Lk,Nspin,Nspin,Norb,Norb,Lreal));Gkreal=zero
    allocate(Akreal(Lk,Lreal));Akreal=zero
    allocate(kgrid(Lk,size(kpath,2)))
    !
    call TB_build_kgrid(kpath,Nkpath,kgrid)
    !
    !
    print*,"Retrieving G(k,w)"
    call start_timer
    do ik=1,Lk
        call periodize_g_scheme([0.d0,kgrid(ik,3)])
        do iw=1,Lreal
            Gkreal(ik,:,:,:,:,iw)=gfreal_periodized(:,:,:,:,iw)
        enddo
        call eta(ik,Lk)
    enddo
    call stop_timer
    !
    print*,"Calculating A(k,w)"
    do ispin=1,Nspin
      do iorb=1,Norb
        Akreal = Akreal - dimag(Gkreal(:,ispin,ispin,iorb,iorb,:))/pi/Nspin/Norb
      enddo
    enddo
    print*,"Printing"
    call splot3d("Akw_real.dat",kgrid(:,3),wr,Akreal) 
    endif
end subroutine get_Akw


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

