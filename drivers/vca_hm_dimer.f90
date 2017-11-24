program vca_hm_dimer
  USE SCIFOR
  USE DMFT_TOOLS
  !
  USE VCA
  !
  implicit none
  integer                            :: Nsys,Nb,Nsys_bath
  integer                            :: ilat,jlat
  integer                            :: i,j
  integer                            :: ix
  logical                            :: converged
  real(8)                            :: wband

  !The local hybridization function:
  real(8),dimension(:),allocatable   :: Bath,Params
  real(8),dimension(:,:),allocatable :: Hsys,Hcluster,Hsys_bath
  character(len=16)                  :: finput
  real(8)                            :: ts,Vh,Eh,tol
  integer                            :: Nc,Lx
  integer                            :: unit
  real(8),dimension(:),allocatable   :: ts_array,omega_array
  !
  real(8),dimension(1)               :: array
  integer,dimension(1)               :: min_loc
  integer                            :: iter
  real(8)                            :: fresult,omega
  integer                            :: nloop
  logical                            :: wloop,wmin

  call parse_cmd_variable(finput,"FINPUT",default='inputVCA.conf')
  call parse_input_variable(ts,"ts",finput,default=1d0)
  call parse_input_variable(Vh,"Vh",finput,default=0d0)
  call parse_input_variable(Eh,"Eh",finput,default=0d0)
  call parse_input_variable(tol,"TOL",finput,default=1.d-10)
  call parse_input_variable(wloop,"wloop",finput,default=.false.,comment="T: includes loop over Vh")
  call parse_input_variable(wmin,"wmin",finput,default=.true.,comment="T: includes global minimization")
  call parse_input_variable(nloop,"NLOOP",finput,default=100)
  call vca_read_input(trim(finput))

  if(Nlat/=1)stop "Nlat != 1"
  if(Nbath/=1)stop "Nbath != 1"
  if(Norb/=1)stop "Norb != 1"
  if(Nspin/=1)stop "Nspin != 1"


  !Add DMFT CTRL Variables:
  call add_ctrl_var(Nlat,"NLAT")
  call add_ctrl_var(Norb,"norb")
  call add_ctrl_var(Nspin,"nspin")
  call add_ctrl_var(beta,"beta")
  call add_ctrl_var(xmu,"xmu")
  call add_ctrl_var(wini,'wini')
  call add_ctrl_var(wfin,'wfin')
  call add_ctrl_var(eps,"eps")


  !>get main system dimensions:
  Nc   = vca_get_cluster_dimension()
  Nsys = vca_get_system_dimension()


  !>build full system lattice tb Hamiltonian:
  allocate(Hsys(Nsys,Nsys))
  Hsys = Htb_chain(Nsys,ts)

  !>build cluster tb Hamiltonian
  allocate(Hcluster(Nc,Nc))
  Hcluster = Htb_chain(Nc,ts)

  !>allocate bath system:
  Nb = vca_get_bath_dimension()
  allocate(Bath(Nb))

  !> Init VCA solver (fill bath if any)
  call vca_init_solver(bath)

  write(*,"(A)")"Done init solver"

  !>build \tilde{H}_sys = H_sys "\plus" bath
  Nsys_bath = vca_get_system_dimension(with_bath=.true.)
  allocate(Hsys_bath(Nsys_bath,Nsys_bath))
  call vca_build_Hsystem(Hsys,bath,Hsys_bath)
  write(*,"(A)")"done building Hsystem"

  allocate(params(Nb+1))
  params=0d0
  params(1) = Hcluster(1,1)
  params(2:) = bath


  if(wloop)then
     allocate(ts_array(Nloop))
     allocate(omega_array(Nloop))
     ts_array = linspace(0.d0,0.4d0,Nloop)
     do i=1,Nloop
        omega_array(i)=solve_vca_HMdimer(ts_array(i))
     enddo
     call splot("sft_Omega_loopVSts.dat",ts_array,omega_array)
     min_loc = minloc(omega_array)
     write(800,*)min_loc,ts_array(min_loc(1)),omega_array(min_loc(1))
     stop
  endif


  ! if(wmin)then
  !    print*,"Guess:",ts
  !    call  brent(solve_vca1d,ts,[0d0,1d0])
  !    print*,"Result ts : ",ts
  !    stop
  ! endif

  ! allocate(omega_array(1))
  ! omega_array(1)=solve_vca1d(ts)
  ! open(10,file="sft_Omega_loopVSts.dat")
  ! write(10,*)ts,omega_array(1)
  ! close(10)

contains



  function solve_vca_HMdimer(array) result(Omega)
    real(8),dimension(1)                   :: array
    real(8)                                :: Omega
    real(8),dimension(Nlat,Nlat)           :: Hcluster
    real(8),dimension(Nsys_bath,Nsys_bath) :: Htile_bath
    real(8),dimension(Nsys_bath,Nsys_bath) :: Vmat    
    !
    print*,"------ INIT ------"
    !> bcast array to Hcluster and Bath: this is strongly problem dependent
    ! Hcluster(1,1) = array(1)
    Bath          = [Eh,array(1)]

    call vca_tile_Hcluster(Hcluster,Bath,Htile_bath)
    !
    Vmat = Hsys_bath - Htile_bath
    !
    call vca_solve(vca_los2nnn_reshape(Hcluster,Nlat,Norb,Nspin),Vmat,bath) 
    call vca_get_sft_potential(omega)
    print*,""
    print*,"------ DONE ------"
    !
  end function solve_vca_HMdimer



  function Htb_chain(N,ts) result(H0)
    integer                  :: N
    real(8)                  :: ts
    character(len=64)        :: file_
    real(8),dimension(N,N)   :: H0
    integer                  :: i,j,k,ilink,link(2),ii,jj
    integer                  :: unit
    !
    file_ = "Tsys_matrix.dat"
    !
    H0 = 0d0
    !
    do i=1,N
       !right hop
       link(1)= i + 1
       if(i==N)link(1)=0
       !left  hop
       link(2)= i - 1
       if(i==1)link(2)=0
       !
       do ilink=1,2
          j = link(ilink)
          if(j<=0)cycle
          H0(i,j) = -ts
       enddo
       !
    enddo
    open(free_unit(unit),file=trim(file_))
    do i=1,N
       write(unit,"(1000000(F5.2,1x))")(H0(i,j),j=1,N)
    enddo
    close(unit)
    !
  end function Htb_chain




















  ! subroutine print_chain_bath_structure(Hmat,file)
  !   integer                   :: Nvec(2)
  !   integer                   :: Nx
  !   integer                   :: Ny
  !   real(8),dimension(:,:)    :: Hmat
  !   character(len=*),optional :: file
  !   character(len=32)         :: file_
  !   integer                   :: unit
  !   integer                   :: Nh
  !   integer                   :: i,j
  !   integer                   :: ix,iy
  !   integer                   :: jx,jy
  !   integer                   :: ilat,jlat
  !   integer                   :: ispin,jspin
  !   integer                   :: iorb,jorb
  !   real(8)                   :: Ri(2),Rj(2),x,y
  !   real(8)                   :: spin_shift(2),orb_shift(Norb)
  !   real(8),parameter         :: xshift=1d-1,dorb=4d-2,pwidth=4d-2
  !   logical                   :: bool
  !   type(rgb_color),parameter :: corb(5)=[red1,blue1,green3,magenta1,black]
  !   character(len=10)         :: chpoint
  !   integer                   :: Xmin,Xmax,Ymin,Ymax
  !   !
  !   file_="lattice"; if(present(file))file_=file
  !   !
  !   Nx   = Ncopies
  !   Ny   = 2
  !   Nlat = Nx*Ny
  !   !
  !   call assert_shape(Hmat,[Nlat,Nsys],"proint_2DLattice_structure","Hmat")
  !   !
  !   !Setup some shifts:
  !   !SPIN:
  !   spin_shift=0d0
  !   if(Nspin>1)then
  !      spin_shift(1)= xshift
  !      spin_shift(2)=-xshift
  !   endif
  !   !ORBITAL:
  !   orb_shift=0d0
  !   if(Norb>1)then
  !      Nh = Norb/2
  !      if(mod(Norb,2)==0)then
  !         do iorb=1,Nh
  !            orb_shift(Nh+iorb)   = iorb*dorb-dorb/2d0
  !            orb_shift(Nh-iorb+1) =-iorb*dorb+dorb/2d0
  !         enddo
  !      else          
  !         do iorb=1,Norb
  !            orb_shift(iorb)      = -Nh*dorb + (iorb-1)*dorb
  !         enddo
  !      endif
  !   endif
  !   !
  !   !PRINT POINTS:
  !   open(free_unit(unit),file=trim(file)//"_points.dat")
  !   do ix=1,Nx
  !      do iy=1,Ny
  !         do ispin=1,Nspin
  !            do iorb=1,Norb
  !               Ri(1) = ix + spin_shift(ispin) 
  !               Ri(2) = iy + orb_shift(iorb)
  !               write(unit,"(3F12.5,I12)")Ri(1),Ri(2),pwidth,rgb(corb(mod(iorb,5)))
  !            enddo
  !         enddo
  !      enddo
  !   enddo
  !   close(unit)
  !   !
  !   !PRINT LINES:
  !   open(free_unit(unit),file=trim(file)//"_lines.dat")
  !   do ix=1,Nx
  !      do iy=1,Ny
  !         ilat = iy + (ix-1)*Ny
  !         !
  !         do jx=ix,Nx
  !            do jy=iy,Ny
  !               jlat = jy + (jx-1)*Ny
  !               !
  !               bool=.false.
  !               do ispin=1,Nspin
  !                  do jspin=1,Nspin
  !                     do iorb=1,Norb
  !                        do jorb=1,Norb
  !                           i = iorb + (ispin-1)*Nspin + (ilat-1)*Nspin*Norb
  !                           j = jorb + (jspin-1)*Nspin + (jlat-1)*Nspin*Norb
  !                           if(abs(Hmat(i,j))/=0d0)bool=.true.
  !                        enddo
  !                     enddo
  !                  enddo
  !               enddo
  !               if(bool)then
  !                  Ri = [ix,iy]
  !                  Rj = [jx,jy]
  !                  write(unit,"(2F12.5)")Ri(1),Ri(2)
  !                  write(unit,"(2F12.5)")Rj(1),Rj(2)
  !                  write(unit,*)""
  !               endif
  !               !
  !            enddo
  !         enddo
  !      enddo
  !   enddo
  !   !
  !   if(.true.)then
  !      ix = 1
  !      iy = 1
  !      write(unit,*)ix-0.25d0,iy
  !      write(unit,*)ix,iy
  !      write(unit,*)""
  !      ix = Nx
  !      write(unit,*)ix+0.25d0,iy
  !      write(unit,*)ix,iy
  !      write(unit,*)""
  !   endif
  !   !
  !   close(unit)
  !   !
  !   Xmin = 0    ; Ymin  = 0
  !   Xmax = Nx+1 ; Ymax  = Ny+1
  !   !
  !   open(free_unit(unit),file=trim(file)//"_structure.gp")
  !   write(unit,*)"set term wxt"
  !   write(unit,*)"#set terminal pngcairo size 350,262 enhanced font 'Verdana,10'"
  !   write(unit,*)"#set out '"//trim(file)//"_structure.png'"
  !   write(unit,*)""
  !   write(unit,*)"#set terminal svg size 350,262 fname 'Verdana, Helvetica, Arial, sans-serif'"
  !   write(unit,*)"#set out '"//trim(file)//"_structure.svg'"
  !   write(unit,*)""
  !   write(unit,*)"#set term postscript eps enhanced color 'Times'"
  !   write(unit,*)"#set output '|ps2pdf - "//trim(file)//"_structure.pdf'"
  !   write(unit,*)""
  !   do iorb=1,Norb
  !      chpoint=str(0.95d0-(iorb-1)*0.05d0)
  !      write(unit,"(A)")str("set label 'Orb "//str(iorb)//"' tc rgb "//str(rgb(corb(iorb)))//&
  !           " at graph 0.9,"//reg(chpoint)//" font 'Times-Italic,11'")
  !   enddo
  !   write(unit,*)"unset key"
  !   write(unit,*)"set bars small"
  !   write(unit,*)"plot ["//str(Xmin)//":"//str(Xmax)//"]["//str(Ymin)//":"//str(Ymax)//"] "&
  !        //"'"//trim(file)//"_lines.dat' u 1:2 w l lw 1 lc rgb 'grey1',"&
  !        //"'"//trim(file)//"_points.dat' u 1:2:3:4 with xerrorbars lc rgb variable ps 0 lw 2 "
  !   close(unit)
  !   call system("chmod +x "//trim(file)//"_structure.gp")
  ! end subroutine print_chain_bath_structure

end program vca_hm_dimer



