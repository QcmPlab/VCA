program vca_test
  USE SCIFOR
  USE DMFT_TOOLS
  USE MPI
  USE VCA
  !
  implicit none
  integer                                         :: Nlso,Nsys
  integer                                         :: ilat,jlat
  integer                                         :: i,j
  integer                                         :: ix,iy
  logical                                         :: converged
  real(8)                                         :: wband

  !The local hybridization function:
  real(8),dimension(:,:),allocatable              :: Tsys,Tref,Vmat,Htb,Mmat,dens
  real(8),allocatable,dimension(:)                :: wm,wr
  complex(8)                                      :: iw
  complex(8),allocatable,dimension(:,:,:,:,:,:,:) :: Gmats,Greal
  character(len=16)                               :: finput
  real(8)                                         :: ts
  integer                                         :: Nx,Ny,Lx,Ly,Rx,Ry
  integer                                         :: unit
  integer                                         :: comm,rank
  logical                                         :: master


  call init_MPI()
  comm = MPI_COMM_WORLD
  call StartMsg_MPI(comm)
  rank = get_Rank_MPI(comm)
  master = get_Master_MPI(comm)

  call parse_cmd_variable(finput,"FINPUT",default='inputVCA.conf')
  call parse_input_variable(ts,"ts",finput,default=1d0)
  call parse_input_variable(Rx,"Rx",finput,default=1,comment="Ratio L/Lc=Rx along X-directions, aka # of copies along X")
  call parse_input_variable(Ry,"Ry",finput,default=1,comment="Ratio L/Lc=Ry along Y-directions, aka # of copies along Y")
  !
  call vca_read_input(trim(finput),comm)


  if(Norb/=1)stop "Norb != 1"
  if(Nspin/=1)stop "Nspin != 1"
  !
  Nx=4
  Ny=4
  Lx   = Rx*Nx
  Ly   = Ry*Ny
  !
  Nlat = Nx*Ny
  Nsys = Lx*Ly
  !
  Nlso = Nlat*Norb*Nspin

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

  !>build full system lattice tb hamiltonian
  Tsys = Htb_square_lattice(Lx,Ly,ts,file="Tsys_matrix.dat")

  !>build cluster tight binding hamiltonian (one could extract it)
  Htb  = Htb_square_lattice(Nx,Ny,ts,file="Htb_matrix.dat")

  !>build Tref: tiling of Htb decoupled clusters;
  call vca_tile_Treference(Htb,[Rx,Ry],Tref)

  !>build Vmat=Tsys-Tref
  Vmat = Tsys - Tref

  !>printing lattice structures.
  call print_2DLattice_Structure(Htb,[Nx,Ny],1,1,file="Htb")
  call print_2DLattice_Structure(Tsys,[Lx,Ly],1,1,file="Tsys")
  call print_2DLattice_Structure(Tref,[Lx,Ly],1,1,file="Tref")
  call print_2DLattice_Structure(Vmat,[Lx,Ly],1,1,file="Vmat")  

  call vca_init_solver(comm)

  call vca_solve(comm,one*vca_lso2nnn_reshape(Htb,Nlat,Nspin,Norb))


  ! allocate(wm(Lmats),wr(Lreal))
  ! wm = pi/beta*(2*arange(1,Lmats)-1)
  ! wr = linspace(wini,wfin,Lreal)

  ! allocate(Gmats(Nlat,Nlat,Nspin,Nspin,Norb,Norb,Lmats))
  ! allocate(Greal(Nlat,Nlat,Nspin,Nspin,Norb,Norb,Lmats))
  ! allocate(dens(Nlat,Nlat))

  ! call vca_get_gimp_matsubara(Gmats)
  ! call vca_get_gimp_realaxis(Greal)
  ! do ilat=1,Nlat
  !    do jlat=1,Nlat
  !       call splot("Gimp_i"//str(ilat,3)//"_j"//str(jlat,3)//"_l11_s1_iw.vca",wm,Gmats(ilat,jlat,1,1,1,1,:))
  !       call splot("Gimp_i"//str(ilat,3)//"_j"//str(jlat,3)//"_l11_s1_realw.vca",wr,Greal(ilat,jlat,1,1,1,1,:))
  !       dens(ilat,jlat) = fft_get_density(Gmats(ilat,jlat,1,1,1,1,:),beta)
  !    enddo
  ! enddo
  ! deallocate(Gmats,Greal,dens)




  call finalize_MPI()





contains




  function Htb_square_lattice(Nrow,Ncol,ts,file) result(H0)
    integer                                :: Nrow
    integer                                :: Ncol
    real(8)                                :: ts
    character(len=*),optional              :: file
    character(len=64)                      :: file_
    real(8),dimension(Nrow*Ncol,Nrow*Ncol) :: H0
    integer                                :: i,jj,row,col,link(4),j
    integer                                :: unit
    !
    file_ = "Htb_square_lattice.dat";if(present(file))file_=file
    !
    H0 = 0.d0
    !
    !
    do row=0,Nrow-1
       do col=0,Ncol-1
          i=col+ 1 + row*Ncol
          !
          !
          !right hop
          link(1)= i + 1     
          if((col+1)==Ncol) link(1)=0  
          !left  hop
          link(3)= i - 1    
          if((col-1)<0)link(3)=0  
          !up    hop
          link(2)= i + Ncol 
          if((row+1)==Nrow)link(2)=0  
          !down  hop
          link(4)= i - Ncol 
          if((row-1)<0)link(4)=0  
          !
          do jj=1,4
             if(link(jj)>0)H0(i,link(jj))=-ts !! ts must be negative.
          enddo
          !
       enddo
    enddo
    open(free_unit(unit),file=trim(file_))
    do i=1,Nrow*Ncol
       write(unit,"(5000(F5.2,1x))")(H0(i,j),j=1,Nrow*Ncol)
    enddo
    close(unit)
  end function Htb_square_lattice





  subroutine vca_tile_Treference(Tcluster,Rvec,Tref,iprint)
    real(8),dimension(:,:) :: Tcluster
    integer,dimension(:)   :: Rvec
    real(8),dimension(:,:) :: Tref
    logical,optional       :: iprint
    integer                :: Nlat,Nratio
    integer                :: i,j,icopy,unit
    !
    Nlat   = size(Tcluster,1)
    Nratio = product(Rvec)
    Nsys   = Nratio*Nlat
    call assert_shape(Tcluster,[Nlat,Nlat],"vca_tile_Treference","Tcluster")
    call assert_shape(Tref,[Nsys,Nsys],"vca_tile_Treference","Tref")
    do icopy=1,Nratio
       do ix=1,Nlat
          do iy=1,Nlat
             i = ix + (icopy-1)*Nlat
             j = iy + (icopy-1)*Nlat
             Tref(i,j) = Tcluster(ix,iy)
          enddo
       enddo
    enddo
    !
    if(present(iprint))then
       if(iprint)then
          open(free_unit(unit),file="Tref_matrix.dat")
          do i=1,Nsys
             write(unit,"(5000(F5.2,1x))")(Tref(i,j),j=1,Nsys)
          enddo
          close(unit)
       endif
    endif
  end subroutine vca_tile_Treference



  ! subroutine get_non_interacting_GF()
  !   complex(8),dimension(Nlso,Nlso) :: zeta
  !   real(8) :: wm(Lmats)
  !   real(8) :: wr(Lreal)
  !   !
  !   !Get GMats(iw)
  !   wm = pi/beta*(2*arange(1,Lmats)-1)
  !   do i=1,Lmats
  !      zeta = (xi*wm(i)+xmu)*eye(Nlso) - Htb
  !      call inv(zeta)
  !      Gmats(:,:,:,:,:,:,i) = lso2nnn_reshape(zeta,Nlat,Nspin,Norb)
  !   enddo
  !   do ilat=1,Nlat
  !      do jlat=1,Nlat
  !         call splot("Gmats_i"//str(ilat,3)//"_j"//str(jlat,3)//"_l11_s1_iw.nint",wm,Gmats(ilat,jlat,1,1,1,1,:))
  !      enddo
  !   enddo
  !   !
  !   wr = linspace(wini,wfin,Lreal)
  !   do i=1,Lreal
  !      zeta = (wr(i)+xi*eps + xmu)*eye(Nlso) - Htb
  !      call inv(zeta)
  !      Greal(:,:,:,:,:,:,i) = lso2nnn_reshape(zeta,Nlat,Nspin,Norb)
  !   enddo
  !   do ilat=1,Nlat
  !      do jlat=1,Nlat
  !         call splot("Greal_i"//str(ilat,3)//"_j"//str(jlat,3)//"_l11_s1_realw.nint",wr,Greal(ilat,jlat,1,1,1,1,:))
  !      enddo
  !   enddo
  !   !
  ! end subroutine get_non_interacting_GF







  subroutine print_2Dlattice_structure(Hmat,Nvec,Nspin,Norb,pbc,file)
    integer                   :: Nvec(2)
    integer                   :: Nx
    integer                   :: Ny
    integer                   :: Nspin
    integer                   :: Norb
    real(8),dimension(:,:)    :: Hmat
    logical,optional          :: pbc(2)
    character(len=*),optional :: file
    character(len=32)         :: file_
    logical                   :: pbc_(2)
    integer                   :: unit
    integer                   :: Nlat,Nh
    integer                   :: i,j
    integer                   :: ix,iy
    integer                   :: jx,jy
    integer                   :: ilat,jlat
    integer                   :: ispin,jspin
    integer                   :: iorb,jorb
    real(8)                   :: Ri(2),Rj(2),x,y
    real(8)                   :: spin_shift(2),orb_shift(Norb)
    real(8),parameter         :: xshift=1d-1,dorb=4d-2,pwidth=4d-2
    logical                   :: bool
    type(rgb_color),parameter :: corb(5)=[red1,blue1,green3,magenta1,black]
    character(len=10)         :: chpoint
    integer                   :: Xmin,Xmax,Ymin,Ymax
    !
    pbc_=.false.   ; if(present(pbc))pbc_=pbc
    file_="lattice"; if(present(file))file_=file
    !
    Nx = Nvec(1)
    Ny = Nvec(2)
    Nlat = Nx*Ny
    !
    call assert_shape(Hmat,[Nlat*Nspin*Norb,Nlat*Nspin*Norb],"proint_2DLattice_structure","Hmat")
    !
    !Setup some shifts:
    !SPIN:
    spin_shift=0d0
    if(Nspin>1)then
       spin_shift(1)= xshift
       spin_shift(2)=-xshift
    endif
    !ORBITAL:
    orb_shift=0d0
    if(Norb>1)then
       Nh = Norb/2
       if(mod(Norb,2)==0)then
          do iorb=1,Nh
             orb_shift(Nh+iorb)   = iorb*dorb-dorb/2d0
             orb_shift(Nh-iorb+1) =-iorb*dorb+dorb/2d0
          enddo
       else          
          do iorb=1,Norb
             orb_shift(iorb)      = -Nh*dorb + (iorb-1)*dorb
          enddo
       endif
    endif
    !
    !PRINT POINTS:
    open(free_unit(unit),file=trim(file)//"_points.dat")
    do ix=1,Nx
       do iy=1,Ny
          do ispin=1,Nspin
             do iorb=1,Norb
                Ri(1) = ix + spin_shift(ispin) 
                Ri(2) = iy + orb_shift(iorb)
                write(unit,"(3F12.5,I12)")Ri(1),Ri(2),pwidth,rgb(corb(mod(iorb,5)))
             enddo
          enddo
       enddo
    enddo
    close(unit)
    !
    !PRINT LINES:
    open(free_unit(unit),file=trim(file)//"_lines.dat")
    do ix=1,Nx
       do iy=1,Ny
          ilat = iy + (ix-1)*Ny
          !
          do jx=ix,Nx
             do jy=iy,Ny
                jlat = jy + (jx-1)*Ny
                !
                bool=.false.
                do ispin=1,Nspin
                   do jspin=1,Nspin
                      do iorb=1,Norb
                         do jorb=1,Norb
                            i = iorb + (ispin-1)*Nspin + (ilat-1)*Nspin*Norb
                            j = jorb + (jspin-1)*Nspin + (jlat-1)*Nspin*Norb
                            if(abs(Hmat(i,j))/=0d0)bool=.true.
                         enddo
                      enddo
                   enddo
                enddo
                if(bool)then
                   Ri = [ix,iy]
                   Rj = [jx,jy]
                   write(unit,"(2F12.5)")Ri(1),Ri(2)
                   write(unit,"(2F12.5)")Rj(1),Rj(2)
                   write(unit,*)""
                endif
                !
             enddo
          enddo
       enddo
    enddo
    !
    if(pbc_(1))then
       ix = 1
       do iy=1,Ny
          write(unit,*)ix-0.25d0,iy
          write(unit,*)ix,iy
          write(unit,*)""
       enddo
       ix = Nx
       do iy=1,Ny
          write(unit,*)ix+0.25d0,iy
          write(unit,*)ix,iy
          write(unit,*)""
       enddo
    endif
    !
    if(pbc_(2))then
       iy = 1
       do ix=1,Nx
          write(unit,*)ix,iy-0.25d0
          write(unit,*)ix,iy
          write(unit,*)""
       enddo
       iy = Ny
       do ix=1,Nx
          write(unit,*)ix,iy+0.25d0
          write(unit,*)ix,iy
          write(unit,*)""
       enddo
    endif
    close(unit)
    !
    Xmin = 0    ; Ymin  = 0
    Xmax = Nx+1 ; Ymax  = Ny+1
    !
    open(free_unit(unit),file=trim(file)//"_structure.gp")
    write(unit,*)"set term wxt"
    write(unit,*)"#set terminal pngcairo size 350,262 enhanced font 'Verdana,10'"
    write(unit,*)"#set out '"//trim(file)//"_structure.png'"
    write(unit,*)""
    write(unit,*)"#set terminal svg size 350,262 fname 'Verdana, Helvetica, Arial, sans-serif'"
    write(unit,*)"#set out '"//trim(file)//"_structure.svg'"
    write(unit,*)""
    write(unit,*)"#set term postscript eps enhanced color 'Times'"
    write(unit,*)"#set output '|ps2pdf - "//trim(file)//"_structure.pdf'"
    write(unit,*)""
    do iorb=1,Norb
       chpoint=str(0.95d0-(iorb-1)*0.05d0)
       write(unit,"(A)")str("set label 'Orb "//str(iorb)//"' tc rgb "//str(rgb(corb(iorb)))//&
            " at graph 0.9,"//reg(chpoint)//" font 'Times-Italic,11'")
    enddo
    write(unit,*)"unset key"
    write(unit,*)"set bars small"
    write(unit,*)"plot ["//str(Xmin)//":"//str(Xmax)//"]["//str(Ymin)//":"//str(Ymax)//"] "&
         //"'"//trim(file)//"_lines.dat' u 1:2 w l lw 1 lc rgb 'grey1',"&
         //"'"//trim(file)//"_points.dat' u 1:2:3:4 with xerrorbars lc rgb variable ps 0 lw 2 "
    close(unit)
    call system("chmod +x "//trim(file)//"_structure.gp")
  end subroutine print_2Dlattice_structure


end program vca_test



