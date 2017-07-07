program vca_test
  USE SCIFOR
  USE DMFT_TOOLS
  !
  USE VCA
  !
  implicit none
  integer                                         :: iloop,Nlso,iw
  integer                                         :: ilat,jlat
  logical                                         :: converged
  real(8)                                         :: wband

  !The local hybridization function:
  real(8),allocatable,dimension(:)                :: wm,wr
  real(8),allocatable                             :: Hloc(:,:,:,:,:,:)
  complex(8),dimension(:,:),allocatable           :: zeta
  complex(8),allocatable,dimension(:,:,:,:,:,:,:) :: Gmats,Greal
  character(len=16)                               :: finput
  real(8)                                         :: ts
  real(8),dimension(:,:),allocatable              :: Htb
  integer :: Nx,Ny

  call parse_cmd_variable(finput,"FINPUT",default='inputVCA.conf')
  call parse_input_variable(ts,"ts",finput,default=1d0)
  call parse_input_variable(Nx,"NX",finput,default=2)
  call parse_input_variable(Ny,"Ny",finput,default=2)
  !
  call vca_read_input(trim(finput))

  !Add DMFT CTRL Variables:
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
  Nlat = Nx*Ny
  ! if(Nlat /= 4)stop "Nlat != 4"
  Nlso = Nlat*Nspin*Norb


  allocate(Gmats(Nlat,Nlat,Nspin,Nspin,Norb,Norb,Lmats))
  allocate(Greal(Nlat,Nlat,Nspin,Nspin,Norb,Norb,Lmats))
  allocate(wm(Lmats),wr(Lreal))
  allocate(zeta(Nlso,Nlso))


  !
  allocate(Hloc(Nlat,Nlat,Nspin,Nspin,Norb,Norb))
  Hloc=zero


  allocate(Htb(Nx,Ny))
  Htb = Htb_square_lattice(Nx,Ny,ts)
  ! Htb = 0d0
  ! Htb(1,:) = [0d0, -ts, 0d0, -ts]
  ! Htb(2,:) = [-ts, 0d0, -ts, 0d0]
  ! Htb(3,:) = [0d0, -ts, 0d0, -ts]
  ! Htb(4,:) = [-ts, 0d0, -ts, 0d0]

  Hloc = lso2nnn_reshape(Htb,Nlat,Nspin,Norb)


  !Get GMats(iw)
  wm = pi/beta*(2*arange(1,Lmats)-1)
  do iw=1,Lmats
     zeta = (xi*wm(iw)+xmu)*eye(Nlso) - Htb
     call inv(zeta)
     Gmats(:,:,:,:,:,:,iw) = lso2nnn_reshape(zeta,Nlat,Nspin,Norb)
  enddo
  do ilat=1,Nlat
     do jlat=1,Nlat
        call splot("Gmats_i"//str(ilat,3)//"_j"//str(jlat,3)//"_l11_s1_iw.vca",wm,Gmats(ilat,jlat,1,1,1,1,:))
     enddo
  enddo

  wr = linspace(wini,wfin,Lreal)
  do iw=1,Lreal
     zeta = (wr(iw)+xi*eps + xmu)*eye(Nlso) - Htb
     call inv(zeta)
     Greal(:,:,:,:,:,:,iw) = lso2nnn_reshape(zeta,Nlat,Nspin,Norb)
  enddo
  do ilat=1,Nlat
     do jlat=1,Nlat
        call splot("Greal_i"//str(ilat,3)//"_j"//str(jlat,3)//"_l11_s1_realw.vca",wr,Greal(ilat,jlat,1,1,1,1,:))
     enddo
  enddo





  call vca_init_solver(Hloc)


  call vca_diag(Hloc) 






contains



  function Htb_square_lattice(Nrow,Ncol,ts) result(H0)
    integer                                :: Nrow
    integer                                :: Ncol
    real(8)                        :: ts
    real(8),dimension(Nrow*Ncol,Nrow*Ncol) :: H0
    integer                                :: i,jj,row,col,link(4),j
    integer                                :: unit
    !
    !
    H0 = 0.d0
    unit=free_unit()
    !
    !+- 2D LATTICE (NROW x NCOL) -+!
    if(Nlat /= Nrow*Ncol) stop "get_lattice_hamiltonian error: Nlat != Nrow*Ncol"
    !THESE ARE STILL GLOBAL VARIABLES...
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
    open(unit,file='Htb_square_lattice.ed')
    do i=1,Nrow*Ncol
       write(unit,"(5000(F5.2,1x))")(H0(i,j),j=1,Nrow*Ncol)
    enddo
    close(unit)
  end function Htb_square_lattice


end program vca_test



