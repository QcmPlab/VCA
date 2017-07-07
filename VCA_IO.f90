MODULE VCA_IO
  USE VCA_VARS_GLOBAL
  USE VCA_AUX_FUNX
  USE SF_LINALG
  USE SF_ARRAYS, only: linspace,arange
  USE SF_IOTOOLS, only: str,reg,free_unit,splot,sread
  implicit none
  private

  !Retrieve imp GF through routines.
  interface vca_get_gimp_matsubara
     module procedure vca_get_gimp_matsubara_1
     module procedure vca_get_gimp_matsubara_2
     module procedure vca_get_gimp_matsubara_3
  end interface vca_get_gimp_matsubara

  interface vca_get_gimp_realaxis
     module procedure vca_get_gimp_realaxis_1
     module procedure vca_get_gimp_realaxis_2
     module procedure vca_get_gimp_realaxis_3
  end interface vca_get_gimp_realaxis


  !Retrieve static common observables  
  interface vca_get_dens
     module procedure vca_get_dens_1
     module procedure vca_get_dens_2
  end interface vca_get_dens

  interface vca_get_mag
     module procedure vca_get_mag_1
     module procedure vca_get_mag_2
  end interface vca_get_mag

  interface vca_get_docc
     module procedure vca_get_docc_1
     module procedure vca_get_docc_2
  end interface vca_get_docc



  public :: vca_get_gimp_matsubara
  public :: vca_get_gimp_realaxis
  public :: vca_get_dens
  public :: vca_get_mag
  public :: vca_get_docc
  public :: vca_print_impG
  public :: vca_read_impG


  !Frequency and time arrays:
  !=========================================================
  real(8),dimension(:),allocatable :: wm,tau,wr,vm
  character(len=64)                :: suffix




contains


  !+------------------------------------------------------------------+
  !PURPOSE  : Print/Read Cluster Greens Functions
  !+------------------------------------------------------------------+
  subroutine vca_print_impG
    integer                                           :: ilat,jlat
    integer                                           :: iorb,jorb
    integer                                           :: ispin
    character(len=20)                                 :: suffix
    !
    if(.not.allocated(wm))allocate(wm(Lmats))
    if(.not.allocated(wr))allocate(wr(Lreal))
    wm     = pi/beta*real(2*arange(1,Lmats)-1,8)
    wr     = linspace(wini,wfin,Lreal)
    !
    !Print the impurity functions:
    do ispin=1,Nspin
       do ilat=1,Nlat
          do jlat=1,Nlat
             do iorb=1,Norb
                do jorb=1,Norb
                   suffix="_i"//str(ilat,3)//"_j"//str(jlat,3)//"_lm"//str(iorb)//str(jorb)//"_s"//str(ispin)
                   call splot("impG"//reg(suffix)//"_iw"//reg(file_suffix)//".ed"   ,wm,impGmats(ilat,jlat,ispin,ispin,iorb,jorb,:))
                   call splot("impG"//reg(suffix)//"_realw"//reg(file_suffix)//".ed",wr,impGreal(ilat,jlat,ispin,ispin,iorb,jorb,:))
                enddo
             enddo
          enddo
       enddo
    enddo
    !
    if(allocated(wm))deallocate(wm)
    if(allocated(wr))deallocate(wr)
    !
  end subroutine vca_print_impG

  subroutine vca_read_impG
    integer                                           :: ilat,jlat
    integer                                           :: iorb,jorb
    integer                                           :: ispin
    character(len=20)                                 :: suffix
    !
    if(.not.allocated(wm))allocate(wm(Lmats))
    if(.not.allocated(wr))allocate(wr(Lreal))
    wm     = pi/beta*real(2*arange(1,Lmats)-1,8)
    wr     = linspace(wini,wfin,Lreal)
    !
    !Print the impurity functions:
    do ispin=1,Nspin
       do ilat=1,Nlat
          do jlat=1,Nlat
             do iorb=1,Norb
                do jorb=1,Norb
                   suffix="_i"//str(ilat,3)//"_j"//str(jlat,3)//"_lm"//str(iorb)//str(jorb)//"_s"//str(ispin)
                   call sread("impG"//reg(suffix)//"_iw"//reg(file_suffix)//".ed"   ,wm,impGmats(ilat,jlat,ispin,ispin,iorb,jorb,:))
                   call sread("impG"//reg(suffix)//"_realw"//reg(file_suffix)//".ed",wr,impGreal(ilat,jlat,ispin,ispin,iorb,jorb,:))
                enddo
             enddo
          enddo
       enddo
    enddo
    !
    if(allocated(wm))deallocate(wm)
    if(allocated(wr))deallocate(wr)
    !
  end subroutine vca_read_impG


  !+-----------------------------------------------------------------------------+!
  ! PURPOSE: Retrieve measured values of the impurity green's functions 
  !+-----------------------------------------------------------------------------+!
  !NORMAL, MATSUBARA GREEN'S FUNCTIONS
  subroutine vca_get_gimp_matsubara_1(Gmats)
    complex(8),dimension(Nlat,Nlat,Nspin,Nspin,Norb,Norb,Lmats),intent(inout) :: Gmats
    Gmats = impGmats
  end subroutine vca_get_gimp_matsubara_1

  subroutine vca_get_gimp_matsubara_2(Gmats)
    complex(8),dimension(Nlat*Nspin*Norb,Nlat*Nspin*Norb,Lmats),intent(inout) :: Gmats
    integer  :: io,jo,iorb,jorb,ispin,jspin,ilat,jlat
    do ilat=1,Nlat
       do jlat=1,Nlat
          do ispin=1,Nspin
             do jspin=1,Nspin
                do iorb=1,Norb
                   do jorb=1,Norb
                      io = iorb + (ispin-1)*Norb + (ilat-1)*Nspin*Norb
                      jo = jorb + (jspin-1)*Norb + (jlat-1)*Nspin*Norb
                      Gmats(io,jo,:) = impGmats(ilat,jlat,ispin,jspin,iorb,jorb,:)
                   enddo
                enddo
             enddo
          enddo
       enddo
    enddo
  end subroutine vca_get_gimp_matsubara_2

  subroutine vca_get_gimp_matsubara_3(Gmats,ilat,jlat,ispin,jspin,iorb,jorb)
    complex(8),dimension(Lmats),intent(inout) :: Gmats
    integer                                   :: ilat,jlat,iorb,jorb,ispin,jspin
    Gmats(:) = impGmats(ilat,jlat,ispin,jspin,iorb,jorb,:)
  end subroutine vca_get_gimp_matsubara_3



  !NORMAL, REALAXIS GREEN'S FUNCTIONS
  subroutine vca_get_gimp_realaxis_1(Greal)
    complex(8),dimension(Nlat,Nlat,Nspin,Nspin,Norb,Norb,Lreal),intent(inout) :: Greal
    Greal = impGreal
  end subroutine vca_get_gimp_realaxis_1

  subroutine vca_get_gimp_realaxis_2(Greal)
    complex(8),dimension(Nlat*Nspin*Norb,Nlat*Nspin*Norb,Lreal),intent(inout) :: Greal
    integer  :: io,jo,iorb,jorb,ispin,jspin,ilat,jlat
    do ilat=1,Nlat
       do jlat=1,Nlat
          do ispin=1,Nspin
             do jspin=1,Nspin
                do iorb=1,Norb
                   do jorb=1,Norb
                      io = iorb + (ispin-1)*Norb + (ilat-1)*Nspin*Norb
                      jo = jorb + (jspin-1)*Norb + (jlat-1)*Nspin*Norb
                      Greal(io,jo,:) = impGreal(ilat,jlat,ispin,jspin,iorb,jorb,:)
                   enddo
                enddo
             enddo
          enddo
       enddo
    enddo
  end subroutine vca_get_gimp_realaxis_2

  subroutine vca_get_gimp_realaxis_3(Greal,ilat,jlat,ispin,jspin,iorb,jorb)
    complex(8),dimension(Lreal),intent(inout) :: Greal
    integer                                   :: ilat,jlat,iorb,jorb,ispin,jspin
    Greal(:) = impGreal(ilat,jlat,ispin,jspin,iorb,jorb,:)
  end subroutine vca_get_gimp_realaxis_3






  !+-----------------------------------------------------------------------------+!
  ! PURPOSE: Retrieve measured values of the local observables
  !+-----------------------------------------------------------------------------+!
  subroutine vca_get_dens_1(dens)
    real(8),dimension(Nlat,Norb) :: dens
    dens = imp_dens
  end subroutine vca_get_dens_1
  !
  subroutine vca_get_mag_1(mag)
    real(8),dimension(Nlat,Norb) :: mag
    mag = imp_dens_up - imp_dens_dw
  end subroutine vca_get_mag_1
  !
  subroutine vca_get_docc_1(docc)
    real(8),dimension(Nlat,Norb) :: docc
    docc = imp_docc
  end subroutine vca_get_docc_1

  subroutine vca_get_dens_2(dens,iorb)
    real(8),dimension(Nlat) :: dens
    integer                 :: iorb
    if(iorb>Norb)stop "imp_get_dens error: orbital index > N_orbital"
    dens = imp_dens(:,iorb)
  end subroutine vca_get_dens_2
  !
  subroutine vca_get_mag_2(mag,iorb)
    real(8),dimension(Nlat) :: mag
    integer                 :: iorb
    if(iorb>Norb)stop "imp_get_mag error: orbital index > N_orbital"
    mag = imp_dens_up(:,iorb) - imp_dens_dw(:,iorb)
  end subroutine vca_get_mag_2
  !
  subroutine vca_get_docc_2(docc,iorb)
    real(8),dimension(Nlat) :: docc
    integer                 :: iorb
    if(iorb>Norb)stop "imp_get_docc error: orbital index > N_orbital"
    docc = imp_docc(:,iorb)
  end subroutine vca_get_docc_2

END MODULE VCA_IO
