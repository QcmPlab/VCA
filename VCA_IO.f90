MODULE VCA_IO
  USE VCA_VARS_GLOBAL
  USE VCA_AUX_FUNX
  !
  USE SF_LINALG
  USE SF_ARRAYS, only: linspace,arange
  USE SF_IOTOOLS, only: str,reg,free_unit,splot,sread
  implicit none
  private

  !Retrieve imp GF through routines.
  interface vca_get_gf_matsubara
     module procedure vca_get_gimp_matsubara_1
     module procedure vca_get_gimp_matsubara_2
     module procedure vca_get_gimp_matsubara_3
  end interface vca_get_gf_matsubara

  interface vca_get_gf_realaxis
     module procedure vca_get_gimp_realaxis_1
     module procedure vca_get_gimp_realaxis_2
     module procedure vca_get_gimp_realaxis_3
  end interface vca_get_gf_realaxis


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



  public :: vca_get_gf_matsubara
  public :: vca_get_gf_realaxis
  public :: vca_get_dens
  public :: vca_get_mag
  public :: vca_get_docc
  public :: vca_get_Nexc
  

  !Frequency and time arrays:
  !=========================================================
  real(8),dimension(:),allocatable :: wm,tau,wr,vm
  character(len=64)                :: suffix




contains




  !+-----------------------------------------------------------------------------+!
  ! PURPOSE: Retrieve measured values of the impurity green's functions 
  !+-----------------------------------------------------------------------------+!
  !NORMAL, MATSUBARA GREEN'S FUNCTIONS
  subroutine vca_get_gimp_matsubara_1(Gmats)
    integer                                                                   :: Nexc
    integer                                                                   :: ilat,jlat
    integer                                                                   :: iorb,jorb
    integer                                                                   :: ispin
    integer                                                                   :: iexc
    integer                                                                   :: i
    real(8)                                                                   :: weight,de
    complex(8),dimension(Nlat,Nlat,Nspin,Nspin,Norb,Norb,Lmats),intent(inout) :: Gmats
    !
    call allocate_grids()
    !
    Gmats = dcmplx(0d0,0d0)
    Nexc = size(Lmatrix,7)
    !
    do ispin=1,Nspin
       do ilat=1,Nlat
          do jlat=1,Nlat
             do iorb=1,Norb
                do jorb=1,Norb
                   !
                   do iexc=1,Nexc
                      weight = cdgQmatrix(ilat,ispin,iorb,iexc)*cQmatrix(jlat,ispin,jorb,iexc)
                      de     = Lmatrix(ilat,jlat,ispin,ispin,iorb,jorb,iexc)
                      Gmats(ilat,jlat,ispin,ispin,iorb,jorb,:) = Gmats(ilat,jlat,ispin,ispin,iorb,jorb,:) + weight/(xi*wm(:)-de)
                   enddo
                   !
                enddo
             enddo
          enddo
       enddo
    enddo
    !
    call deallocate_grids()
  end subroutine vca_get_gimp_matsubara_1

  subroutine vca_get_gimp_matsubara_2(Gmats)
    integer                                                                   :: Nexc
    integer                                                                   :: ilat,jlat
    integer                                                                   :: iorb,jorb
    integer                                                                   :: ispin
    integer                                                                   :: iexc
    integer                                                                   :: i,io,jo
    real(8)                                                                   :: weight,de
    complex(8),dimension(Nlat*Nspin*Norb,Nlat*Nspin*Norb,Lmats),intent(inout) :: Gmats
    !
    call allocate_grids()
    !
    Gmats = dcmplx(0d0,0d0)
    Nexc = size(Lmatrix,7)
    !
    do ispin=1,Nspin
       do ilat=1,Nlat
          do jlat=1,Nlat
             do iorb=1,Norb
                do jorb=1,Norb
                   !
                   io = iorb + (ispin-1)*Norb + (ilat-1)*Nspin*Norb
                   jo = jorb + (ispin-1)*Norb + (jlat-1)*Nspin*Norb
                   do iexc=1,Nexc
                      weight = cdgQmatrix(ilat,ispin,iorb,iexc)*cQmatrix(jlat,ispin,jorb,iexc)
                      de     = Lmatrix(ilat,jlat,ispin,ispin,iorb,jorb,iexc)
                      Gmats(io,jo,:) = Gmats(io,jo,:) + weight/(xi*wm(:)-de)
                   enddo
                   !
                enddo
             enddo
          enddo
       enddo
    enddo
  end subroutine vca_get_gimp_matsubara_2

  subroutine vca_get_gimp_matsubara_3(Gmats,ilat,jlat,ispin,jspin,iorb,jorb)
    integer                                   :: Nexc
    integer                                   :: ilat,jlat
    integer                                   :: iorb,jorb
    integer                                   :: ispin,jspin
    integer                                   :: iexc
    integer                                   :: i
    real(8)                                   :: weight,de
    complex(8),dimension(Lmats),intent(inout) :: Gmats
    !
    call allocate_grids()
    !
    Gmats = dcmplx(0d0,0d0)
    Nexc = size(Lmatrix,7)
    !
    do iexc=1,Nexc
       weight = cdgQmatrix(ilat,ispin,iorb,iexc)*cQmatrix(jlat,ispin,jorb,iexc)
       de     = Lmatrix(ilat,jlat,ispin,ispin,iorb,jorb,iexc)
       Gmats(:) = Gmats(:) + weight/(xi*wm(:)-de)
    enddo
    !
  end subroutine vca_get_gimp_matsubara_3




  !NORMAL, REALAXIS GREEN'S FUNCTIONS
  subroutine vca_get_gimp_realaxis_1(Greal)
    integer                                                                   :: Nexc
    integer                                                                   :: ilat,jlat
    integer                                                                   :: iorb,jorb
    integer                                                                   :: ispin
    integer                                                                   :: iexc
    integer                                                                   :: i
    real(8)                                                                   :: weight,de
    complex(8),dimension(Nlat,Nlat,Nspin,Nspin,Norb,Norb,Lmats),intent(inout) :: Greal
    !
    call allocate_grids()
    !
    Greal = dcmplx(0d0,0d0)
    Nexc = size(Lmatrix,7)
    !
    do ispin=1,Nspin
       do ilat=1,Nlat
          do jlat=1,Nlat
             do iorb=1,Norb
                do jorb=1,Norb
                   !
                   do iexc=1,Nexc
                      weight = cdgQmatrix(ilat,ispin,iorb,iexc)*cQmatrix(jlat,ispin,jorb,iexc)
                      de     = Lmatrix(ilat,jlat,ispin,ispin,iorb,jorb,iexc)
                      Greal(ilat,jlat,ispin,ispin,iorb,jorb,:) = Greal(ilat,jlat,ispin,ispin,iorb,jorb,:) + weight/(wr(:)+xi*eps-de)
                   enddo
                   !
                enddo
             enddo
          enddo
       enddo
    enddo
    !
    call deallocate_grids()
  end subroutine vca_get_gimp_realaxis_1

  subroutine vca_get_gimp_realaxis_2(Greal)
    integer                                                                   :: Nexc
    integer                                                                   :: ilat,jlat
    integer                                                                   :: iorb,jorb
    integer                                                                   :: ispin
    integer                                                                   :: iexc
    integer                                                                   :: i,io,jo
    real(8)                                                                   :: weight,de
    complex(8),dimension(Nlat*Nspin*Norb,Nlat*Nspin*Norb,Lmats),intent(inout) :: Greal
    !
    call allocate_grids()
    !
    Greal = dcmplx(0d0,0d0)
    Nexc = size(Lmatrix,7)
    !
    do ispin=1,Nspin
       do ilat=1,Nlat
          do jlat=1,Nlat
             do iorb=1,Norb
                do jorb=1,Norb
                   !
                   io = iorb + (ispin-1)*Norb + (ilat-1)*Nspin*Norb
                   jo = jorb + (ispin-1)*Norb + (jlat-1)*Nspin*Norb
                   do iexc=1,Nexc
                      weight = cdgQmatrix(ilat,ispin,iorb,iexc)*cQmatrix(jlat,ispin,jorb,iexc)
                      de     = Lmatrix(ilat,jlat,ispin,ispin,iorb,jorb,iexc)
                      Greal(io,jo,:) = Greal(io,jo,:) + weight/(xi*wm(:)-de)
                   enddo
                   !
                enddo
             enddo
          enddo
       enddo
    enddo
  end subroutine vca_get_gimp_realaxis_2

  subroutine vca_get_gimp_realaxis_3(Greal,ilat,jlat,ispin,jspin,iorb,jorb)
    integer                                   :: Nexc
    integer                                   :: ilat,jlat
    integer                                   :: iorb,jorb
    integer                                   :: ispin,jspin
    integer                                   :: iexc
    integer                                   :: i
    real(8)                                   :: weight,de
    complex(8),dimension(Lmats),intent(inout) :: Greal
    !
    call allocate_grids()
    !
    Greal = dcmplx(0d0,0d0)
    Nexc = size(Lmatrix,7)
    !
    do iexc=1,Nexc
       weight = cdgQmatrix(ilat,ispin,iorb,iexc)*cQmatrix(jlat,ispin,jorb,iexc)
       de     = Lmatrix(ilat,jlat,ispin,ispin,iorb,jorb,iexc)
       Greal(:) = Greal(:) + weight/(xi*wm(:)-de)
    enddo
    !
  end subroutine vca_get_gimp_realaxis_3





  subroutine vca_get_Nexc(n)
    integer :: n
    N = size(Lmatrix,7)
  end subroutine vca_get_Nexc




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



  !+------------------------------------------------------------------+
  !PURPOSE  : Allocate arrays and setup frequencies and times
  !+------------------------------------------------------------------+
  subroutine allocate_grids
    integer :: i
    if(.not.allocated(wm))allocate(wm(Lmats))
    if(.not.allocated(vm))allocate(vm(0:Lmats))          !bosonic frequencies
    if(.not.allocated(wr))allocate(wr(Lreal))
    if(.not.allocated(tau))allocate(tau(0:Ltau))
    wm     = pi/beta*real(2*arange(1,Lmats)-1,8)
    do i=0,Lmats
       vm(i) = pi/beta*2.d0*dble(i)
    enddo
    wr     = linspace(wini,wfin,Lreal)
    tau(0:)= linspace(0.d0,beta,Ltau+1)
  end subroutine allocate_grids


  subroutine deallocate_grids
    if(allocated(wm))deallocate(wm)
    if(allocated(vm))deallocate(vm)
    if(allocated(tau))deallocate(tau)
    if(allocated(wr))deallocate(wr)
  end subroutine deallocate_grids



END MODULE VCA_IO
