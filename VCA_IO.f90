MODULE VCA_IO
  USE VCA_VARS_GLOBAL
  USE VCA_AUX_FUNX
  !
  ! USE SF_LINALG
  USE SF_ARRAYS, only: linspace,arange
  USE SF_IOTOOLS, only: str,reg
  USE SF_MISC,only : assert_shape
  implicit none
  private


  public :: vca_get_Gcluster_matsubara
  public :: vca_get_Gcluster_realaxis
  !
  public :: vca_get_Gsystem_matsubara
  public :: vca_get_Gsystem_realaxis
  !
  public :: vca_get_dens
  public :: vca_get_mag
  public :: vca_get_docc
  public :: vca_get_sft_potential



  !Frequency and time arrays:
  !=========================================================
  real(8),dimension(:),allocatable :: wm,tau,wr,vm
  character(len=64)                :: suffix




contains




  !+-----------------------------------------------------------------------------+!
  ! PURPOSE: Build the cluster GF from cluster Qmatrix 
  !+-----------------------------------------------------------------------------+!
  subroutine vca_get_Gcluster_matsubara(Gmats)
    complex(8),dimension(Nlat,Nlat,Norb,Norb,Nspin,Nspin,Lmats),intent(inout) :: Gmats
    call Qmatrix_to_matsubara_gf(Gmats,Qcluster)
  end subroutine vca_get_Gcluster_matsubara

  subroutine vca_get_Gcluster_realaxis(Greal)
    complex(8),dimension(Nlat,Nlat,Norb,Norb,Nspin,Nspin,Lreal),intent(inout) :: Greal
    call Qmatrix_to_realaxis_gf(Greal,Qcluster)
  end subroutine vca_get_Gcluster_realaxis



  !+-----------------------------------------------------------------------------+!
  ! PURPOSE: Build the system GF from system Qmatrix 
  !+-----------------------------------------------------------------------------+!
  subroutine vca_get_Gsystem_matsubara(Gmats)
    complex(8),dimension(:,:,:,:,:,:,:),intent(inout) :: Gmats
    integer                                           :: Nsys
    Nsys = size(Gmats,1)
    call assert_shape(Gmats,[Nsys,Nsys,Norb,Norb,Nspin,Nspin,Lmats],"vca_get_Gsystem_matsubara","Gmats")
    call Qmatrix_to_matsubara_gf(Gmats,Qsystem)
  end subroutine vca_get_Gsystem_matsubara

  subroutine vca_get_Gsystem_realaxis(Greal)
    complex(8),dimension(:,:,:,:,:,:,:),intent(inout) :: Greal
    integer                                           :: Nsys
    Nsys = size(Greal,1)
    call assert_shape(Greal,[Nsys,Nsys,Norb,Norb,Nspin,Nspin,Lmats],"vca_get_Gsystem_realaxis","Greal")
    call Qmatrix_to_realaxis_gf(Greal,Qsystem)
  end subroutine vca_get_Gsystem_realaxis











  !+-----------------------------------------------------------------------------+!
  ! PURPOSE: Retrieve measured values of the impurity green's functions 
  !+-----------------------------------------------------------------------------+!
  !NORMAL, MATSUBARA GREEN'S FUNCTIONS
  subroutine Qmatrix_to_matsubara_gf(Gmats,Matrix)
    complex(8),dimension(Nlat,Nlat,Norb,Norb,Nspin,Nspin,Lmats),intent(inout) :: Gmats
    type(Qmatrix)                                                             :: Matrix
    integer                                                                   :: ispin
    integer                                                                   :: ilat,jlat
    integer                                                                   :: iorb,jorb
    integer                                                                   :: iexc,Nexc
    integer                                                                   :: i,is,js
    real(8)                                                                   :: weight,de
    !
    if(.not.Matrix%allocated)stop "Qmatrix_to_matsubara_gf ERROR: Matrix not allocated"
    !
    call allocate_grids()
    !
    Gmats = dcmplx(0d0,0d0)
    !
    Nexc = Matrix%Nexc
    !
    do ilat=1,Nlat
       do jlat=1,Nlat
          do iorb=1,Norb
             do jorb=1,Norb
                do ispin=1,Nspin
                   !
                   is = index_stride_los(ilat,iorb,ispin)
                   js = index_stride_los(jlat,jorb,ispin)
                   !
                   do iexc=1,Nexc
                      weight = Matrix%c(is,iexc)*Matrix%cdg(iexc,js)
                      de     = Matrix%poles(iexc)
                      Gmats(ilat,jlat,iorb,jorb,ispin,ispin,:) = Gmats(ilat,jlat,iorb,jorb,ispin,ispin,:) + weight/(xi*wm(:)-de)
                   enddo
                   !
                enddo
             enddo
          enddo
       enddo
    enddo
    !
    call deallocate_grids()
  end subroutine Qmatrix_to_matsubara_gf

  !NORMAL, REALAXIS GREEN'S FUNCTIONS
  subroutine Qmatrix_to_realaxis_gf(Greal,Matrix)
    complex(8),dimension(Nlat,Nlat,Norb,Norb,Nspin,Nspin,Lreal),intent(inout) :: Greal
    type(Qmatrix)                                                             :: Matrix
    integer                                                                   :: ispin
    integer                                                                   :: ilat,jlat
    integer                                                                   :: iorb,jorb
    integer                                                                   :: iexc,Nexc
    integer                                                                   :: i,is,js
    real(8)                                                                   :: weight,de
    !
    if(.not.Matrix%allocated)stop "Qmatrix_to_realaxis_gf ERROR: Matrix not allocated"        
    !
    call allocate_grids()
    !
    Greal = dcmplx(0d0,0d0)
    !
    Nexc = Matrix%Nexc
    !
    do ilat=1,Nlat
       do jlat=1,Nlat
          do iorb=1,Norb
             do jorb=1,Norb
                do ispin=1,Nspin
                   !
                   is = index_stride_los(ilat,iorb,ispin)
                   js = index_stride_los(jlat,jorb,ispin)
                   !
                   do iexc=1,Nexc
                      weight = Matrix%c(is,iexc)*Matrix%cdg(iexc,js)
                      de     = Matrix%poles(iexc)
                      Greal(ilat,jlat,iorb,jorb,ispin,ispin,:) = Greal(ilat,jlat,iorb,jorb,ispin,ispin,:) + weight/(wr(:)+xi*eps-de)
                   enddo
                   !
                enddo
             enddo
          enddo
       enddo
    enddo
    !
    call deallocate_grids()
  end subroutine Qmatrix_to_realaxis_gf






  !+-----------------------------------------------------------------------------+!
  ! PURPOSE: Retrieve measured values of the local observables
  !+-----------------------------------------------------------------------------+!
  subroutine vca_get_dens(dens)
    real(8),dimension(Nlat,Norb) :: dens
    dens = imp_dens
  end subroutine vca_get_dens
  !
  subroutine vca_get_mag(mag)
    real(8),dimension(Nlat,Norb) :: mag
    mag = imp_dens_up - imp_dens_dw
  end subroutine vca_get_mag
  !
  subroutine vca_get_docc(docc)
    real(8),dimension(Nlat,Norb) :: docc
    docc = imp_docc
  end subroutine vca_get_docc
  !
  subroutine vca_get_sft_potential(potential)
    real(8)  :: potential
    potential = sft_potential
  end subroutine vca_get_sft_potential


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
