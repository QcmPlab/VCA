MODULE VCA_IO
  USE VCA_VARS_GLOBAL
  USE VCA_AUX_FUNX
  !
  ! USE SF_LINALG
  USE SF_ARRAYS, only: linspace,arange
  USE SF_IOTOOLS, only: str,reg,splot
  USE SF_MISC,only : assert_shape
  implicit none
  private



  !Retrieve self-energy through routines:
  interface vca_get_sigma_matsubara
     module procedure vca_get_sigma_matsubara_full
     module procedure vca_get_sigma_matsubara_ij
  end interface vca_get_sigma_matsubara

  interface vca_get_sigma_realaxis
     module procedure vca_get_sigma_real_full
     module procedure vca_get_sigma_real_ij
  end interface vca_get_sigma_realaxis


  !Retrieve imp GF through routines.
  interface vca_get_gimp_matsubara
     module procedure vca_get_gimp_matsubara_full
     module procedure vca_get_gimp_matsubara_ij
  end interface vca_get_gimp_matsubara

  interface vca_get_gimp_realaxis
     module procedure vca_get_gimp_real_full
     module procedure vca_get_gimp_real_ij
  end interface vca_get_gimp_realaxis


  !
  public :: vca_get_sigma_matsubara
  public :: vca_get_sigma_realaxis
  !
  public :: vca_get_gimp_matsubara
  public :: vca_get_gimp_realaxis
  !
  public :: vca_print_impSigma
  public :: vca_print_impG
  public :: vca_print_impG0
  !
  public :: vca_get_dens
  public :: vca_get_mag
  public :: vca_get_docc
  public :: vca_get_sft_potential






contains


  !+------------------------------------------------------------------+
  !PURPOSE  : Print impurity Functions case:
  ! - impSigma
  ! - impG
  ! - impG0
  !+------------------------------------------------------------------+
  subroutine vca_print_impSigma
    integer                                           :: i
    character(len=64)                                 :: suffix
    integer :: ilat,jlat,iorb,jorb,ispin,jspin
    !
    call vca_allocate_time_freq_arrays()
    !Print the impurity Sigma:
    do ilat=1,Nlat
       do jlat=1,Nlat
          do iorb=1,Norb
             do ispin=1,Nspin
                suffix="_Isite"//str(ilat,4)//"_Jsite"//str(jlat,4)//"_l"//str(iorb)//str(iorb)//"_s"//str(ispin)
                call splot("impSigma"//reg(suffix)//"_iw"//reg(file_suffix)//".vca"   ,wm,impSmats(ilat,jlat,ispin,ispin,iorb,iorb,:))
                call splot("impSigma"//reg(suffix)//"_realw"//reg(file_suffix)//".vca",wr,impSreal(ilat,jlat,ispin,ispin,iorb,iorb,:))
             enddo
          enddo
       enddo
    enddo
    call vca_deallocate_time_freq_arrays()
  end subroutine vca_print_impSigma

  subroutine vca_print_impG
    integer                                           :: i
    character(len=64)                                 :: suffix
    integer :: ilat,jlat,iorb,jorb,ispin,jspin
    !
    call vca_allocate_time_freq_arrays()
    !Print the impurity GF:
    do ilat=1,Nlat
       do jlat=1,Nlat
          do iorb=1,Norb
             do ispin=1,Nspin
                suffix="_Isite"//str(ilat,4)//"_Jsite"//str(jlat,4)//"_l"//str(iorb)//str(iorb)//"_s"//str(ispin)
                call splot("impG"//reg(suffix)//"_iw"//reg(file_suffix)//".vca"   ,wm,impGmats(ilat,jlat,ispin,ispin,iorb,iorb,:))
                call splot("impG"//reg(suffix)//"_realw"//reg(file_suffix)//".vca",wr,impGreal(ilat,jlat,ispin,ispin,iorb,iorb,:))
             enddo
          enddo
       enddo
    enddo
    call vca_deallocate_time_freq_arrays()
    !
  end subroutine vca_print_impG

  subroutine vca_print_impG0
    integer                                           :: i
    character(len=64)                                 :: suffix
    integer :: ilat,jlat,iorb,jorb,ispin,jspin
    !
    call vca_allocate_time_freq_arrays()
    !Print the impurity GF:
    do ilat=1,Nlat
       do jlat=1,Nlat
          do iorb=1,Norb
             do ispin=1,Nspin
                suffix="_Isite"//str(ilat,4)//"_Jsite"//str(jlat,4)//"_l"//str(iorb)//str(iorb)//"_s"//str(ispin)
                call splot("impG0"//reg(suffix)//"_iw"//reg(file_suffix)//".vca"   ,wm,impG0mats(ilat,jlat,ispin,ispin,iorb,iorb,:))
                call splot("impG0"//reg(suffix)//"_realw"//reg(file_suffix)//".vca",wr,impG0real(ilat,jlat,ispin,ispin,iorb,iorb,:))
             enddo
          enddo
       enddo
    enddo
    call vca_deallocate_time_freq_arrays()
    !
  end subroutine vca_print_impG0






  !+------------------------------------------------------------------+
  !PURPOSE  : Retrieve impurity Functions:
  ! - impG (Mats,real)
  ! - impSigma  (Mats,real)
  !+------------------------------------------------------------------+
  subroutine vca_get_gimp_matsubara_full(Gmats)
    complex(8),dimension(Nlat,Nlat,Nspin,Nspin,Norb,Norb,Lmats),intent(inout) :: Gmats
    Gmats = impGmats
  end subroutine vca_get_gimp_matsubara_full

  subroutine vca_get_gimp_matsubara_ij(Gmats,ilat,jlat)
    integer                                                         :: ilat,jlat
    complex(8),dimension(Nspin,Nspin,Norb,Norb,Lmats),intent(inout) :: Gmats
    Gmats = impGmats(ilat,jlat,:,:,:,:,:)
  end subroutine vca_get_gimp_matsubara_ij

  subroutine vca_get_sigma_matsubara_full(Smats)
    complex(8),dimension(Nlat,Nlat,Nspin,Nspin,Norb,Norb,Lmats),intent(inout) :: Smats
    Smats = impSmats
  end subroutine vca_get_sigma_matsubara_full

  subroutine vca_get_sigma_matsubara_ij(Smats,ilat,jlat)
    integer                                                         :: ilat,jlat
    complex(8),dimension(Nspin,Nspin,Norb,Norb,Lmats),intent(inout) :: Smats
    Smats = impSmats(ilat,jlat,:,:,:,:,:)
  end subroutine vca_get_sigma_matsubara_ij

  subroutine vca_get_gimp_real_full(Greal)
    complex(8),dimension(Nlat,Nlat,Nspin,Nspin,Norb,Norb,Lreal),intent(inout) :: Greal
    Greal = impGreal
  end subroutine vca_get_gimp_real_full

  subroutine vca_get_gimp_real_ij(Greal,ilat,jlat)
    integer                                                         :: ilat,jlat
    complex(8),dimension(Nspin,Nspin,Norb,Norb,Lreal),intent(inout) :: Greal
    Greal = impGreal(ilat,jlat,:,:,:,:,:)
  end subroutine vca_get_gimp_real_ij

  subroutine vca_get_sigma_real_full(Sreal)
    complex(8),dimension(Nlat,Nlat,Nspin,Nspin,Norb,Norb,Lreal),intent(inout) :: Sreal
    Sreal = impSreal
  end subroutine vca_get_sigma_real_full

  subroutine vca_get_sigma_real_ij(Sreal,ilat,jlat)
    integer                                                         :: ilat,jlat
    complex(8),dimension(Nspin,Nspin,Norb,Norb,Lreal),intent(inout) :: Sreal
    Sreal = impSreal(ilat,jlat,:,:,:,:,:)
  end subroutine vca_get_sigma_real_ij







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



END MODULE VCA_IO





! !+-----------------------------------------------------------------------------+!
! ! PURPOSE: Build the cluster GF from cluster Qmatrix 
! !+-----------------------------------------------------------------------------+!
! subroutine vca_get_Gcluster_matsubara(Gmats)
!   complex(8),dimension(Nlat,Nlat,Norb,Norb,Nspin,Nspin,Lmats),intent(inout) :: Gmats
!   call Qmatrix_to_matsubara_gf(Gmats,Qcluster)
! end subroutine vca_get_Gcluster_matsubara

! subroutine vca_get_Gcluster_realaxis(Greal)
!   complex(8),dimension(Nlat,Nlat,Norb,Norb,Nspin,Nspin,Lreal),intent(inout) :: Greal
!   call Qmatrix_to_realaxis_gf(Greal,Qcluster)
! end subroutine vca_get_Gcluster_realaxis



! !+-----------------------------------------------------------------------------+!
! ! PURPOSE: Build the system GF from system Qmatrix 
! !+-----------------------------------------------------------------------------+!
! subroutine vca_get_Gsystem_matsubara(Gmats)
!   complex(8),dimension(:,:,:,:,:,:,:),intent(inout) :: Gmats
!   integer                                           :: Nsys
!   Nsys = size(Gmats,1)
!   call assert_shape(Gmats,[Nsys,Nsys,Norb,Norb,Nspin,Nspin,Lmats],"vca_get_Gsystem_matsubara","Gmats")
!   call Qmatrix_to_matsubara_gf(Gmats,Qsystem)
! end subroutine vca_get_Gsystem_matsubara

! subroutine vca_get_Gsystem_realaxis(Greal)
!   complex(8),dimension(:,:,:,:,:,:,:),intent(inout) :: Greal
!   integer                                           :: Nsys
!   Nsys = size(Greal,1)
!   call assert_shape(Greal,[Nsys,Nsys,Norb,Norb,Nspin,Nspin,Lmats],"vca_get_Gsystem_realaxis","Greal")
!   call Qmatrix_to_realaxis_gf(Greal,Qsystem)
! end subroutine vca_get_Gsystem_realaxis











! !+-----------------------------------------------------------------------------+!
! ! PURPOSE: Retrieve measured values of the impurity green's functions 
! !+-----------------------------------------------------------------------------+!
! !NORMAL, MATSUBARA GREEN'S FUNCTIONS
! subroutine Qmatrix_to_matsubara_gf(Gmats,Matrix)
!   complex(8),dimension(Nlat,Nlat,Norb,Norb,Nspin,Nspin,Lmats),intent(inout) :: Gmats
!   type(Qmatrix)                                                             :: Matrix
!   integer                                                                   :: ispin
!   integer                                                                   :: ilat,jlat
!   integer                                                                   :: iorb,jorb
!   integer                                                                   :: iexc,Nexc
!   integer                                                                   :: i,is,js
!   real(8)                                                                   :: weight,de
!   !
!   if(.not.Matrix%allocated)stop "Qmatrix_to_matsubara_gf ERROR: Matrix not allocated"
!   !
!   call allocate_grids()
!   !
!   Gmats = dcmplx(0d0,0d0)
!   !
!   Nexc = Matrix%Nexc
!   !
!   do ilat=1,Nlat
!      do jlat=1,Nlat
!         do iorb=1,Norb
!            do jorb=1,Norb
!               do ispin=1,Nspin
!                  !
!                  is = index_stride_los(ilat,iorb,ispin)
!                  js = index_stride_los(jlat,jorb,ispin)
!                  !
!                  do iexc=1,Nexc
!                     weight = Matrix%c(is,iexc)*Matrix%cdg(iexc,js)
!                     de     = Matrix%poles(iexc)
!                     Gmats(ilat,jlat,iorb,jorb,ispin,ispin,:) = Gmats(ilat,jlat,iorb,jorb,ispin,ispin,:) + weight/(xi*wm(:)-de)
!                  enddo
!                  !
!               enddo
!            enddo
!         enddo
!      enddo
!   enddo
!   !
!   call deallocate_grids()
! end subroutine Qmatrix_to_matsubara_gf

! !NORMAL, REALAXIS GREEN'S FUNCTIONS
! subroutine Qmatrix_to_realaxis_gf(Greal,Matrix)
!   complex(8),dimension(Nlat,Nlat,Norb,Norb,Nspin,Nspin,Lreal),intent(inout) :: Greal
!   type(Qmatrix)                                                             :: Matrix
!   integer                                                                   :: ispin
!   integer                                                                   :: ilat,jlat
!   integer                                                                   :: iorb,jorb
!   integer                                                                   :: iexc,Nexc
!   integer                                                                   :: i,is,js
!   real(8)                                                                   :: weight,de
!   !
!   if(.not.Matrix%allocated)stop "Qmatrix_to_realaxis_gf ERROR: Matrix not allocated"        
!   !
!   call allocate_grids()
!   !
!   Greal = dcmplx(0d0,0d0)
!   !
!   Nexc = Matrix%Nexc
!   !
!   do ilat=1,Nlat
!      do jlat=1,Nlat
!         do iorb=1,Norb
!            do jorb=1,Norb
!               do ispin=1,Nspin
!                  !
!                  is = index_stride_los(ilat,iorb,ispin)
!                  js = index_stride_los(jlat,jorb,ispin)
!                  !
!                  do iexc=1,Nexc
!                     weight = Matrix%c(is,iexc)*Matrix%cdg(iexc,js)
!                     de     = Matrix%poles(iexc)
!                     Greal(ilat,jlat,iorb,jorb,ispin,ispin,:) = Greal(ilat,jlat,iorb,jorb,ispin,ispin,:) + weight/(wr(:)+xi*eps-de)
!                  enddo
!                  !
!               enddo
!            enddo
!         enddo
!      enddo
!   enddo
!   !
!   call deallocate_grids()
! end subroutine Qmatrix_to_realaxis_gf
