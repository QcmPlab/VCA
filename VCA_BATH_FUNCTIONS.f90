MODULE VCA_BATH_FUNCTIONS
  USE SF_CONSTANTS, only: zero
  USE SF_IOTOOLS, only:free_unit,reg,file_length,txtfy
  USE SF_LINALG, only: eye,inv,zeye,inv_her
  USE VCA_INPUT_VARS
  USE VCA_VARS_GLOBAL
  USE VCA_BATH_SETUP
  USE VCA_AUX_FUNX
  implicit none

  private

  interface delta_bath_freq
    module procedure delta_bath_freq_main
  endinterface delta_bath_freq


  !##################################################################
  !
  !\DELTA HYBRIDIZATION FUNCTION MATSUBARA
  !
  !##################################################################
  !NORMAL
  interface delta_bath_mats
     module procedure delta_bath_mats_main
     module procedure delta_bath_mats_main_
  end interface delta_bath_mats


  !##################################################################
  !
  !\DELTA HYBRIDIZATION FUNCTION REAL
  !
  !##################################################################
  !NORMAL
  interface delta_bath_real
     module procedure delta_bath_real_main
     module procedure delta_bath_real_main_
  end interface delta_bath_real
  !



  !##################################################################
  !
  !NON-INTERACTING GREEN'S FUNCTION MATSUBARA
  !
  !##################################################################
  !NORMAL
  interface g0and_bath_mats
     module procedure g0and_bath_mats_main
     module procedure g0and_bath_mats_main_
  end interface g0and_bath_mats
  !



  !##################################################################
  !
  !INVERSE NON-INTERACTING GREEN'S FUNCTION MATSUBARA
  !
  !##################################################################
  !NORMAL
  interface invg0_bath_mats
     module procedure invg0_bath_mats_main
     module procedure invg0_bath_mats_main_
  end interface invg0_bath_mats



  !##################################################################
  !
  !NON-INTERACTING GREEN'S FUNCTION REAL-AXIS
  !
  !##################################################################
  !NORMAL
  interface g0and_bath_real
     module procedure g0and_bath_real_main
     module procedure g0and_bath_real_main_
  end interface g0and_bath_real
  !



  !##################################################################
  !
  !INVERSE NON-INTERACTING GREEN'S FUNCTION REAL-AXIS
  !
  !##################################################################
  !NORMAL
  interface invg0_bath_real
     module procedure invg0_bath_real_main
     module procedure invg0_bath_real_main_
  end interface invg0_bath_real
  !
  public :: delta_bath_freq

  public :: delta_bath_mats
  public :: g0and_bath_mats
  public :: invg0_bath_mats
  !
  public :: delta_bath_real
  public :: g0and_bath_real
  public :: invg0_bath_real


  integer :: ilat,jlat,iorb,jorb,ispin,jspin


contains




  !##################################################################
  !
  !     DELTA FUNCTIONS
  !     G0 FUNCTIONS
  !     G0^{-1} FUNCTIONS
  !     at a point x from type(effective_bath) :: vca_bath
  !
  !##################################################################
  !+-----------------------------------------------------------------------------+!
  !PURPOSE:  Delta functions for specific frequency:
  !+-----------------------------------------------------------------------------+!
  !NORMAL:
  function delta_bath_freq_main(x,vca_bath_) result(Delta)
    complex(8),intent(in)                                                             :: x
    type(effective_bath)                                                              :: vca_bath_
    complex(8),dimension(Nlat,Nlat,Nspin,Nspin,Norb,Norb)                             :: Delta
    integer                                                                           :: i,ih,L
    integer                                                                           :: io,jo
    complex(8),dimension(Nlat_bath*Nspin*Norb_bath,Nlat_bath*Nspin*Norb_bath)         :: hbath
    complex(8),dimension(Nlat*Nspin*Norb,Nlat_bath*Nspin*Norb_bath)                   :: vbath_ver
    complex(8),dimension(Nlat_bath*Nspin*Norb_bath,Nlat*Nspin*Norb)                   :: vbath_hor
    !
    Delta=zero
    !
    !
    if(.not.vca_bath_%status)return
    !
    hbath = (x+XMU)*zeye(Nlat_bath*Nspin*Norb_bath) - vca_nnn2lso_reshape(vca_bath_%h,Nlat_bath,Nspin,Norb_bath)
    vbath_ver = vca_rectangular_n2j_reshape(vca_bath_%v,Nlat,Nlat_bath,Nspin,Nspin,Norb,Norb_bath)
    vbath_hor = conjg(transpose(vbath_ver))
    !
    call inv(hbath)
    !
    Delta = vca_lso2nnn_reshape(matmul(vbath_ver,matmul(hbath,vbath_hor)),Nlat,Nspin,Norb)
    !
  end function delta_bath_freq_main

  !+-----------------------------------------------------------------------------+!
  !PURPOSE:  Delta functions on the Matsubara axis:
  !+-----------------------------------------------------------------------------+!
  !NORMAL:
  function delta_bath_mats_main(x,vca_bath_) result(Delta)
    complex(8),dimension(:),intent(in)                                                 :: x
    type(effective_bath)                                                               :: vca_bath_
    complex(8),dimension(Nlat,Nlat,Nspin,Nspin,Norb,Norb,size(x))                      :: Delta
    integer                                                                            :: i,ih,L
    integer                                                                            :: io,jo
    complex(8),dimension(Nlat_bath*Nspin*Norb_bath,Nlat_bath*Nspin*Norb_bath)          :: hbath
    complex(8),dimension(Nlat*Nspin*Norb,Nlat_bath*Nspin*Norb_bath)                    :: vbath_ver
    complex(8),dimension(Nlat_bath*Nspin*Norb_bath,Nlat*Nspin*Norb)                    :: vbath_hor
    !
    Delta=zero
    !
    L = size(x)
    !
    if(.not.vca_bath_%status)return
    vbath_ver = vca_rectangular_n2j_reshape(vca_bath_%v,Nlat,Nlat_bath,Nspin,Nspin,Norb,Norb_bath)
    vbath_hor = conjg(transpose(vbath_ver))
    !
    do i=1,L
      hbath = (x(i)+XMU)*zeye(Nlat_bath*Nspin*Norb_bath) - vca_nnn2lso_reshape(vca_bath_%h,Nlat_bath,Nspin,Norb_bath)
      call inv(hbath)
      Delta(:,:,:,:,:,:,i) = vca_lso2nnn_reshape(matmul(vbath_ver,matmul(hbath,vbath_hor)),Nlat,Nspin,Norb)
    enddo
  end function delta_bath_mats_main



  !+-----------------------------------------------------------------------------+!
  !PURPOSE:  Delta functions on the Real axis:
  !+-----------------------------------------------------------------------------+!
  function delta_bath_real_main(x,vca_bath_) result(Delta)
    complex(8),dimension(:),intent(in)                                                 :: x
    type(effective_bath)                                                               :: vca_bath_
    complex(8),dimension(Nlat,Nlat,Nspin,Nspin,Norb,Norb,size(x))                      :: Delta
    integer                                                                            :: i,ih,L
    integer                                                                            :: io,jo
    complex(8),dimension(Nlat_bath*Nspin*Norb_bath,Nlat_bath*Nspin*Norb_bath)          :: hbath
    complex(8),dimension(Nlat*Nspin*Norb,Nlat_bath*Nspin*Norb_bath)                    :: vbath_ver
    complex(8),dimension(Nlat_bath*Nspin*Norb_bath,Nlat*Nspin*Norb)                    :: vbath_hor
    !
    Delta=zero
    !
    L = size(x)
    !
    if(.not.vca_bath_%status)return
    vbath_ver = vca_rectangular_n2j_reshape(vca_bath_%v,Nlat,Nlat_bath,Nspin,Nspin,Norb,Norb_bath)
    vbath_hor = conjg(transpose(vbath_ver))
    !
    do i=1,L
      hbath = (x(i)+XMU)*zeye(Nlat_bath*Nspin*Norb_bath) - vca_nnn2lso_reshape(vca_bath_%h,Nlat_bath,Nspin,Norb_bath)
      call inv(hbath)
      Delta(:,:,:,:,:,:,i) = vca_lso2nnn_reshape(matmul(vbath_ver,matmul(hbath,vbath_hor)),Nlat,Nspin,Norb)
    enddo
    !
  end function delta_bath_real_main












  !+-----------------------------------------------------------------------------+!
  !PURPOSE:  G0 non-interacting Green's functions on the Matsubara axis
  !+-----------------------------------------------------------------------------+!
  !NORMAL:
  function g0and_bath_mats_main(x,vca_bath_) result(G0and)
    complex(8),dimension(:),intent(in)                            :: x
    type(effective_bath)                                          :: vca_bath_
    complex(8),dimension(Nlat,Nlat,Nspin,Nspin,Norb,Norb,size(x)) :: G0and,Delta
    complex(8),dimension(Nlat*Nspin*Norb,Nlat*Nspin*Norb)         :: tmp
    integer                                                       :: iorb,jorb,ispin,jspin,io,jo,Nso,i,L
    real(8),dimension(size(x))                                    :: det
    complex(8),dimension(size(x))                                 :: fg,ff
    complex(8),dimension(:,:),allocatable                         :: fgorb,zeta
    !
    G0and = zero
    !
    L=size(x)
    !
    !
    Delta = delta_bath_mats(x,vca_bath_)
    do i=1,L
      tmp=(x(i) + xmu)*Zeye(Nlat*Nspin*Norb) - vca_nnn2lso_reshape(impHloc - Delta(:,:,:,:,:,:,i),Nlat,Nspin,Norb)
      call inv(tmp)
      G0and(:,:,:,:,:,:,i) = vca_lso2nnn_reshape(tmp,Nlat,Nspin,Norb)
    enddo
    !
  end function g0and_bath_mats_main





  !+-----------------------------------------------------------------------------+!
  !PURPOSE:  G0 non-interacting Green's functions on the real-axis
  !+-----------------------------------------------------------------------------+!
  !NORMAL:
  function g0and_bath_real_main(x,vca_bath_) result(G0and)
    complex(8),dimension(:),intent(in)                            :: x
    type(effective_bath)                                          :: vca_bath_
    complex(8),dimension(Nlat,Nlat,Nspin,Nspin,Norb,Norb,size(x)) :: G0and,Delta
    complex(8),dimension(Nlat*Nspin*Norb,Nlat*Nspin*Norb)         :: tmp
    integer                                                       :: iorb,jorb,ispin,jspin,io,jo,Nso,i,L
    complex(8),dimension(size(x))                                 :: det,fg,ff
    complex(8),dimension(:,:),allocatable                         :: fgorb,zeta
    !
    G0and = zero
    !
    L = size(x)
    !
    Delta = delta_bath_real(x,vca_bath_)
    do i=1,L
      tmp=(x(i) + xmu)*Zeye(Nlat*Nspin*Norb) - vca_nnn2lso_reshape(impHloc - Delta(:,:,:,:,:,:,i),Nlat,Nspin,Norb)
      call inv(tmp)
      G0and(:,:,:,:,:,:,i) = vca_lso2nnn_reshape(tmp,Nlat,Nspin,Norb)
    enddo
    !
  end function g0and_bath_real_main




  !+-----------------------------------------------------------------------------+!
  !PURPOSE:  G0^{-1} non-interacting Green's functions on the Matsubara axis:
  !+-----------------------------------------------------------------------------+!
  !NORMAL:
  function invg0_bath_mats_main(x,vca_bath_) result(G0and)
    complex(8),dimension(:),intent(in)                            :: x
    type(effective_bath)                                          :: vca_bath_
    complex(8),dimension(Nlat,Nlat,Nspin,Nspin,Norb,Norb,size(x)) :: G0and,Delta
    integer                                                       :: i,iorb,jorb,ispin,jspin,io,jo,Nso,L
    complex(8),dimension(:,:),allocatable                         :: zeta!,fgorb
    !
    G0and = zero
    !
    L=size(x)
    !
    Delta = delta_bath_mats(x,vca_bath_)
    do i=1,L
      G0and(:,:,:,:,:,:,i) =  vca_lso2nnn_reshape((x(i) + xmu)*Zeye(Nlat*Nspin*Norb) - vca_nnn2lso_reshape(impHloc - Delta(:,:,:,:,:,:,i),Nlat,Nspin,Norb),Nlat,Nspin,Norb)
    enddo
  end function invg0_bath_mats_main


  !+-----------------------------------------------------------------------------+!
  !PURPOSE:  G0 non-interacting Green's functions on the real-axis:
  !+-----------------------------------------------------------------------------+!
  function invg0_bath_real_main(x,vca_bath_) result(G0and)
    complex(8),dimension(:),intent(in)                            :: x
    type(effective_bath)                                          :: vca_bath_
    complex(8),dimension(Nlat,Nlat,Nspin,Nspin,Norb,Norb,size(x)) :: G0and,Delta
    integer                                                       :: iorb,jorb,ispin,jspin,io,jo,Nso,i,L
    complex(8),dimension(:,:),allocatable                         :: zeta!,fgorb
    !
    G0and = zero
    !
    L = size(x)
    !
    Delta = delta_bath_real(x,vca_bath_)
    do i=1,L
      G0and(:,:,:,:,:,:,i) =  vca_lso2nnn_reshape((x(i) + xmu)*Zeye(Nlat*Nspin*Norb) - vca_nnn2lso_reshape(impHloc - Delta(:,:,:,:,:,:,i),Nlat,Nspin,Norb),Nlat,Nspin,Norb)
    enddo
    !    !
  end function invg0_bath_real_main









  !##################################################################
  !
  !     DELTA FUNCTIONS
  !     G0 FUNCTIONS
  !     G0^{-1} FUNCTIONS
  !     at a point x from real(8),dimension(:) :: bath_array
  !
  !##################################################################
  function delta_bath_mats_main_(x,h_in,v_in) result(Delta)
    complex(8),dimension(:),intent(in)                            :: x
    type(effective_bath)                                          :: vca_bath_
    complex(8),dimension(Nlat,Nlat,Nspin,Nspin,Norb,Norb,size(x)) :: Delta
    complex(8),dimension(:,:,:,:,:,:)   :: h_in, v_in
    logical                                                       :: check
    call vca_allocate_bath(vca_bath_)
    call vca_set_bath(h_in,v_in,vca_bath_)
    Delta = delta_bath_mats_main(x,vca_bath_)
    call vca_deallocate_bath(vca_bath_)
  end function delta_bath_mats_main_

  function delta_bath_real_main_(x,h_in,v_in) result(Delta)
    complex(8),dimension(:),intent(in)                            :: x
    type(effective_bath)                                          :: vca_bath_
    complex(8),dimension(Nlat,Nlat,Nspin,Nspin,Norb,Norb,size(x)) :: Delta
    complex(8),dimension(:,:,:,:,:,:)   :: h_in, v_in
    logical                                                       :: check
    call vca_allocate_bath(vca_bath_)
    call vca_set_bath(h_in,v_in,vca_bath_)
    Delta = delta_bath_real_main(x,vca_bath_)
    call vca_deallocate_bath(vca_bath_)
  end function delta_bath_real_main_

  function g0and_bath_mats_main_(x,h_in,v_in) result(G0and)
    complex(8),dimension(:),intent(in)                  :: x
    type(effective_bath)                                :: vca_bath_
    complex(8),dimension(Nlat,Nlat,Nspin,Nspin,Norb,Norb,size(x)) :: G0and
    complex(8),dimension(:,:,:,:,:,:)   :: h_in, v_in

    call vca_allocate_bath(vca_bath_)
    call vca_set_bath(h_in,v_in,vca_bath_)
    G0and = g0and_bath_mats_main(x,vca_bath_)
    call vca_deallocate_bath(vca_bath_)
  end function g0and_bath_mats_main_

  function g0and_bath_real_main_(x,h_in,v_in) result(G0and)
    complex(8),dimension(:),intent(in)                  :: x
    type(effective_bath)                                :: vca_bath_
    complex(8),dimension(Nlat,Nlat,Nspin,Nspin,Norb,Norb,size(x)) :: G0and
    complex(8),dimension(:,:,:,:,:,:)                   :: h_in, v_in

    call vca_allocate_bath(vca_bath_)
    call vca_set_bath(h_in,v_in,vca_bath_)
    G0and = g0and_bath_real_main(x,vca_bath_)
    call vca_deallocate_bath(vca_bath_)
  end function g0and_bath_real_main_

  function invg0_bath_mats_main_(x,h_in,v_in) result(G0and)
    complex(8),dimension(:),intent(in)                  :: x
    type(effective_bath)                                :: vca_bath_
    complex(8),dimension(Nlat,Nlat,Nspin,Nspin,Norb,Norb,size(x)) :: G0and
    complex(8),dimension(:,:,:,:,:,:)                   :: h_in, v_in

    call vca_allocate_bath(vca_bath_)
    call vca_set_bath(h_in,v_in,vca_bath_)
    G0and = invg0_bath_mats_main(x,vca_bath_)
    call vca_deallocate_bath(vca_bath_)
  end function invg0_bath_mats_main_

  function invg0_bath_real_main_(x,h_in,v_in) result(G0and)
    complex(8),dimension(:),intent(in)                  :: x
    type(effective_bath)                                :: vca_bath_
    complex(8),dimension(Nlat,Nlat,Nspin,Nspin,Norb,Norb,size(x)) :: G0and
   complex(8),dimension(:,:,:,:,:,:)                    :: h_in, v_in

    call vca_allocate_bath(vca_bath_)
    call vca_set_bath(h_in,v_in,vca_bath_)
    G0and = invg0_bath_real_main(x,vca_bath_)
    call vca_deallocate_bath(vca_bath_)
  end function invg0_bath_real_main_








































END MODULE VCA_BATH_FUNCTIONS
