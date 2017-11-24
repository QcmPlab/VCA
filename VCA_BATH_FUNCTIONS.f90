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



  public :: delta_bath_mats
  public :: delta_bath_real


  integer :: ilat,jlat,iorb,jorb,ispin,jspin


contains




  !##################################################################
  !
  !     DELTA FUNCTIONS
  !     G0 FUNCTIONS
  !     G0^{-1} FUNCTIONS
  !
  !##################################################################

  !+-------------------------------------------------------------------+
  !PURPOSE  : compute the hybridization function at a point x from
  ! type(effective_bath) :: vca_bath
  ! OR
  ! real(8),dimension(:) :: bath_array
  !+-------------------------------------------------------------------+
  !+-----------------------------------------------------------------------------+!
  !PURPOSE:  Delta and Fdelta functions on the Matsubara axis:
  ! _1 : input type(effective_bath) vca_bath
  ! _2 : input array bath
  ! Delta_ : normal
  !+-----------------------------------------------------------------------------+!
  !NORMAL:
  function delta_bath_mats_main(x,vca_bath_) result(Delta)
    complex(8),dimension(:),intent(in)                            :: x
    type(effective_bath)                                          :: vca_bath_
    complex(8),dimension(Nlat,Nlat,Nspin,Nspin,Norb,Norb,size(x)) :: Delta
    integer                                                       :: i,ih,L
    integer                                                       :: ibath
    integer                                                       :: io,jo
    real(8),dimension(Nbath)                                      :: eps,vps
    real(8),dimension(Nlat,Norb,Nbath)                            :: vops
    !
    Delta=zero
    !
    L = size(x)
    !
    select case(bath_type)
    case default                !normal: only _{aa} are allowed (no inter-orbital local mixing)
       !\Delta_{aa} = \sum_k [ V_{a}(k) * V_{a}(k)/(iw_n - E_{a}(k)) ]
       do ilat=1,Nlat
          do iorb=1,Norb
             do ispin=1,Nspin
                eps = vca_bath_%e(ilat,iorb,ispin,1:Nbath)
                vps = vca_bath_%v(ilat,iorb,ispin,1:Nbath)
                do i=1,L
                   Delta(ilat,ilat,ispin,ispin,iorb,iorb,i) = sum( vps(:)*vps(:)/(x(i) - eps(:)) )
                enddo
             enddo
          enddo
       enddo
       !
       !
    case ("hybrid")             !hybrid: all _{ab} components allowed (inter-orbital local mixing present)
       !\Delta_{ab} = \sum_k [ V_{a}(k) * V_{b}(k)/(iw_n - E(k)) ]
       do ispin=1,Nspin
          eps  = vca_bath_%e(1     ,ispin,1     ,1:Nbath)
          vops = vca_bath_%v(1:Nlat,ispin,1:Norb,1:Nbath)
          do ilat=1,Nlat
             do jlat=1,Nlat
                do iorb=1,Norb
                   do jorb=1,Norb
                      do i=1,L
                         Delta(ilat,jlat,iorb,jorb,ispin,ispin,i) = sum( vops(ilat,iorb,:)*vops(jlat,jorb,:)/(x(i) - eps(:)) )
                      enddo
                   enddo
                enddo
             enddo
          enddo
       enddo
       !
    end select
  end function delta_bath_mats_main


  function delta_bath_mats_main_(x,bath_) result(Delta)
    complex(8),dimension(:),intent(in)                            :: x
    type(effective_bath)                                          :: vca_bath_
    complex(8),dimension(Nlat,Nlat,Nspin,Nspin,Norb,Norb,size(x)) :: Delta
    real(8),dimension(:)                                          :: bath_
    logical                                                       :: check
    check= check_bath_dimension(bath_)
    if(.not.check)stop "delta_bath_mats_main_ error: wrong bath dimensions"
    call vca_allocate_bath(vca_bath_)
    call vca_set_bath(bath_,vca_bath_)
    Delta = delta_bath_mats_main(x,vca_bath_)
    call vca_deallocate_bath(vca_bath_)
  end function delta_bath_mats_main_












  !+-----------------------------------------------------------------------------+!
  !PURPOSE:  Delta and Fdelta functions on the Real axis:
  ! _1 : input type(effective_bath) vca_bath
  ! _2 : input array bath
  ! Delta_ : normal
  ! Fdelta_: anomalous
  !+-----------------------------------------------------------------------------+!
  function delta_bath_real_main(x,vca_bath_) result(Delta)
    complex(8),dimension(:),intent(in)                            :: x
    type(effective_bath)                                          :: vca_bath_
    complex(8),dimension(Nlat,Nlat,Nspin,Nspin,Norb,Norb,size(x)) :: Delta
    integer                                                       :: i,ih,L
    integer                                                       :: ibath
    integer                                                       :: io,jo
    real(8),dimension(Nbath)                                      :: eps,vps
    real(8),dimension(Nlat,Norb,Nbath)                            :: vops
    !
    !
    Delta=zero
    !
    L = size(x)
    !
    select case(bath_type)
    case default
       !\Delta_{aa} = \sum_k [ V_{a}(k) * V_{a}(k)/(w+i\h - E_{a}(k)) ]
       do ilat=1,Nlat
          do iorb=1,Norb
             do ispin=1,Nspin
                eps = vca_bath_%e(ilat,iorb,ispin,1:Nbath)
                vps = vca_bath_%v(ilat,iorb,ispin,1:Nbath)
                do i=1,L
                   Delta(ilat,ilat,ispin,ispin,iorb,iorb,i) = sum( vps(:)*vps(:)/(x(i) - eps(:)) )
                enddo
             enddo
          enddo
       enddo
       !
    case ("hybrid")             !hybrid: all _{ab} components allowed (inter-orbital local mixing present)
       !\Delta_{ab} = \sum_k [ V_{a}(k) * V_{b}(k)/(iw_n - E(k)) ]
       do ispin=1,Nspin
          eps  = vca_bath_%e(1     ,ispin,1     ,1:Nbath)
          vops = vca_bath_%v(1:Nlat,ispin,1:Norb,1:Nbath)
          do ilat=1,Nlat
             do jlat=1,Nlat
                do iorb=1,Norb
                   do jorb=1,Norb
                      do i=1,L
                         Delta(ilat,jlat,iorb,jorb,ispin,ispin,i) = sum( vops(ilat,iorb,:)*vops(jlat,jorb,:)/(x(i) - eps(:)) )
                      enddo
                   enddo
                enddo
             enddo
          enddo
       enddo
    end select
  end function delta_bath_real_main


  function delta_bath_real_main_(x,bath_) result(Delta)
    complex(8),dimension(:),intent(in)                  :: x
    type(effective_bath)                                :: vca_bath_
    complex(8),dimension(Nlat,Nlat,Nspin,Nspin,Norb,Norb,size(x)) :: Delta
    real(8),dimension(:)                                :: bath_
    logical                                             :: check
    check= check_bath_dimension(bath_)
    if(.not.check)stop "delta_bath_real_main_ error: wrong bath dimensions"
    call vca_allocate_bath(vca_bath_)
    call vca_set_bath(bath_,vca_bath_)
    Delta = delta_bath_real_main(x,vca_bath_)
    call vca_deallocate_bath(vca_bath_)
  end function delta_bath_real_main_





















































END MODULE VCA_BATH_FUNCTIONS
