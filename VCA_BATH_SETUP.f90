MODULE VCA_BATH_SETUP
  USE SF_CONSTANTS, only: zero
  USE SF_IOTOOLS, only:free_unit,reg,file_length,txtfy
  USE SF_LINALG, only: eye,inv
  USE SF_MISC, only: assert_shape
  USE VCA_INPUT_VARS
  USE VCA_VARS_GLOBAL
  USE VCA_AUX_FUNX
  implicit none

  private



  !##################################################################
  !
  !     USER BATH ROUTINES:
  !
  !##################################################################
  public :: get_bath_dimension
  public :: check_bath_dimension



  !##################################################################
  !
  !     VCA BATH ROUTINES:
  !
  !##################################################################
  !VCA BATH procedures:
  public :: allocate_vca_bath               !INTERNAL (for effective_bath)
  public :: deallocate_vca_bath             !INTERNAL (for effective_bath)
  public :: init_vca_bath                   !INTERNAL (for effective_bath)
  public :: write_vca_bath                  !INTERNAL (for effective_bath)
  public :: save_vca_bath                   !INTERNAL (for effective_bath)
  public :: set_vca_bath                    !INTERNAL (for effective_bath)
  public :: get_vca_bath                    !INTERNAL (for effective_bath)



contains

  


  !##################################################################
  !
  !     USER BATH ROUTINES:
  !
  !##################################################################
  !+-------------------------------------------------------------------+
  !PURPOSE  : Inquire the correct bath size to allocate the 
  ! the bath array in the calling program.
  !+-------------------------------------------------------------------+
  function get_bath_dimension(ispin_) result(bath_size)
    integer,optional               :: ispin_
    integer                        :: bath_size,ndx,ispin,iorb,jspin,jorb,io,jo
    !
    select case(bath_type)
    case default
       !e:[Nlat][Norb][Nspin][Nbath] + v:[Nlat][Norb][Nspin][Nbath]
       bath_size = Nlat*Norb*Nbath + Nlat*Norb*Nbath
       if(.not.present(ispin_))bath_size=Nspin*bath_size
    case('hybrid')
       !e:[1][1][Nspin][Nbath] + v [Nlat][Norb][Nspin][Nbath]
       bath_size = Nbath + Nlat*Norb*Nbath
       if(.not.present(ispin_))bath_size=Nspin*bath_size
    end select
  end function get_bath_dimension










  !##################################################################
  !
  !     VCA BATH ROUTINES:
  !
  !##################################################################
  !+-------------------------------------------------------------------+
  !PURPOSE  : Allocate the ED bath
  !+-------------------------------------------------------------------+
  subroutine allocate_vca_bath(vca_bath_)
    type(effective_bath) :: vca_bath_
    !
    if(vca_bath_%status)call deallocate_vca_bath(vca_bath_)
    !
    if(Nbath==0)stop "Allocate_vca_bath ERROR: Nbath==0"
    !
    select case(bath_type)
    case default
       !
       allocate(vca_bath_%e(Nlat,Norb,Nspin,Nbath))  !local energies of the bath per site,orb
       allocate(vca_bath_%v(Nlat,Norb,Nspin,Nbath))  !same-spin hybridization per site,orb
       !
    case('hybrid')
       !
       allocate(vca_bath_%e(1,1,Nspin,Nbath))        !local energies of the bath stand-alone
       allocate(vca_bath_%v(Nlat,Norb,Nspin,Nbath))  !same-spin hybridization, connects site,orb to bath sites
       !
    end select
    vca_bath_%status=.true.
  end subroutine allocate_vca_bath



  

  !+-------------------------------------------------------------------+
  !PURPOSE  : Deallocate the ED bath
  !+-------------------------------------------------------------------+
  subroutine deallocate_vca_bath(vca_bath_)
    type(effective_bath) :: vca_bath_
    if(allocated(vca_bath_%e))   deallocate(vca_bath_%e)
    if(allocated(vca_bath_%v))   deallocate(vca_bath_%v)
    vca_bath_%status=.false.
  end subroutine deallocate_vca_bath

  



  !+------------------------------------------------------------------+
  !PURPOSE  : Initialize the VCA loop, builindg H parameters and/or 
  !reading previous (converged) solution
  !+------------------------------------------------------------------+
  subroutine init_vca_bath(vca_bath_)
    type(effective_bath) :: vca_bath_
    integer              :: i,unit,flen,Nh
    integer              :: io,jo,iorb,ispin,jorb,jspin,ilat
    logical              :: IOfile
    real(8)              :: de,noise_tot
    character(len=21)    :: space
    real,dimension(Nbath):: ran
    !
    if(.not.vca_bath_%status)stop "init_vca_bath error: bath not allocated"
    !
    if(Nbath==0)stop "Allocate_vca_bath ERROR: Nbath==0"
    !
    !Get energies:
    call random_number(ran)
    forall(ilat=1:size(vca_bath_%e,1),iorb=1:size(vca_bath_%e,2),ispin=1:Nspin)&
         vca_bath_%e(ilat,iorb,ispin,:) = impHloc(ilat,ilat,iorb,iorb,ispin,ispin) + ran/10d0
    !
    !Get spin-keep yhbridizations
    do i=1,Nbath
       vca_bath_%v(:,:,:,i)=max(0.1d0,1.d0/sqrt(dble(Nbath)))
    enddo
    !
    !Read from file if exist:
    !
    inquire(file="bath.restart",exist=IOfile)
    if(IOfile)then
       write(LOGfile,"(A)")"Reading bath from file bath.restart"
       unit = free_unit()
       !
       open(unit,file="bath.restart")
       !
       read(unit,*)          !read the header:
       select case(bath_type)
       case default
          do ilat=1,Nlat
             do i=1,Nbath
                read(unit,*)((&
                     vca_bath_%e(ilat,iorb,ispin,i),&
                     vca_bath_%v(ilat,iorb,ispin,i),&
                     iorb=1,Norb),ispin=1,Nspin)
             enddo
          enddo
          !
       case ('hybrid')
          do ilat=1,Nlat
             do i=1,Nbath
                read(unit,*)(&
                     vca_bath_%e(1,1,ispin,i),&
                     (&
                     vca_bath_%v(ilat,iorb,ispin,i),&
                     iorb=1,Norb),ispin=1,Nspin)
             enddo
          enddo
       end select
       close(unit)
    endif
  end subroutine init_vca_bath




  !+-------------------------------------------------------------------+
  !PURPOSE  : write out the bath to a given unit with 
  ! the following column formatting: 
  ! [(Ek_iorb,Vk_iorb)_iorb=1,Norb]_ispin=1,Nspin
  !+-------------------------------------------------------------------+
  subroutine write_vca_bath(vca_bath_,unit)
    type(effective_bath) :: vca_bath_
    integer,optional     :: unit
    integer              :: unit_
    integer              :: i
    integer              :: io,jo,iorb,ispin,ilat
    complex(8)           :: hybr_aux
    complex(8)           :: hrep_aux(Nspin*Norb,Nspin*Norb)
    !
    unit_=LOGfile;if(present(unit))unit_=unit
    !
    if(.not.vca_bath_%status)stop "write_vca_bath error: bath not allocated"
    !
    select case(bath_type)
    case default
       write(unit_,"(90(A21,1X))")((&
            "#Ek_l"//str(iorb)//"_s"//str(ispin),"Vk_l"//str(iorb)//"_s"//str(ispin),&
            iorb=1,Norb),ispin=1,Nspin)
       do ilat=1,Nlat
          do i=1,Nbath
             write(unit_,"(90(F21.12,1X))")((&
                  vca_bath_%e(ilat,iorb,ispin,i),&
                  vca_bath_%v(ilat,iorb,ispin,i),&
                  iorb=1,Norb),ispin=1,Nspin)
          enddo
       enddo
       !
    case('hybrid')
       !
       write(unit_,"(90(A21,1X))")(&
            "#Ek_s"//reg(txtfy(ispin)),("Vk_l"//reg(txtfy(iorb))//"_s"//reg(txtfy(ispin)),&
            iorb=1,Norb),ispin=1,Nspin)
       do ilat=1,Nlat
          do i=1,Nbath
             write(unit_,"(90(F21.12,1X))")(&
                  vca_bath_%e(ilat,1,ispin,i),&
                  (vca_bath_%v(ilat,iorb,ispin,i),&
                  iorb=1,Norb),ispin=1,Nspin)
          enddo
       enddo
    end select
  end subroutine write_vca_bath






  !+-------------------------------------------------------------------+
  !PURPOSE  : save the bath to a given file using the write bath
  ! procedure and formatting: 
  !+-------------------------------------------------------------------+
  subroutine save_vca_bath(vca_bath_,file,used)
    type(effective_bath)      :: vca_bath_
    character(len=*),optional :: file
    character(len=256)        :: file_
    logical,optional          :: used
    logical                   :: used_
    character(len=16)         :: extension
    integer                   :: unit_
    !
    if(.not.vca_bath_%status)stop "save_vca_bath error: bath is not allocated"
    !
    used_=.false.;if(present(used))used_=used
    extension=".restart";if(used_)extension=".used"
    !
    file_="bath"//reg(extension)
    if(present(file))file_=reg(file)
    unit_=free_unit()
    open(unit_,file=reg(file_))
    call write_vca_bath(vca_bath_,unit_)
    close(unit_)
  end subroutine save_vca_bath




  !+-------------------------------------------------------------------+
  !PURPOSE  : set the bath components from a given user provided 
  ! bath-array 
  !+-------------------------------------------------------------------+
  subroutine set_vca_bath(bath_,vca_bath_)
    real(8),dimension(:)   :: bath_
    type(effective_bath)   :: vca_bath_
    integer                :: stride,io,jo,i
    integer                :: iorb,ispin,jorb,jspin,ibath,ilat
    logical                :: check
    complex(8)             :: hrep_aux(Nspin*Norb,Nspin*Norb)
    complex(8)             :: U(Nspin*Norb,Nspin*Norb)
    complex(8)             :: Udag(Nspin*Norb,Nspin*Norb)
    real(8)                :: element_R,element_I,eps_k,lambda_k
    !
    if(.not.vca_bath_%status)stop "set_vca_bath error: bath not allocated"
    !
    check = check_bath_dimension(bath_)
    if(.not.check)stop "set_vca_bath error: wrong bath dimensions"
    !
    select case(bath_type)
    case default
       stride = 0
       do ilat=1,Nlat
          do iorb=1,Norb
             do ispin=1,Nspin
                do i=1,Nbath
                   io = stride + i + index_stride_los(ilat,iorb,ispin)
                   vca_bath_%e(ilat,iorb,ispin,i) = bath_(io)
                enddo
             enddo
          enddo
       enddo
       stride = Nlat*Nspin*Norb*Nbath
       do ilat=1,Nlat
          do iorb=1,Norb
             do ispin=1,Nspin
                do i=1,Nbath
                   io = stride + i + index_stride_los(ilat,iorb,ispin)
                   vca_bath_%v(ilat,iorb,ispin,i) = bath_(io)
                enddo
             enddo
          enddo
       enddo
       !
    case ('hybrid')
       stride = 0
       do ispin=1,Nspin
          do i=1,Nbath
             io = stride + i + (ispin-1)*Nbath
             vca_bath_%e(1,1,ispin,i) = bath_(io)
          enddo
       enddo
       stride = Nspin*Nbath
       do ilat=1,Nlat
          do iorb=1,Norb
             do ispin=1,Nspin
                do i=1,Nbath
                   io = stride + i + index_stride_los(ilat,iorb,ispin)
                   vca_bath_%v(ilat,iorb,ispin,i) = bath_(io)
                enddo
             enddo
          enddo
       enddo
    end select
  end subroutine set_vca_bath



  !+-------------------------------------------------------------------+
  !PURPOSE  : copy the bath components back to a 1-dim array 
  !+-------------------------------------------------------------------+
  subroutine get_vca_bath(vca_bath_,bath_)
    type(effective_bath)   :: vca_bath_
    real(8),dimension(:)   :: bath_
    complex(8)             :: hrep_aux(Nspin*Norb,Nspin*Norb)
    complex(8)             :: U(Nspin*Norb,Nspin*Norb)
    complex(8)             :: Udag(Nspin*Norb,Nspin*Norb)
    integer                :: stride,io,jo,i,ilat
    integer                :: iorb,ispin,jorb,jspin,ibath
    logical                :: check
    !
    if(.not.vca_bath_%status)stop "get_vca_bath error: bath not allocated"
    !
    check=check_bath_dimension(bath_)
    if(.not.check)stop "get_vca_bath error: wrong bath dimensions"
    !
    select case(bath_type)
    case default
       stride = 0
       do ilat=1,Nlat
          do iorb=1,Norb
             do ispin=1,Nspin
                do i=1,Nbath
                   io = stride + i + index_stride_los(ilat,iorb,ispin)
                   bath_(io) = vca_bath_%e(ilat,iorb,ispin,i)
                enddo
             enddo
          enddo
       enddo
       stride = Nlat*Nspin*Norb*Nbath
       do ilat=1,Nlat
          do iorb=1,Norb
             do ispin=1,Nspin
                do i=1,Nbath
                   io = stride + i + index_stride_los(ilat,iorb,ispin)
                   bath_(io) = vca_bath_%v(ilat,iorb,ispin,i)
                enddo
             enddo
          enddo
       enddo
       !
    case ('hybrid')
       stride = 0
       do ispin=1,Nspin
          do i=1,Nbath
             io = stride + i + (ispin-1)*Nbath
             bath_(io) = vca_bath_%e(1,1,ispin,i)
          enddo
       enddo
       stride = Nspin*Nbath
       do ilat=1,Nlat
          do iorb=1,Norb
             do ispin=1,Nspin
                do i=1,Nbath
                   io = stride + i + index_stride_los(ilat,iorb,ispin)
                   bath_(io) = vca_bath_%v(ilat,iorb,ispin,i)
                enddo
             enddo
          enddo
       enddo
       !
    end select
  end subroutine get_vca_bath







  !##################################################################
  !
  !     USER BATH CHECKS:
  !
  !##################################################################
  !+-------------------------------------------------------------------+
  !PURPOSE  : Check if the dimension of the bath array are consistent
  !+-------------------------------------------------------------------+
  function check_bath_dimension(bath_) result(bool)
    real(8),dimension(:)           :: bath_
    integer                        :: Ntrue
    logical                        :: bool
    Ntrue = get_bath_dimension()
    bool  = ( size(bath_) == Ntrue )
  end function check_bath_dimension






END MODULE VCA_BATH_SETUP
