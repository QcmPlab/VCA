MODULE VCA_AUX_FUNX
  USE VCA_INPUT_VARS
  USE VCA_VARS_GLOBAL
  USE VCA_BATH_SETUP
  USE SF_TIMER
  USE SF_LINALG
  USE SF_MISC, only: assert_shape
  USE SF_IOTOOLS, only:free_unit,reg
  implicit none
  private



  interface vca_lso2nnn_reshape
     module procedure :: d_nlso2nnn_scalar
     module procedure :: c_nlso2nnn_scalar
  end interface vca_lso2nnn_reshape


  interface vca_nnn2lso_reshape
     module procedure :: d_nnn2nlso_scalar
     module procedure :: c_nnn2nlso_scalar
  end interface vca_nnn2lso_reshape


  interface vca_set_Hcluster
     module procedure :: set_Hcluster_nnn
     module procedure :: set_Hcluster_lso
  end interface vca_set_Hcluster

  interface vca_set_Hk
     module procedure :: set_Hk_nnn
     module procedure :: set_Hk_lso
  end interface vca_set_Hk

  interface vca_print_Hcluster
     module procedure :: print_Hcluster_nnn
     module procedure :: print_Hcluster_lso
  end interface vca_print_Hcluster


  public :: vca_get_cluster_dimension
  !
  public :: vca_set_Hcluster
  public :: vca_set_Hk
  public :: vca_print_Hcluster
  !
  public :: vca_lso2nnn_reshape
  public :: vca_nnn2lso_reshape
  !
  !
  public :: search_chemical_potential





contains

  !##################################################################
  !                   DIMENSION PROCEDURES
  !##################################################################
  function vca_get_cluster_dimension(with_bath) result(Ncluster)
    logical,optional :: with_bath
    logical          :: bool
    integer          :: Ns
    integer          :: Ncluster
    !
    bool = .false. ; if(present(with_bath))bool = with_bath
    !
    !Count how many levels are there in the cluster:
    Ns = Nlat*Norb
    if(bool)then
       ! select case(bath_type)
       ! case default
       Ns = (Nbath+1)*Nlat*Norb !Norb per site plus Nbath per orb per site
       ! case ('hybrid')
       !    Ns = Nbath+Nlat*Norb     !Norb per site plus shared Nbath sites
       ! end select
    endif
    !
    !Count the spin:
    Ncluster  = Nspin*Ns
    !
  end function vca_get_cluster_dimension




  !##################################################################
  !                   HCLUSTER ROUTINES
  !##################################################################
  !+------------------------------------------------------------------+
  !PURPOSE  : 
  !+------------------------------------------------------------------+
  subroutine print_Hcluster_nnn(hloc,file)
    complex(8),dimension(:,:,:,:,:,:) :: hloc
    character(len=*),optional         :: file
    integer                           :: ilat,jlat
    integer                           :: iorb,jorb
    integer                           :: ispin,jspin
    integer                           :: unit
    character(len=32)                 :: fmt
    !
    call assert_shape(Hloc,[Nlat,Nlat,Nspin,Nspin,Norb,Norb],"print_Hcluster_nnn","Hloc")
    !
    unit=LOGfile;
    !
    if(present(file))then
       open(free_unit(unit),file=reg(file))
       write(LOGfile,"(A)")"print_Hloc to file :"//reg(file)
    endif
    write(fmt,"(A,I0,A)")"(",Nlat*Nspin*Norb,"F9.2)"
    do ilat=1,Nlat
       do ispin=1,Nspin
          do iorb=1,Norb
             write(unit,fmt)(((dreal(Hloc(ilat,jlat,ispin,jspin,iorb,jorb)),jlat=1,Nlat),jspin=1,Nspin),jorb=1,Norb)
          enddo
       enddo
    enddo
    write(unit,*)""
    do ilat=1,Nlat
       do ispin=1,Nspin
          do iorb=1,Norb
             write(unit,fmt)(((dimag(Hloc(ilat,jlat,ispin,jspin,iorb,jorb)),jlat=1,Nlat),jspin=1,Nspin),jorb=1,Norb)
          enddo
       enddo
    enddo
    write(unit,*)""
    if(present(file))close(unit)
  end subroutine print_Hcluster_nnn
  !
  subroutine print_Hcluster_lso(hloc,file)
    complex(8),dimension(:,:) :: hloc
    character(len=*),optional :: file
    integer                   :: unit,is,js
    character(len=32)         :: fmt
    !
    call assert_shape(Hloc,[Nlat*Nspin*Norb,Nlat*Nspin*Norb],"print_Hcluster_lso","Hloc")
    !
    unit=LOGfile;
    !
    if(present(file))then
       open(free_unit(unit),file=reg(file))
       write(LOGfile,"(A)")"print_Hloc to file :"//reg(file)
    endif
    write(fmt,"(A,I0,A)")"(",Nlat*Nspin*Norb,"A)"
    do is=1,Nlat*Nspin*Norb
       write(unit,fmt)(str(Hloc(is,js),3)//" ",js=1,Nlat*Nspin*Norb)
    enddo
    write(unit,*)""
    if(present(file))close(unit)
  end subroutine print_Hcluster_lso


  !+------------------------------------------------------------------+
  !PURPOSE  : 
  !+------------------------------------------------------------------+
  subroutine set_Hcluster_nnn(hloc)
    complex(8),dimension(:,:,:,:,:,:) :: Hloc ![Nlat][Nlat][Nspin][Nspin][Norb][Norb]
    call assert_shape(Hloc,[Nlat,Nlat,Nspin,Nspin,Norb,Norb],"set_Hcluster_nnn","Hloc")
    !
    impHloc = Hloc
    !
    write(LOGfile,"(A)")"Set Hcluster: done"
    if(verbose>2)call vca_print_Hcluster(impHloc)
  end subroutine set_Hcluster_nnn

  subroutine set_Hcluster_lso(hloc)
    complex(8),dimension(:,:) :: Hloc ![Nlat*Nspin*Norb][Nlat*Nspin*Norb]
    call assert_shape(Hloc,[Nlat*Nspin*Norb,Nlat*Nspin*Norb],"set_Hcluster_lso","Hloc")
    !
    impHloc = vca_lso2nnn_reshape(Hloc,Nlat,Nspin,Norb)
    !
    write(LOGfile,"(A)")"Set Hcluster: done"
    if(verbose>2)call vca_print_Hcluster(impHloc)
  end subroutine set_Hcluster_lso


  subroutine set_Hk_nnn(hloc)
    complex(8),dimension(:,:,:,:,:,:,:) :: Hloc ![Nlat][Nlat][Nspin][Nspin][Norb][Norb]
    call assert_shape(Hloc,[Nlat,Nlat,Nspin,Nspin,Norb,Norb,Nkpts**ndim],"set_Hk_nnn","Hloc") 
    !
    impHk = Hloc
    !
    write(LOGfile,"(A)")"Set Hk: done"
    !if(verbose>2)call vca_print_Hcluster(impHk)
  end subroutine set_Hk_nnn

  subroutine set_Hk_lso(hloc)
    complex(8),dimension(:,:,:) :: Hloc ![Nlat*Nspin*Norb][Nlat*Nspin*Norb][Nkpts**Ndim]
    integer                     :: i
    call assert_shape(Hloc,[Nlat*Nspin*Norb,Nlat*Nspin*Norb,Nkpts**ndim],"set_Hk_lso","Hloc") 
    !
    do i=1,Nkpts**Ndim 
        impHk(:,:,:,:,:,:,i) = vca_lso2nnn_reshape(Hloc(:,:,i),Nlat,Nspin,Norb)
    enddo
    !
    write(LOGfile,"(A)")"Set Hk: done"
    !if(verbose>2)call vca_print_Hcluster(impHk)
  end subroutine set_Hk_lso








  !+-----------------------------------------------------------------------------+!
  !PURPOSE: 
  ! reshape a matrix from/to the
  ! [Nlso][Nlso] and
  ! [Nlat][Nlat][Norb][Norb][Nspin][Nspin]
  ! shapes.
  ! - dble & cmplx
  ! 0-Reshape a scalar array, dble & cmplx
  ! 1-Reshape a function array (:,:,:,:,:,:,1:L)
  !+-----------------------------------------------------------------------------+!
  function d_nlso2nnn_scalar(Hlso,Nlat,Nspin,Norb) result(Hnnn)
    integer                                            :: Nlat,Nspin,Norb
    real(8),dimension(Nlat*Nspin*Norb,Nlat*Nspin*Norb) :: Hlso
    real(8),dimension(Nlat,Nlat,Nspin,Nspin,Norb,Norb) :: Hnnn
    integer                                            :: ilat,jlat
    integer                                            :: iorb,jorb
    integer                                            :: ispin,jspin
    integer                                            :: is,js
    Hnnn=zero
    do ilat=1,Nlat
       do jlat=1,Nlat
          do ispin=1,Nspin
             do jspin=1,Nspin
                do iorb=1,Norb
                   do jorb=1,Norb
                      is = index_stride_lso(ilat,ispin,iorb)
                      js = index_stride_lso(jlat,jspin,jorb)
                      Hnnn(ilat,jlat,ispin,jspin,iorb,jorb) = Hlso(is,js)
                   enddo
                enddo
             enddo
          enddo
       enddo
    enddo
  end function d_nlso2nnn_scalar
  !
  function d_nnn2nlso_scalar(Hnnn,Nlat,Nspin,Norb) result(Hlso)
    integer                                            :: Nlat,Nspin,Norb
    real(8),dimension(Nlat,Nlat,Nspin,Nspin,Norb,Norb) :: Hnnn
    real(8),dimension(Nlat*Nspin*Norb,Nlat*Nspin*Norb) :: Hlso
    integer                                            :: ilat,jlat
    integer                                            :: iorb,jorb
    integer                                            :: ispin,jspin
    integer                                            :: is,js
    Hlso=zero
    do ilat=1,Nlat
       do jlat=1,Nlat
          do ispin=1,Nspin
             do jspin=1,Nspin
                do iorb=1,Norb
                   do jorb=1,Norb
                      is = index_stride_lso(ilat,ispin,iorb)
                      js = index_stride_lso(jlat,jspin,jorb)
                      Hlso(is,js) = Hnnn(ilat,jlat,ispin,jspin,iorb,jorb)
                   enddo
                enddo
             enddo
          enddo
       enddo
    enddo
  end function d_nnn2nlso_scalar
  !
  function c_nlso2nnn_scalar(Hlso,Nlat,Nspin,Norb) result(Hnnn)
    integer                                               :: Nlat,Nspin,Norb
    complex(8),dimension(Nlat*Nspin*Norb,Nlat*Nspin*Norb) :: Hlso
    complex(8),dimension(Nlat,Nlat,Nspin,Nspin,Norb,Norb) :: Hnnn
    integer                                               :: ilat,jlat
    integer                                               :: iorb,jorb
    integer                                               :: ispin,jspin
    integer                                               :: is,js
    Hnnn=zero
    do ilat=1,Nlat
       do jlat=1,Nlat
          do ispin=1,Nspin
             do jspin=1,Nspin
                do iorb=1,Norb
                   do jorb=1,Norb
                      is = index_stride_lso(ilat,ispin,iorb)
                      js = index_stride_lso(jlat,jspin,jorb)
                      Hnnn(ilat,jlat,ispin,jspin,iorb,jorb) = Hlso(is,js)
                   enddo
                enddo
             enddo
          enddo
       enddo
    enddo
  end function c_nlso2nnn_scalar
  !
  function c_nnn2nlso_scalar(Hnnn,Nlat,Nspin,Norb) result(Hlso)
    integer                                               :: Nlat,Nspin,Norb
    complex(8),dimension(Nlat,Nlat,Nspin,Nspin,Norb,Norb) :: Hnnn
    complex(8),dimension(Nlat*Nspin*Norb,Nlat*Nspin*Norb) :: Hlso
    integer                                               :: ilat,jlat
    integer                                               :: iorb,jorb
    integer                                               :: ispin,jspin
    integer                                               :: is,js
    Hlso=zero
    do ilat=1,Nlat
       do jlat=1,Nlat
          do ispin=1,Nspin
             do jspin=1,Nspin
                do iorb=1,Norb
                   do jorb=1,Norb
                      is = index_stride_lso(ilat,ispin,iorb)
                      js = index_stride_lso(jlat,jspin,jorb)
                      Hlso(is,js) = Hnnn(ilat,jlat,ispin,jspin,iorb,jorb)
                   enddo
                enddo
             enddo
          enddo
       enddo
    enddo
  end function c_nnn2nlso_scalar


















  !##################################################################
  !##################################################################
  ! ROUTINES TO SEARCH CHEMICAL POTENTIAL UP TO SOME ACCURACY
  ! can be used to fix any other *var so that  *ntmp == nread
  !##################################################################
  !##################################################################

  !+------------------------------------------------------------------+
  !PURPOSE  : 
  !+------------------------------------------------------------------+
  subroutine search_chemical_potential(var,ntmp,converged)
    real(8),intent(inout) :: var
    real(8),intent(in)    :: ntmp
    logical,intent(inout) :: converged
    logical               :: bool
    real(8)               :: ndiff
    integer,save          :: count=0,totcount=0,i
    integer,save          :: nindex=0
    integer               :: nindex_old(3)
    real(8)               :: ndelta_old,nratio
    integer,save          :: nth_magnitude=-2,nth_magnitude_old=-2
    real(8),save          :: nth=1.d-2
    logical,save          :: ireduce=.true.
    integer               :: unit
    !
    ndiff=ntmp-nread
    nratio = 0.5d0;!nratio = 1.d0/(6.d0/11.d0*pi)
    !
    !check actual value of the density *ntmp* with respect to goal value *nread*
    count=count+1
    totcount=totcount+1
    if(count>2)then
       do i=1,2
          nindex_old(i+1)=nindex_old(i)
       enddo
    endif
    nindex_old(1)=nindex
    !
    if(ndiff >= nth)then
       nindex=-1
    elseif(ndiff <= -nth)then
       nindex=1
    else
       nindex=0
    endif
    !
    ndelta_old=ndelta
    bool=nindex/=0.AND.( (nindex+nindex_old(1)==0).OR.(nindex+sum(nindex_old(:))==0) )
    !if(nindex_old(1)+nindex==0.AND.nindex/=0)then !avoid loop forth and back
    if(bool)then
       ndelta=ndelta_old*nratio !decreasing the step
    else
       ndelta=ndelta_old
    endif
    !
    if(ndelta_old<1.d-9)then
       ndelta_old=0.d0
       nindex=0
    endif
    !update chemical potential
    var=var+dble(nindex)*ndelta
    !xmu=xmu+dble(nindex)*ndelta
    !
    !Print information
    write(LOGfile,"(A,f16.9,A,f15.9)")"n    = ",ntmp," /",nread
    if(nindex>0)then
       write(LOGfile,"(A,es16.9,A)")"shift= ",nindex*ndelta," ==>"
    elseif(nindex<0)then
       write(LOGfile,"(A,es16.9,A)")"shift= ",nindex*ndelta," <=="
    else
       write(LOGfile,"(A,es16.9,A)")"shift= ",nindex*ndelta," == "
    endif
    write(LOGfile,"(A,f15.9)")"var  = ",var
    write(LOGfile,"(A,ES16.9,A,ES16.9)")"dn   = ",ndiff,"/",nth
    unit=free_unit()
    open(unit,file="search_mu_iteration"//reg(file_suffix)//".ed",position="append")
    write(unit,*)var,ntmp,ndiff
    close(unit)
    !
    !check convergence within actual threshold
    !if reduce is activetd
    !if density is in the actual threshold
    !if DMFT is converged
    !if threshold is larger than nerror (i.e. this is not last loop)
    bool=ireduce.AND.(abs(ndiff)<nth).AND.converged.AND.(nth>nerr)
    if(bool)then
       nth_magnitude_old=nth_magnitude        !save old threshold magnitude
       nth_magnitude=nth_magnitude_old-1      !decrease threshold magnitude || floor(log10(abs(ntmp-nread)))
       nth=max(nerr,10.d0**(nth_magnitude))   !set the new threshold 
       count=0                                !reset the counter
       converged=.false.                      !reset convergence
       ndelta=ndelta_old*nratio                  !reduce the delta step
       !
    endif
    !
    !if density is not converged set convergence to .false.
    if(abs(ntmp-nread)>nth)converged=.false.
    !
    !check convergence for this threshold
    !!---if smallest threshold-- NO MORE
    !if reduce is active (you reduced the treshold at least once)
    !if # iterations > max number
    !if not yet converged
    !set threshold back to the previous larger one.
    !bool=(nth==nerr).AND.ireduce.AND.(count>niter).AND.(.not.converged)
    bool=ireduce.AND.(count>niter).AND.(.not.converged)
    if(bool)then
       ireduce=.false.
       nth=10.d0**(nth_magnitude_old)
    endif
    !
    write(LOGfile,"(A,I5)")"count= ",count
    write(LOGfile,"(A,L2)")"Converged=",converged
    print*,""
    !
  end subroutine search_chemical_potential















END MODULE VCA_AUX_FUNX




