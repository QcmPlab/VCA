MODULE VCA_AUX_FUNX
  USE VCA_VARS_GLOBAL
  USE VCA_BATH_SETUP
  USE SF_MISC, only: assert_shape
  USE SF_IOTOOLS, only:free_unit,reg
  implicit none
  private



  interface vca_los2nnn_reshape
     module procedure :: d_nlos2nnn_scalar
     module procedure :: c_nlos2nnn_scalar
  end interface vca_los2nnn_reshape


  interface vca_nnn2los_reshape
     module procedure :: d_nnn2nlos_scalar
     module procedure :: c_nnn2nlos_scalar
  end interface vca_nnn2los_reshape


  interface vca_set_Hcluster
     module procedure :: set_Hcluster_nnn
     module procedure :: set_Hcluster_los
  end interface vca_set_Hcluster


  interface vca_print_Hcluster
     module procedure :: print_Hcluster_nnn
     module procedure :: print_Hcluster_los
  end interface vca_print_Hcluster


  public :: vca_loS2nnn_reshape
  public :: vca_nnn2loS_reshape
  !
  public :: vca_set_Hcluster
  !
  public :: vca_print_Hcluster
  !
  public :: search_chemical_potential





contains



  !##################################################################
  !                   HCLUSTER ROUTINES
  !##################################################################
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

  subroutine set_Hcluster_los(hloc)
    complex(8),dimension(:,:) :: Hloc ![Nlat*Nspin*Norb][Nlat*Nspin*Norb]
    call assert_shape(Hloc,[Nlat*Nspin*Norb,Nlat*Nspin*Norb],"set_Hcluster_los","Hloc")
    !
    impHloc = vca_los2nnn_reshape(Hloc,Nlat,Nspin,Norb)
    !
    write(LOGfile,"(A)")"Set Hcluster: done"
    if(verbose>2)call vca_print_Hcluster(impHloc)
  end subroutine set_Hcluster_los





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
  subroutine print_Hcluster_los(hloc,file)
    complex(8),dimension(:,:) :: hloc
    character(len=*),optional :: file
    integer                   :: ilat,jlat
    integer                   :: iorb,jorb
    integer                   :: ispin,jspin
    integer                   :: unit,is,js
    character(len=32)         :: fmt
    !
    call assert_shape(Hloc,[Nlat*Nspin*Norb,Nlat*Nspin*Norb],"print_Hcluster_los","Hloc")
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
  end subroutine print_Hcluster_los












  !+-----------------------------------------------------------------------------+!
  !PURPOSE: 
  ! reshape a matrix from/to the
  ! [Nlos][Nlos] and
  ! [Nlat][Nlat][Norb][Norb][Nspin][Nspin]
  ! shapes.
  ! - dble & cmplx
  ! 0-Reshape a scalar array, dble & cmplx
  ! 1-Reshape a function array (:,:,:,:,:,:,1:L)
  !+-----------------------------------------------------------------------------+!
  function d_nlos2nnn_scalar(Hlso,Nlat,Nspin,Norb) result(Hnnn)
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
                      is = index_stride_los(ilat,ispin,iorb)
                      js = index_stride_los(jlat,jspin,jorb)
                      Hnnn(ilat,jlat,ispin,jspin,iorb,jorb) = Hlso(is,js)
                   enddo
                enddo
             enddo
          enddo
       enddo
    enddo
  end function d_nlos2nnn_scalar
  !
  function d_nnn2nlos_scalar(Hnnn,Nlat,Nspin,Norb) result(Hlso)
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
                      is = index_stride_los(ilat,ispin,iorb)
                      js = index_stride_los(jlat,jspin,jorb)
                      Hlso(is,js) = Hnnn(ilat,jlat,ispin,jspin,iorb,jorb)
                   enddo
                enddo
             enddo
          enddo
       enddo
    enddo
  end function d_nnn2nlos_scalar
  !
  function c_nlos2nnn_scalar(Hlso,Nlat,Nspin,Norb) result(Hnnn)
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
                      is = index_stride_los(ilat,ispin,iorb)
                      js = index_stride_los(jlat,jspin,jorb)
                      Hnnn(ilat,jlat,ispin,jspin,iorb,jorb) = Hlso(is,js)
                   enddo
                enddo
             enddo
          enddo
       enddo
    enddo
  end function c_nlos2nnn_scalar
  !
  function c_nnn2nlos_scalar(Hnnn,Nlat,Nspin,Norb) result(Hlso)
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
                      is = index_stride_los(ilat,ispin,iorb)
                      js = index_stride_los(jlat,jspin,jorb)
                      Hlso(is,js) = Hnnn(ilat,jlat,ispin,jspin,iorb,jorb)
                   enddo
                enddo
             enddo
          enddo
       enddo
    enddo
  end function c_nnn2nlos_scalar


















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





! interface vca_tile_Hcluster
!    module procedure :: vca_tile_Hcluster_normal
!    module procedure :: vca_tile_Hcluster_bath
! end interface vca_tile_Hcluster


! interface vca_build_Hsystem
!    module procedure :: vca_build_Hsystem_normal
!    module procedure :: vca_build_Hsystem_bath
! end interface vca_build_Hsystem


! public :: vca_get_system_dimension
! public :: vca_get_cluster_dimension
! !
! public :: vca_tile_Hcluster
! public :: vca_build_Hsystem
!  

! !##################################################################
! !                   DIMENSION PROCEDURES
! !##################################################################
! function vca_get_cluster_dimension(with_bath) result(Ncluster)
!   logical,optional :: with_bath
!   logical          :: bool
!   integer          :: Ns
!   integer          :: Ncluster
!   !
!   bool = .false. ; if(present(with_bath))bool = with_bath
!   !
!   !Count how many levels are there in the cluster:
!   Ns = Nlat*Norb
!   if(bool)then
!      ! select case(bath_type)
!      ! case default
!      Ns = (Nbath+1)*Nlat*Norb !Norb per site plus Nbath per orb per site
!      ! case ('hybrid')
!      !    Ns = Nbath+Nlat*Norb     !Norb per site plus shared Nbath sites
!      ! end select
!   endif
!   !
!   !Count the spin:
!   Ncluster  = Nspin*Ns
!   !
! end function vca_get_cluster_dimension
! !
! function vca_get_system_dimension(with_bath) result(Nsys)
!   logical,optional :: with_bath
!   logical          :: bool
!   integer          :: Ns, Nlevels
!   integer          :: Nsys
!   !
!   bool = .false. ; if(present(with_bath))bool = with_bath
!   !
!   !Count how many levels are there in the cluster:
!   Ns  = vca_get_cluster_dimension(bool)
!   !
!   !Count the copies:
!   Nsys = Ncopies*Ns
!   !
! end function vca_get_system_dimension





! !##################################################################
! !                  CLUSTER TILING PROCEDURES
! !##################################################################
! subroutine vca_tile_Hcluster_normal(Hcluster,Htile)
!   real(8),dimension(:,:) :: Hcluster
!   real(8),dimension(:,:) :: Htile
!   integer                :: Nsys
!   integer                :: Nc,unit
!   integer                :: i,j,icopy,ic,jc
!   !
!   Nc   = vca_get_cluster_dimension()
!   Nsys = vca_get_system_dimension()
!   !
!   call assert_shape(Hcluster,[Nlat*Nspin*Norb,Nc],"vca_tile_Hcluster","Hcluster")
!   call assert_shape(Htile,[size(Htile,1),Nsys],"vca_tile_Hcluster","Htile")
!   Htile=0d0
!   do icopy=1,Ncopies          !loop over the number of copies of the cluster:
!      !
!      do ic=1,Nc           !loop over the orbital-site-spin index of the cluster 
!         do jc=1,Nc
!            i = ic + (icopy-1)*Nc
!            j = jc + (icopy-1)*Nc
!            Htile(i,j) = Hcluster(ic,jc)
!         enddo
!      enddo
!   enddo
!   !
!   open(free_unit(unit),file="Htile_matrix.dat")
!   do i=1,Nsys
!      write(unit,"(1000000(F5.2,1x))")(Htile(i,j),j=1,Nsys)
!   enddo
!   close(unit)
! end subroutine vca_tile_Hcluster_normal
! !
! subroutine vca_tile_Hcluster_bath(Hcluster,bath,Htile_bath)
!   real(8),dimension(:,:)               :: Hcluster
!   real(8),dimension(:)                 :: bath
!   real(8),dimension(:,:)               :: Htile_bath
!   integer                              :: Nsys,Nsys_bath
!   integer                              :: Nc,Nc_bath,unit
!   integer                              :: i,j,icopy,ic,jc,ilat,iorb,ispin,ii,jj,ibath
!   type(effective_bath)                 :: vca_bath_
!   !
!   Nc        = vca_get_cluster_dimension(.false.)
!   Nc_bath   = vca_get_cluster_dimension(.true.)
!   Nsys      = vca_get_system_dimension(.false.)
!   Nsys_bath = vca_get_system_dimension(.true.)    
!   !
!   call assert_shape(Hcluster,[Nlat*Nspin*Norb,Nc],"vca_tile_Hcluster_bath","Hcluster")
!   call assert_shape(Htile_bath,[size(Htile_bath,1),Nsys_bath],"vca_tile_Hcluster_bath","Htile_bath")
!   !
!   call vca_allocate_bath(vca_bath_)
!   call vca_set_bath(bath,vca_bath_)
!   !
!   Htile_bath=0d0
!   do icopy=1,Ncopies
!      !
!      do ilat=1,Nlat              !# of cluster sites
!         do iorb=1,Norb           !# of orbital per site
!            do ispin=1,Nspin      !# of spin per site
!               ic = index_stride_los(ilat,iorb,ispin)
!               !
!               i  = ic + (icopy-1)*Nc
!               ii = i  + (i-1)*Nbath
!               do ibath=1,Nbath
!                  Htile_bath(ii+ibath,ii+ibath) = vca_bath_%e(ilat,iorb,ispin,ibath)
!                  Htile_bath(ii,ii+ibath)       = vca_bath_%v(ilat,iorb,ispin,ibath)
!                  Htile_bath(ii+ibath,ii)       = vca_bath_%v(ilat,iorb,ispin,ibath)
!               enddo
!               !
!               do jc=1,Nc
!                  j  = jc + (icopy-1)*Nc
!                  jj = j  + (j-1)*Nbath
!                  Htile_bath(ii,jj) = Hcluster(ic,jc)
!               enddo
!               !
!            enddo
!         enddo
!      enddo
!   enddo
!   !
!   call vca_deallocate_bath(vca_bath_)
!   !
!   open(free_unit(unit),file="Htile_bath_matrix.dat")
!   do i=1,Nsys_bath
!      write(unit,"(1000000(F5.2,1x))")(Htile_bath(i,j),j=1,Nsys_bath)
!   enddo
!   close(unit)
! end subroutine vca_tile_Hcluster_bath









! !##################################################################
! !                  SYSTEM BUILDING PROCEDURES
! !##################################################################
! subroutine vca_build_Hsystem_normal(Hsys,Hsys_bath)
!   real(8),dimension(:,:),intent(in)    :: Hsys
!   real(8),dimension(:,:),intent(out)   :: Hsys_bath
!   integer                              :: Nsys
!   !
!   Nsys  = vca_get_system_dimension(with_bath=.false.)
!   !
!   call assert_shape(Hsys,[Nsys,Ncopies*Nlat*Norb*Nspin],"vca_build_Hsystem","Hsys")
!   call assert_shape(Hsys_bath,[Nsys,Nsys],"vca_build_Hsystem","Hsys_bath")
!   Hsys_bath = Hsys
!   !
! end subroutine vca_build_Hsystem_normal

! subroutine vca_build_Hsystem_bath(Hsys,Bath,Hsys_bath)
!   real(8),dimension(:,:),intent(in)  :: Hsys
!   real(8),dimension(:),intent(in)    :: Bath
!   real(8),dimension(:,:),intent(out) :: Hsys_bath
!   integer                            :: Nlos
!   integer                            :: ibath,unit
!   integer                            :: icopy,jcopy
!   integer                            :: i,ii
!   integer                            :: j,jj
!   integer                            :: io,jo,ilat,iorb,ispin    
!   integer                            :: Nsys,Nsys_bath,Nb
!   type(effective_bath)               :: vca_bath_
!   !
!   write(LOGfile,"(A)")"Enter Build Hsystem"
!   !
!   Nsys          = vca_get_system_dimension(with_bath=.false.)
!   Nsys_bath     = vca_get_system_dimension(with_bath=.true.)
!   Nlos          = vca_get_cluster_dimension(with_bath=.false.)
!   !
!   call assert_shape(Hsys,[Nsys,Ncopies*Nlat*Norb*Nspin],"vca_build_Hsystem_bath","Hsys")
!   call assert_shape(Hsys_bath,[Nsys_bath,Nsys_bath],"vca_build_Hsystem_bath","Hsys_bath")
!   !
!   call vca_allocate_bath(vca_bath_)
!   call vca_set_bath(bath,vca_bath_)
!   !
!   Hsys_bath=0d0
!   do icopy=1,Ncopies
!      do ilat=1,Nlat              !# of cluster sites
!         do iorb=1,Norb           !# of orbital per site
!            do ispin=1,Nspin      !# of spin per site
!               io = index_stride_los(ilat,iorb,ispin)
!               i  = io + (icopy-1)*Nlat*Norb*Nspin
!               ii = i  + (i-1)*Nbath
!               do ibath=1,Nbath
!                  Hsys_bath(ii+ibath,ii+ibath) = vca_bath_%e(ilat,iorb,ispin,ibath)
!                  Hsys_bath(ii,ii+ibath)       = vca_bath_%v(ilat,iorb,ispin,ibath)
!                  Hsys_bath(ii+ibath,ii)       = vca_bath_%v(ilat,iorb,ispin,ibath)
!               enddo
!               !
!               do jcopy=1,Ncopies
!                  do jo=1,Nlos
!                     j  = jo + (jcopy-1)*Nlat*Norb*Nspin
!                     jj = j  + (j-1)*Nbath
!                     Hsys_bath(ii,jj) = Hsys(i,j)
!                     Hsys_bath(jj,ii) = Hsys(j,i)
!                  enddo
!               enddo
!               !
!            enddo
!         enddo
!      enddo
!   enddo
!   !
!   call vca_deallocate_bath(vca_bath_)
!   !
!   open(free_unit(unit),file="Hsys_bath_matrix.dat")
!   do i=1,Nsys_bath
!      write(unit,"(1000000(F5.2,1x))")(Hsys_bath(i,j),j=1,Nsys_bath)
!   enddo
!   close(unit)
! end subroutine vca_build_Hsystem_bath
