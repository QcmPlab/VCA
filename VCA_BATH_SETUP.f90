MODULE VCA_BATH_SETUP
   USE SF_CONSTANTS, only: zero
   USE SF_IOTOOLS, only:free_unit,reg,file_length,txtfy
   USE SF_LINALG, only: eye,inv
   USE SF_MISC, only: assert_shape
   USE VCA_INPUT_VARS
   USE VCA_VARS_GLOBAL
   implicit none


   private



  !##################################################################
   !
   !     VCA BATH ROUTINES:
   !
   !##################################################################
   !VCA BATH procedures:
   public :: vca_allocate_bath               !INTERNAL (for effective_bath)
   public :: vca_deallocate_bath             !INTERNAL (for effective_bath)
   public :: vca_init_bath                   !INTERNAL (for effective_bath)
   public :: vca_set_bath                    !INTERNAL (for effective_bath)


contains




   !##################################################################
   !
   !     USER BATH ROUTINES:
   !
   !##################################################################



   !##################################################################
   !
   !     VCA BATH ROUTINES:
   !
   !##################################################################
   !+-------------------------------------------------------------------+
   !PURPOSE  : Allocate the VCA bath
   !+-------------------------------------------------------------------+
   subroutine vca_allocate_bath(vca_bath_)
      type(effective_bath) :: vca_bath_
      !
      if(vca_bath_%status)call vca_deallocate_bath(vca_bath_)
      !
      !
      ! select case(bath_type)
      ! case default!normal_normal
      !
      allocate(vca_bath_%h(Nlat_bath,Nlat_bath,Nspin,Nspin,Norb_bath,Norb_bath))
      allocate(vca_bath_%v(Nlat,Nlat_bath,Nspin,Nspin,Norb,Norb_bath))  !vertical rectangular
      !
      vca_bath_%h=0d0
      vca_bath_%v=0d0
      vca_bath_%status=.true.
   end subroutine vca_allocate_bath





   !+-------------------------------------------------------------------+
   !PURPOSE  : Deallocate the VCA bath
   !+-------------------------------------------------------------------+
   subroutine vca_deallocate_bath(vca_bath_)
      type(effective_bath) :: vca_bath_
      if(allocated(vca_bath_%h))   deallocate(vca_bath_%h)
      if(allocated(vca_bath_%v))   deallocate(vca_bath_%v)
      vca_bath_%status=.false.
   end subroutine vca_deallocate_bath





   !+------------------------------------------------------------------+
   !PURPOSE  : Initialize the VCA loop, builindg H parameters and/or 
   !reading previous (converged) solution
   !+------------------------------------------------------------------+
   subroutine vca_init_bath(vca_bath_)
      type(effective_bath) :: vca_bath_
      integer              :: i,unit,flen,Nh
      integer              :: io,jo,ispin,iorb,jorb,jspin,ilat
      logical              :: IOfile
      real(8)              :: de,hwband
      character(len=21)    :: space
      !
      hwband=1d0
      !
      if(.not.vca_bath_%status)stop "vca_init_bath error: bath not allocated"
      !
      !
      vca_bath_%h = 1d0
      vca_bath_%v = 0d0
      !
      !
   end subroutine vca_init_bath


   !+-------------------------------------------------------------------+
   !PURPOSE  : set the bath components from a given user provided 
   ! bath-array 
   !+-------------------------------------------------------------------+
   subroutine vca_set_bath(h_in,v_in,vca_bath_)
      complex(8),dimension(:,:,:,:,:,:)   :: h_in, v_in
      type(effective_bath)                :: vca_bath_
      integer                             :: ilat,ispin,iorb
      !
      call assert_shape(h_in,[Nlat_bath,Nlat_bath,Nspin,Nspin,Norb_bath,Norb_bath],"print_Hcluster_nnn","Hloc")
      call assert_shape(v_in,[Nlat,Nlat_bath,Nspin,Nspin,Norb,Norb_bath],"print_Hcluster_nnn","Hloc")
      !
      if(.not.vca_bath_%status)stop "set_vca_bath error: bath not allocated"
      !
      !
      !
      vca_bath_%h = h_in
      vca_bath_%v = v_in
      !
      if(hfmode)then
         do ilat=1,Nlat_bath
            do iorb=1,Norb_bath
               do ispin=1,Nspin
                     if( .not. hfshift) vca_bath_%h(ilat,ilat,ispin,ispin,iorb,iorb) = &
                          vca_bath_%h(ilat,ilat,ispin,ispin,iorb,iorb)-(0.5d0*Uloc(iorb)+0.5d0*Ust*(Norb-1)+0.5d0*(Ust-Jh)*(Norb-1))
               enddo
            enddo
         enddo
      endif
      !
   end subroutine vca_set_bath



END MODULE VCA_BATH_SETUP
