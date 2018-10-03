MODULE VCA_OBSERVABLES
  USE VCA_INPUT_VARS
  USE VCA_VARS_GLOBAL
  USE VCA_SETUP
  USE VCA_DIAG
  USE VCA_AUX_FUNX
  USE VCA_EIGENSPACE
  USE VCA_HAMILTONIAN
  !
  USE SF_CONSTANTS, only:zero,pi,xi
  USE SF_IOTOOLS, only:free_unit,reg,txtfy
  USE SF_ARRAYS, only: arange
  USE SF_LINALG
  !
  implicit none
  private
  !
  public :: observables_cluster



  logical,save                         :: iolegend=.true.
  real(8),dimension(:,:),allocatable   :: dens,dens_up,dens_dw
  real(8),dimension(:,:),allocatable   :: docc
  real(8),dimension(:,:),allocatable   :: magz
  real(8),dimension(:,:,:),allocatable :: sz2
  real(8),dimension(:,:,:),allocatable :: zimp,simp
  real(8)                              :: Egs
  real(8)                              :: s2tot
  real(8)                              :: Ei
  !
  integer                            :: iorb,jorb,iorb1,jorb1
  integer                            :: ispin,jspin
  integer                            :: isite,jsite
  integer                            :: ibath
  integer                            :: r,m,k,k1,k2
  integer                            :: iup,idw
  integer                            :: jup,jdw
  integer                            :: mup,mdw
  real(8)                            :: sgn,sgn1,sgn2,sg1,sg2
  real(8)                            :: gs_weight
  !
  real(8)                            :: peso
  real(8)                            :: norm
  !
  integer                            :: i,j,ii
  integer                            :: isector,jsector
  integer                            :: idim,idimUP,idimDW
  !
  real(8),dimension(:),pointer       :: gscvec
  logical                            :: Jcondition



contains 

!+-------------------------------------------------------------------+
 !PURPOSE  : Evaluate and print out many interesting physical qties
!+-------------------------------------------------------------------+

  subroutine observables_cluster()
    integer :: iorb
    write(LOGfile,"(A)")"Get local observables:"
    !select case(vca_method)
    !case ("full")
    !   call vca_diag_observables
    !case ("lanc")
       call vca_lanc_observables
    !case default
    !   stop "observables_cluster error: vca_method != ['full','lanc']"
    !end select
  end subroutine observables_cluster

  !+-------------------------------------------------------------------+
  !PURPOSE  : Evaluate and print out many interesting physical qties
  !+-------------------------------------------------------------------+
  subroutine vca_lanc_observables()
    integer,dimension(Ns)               :: IbUp,IbDw  ![Ns]
    integer,dimension(Ns_Ud,Ns_Orb)     :: Nups,Ndws  ![1,Ns]-[Norb,1+Nbath]
    integer,dimension(2*Ns_Ud)          :: Indices,Jndices
    integer,dimension(Ns_Ud)            :: iDimUps,iDimDws
    integer                             :: i,j,ii
    integer                             :: istate
    integer                             :: isector,jsector
    integer                             :: idim,jdim,ilat
    integer                             :: isz,jsz
    integer                             :: iorb,jorb,ispin,jspin,isite,jsite,ibath
    integer                             :: numstates
    integer                             :: r,m,k
    real(8)                             :: sgn,sgn1,sgn2
    real(8)                             :: weight
    real(8)                             :: Ei
    real(8)                             :: peso
    real(8)                             :: norm
    real(8),dimension(Nlat,Norb)        :: nup,ndw,Sz,nt
    real(8),dimension(:),pointer        :: gscvec
    type(sector_map)                    :: Hi(2*Ns_Ud)
    !
    !LOCAL OBSERVABLES:
    allocate(dens(Nlat,Norb),dens_up(Nlat,Norb),dens_dw(Nlat,Norb))
    allocate(docc(Nlat,Norb))
    allocate(magz(Nlat,Norb),sz2(Nlat,Norb,Norb))
    allocate(simp(Nlat,Norb,Nspin),zimp(Nlat,Norb,Nspin))
    !
    dens    = 0.d0
    dens_up = 0.d0
    dens_dw = 0.d0
    docc    = 0.d0
    magz    = 0.d0
    sz2     = 0.d0
    simp    = 0.d0
    zimp    = 0.d0
    Egs     = state_list%emin
    !
    do istate=1,state_list%size
       isector = es_return_sector(state_list,istate)
       Ei      = es_return_energy(state_list,istate)
       !
#ifdef _MPI
       if(MpiStatus)then
          gscvec => es_return_cvector(MpiComm,state_list,istate)
       else
          gscvec => es_return_cvector(state_list,istate)
       endif
#else
       gscvec => es_return_cvector(state_list,istate)
#endif
       !
       idim    = getdim(isector)
       call get_DimUp(isector,iDimUps)
       call get_DimDw(isector,iDimDws)
       iDimUp = product(iDimUps)
       iDimDw = product(iDimDws)
       !norm=sqrt(dot_product(gscvec,gscvec))
       !
       !if(abs(norm-1.d0)>1.d-9)stop "GS is not normalized"
       !
       peso = 1.d0 ; if(finiteT)peso=exp(-beta*(Ei-Egs))
       peso = peso/zeta_function
       !
       if(MpiMaster)then 
         call build_sector(isector,Hi)
         do i=1,idim
            call state2indices(i,[iDimUps,iDimDws],Indices)
            do ii=1,Ns_Ud
              mup = HI(ii)%map(Indices(ii))
              mdw = HI(ii+Ns_Ud)%map(Indices(ii+Ns_ud))
              Nups(ii,:) = Bdecomp(mup,Ns_Orb)
              Ndws(ii,:) = Bdecomp(mdw,Ns_Orb)
           enddo
           IbUp = Breorder(Nups)
           IbDw = Breorder(Ndws)
           !
           gs_weight=peso*abs(gscvec(i))**2
           !
           !Get operators:
           do ilat=1,Nlat
             do iorb=1,Norb
                nup(ilat,iorb)= ibup(imp_state_index(ilat,iorb,1))
                ndw(ilat,iorb)= ibdw(imp_state_index(ilat,iorb,1))
                sz(ilat,iorb) = (nup(ilat,iorb) - ndw(ilat,iorb))/2d0
                nt(ilat,iorb) =  nup(ilat,iorb) + ndw(ilat,iorb)
             enddo
           enddo
            !
            !Evaluate averages of observables:
           do ilat=1,Nlat
             do iorb=1,Norb
                dens(ilat,iorb)     = dens(ilat,iorb)      +  nt(ilat,iorb)*gs_weight
                dens_up(ilat,iorb)  = dens_up(ilat,iorb)   +  nup(ilat,iorb)*gs_weight
                dens_dw(ilat,iorb)  = dens_dw(ilat,iorb)   +  ndw(ilat,iorb)*gs_weight
                docc(ilat,iorb)     = docc(ilat,iorb)      +  nup(ilat,iorb)*ndw(ilat,iorb)*gs_weight
                magz(ilat,iorb)     = magz(ilat,iorb)      +  (nup(ilat,iorb)-ndw(ilat,iorb))*gs_weight
                sz2(ilat,iorb,iorb) = sz2(ilat,iorb,iorb)  +  (sz(ilat,iorb)*sz(ilat,iorb))*gs_weight
                do jorb=iorb+1,Norb
                   sz2(ilat,iorb,jorb) = sz2(ilat,iorb,jorb)  +  (sz(ilat,iorb)*sz(ilat,jorb))*gs_weight
                   sz2(ilat,jorb,iorb) = sz2(ilat,jorb,iorb)  +  (sz(ilat,jorb)*sz(ilat,iorb))*gs_weight
                enddo
             enddo
             s2tot = s2tot  + (sum(sz))**2*gs_weight
           enddo
         enddo
         call delete_sector(isector,HI)
         if(associated(gscvec))nullify(gscvec)
       endif
    enddo
    if(MPIMASTER)then
        call get_szr
        if(iolegend)call write_legend
        call write_observables()
        !
        write(LOGfile,"(A,10f18.12,f18.12,A)")"dens"//reg(file_suffix)//"=",(sum(dens(:,iorb))/Nlat,iorb=1,Norb),sum(dens)/Nlat
        write(LOGfile,"(A,10f18.12,A)")"docc"//reg(file_suffix)//"=",(sum(docc(:,iorb))/Nlat,iorb=1,Norb)
        if(Nspin==2)write(LOGfile,"(A,10f18.12,A)") "mag "//reg(file_suffix)//"=",(sum(magz(:,iorb))/Nlat,iorb=1,Norb)
        !
        imp_dens_up = dens_up
        imp_dens_dw = dens_dw
        imp_dens    = dens
        imp_docc    = docc
    endif
#ifdef _MPI
    if(MpiStatus)then
       call Bcast_MPI(MpiComm,dens_up)
       call Bcast_MPI(MpiComm,dens_dw)
       call Bcast_MPI(MpiComm,dens)
       call Bcast_MPI(MpiComm,docc)
    endif
#endif
    deallocate(dens,docc,dens_up,dens_dw,magz,sz2)
    deallocate(simp,zimp)
  end subroutine vca_lanc_observables








  !####################################################################
  !                    COMPUTATIONAL ROUTINES
  !####################################################################
  !+-------------------------------------------------------------------+
  !PURPOSE  : get scattering rate and renormalization constant Z
  !+-------------------------------------------------------------------+
  subroutine get_szr()
    integer                  :: ilat,ispin,iorb
    real(8)                  :: wm1,wm2
    wm1 = pi/beta ; wm2=3d0*pi/beta
    do ilat=1,Nlat
       do ispin=1,Nspin
          do iorb=1,Norb
             simp(ilat,iorb,ispin) = dimag(impSmats(ilat,ilat,ispin,ispin,iorb,iorb,1)) - &
                  wm1*(dimag(impSmats(ilat,ilat,ispin,ispin,iorb,iorb,2))-&
                  dimag(impSmats(ilat,ilat,ispin,ispin,iorb,iorb,1)))/(wm2-wm1)
             zimp(ilat,iorb,ispin)   = 1.d0/( 1.d0 + abs( dimag(impSmats(ilat,ilat,ispin,ispin,iorb,iorb,1))/wm1 ))
          enddo
       enddo
    enddo
  end subroutine get_szr



  !+-------------------------------------------------------------------+
  !PURPOSE  : write legend, i.e. info about columns 
  !+-------------------------------------------------------------------+
  subroutine write_legend()
    integer :: unit,iorb,jorb,ispin
    unit = free_unit()
    open(unit,file="observables_info.vca")
    write(unit,"(A1,90(A10,6X))")"#",&
         (reg(txtfy(iorb))//"dens_"//reg(txtfy(iorb)),iorb=1,Norb),&
         (reg(txtfy(Norb+iorb))//"docc_"//reg(txtfy(iorb)),iorb=1,Norb),&
         (reg(txtfy(2*Norb+iorb))//"nup_"//reg(txtfy(iorb)),iorb=1,Norb),&
         (reg(txtfy(3*Norb+iorb))//"ndw_"//reg(txtfy(iorb)),iorb=1,Norb),&
         (reg(txtfy(4*Norb+iorb))//"mag_"//reg(txtfy(iorb)),iorb=1,Norb),&
         reg(txtfy(5*Norb+1))//"egs",&
         ((reg(txtfy(5*Norb+1+(iorb-1)*Norb+jorb))//"sz2_"//reg(txtfy(iorb))//reg(txtfy(jorb)),jorb=1,Norb),iorb=1,Norb),&
         ((reg(txtfy((6+2*Norb)*Norb+2+(ispin-1)*Nspin+iorb))//"z_"//reg(txtfy(iorb))//"s"//reg(txtfy(ispin)),iorb=1,Norb),ispin=1,Nspin),&
         ((reg(txtfy((7+2*Norb)*Norb+2+Nspin+(ispin-1)*Nspin+iorb))//"sig_"//reg(txtfy(iorb))//"s"//reg(txtfy(ispin)),iorb=1,Norb),ispin=1,Nspin)
    close(unit)
    !
    unit = free_unit()
    open(unit,file="parameters_info.vca")
    write(unit,"(A1,90(A14,1X))")"#","1xmu","2beta",&
         (reg(txtfy(2+iorb))//"U_"//reg(txtfy(iorb)),iorb=1,Norb),&
         reg(txtfy(2+Norb+1))//"U'",reg(txtfy(2+Norb+2))//"Jh"
    close(unit)
    !
    iolegend=.false.
  end subroutine write_legend


  !+-------------------------------------------------------------------+
  !PURPOSE  : write observables to file
  !+-------------------------------------------------------------------+
  subroutine write_observables()
    integer :: unit
    integer :: ilat,iorb,jorb,ispin
    unit = free_unit()     
    open(unit,file="parameters_last"//reg(file_suffix)//".vca")
    write(unit,"(90F15.9)")xmu,beta,(uloc(iorb),iorb=1,Norb),Ust,Jh,Jx,Jp
    close(unit)
    !
    do ilat=1,Nlat
       unit = free_unit()
       open(unit,file="observables_last"//reg(file_suffix)//"_site"//str(ilat,3)//".vca")
       write(unit,"(90(F15.9,1X))")&
            (dens(ilat,iorb),iorb=1,Norb),&
            (docc(ilat,iorb),iorb=1,Norb),&
            (dens_up(ilat,iorb),iorb=1,Norb),&
            (dens_dw(ilat,iorb),iorb=1,Norb),&
            (magz(ilat,iorb),iorb=1,Norb),&
            egs,&
            ((sz2(ilat,iorb,jorb),jorb=1,Norb),iorb=1,Norb),&
            ((zimp(ilat,iorb,ispin),iorb=1,Norb),ispin=1,Nspin),&
            ((simp(ilat,iorb,ispin),iorb=1,Norb),ispin=1,Nspin)
       close(unit)
    enddo
  end subroutine write_observables


END MODULE VCA_OBSERVABLES









