MODULE VCA_OBSERVABLES
  USE VCA_VARS_GLOBAL
  USE VCA_SETUP
  USE VCA_DIAG
  USE VCA_AUX_FUNX
  USE VCA_EIGENSPACE
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
  real(8),dimension(:,:,:),allocatable   :: zimp,simp
  real(8)                              :: Egs



contains 

  subroutine observables_cluster()
    integer :: iorb
    write(LOGfile,"(A)")"Get local observables:"
    select case(vca_method)
    case ("full")
       call vca_diag_observables
    case ("lanc")
       call vca_lanc_observables
    case default
       stop "observables_cluster error: vca_method != ['full','lanc']"
    end select
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
    !
    deallocate(dens,docc,dens_up,dens_dw,magz,sz2)
    deallocate(simp,zimp)
  end subroutine observables_cluster






  !+-------------------------------------------------------------------+
  !PURPOSE  : Evaluate and print out many interesting physical qties
  !+-------------------------------------------------------------------+
  subroutine vca_diag_observables()
    integer,dimension(Nlevels)      :: ib
    integer                         :: i,j,ilat,jlat
    integer                         :: izero,istate
    integer                         :: isector,jsector
    integer                         :: idim,jdim
    integer                         :: isz,jsz
    integer                         :: iorb,jorb,ispin,jspin,isite,jsite
    integer                         :: numstates
    integer                         :: r,m,k
    real(8)                         :: sgn,sgn1,sgn2
    real(8)                         :: boltzman_weight
    real(8)                         :: state_weight
    real(8)                         :: weight
    real(8)                         :: Ei
    real(8)                         :: norm
    real(8),dimension(Nlat,Norb)    :: nup,ndw,Sz,nt
    type(sector_map)                :: H,HJ
    complex(8),dimension(:),pointer :: evec
    !
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
    !
    do isector=1,Nsectors
       idim    = getdim(isector)
       call build_sector(isector,H)
       !
       do istate=1,idim
          Ei=espace(isector)%e(istate)
          boltzman_weight=exp(-beta*Ei)/zeta_function
          if(boltzman_weight < cutoff)cycle
          !
          evec => espace(isector)%M(:,istate)
          !
          do i=1,idim
             m=H%map(i)
             ib = bdecomp(m,2*Ns)
             !
             state_weight=evec(i)*evec(i)
             weight = boltzman_weight*state_weight
             !
             !Get operators:
             do ilat=1,Nlat
                do iorb=1,Norb
                   nup(ilat,iorb)= ib(imp_state_index(ilat,iorb,1))
                   ndw(ilat,iorb)= ib(imp_state_index(ilat,iorb,2))
                   sz(ilat,iorb) = (nup(ilat,iorb) - ndw(ilat,iorb))/2d0
                   nt(ilat,iorb) =  nup(ilat,iorb) + ndw(ilat,iorb)
                enddo
             enddo
             !
             !Evaluate averages of observables:
             do ilat=1,Nlat
                do iorb=1,Norb
                   dens(ilat,iorb)     = dens(ilat,iorb)      +  nt(ilat,iorb)*weight
                   dens_up(ilat,iorb)  = dens_up(ilat,iorb)   +  nup(ilat,iorb)*weight
                   dens_dw(ilat,iorb)  = dens_dw(ilat,iorb)   +  ndw(ilat,iorb)*weight
                   docc(ilat,iorb)     = docc(ilat,iorb)      +  nup(ilat,iorb)*ndw(ilat,iorb)*weight
                   magz(ilat,iorb)     = magz(ilat,iorb)      +  (nup(ilat,iorb)-ndw(ilat,iorb))*weight
                   sz2(ilat,iorb,iorb) = sz2(ilat,iorb,iorb)  +  (sz(ilat,iorb)*sz(ilat,iorb))*weight
                   do jorb=iorb+1,Norb
                      sz2(ilat,iorb,jorb) = sz2(ilat,iorb,jorb)  +  (sz(ilat,iorb)*sz(ilat,jorb))*weight
                      sz2(ilat,jorb,iorb) = sz2(ilat,jorb,iorb)  +  (sz(ilat,jorb)*sz(ilat,iorb))*weight
                   enddo
                enddo
             enddo
             !
          enddo
       enddo
       deallocate(H%map)
    enddo
    !
  end subroutine vca_diag_observables












  !+-------------------------------------------------------------------+
  !PURPOSE  : Evaluate and print out many interesting physical qties
  !+-------------------------------------------------------------------+
  subroutine vca_lanc_observables()
    integer,dimension(Nlevels)      :: ib
    integer                         :: i,j
    integer                         :: izero
    integer                         :: isector,jsector
    integer                         :: idim,jdim,ilat
    integer                         :: isz,jsz
    integer                         :: iorb,jorb,ispin,jspin,isite,jsite,ibath
    integer                         :: numstates
    integer                         :: r,m,k
    real(8)                         :: sgn,sgn1,sgn2
    real(8)                         :: weight
    real(8)                         :: Ei
    real(8)                         :: peso
    real(8)                         :: norm
    real(8),dimension(Nlat,Norb)    :: nup,ndw,Sz,nt
    complex(8),dimension(:),pointer :: gscvec
    type(sector_map)                :: H,HJ
    complex(8),allocatable          :: vvinit(:)
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
    numstates=state_list%size
    do izero=1,numstates
       isector = es_return_sector(state_list,izero)
       Ei      = es_return_energy(state_list,izero)
       idim    = getdim(isector)
       !
       gscvec  => es_return_cvector(state_list,izero)
       norm=sqrt(dot_product(gscvec,gscvec))
       if(abs(norm-1.d0)>1.d-9)stop "GS is not normalized"
       !
       peso = 1.d0 ; if(finiteT)peso=exp(-beta*(Ei-Egs))
       peso = peso/zeta_function
       !
       call build_sector(isector,H)
       !
       !pdens=0d0
       do i=1,idim
          m=H%map(i)
          ib = bdecomp(m,2*Ns)
          !
          weight=peso*abs(gscvec(i))**2
          !
          !Get operators:
          !Get operators:
          do ilat=1,Nlat
             do iorb=1,Norb
                nup(ilat,iorb)= ib(imp_state_index(ilat,iorb,1))
                ndw(ilat,iorb)= ib(imp_state_index(ilat,iorb,2))
                sz(ilat,iorb) = (nup(ilat,iorb) - ndw(ilat,iorb))/2d0
                nt(ilat,iorb) =  nup(ilat,iorb) + ndw(ilat,iorb)
             enddo
          enddo

          !
          !Evaluate averages of observables:
          do ilat=1,Nlat
             do iorb=1,Norb
                dens(ilat,iorb)     = dens(ilat,iorb)      +  nt(ilat,iorb)*weight
                dens_up(ilat,iorb)  = dens_up(ilat,iorb)   +  nup(ilat,iorb)*weight
                dens_dw(ilat,iorb)  = dens_dw(ilat,iorb)   +  ndw(ilat,iorb)*weight
                docc(ilat,iorb)     = docc(ilat,iorb)      +  nup(ilat,iorb)*ndw(ilat,iorb)*weight
                magz(ilat,iorb)     = magz(ilat,iorb)      +  (nup(ilat,iorb)-ndw(ilat,iorb))*weight
                sz2(ilat,iorb,iorb) = sz2(ilat,iorb,iorb)  +  (sz(ilat,iorb)*sz(ilat,iorb))*weight
                do jorb=iorb+1,Norb
                   sz2(ilat,iorb,jorb) = sz2(ilat,iorb,jorb)  +  (sz(ilat,iorb)*sz(ilat,jorb))*weight
                   sz2(ilat,jorb,iorb) = sz2(ilat,jorb,iorb)  +  (sz(ilat,jorb)*sz(ilat,iorb))*weight
                enddo
             enddo
          enddo

       enddo
       if(associated(gscvec))nullify(gscvec)
       deallocate(H%map)
    enddo
    !
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









