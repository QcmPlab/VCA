MODULE VCA_OBSERVABLES
  USE VCA_VARS_GLOBAL
  USE VCA_SETUP
  USE VCA_DIAG
  USE VCA_AUX_FUNX
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
  real(8)                              :: Egs



contains 





  !+-------------------------------------------------------------------+
  !PURPOSE  : Evaluate and print out many interesting physical qties
  !+-------------------------------------------------------------------+
  subroutine observables_cluster()
    integer,dimension(Nlevels)   :: ib
    integer                      :: i,j,ilat,jlat
    integer                      :: izero,istate
    integer                      :: isector,jsector
    integer                      :: idim,jdim
    integer                      :: isz,jsz
    integer                      :: iorb,jorb,ispin,jspin,isite,jsite
    integer                      :: numstates
    integer                      :: r,m,k
    real(8)                      :: sgn,sgn1,sgn2
    real(8)                      :: boltzman_weight
    real(8)                      :: state_weight
    real(8)                      :: weight
    real(8)                      :: Ei
    real(8)                      :: norm
    real(8),dimension(Nlat,Norb) :: nup,ndw,Sz,nt
    type(sector_map)             :: H,HJ
    real(8),allocatable          :: vvinit(:)
    real(8),dimension(:),pointer :: evec
    !
    !
    !LOCAL OBSERVABLES:
    allocate(dens(Nlat,Norb),dens_up(Nlat,Norb),dens_dw(Nlat,Norb))
    allocate(docc(Nlat,Norb))
    allocate(magz(Nlat,Norb),sz2(Nlat,Norb,Norb))
    !
    dens    = 0.d0
    dens_up = 0.d0
    dens_dw = 0.d0
    docc    = 0.d0
    magz    = 0.d0
    sz2     = 0.d0
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
                   nup(ilat,iorb)= ib(state_index(ilat,iorb,1))
                   ndw(ilat,iorb)= ib(state_index(ilat,iorb,2))
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
  end subroutine observables_cluster






  !####################################################################
  !                    COMPUTATIONAL ROUTINES
  !####################################################################



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
         ((reg(txtfy(5*Norb+1+(iorb-1)*Norb+jorb))//"sz2_"//reg(txtfy(iorb))//reg(txtfy(jorb)),jorb=1,Norb),iorb=1,Norb)
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
            ((sz2(ilat,iorb,jorb),jorb=1,Norb),iorb=1,Norb)
       close(unit)
    enddo
  end subroutine write_observables


end MODULE VCA_OBSERVABLES









