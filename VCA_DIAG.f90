module VCA_DIAG
  USE VCA_VARS_GLOBAL
  USE VCA_SETUP
  !
  USE SF_CONSTANTS
  USE SF_LINALG, only: eigh
  USE SF_TIMER,  only: start_timer,stop_timer,eta
  USE SF_IOTOOLS, only:reg,free_unit
  USE SF_MISC, only: assert_shape
  !
  implicit none
  private


  integer                      :: Hsector=0
  logical                      :: Hstatus=.false.
  type(sector_map)             :: H,Hup,Hdw


  !>build sparse hamiltonian of the sector
  public  :: buildH_c

  !>auxiliary routines
  public  :: setup_Hv_sector
  public  :: delete_Hv_sector

  !>diagonalize cluster Hamiltonian
  public  :: diagonalize_cluster



contains




  !> Set the actual Sector
  subroutine setup_Hv_sector(isector)
    integer                   :: isector
    Hsector=isector
    Hstatus=.true.
    call build_sector(isector,H)
  end subroutine setup_Hv_sector

  !> Reset the sector 
  subroutine delete_Hv_sector()
    call delete_sector(Hsector,H)
    Hsector=0
    Hstatus=.false.
  end subroutine delete_Hv_sector



  !+-------------------------------------------------------------------+
  !PURPOSE  : diagonalize the Hamiltonian in each sector and find the 
  ! spectrum DOUBLE COMPLEX
  !+------------------------------------------------------------------+
  subroutine diagonalize_cluster
    integer                     :: nup,ndw,isector,dim
    integer                     :: sz,nt
    integer                     :: i,j,unit
    real(8),dimension(Nsectors) :: e0 
    real(8)                     :: egs
    logical                     :: Tflag
    !
    e0=1000.d0
    write(LOGfile,"(A)")"Diagonalize Cluster Hc+Hint:"
    call start_timer()
    !
    !
    sector: do isector=1,Nsectors
       !
       Dim      = getdim(isector)
       !
       if(verbose==3)then
          nup  = getnup(isector)
          ndw  = getndw(isector)
          write(LOGfile,"(A,I4,A6,I2,A6,I2,A6,I15)")"Solving sector:",isector,", nup:",nup,", ndw:",ndw,", dim=",getdim(isector)
       elseif(verbose==1.OR.verbose==2)then
          call eta(isector,Nsectors,LOGfile)
       endif
       !
       call setup_Hv_sector(isector)
       call buildH_c(espace(isector)%M)
       call delete_Hv_sector()
       call eigh(espace(isector)%M,espace(isector)%e,'V','U')
       if(dim==1)espace(isector)%M=one
       !
       e0(isector)=minval(espace(isector)%e)
       !
    enddo sector
    !
    call stop_timer(LOGfile)
    !
    !Get the ground state energy and rescale energies
    egs=minval(e0)
    forall(isector=1:Nsectors)espace(isector)%e = espace(isector)%e - egs
    !
    !Get the partition function Z
    zeta_function=0.d0;zeta_function=0.d0
    do isector=1,Nsectors
       dim=getdim(isector)
       do i=1,dim
          zeta_function=zeta_function+exp(-beta*espace(isector)%e(i))
       enddo
    enddo
    !
    write(LOGfile,"(A)")"DIAG summary:"
    open(free_unit(unit),file='egs'//reg(file_suffix)//".vca",position='append')
    do isector=1,Nsectors
       if(e0(isector)/=0d0)cycle
       nup  = getnup(isector)
       ndw  = getndw(isector)
       write(LOGfile,"(A,F20.12,2I4)")'Egs =',e0(isector),nup,ndw
       write(unit,"(F20.12,2I4)")e0(isector),nup,ndw
    enddo
    write(LOGfile,"(A,F20.12)")'Z   =',zeta_function
    close(unit)
    return
  end subroutine diagonalize_cluster





  !> Build the Cluster Hamiltonian (into Hmat)
  subroutine buildH_c(Hmat)
    complex(8),dimension(:,:)    :: Hmat
    integer                      :: isector
    integer,dimension(Nlevels)   :: ib
    real(8),dimension(Nlat,Norb) :: nup,ndw
    integer                      :: dim,dimUp,dimDw
    integer                      :: i
    integer                      :: m
    integer                      :: ishift,ishift_up,ishift_dw
    integer                      :: j,ms,impi
    integer                      :: ilat,jlat,iorb,jorb,ispin,jspin,is,js
    integer                      :: i_up,i_dw,j_up,j_dw
    integer                      :: kp,k1,k2,k3,k4
    integer                      :: alfa,beta
    real(8)                      :: sg1,sg2,sg3,sg4
    complex(8)                   :: htmp,htmpup,htmpdw
    logical                      :: Jcondition
    integer                      :: first_state,last_state
    integer                      :: first_state_up,last_state_up
    integer                      :: first_state_dw,last_state_dw
    !
    if(.not.Hstatus)stop "buildH_c ERROR: Hsector NOT set"
    isector=Hsector
    !
    dim=getdim(isector)
    !
    call assert_shape(Hmat,[dim,dim],"buildH_c","Hmat")
    !
    Hmat=zero
    !
    !-----------------------------------------------!
    states: do i=1,Dim
       m = H%map(i)
       impi = i
       ib = bdecomp(m,2*Ns)
       !
       do ilat=1,Nlat
          do iorb=1,Norb
             nup(ilat,iorb)=dble(ib(state_index_up(ilat,iorb)))
             ndw(ilat,iorb)=dble(ib(state_index_dw(ilat,iorb)))
          enddo
       enddo
       !
       !CLUSTER  HAMILTONIAN
       !Diagonal Elements, i.e. local part
       htmp = zero
       htmp = htmp - xmu*(sum(nup)+sum(ndw))
       !
       do ilat=1,Nlat
          do iorb=1,Norb
             htmp = htmp + impHloc(ilat,ilat,1,1,iorb,iorb)*nup(ilat,iorb)
             htmp = htmp + impHloc(ilat,ilat,Nspin,Nspin,iorb,iorb)*ndw(ilat,iorb)
          enddo
       enddo
       !
       Hmat(i,i)=Hmat(i,i) + htmp
       !
       !
       !Off-diagonal elements, i.e. non-local part
       !this loop considers only the site-orbitals off-diagonal terms
       !because iorb=jorb can not have simultaneously
       !occupation 0 and 1, as required by this if Jcondition:  
       do ilat=1,Nlat
          do jlat=1,Nlat
             do iorb=1,Norb
                do jorb=1,Norb
                   !UP
                   ispin=1
                   jspin=1
                   !
                   is = state_index_up(ilat,iorb)
                   js = state_index_up(jlat,jorb)
                   !
                   Jcondition = (impHloc(ilat,jlat,ispin,jspin,iorb,jorb)/=zero).AND.(ib(js)==1).AND.(ib(is)==0)
                   if (Jcondition) then
                      call c(js,m,k1,sg1)
                      call cdg(is,k1,k2,sg2)
                      j = binary_search(H%map,k2)
                      htmp = impHloc(ilat,jlat,ispin,jspin,iorb,jorb)*sg1*sg2
                      !
                      Hmat(i,j) = Hmat(i,j) + htmp
                      !
                   endif
                   !DW
                   ispin=Nspin
                   jspin=Nspin
                   !
                   is = state_index_dw(ilat,iorb)
                   js = state_index_dw(jlat,jorb)
                   !
                   Jcondition = (impHloc(ilat,jlat,ispin,jspin,iorb,jorb)/=zero).AND.(ib(js)==1).AND.(ib(is)==0)
                   if (Jcondition) then
                      call c(js,m,k1,sg1)
                      call cdg(is,k1,k2,sg2)
                      j = binary_search(H%map,k2)
                      htmp = impHloc(ilat,jlat,Nspin,Nspin,iorb,jorb)*sg1*sg2
                      !
                      Hmat(i,j) = Hmat(i,j) + htmp
                      !
                   endif
                enddo
             enddo
          enddo
       enddo
       !
       !LOCAL INTERACTION
       !density-density interaction: same orbital, opposite spins:
       ! = \sum_\a U_\a*(n_{\a,up}*n_{\a,dw})
       htmp = zero
       do ilat=1,Nlat
          do iorb=1,Norb
             htmp = htmp + Uloc(iorb)*nup(ilat,iorb)*ndw(ilat,iorb)
          enddo
       enddo
       if(Norb>1)then
          !density-density interaction: different orbitals, opposite spins:
          ! =   U'   *     sum_{i/=j} [ n_{i,up}*n_{j,dw} + n_{j,up}*n_{i,dw} ]
          ! =  (Uloc-2*Jh)*sum_{i/=j} [ n_{i,up}*n_{j,dw} + n_{j,up}*n_{i,dw} ]
          do ilat=1,Nlat
             do iorb=1,Norb
                do jorb=iorb+1,Norb
                   htmp = htmp + Ust*(nup(ilat,iorb)*ndw(ilat,jorb) + nup(ilat,jorb)*ndw(ilat,iorb))
                enddo
             enddo
          enddo
          !density-density interaction: different orbitals, parallel spins
          ! = \sum_{i<j}    U''     *[ n_{i,up}*n_{j,up} + n_{i,dw}*n_{j,dw} ]
          ! = \sum_{i<j} (Uloc-3*Jh)*[ n_{i,up}*n_{j,up} + n_{i,dw}*n_{j,dw} ]
          do ilat=1,Nlat
             do iorb=1,Norb
                do jorb=iorb+1,Norb
                   htmp = htmp + (Ust-Jh)*(nup(ilat,iorb)*nup(ilat,jorb) + ndw(ilat,iorb)*ndw(ilat,jorb))
                enddo
             enddo
          enddo
       endif
       !if using the Hartree-shifted chemical potential: mu=0 for half-filling
       !sum up the contributions of hartree terms:
       if(hfmode)then
          do ilat=1,Nlat
             do iorb=1,Norb
                htmp = htmp - 0.5d0*Uloc(iorb)*(nup(ilat,iorb)+ndw(ilat,iorb)) + 0.25d0*Uloc(iorb)
             enddo
          enddo
          if(Norb>1)then
             do ilat=1,Nlat
                do iorb=1,Norb
                   do jorb=iorb+1,Norb
                      htmp=htmp-0.5d0*Ust*(nup(ilat,iorb)+ndw(ilat,iorb)+nup(ilat,jorb)+ndw(ilat,jorb))+0.25d0*Ust
                      htmp=htmp-0.5d0*(Ust-Jh)*(nup(ilat,iorb)+ndw(ilat,iorb)+nup(ilat,jorb)+ndw(ilat,jorb))+0.25d0*(Ust-Jh)
                   enddo
                enddo
             enddo
          endif
       endif
       !
       Hmat(i,i) = Hmat(i,i) + htmp
       !
       !
       !
       ! SPIN-EXCHANGE (S-E) and PAIR-HOPPING TERMS
       !    S-E: J c^+_iorb_up c^+_jorb_dw c_iorb_dw c_jorb_up  (i.ne.j) 
       !    S-E: J c^+_{iorb} c^+_{jorb+Ns} c_{iorb+Ns} c_{jorb}
       if(Norb>1.AND.Jhflag)then
          do ilat=1,Nlat
             do iorb=1,Norb
                do jorb=1,Norb
                   i_up = state_index_up(ilat,iorb)
                   i_dw = state_index_dw(ilat,iorb)
                   j_up = state_index_up(ilat,jorb)
                   j_dw = state_index_dw(ilat,jorb)
                   Jcondition=(&
                        (iorb/=jorb).AND.&
                        (ib(j_up)==1).AND.&
                        (ib(i_dw)==1).AND.&
                        (ib(j_dw)==0).AND.&
                        (ib(i_up)==0))
                   if(Jcondition)then
                      call c(j_up,m,k1,sg1)
                      call c(i_dw,k1,k2,sg2)
                      call cdg(j_dw,k2,k3,sg3)
                      call cdg(i_up,k3,k4,sg4)
                      j=binary_search(H%map,k4)
                      htmp = one*Jx*sg1*sg2*sg3*sg4
                      !
                      if(j/=0)Hmat(i,j) = Hmat(i,j) + htmp
                      !
                   endif
                enddo
             enddo
          enddo
       endif
       !
       ! PAIR-HOPPING (P-H) TERMS
       !    P-H: J c^+_iorb_up c^+_iorb_dw   c_jorb_dw   c_jorb_up  (i.ne.j) 
       !    P-H: J c^+_{iorb}  c^+_{iorb+Ns} c_{jorb+Ns} c_{jorb}
       if(Norb>1.AND.Jhflag)then
          do ilat=1,Nlat
             do iorb=1,Norb
                do jorb=1,Norb
                   i_up = state_index_up(ilat,iorb)
                   i_dw = state_index_dw(ilat,iorb)
                   j_up = state_index_up(ilat,jorb)
                   j_dw = state_index_dw(ilat,jorb)
                   Jcondition=(&
                        (iorb/=jorb).AND.&
                        (ib(j_up)==1).AND.&
                        (ib(j_dw)==1).AND.&
                        (ib(i_dw)==0).AND.&
                        (ib(i_up)==0))
                   if(Jcondition)then
                      call c(j_up,m,k1,sg1)
                      call c(j_dw,k1,k2,sg2)
                      call cdg(i_dw,k2,k3,sg3)
                      call cdg(i_up,k3,k4,sg4)
                      j=binary_search(H%map,k4)
                      htmp = one*Jp*sg1*sg2*sg3*sg4
                      !
                      if(j/=0)Hmat(i,j) = Hmat(i,j) + htmp
                      !
                   endif
                enddo
             enddo
          enddo
       endif
    enddo states
  end subroutine buildH_c










END MODULE VCA_DIAG









