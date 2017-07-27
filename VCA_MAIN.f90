module VCA_MAIN
  USE VCA_VARS_GLOBAL
  USE VCA_AUX_FUNX
  USE VCA_SETUP
  USE VCA_DIAG
  USE VCA_GREENS_FUNCTIONS
  USE VCA_OBSERVABLES
  !
  USE SF_LINALG,  only: eigh,diag
  USE SF_IOTOOLS, only: str
  USE SF_TIMER,only: start_timer,stop_timer
  USE SF_MISC, only: assert_shape
  implicit none
  private

  !>INIT VCA SOLVER
  public :: vca_init_solver


  !> VCA DIAG CLUSTER
  public :: vca_diag

  !> VCA GF POLES UPDATE
  public :: vca_update_poles

contains



  !+-----------------------------------------------------------------------------+!
  ! PURPOSE: allocate and initialize one or multiple baths -+!
  !+-----------------------------------------------------------------------------+!
  subroutine vca_init_solver(Hloc)
    real(8),intent(in)              :: Hloc(:,:,:,:,:,:) !Hloc(Nlat,Nlat,Nspin,Nspin,Norb,Norb)
    logical,save                    :: isetup=.true.
    !
    write(LOGfile,"(A)")"INIT SOLVER FOR "//trim(file_suffix)
    !
    Nlat = size(Hloc,1)
    call assert_shape(Hloc,[Nlat,Nlat,Nspin,Nspin,Norb,Norb],"vca_init_solver","Hloc")
    !
    !Init Structure & memory
    if(isetup)call init_cluster_structure()
    !
    !Init Hcluster
    call set_Hcluster(Hloc)
    !
    if(isetup)call setup_pointers_normal
    isetup=.false.
    !
  end subroutine vca_init_solver


  !+-----------------------------------------------------------------------------+!
  !PURPOSE: Diag the cluster, reference system
  !+-----------------------------------------------------------------------------+!
  subroutine vca_diag(Hloc)
    real(8),optional,intent(in) :: Hloc(:,:,:,:,:,:) !Hloc(Nlat,Nlat,Nspin,Nspin,Norb,Norb)
    !
    Nlat = size(Hloc,1)
    call assert_shape(Hloc,[Nlat,Nlat,Nspin,Nspin,Norb,Norb],"vca_diag","Hloc")
    !
    if(present(Hloc))call set_Hcluster(Hloc)
    !
    call setup_eigenspace()
    !
    !SOLVE THE QUANTUM IMPURITY PROBLEM:
    call diagonalize_cluster()         !find target states by digonalization of Hamiltonian
    call observables_cluster()         !obtain impurity observables as thermal averages.  
    call build_QL_cluster()            !build the one-particle Q and Lambda matrices, storing spectral sums
    !
    call delete_eigenspace()
  end subroutine vca_diag






  subroutine vca_update_poles(Vmat)
    real(8),dimension(Nlat*Nspin*Norb,Nlat*Nspin*Norb) :: Vmat
    integer                                            :: ilat,jlat
    integer                                            :: iorb,jorb
    integer                                            :: ispin
    integer                                            :: iexc,jexc
    integer                                            :: i,j
    real(8),dimension(Nspin,Nexc,Nexc)                 :: Mmat
    real(8),dimension(Nspin,Nexc)                      :: tmp_Lmat
    real(8),dimension(Nlat,Nspin,Norb,Nexc)            :: tmp_cQ,tmp_cdgQ
    do ispin=1,Nspin
       Mmat(ispin,:,:) = diag(Lmatrix(ispin,:))
    enddo

    do iexc=1,Nexc
       do jexc=1,Nexc
          !
          do ispin=1,Nspin
             do ilat=1,Nlat
                do jlat=1,Nlat
                   do iorb=1,Norb
                      do jorb=1,Norb
                         i = iorb + (ispin-1)*Norb + (ilat-1)*Nspin*Norb
                         j = jorb + (ispin-1)*Norb + (jlat-1)*Nspin*Norb
                         Mmat(ispin,iexc,jexc) = Mmat(ispin,iexc,jexc) + cdgQmatrix(ilat,ispin,iorb,iexc)*Vmat(i,j)*cQmatrix(jlat,ispin,jorb,iexc)
                      enddo
                   enddo
                   !
                enddo
             enddo
          enddo
          !
       enddo
    enddo

    do ispin=1,Nspin
       call eigh(Mmat(ispin,:,:),tmp_Lmat(ispin,:))
    enddo

    tmp_cQ=0d0
    tmp_cdgQ=0d0
    do ispin=1,Nspin
       do ilat=1,Nlat
          do iorb=1,Norb
             do iexc=1,Nexc
                do jexc=1,Nexc
                   tmp_cQ(ilat,ispin,iorb,iexc)   = tmp_cQ(ilat,ispin,iorb,iexc)   + cQmatrix(ilat,ispin,iorb,jexc)*Mmat(ispin,jexc,iexc)
                   tmp_cdgQ(ilat,ispin,iorb,iexc) = tmp_cdgQ(ilat,ispin,iorb,iexc) + Mmat(ispin,iexc,jexc)*cdgQmatrix(ilat,ispin,iorb,jexc)
                enddo
             enddo
          enddo
       enddo
    enddo

    cQmatrix   = tmp_cQ
    cdgQmatrix = tmp_cdgQ
    Lmatrix    = tmp_Lmat
  end subroutine vca_update_poles





  

end module VCA_MAIN
