! > SPARSE MAT-VEC DIRECT ON-THE-FLY PRODUCT
MODULE VCA_HAMILTONIAN_DIRECT_HxV
  USE VCA_HAMILTONIAN_COMMON
  implicit none
  private

  integer                              :: iiup,iidw
  integer                              :: iud,jj
  integer                              :: i,iup,idw
  integer                              :: j,jup,jdw
  integer                              :: m,mup,mdw
  integer                              :: ishift
  integer                              :: isector,jsector
  integer                              :: ms
  integer                              :: impi
  integer                              :: ilat,jlat,iorb,jorb,ispin,jspin,is,js
  integer                              :: kp,k1,k2,k3,k4
  integer                              :: ialfa,ibeta
  real(8)                              :: sg1,sg2,sg3,sg4
  real(8)                              :: htmp,htmpup,htmpdw
  logical                              :: Jcondition
  integer                              :: Nfoo
  !real(8),dimension(:,:,:),allocatable :: diag_hybr ![Nspin,Norb,Nbath]
  !real(8),dimension(:,:,:),allocatable :: bath_diag ![Nspin,Norb/1,Nbath]




  !>Sparse Mat-Vec direct on-the-fly product 
  public  :: directMatVec_main
  !public  :: directMatVec_orbs

contains


  subroutine directMatVec_main(Nloc,vin,Hv)
    integer                             :: Nloc
    real(8),dimension(Nloc)             :: vin
    real(8),dimension(Nloc)             :: Hv
    real(8),dimension(:),allocatable    :: vt,Hvt
    integer,dimension(Ns)               :: ibup,ibdw
    integer,dimension(2*Ns_Ud)          :: Indices,Jndices ![2-2*Norb] CONTROLLA
    integer,dimension(Ns_Ud,Ns_Orb)     :: Nups,Ndws       ![1,Ns]-[Norb,1+Nbath] CONTROLLA
    integer,dimension(Nlat,Norb)        :: Nup,Ndw !CONTROLLA
    !
    if(.not.Hstatus)stop "directMatVec_cc ERROR: Hsector NOT set"
    isector=Hsector
    !
    if(Nloc/=getdim(isector))stop "directMatVec_cc ERROR: Nloc != dim(isector)"
    !
    !Get diagonal hybridization, bath energy
    !include "ED_HAMILTONIAN/diag_hybr_bath.f90"
    !
    !
    Hv=0d0
    !
    !-----------------------------------------------!
    !LOCAL HAMILTONIAN PART: H_loc*vin = vout
    include "VCA_HAMILTONIAN/direct/HxV_local.f90" 
    !
    !UP HAMILTONIAN TERMS
    include "VCA_HAMILTONIAN/direct/HxV_up.f90"
    !    
    !DW HAMILTONIAN TERMS
    include "VCA_HAMILTONIAN/direct/HxV_dw.f90"
    !-----------------------------------------------!
    !
    !deallocate(diag_hybr,bath_diag)
    return
  end subroutine directMatVec_main



END MODULE VCA_HAMILTONIAN_DIRECT_HxV
