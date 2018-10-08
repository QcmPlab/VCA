MODULE VCA_OMEGA
  !
  USE VCA_GF_SHARED
  USE VCA_GF_NORMAL
  USE SCIFOR, only: det
  USE DMFT_TOOLS, only:TB_build_kgrid
  !
  implicit none
  private 

  public                                            :: sum_kmesh
  public                                            :: generate_hcluster
  integer                                           :: lmats_old,lreal_old
  complex(8),allocatable,dimension(:,:,:,:,:,:)     :: t_k, t_prime
  complex(8),allocatable,dimension(:,:)             :: tmp_mat
  real(8)                                           :: mu,t



contains
    
  subroutine generate_hcluster()
    integer                         :: ii,ispin,iorb
    !
    allocate(t_prime(Nlat,Nlat,Nspin,Nspin,Norb,Norb))
    t_prime=zero
    !
    do ispin=1,Nspin !!FOR NOW NSPIN=1
      do iorb=1,Norb !!FOR NOW NORB=1
        do ii=1,Nlat
          t_prime(ii,ii,1,1,1,1)=1.d0
          if(ii+1<=Nlat .and. ii-1 > 0) then
            t_prime(ii,ii+1,1,1,1,1)= t_prime(ii,ii+1,1,1,1,1) + t
            t_prime(ii,ii-1,1,1,1,1)= t_prime(ii,ii+1,1,1,1,1) + t
          endif
          t_prime(1,Nlat,1,1,1,1)= t_prime(1,Nlat,1,1,1,1) + t
          t_prime(Nlat,1,1,1,1,1)= t_prime(Nlat,1,1,1,1,1) + t
        enddo
      enddo
    enddo
  end subroutine generate_hcluster

  !+------------------------------------------------------------------+
  !PURPOSE  : Function that does the k-sum in the RBZ of the VCA variational function(al)
  !+------------------------------------------------------------------+
  function sum_kmesh(omega,params) result(ASD)
    integer                                                  :: ii,jj,kk
    real(8),intent(in)                                       :: omega
    real(8)                                                  :: asd
    real(8),dimension(:),intent(in),optional                 :: params
    real(8),dimension(100,2)                                 :: kgrid
    !
    ASD=0.d0
    mu=0.5d0
    t=1.d0
    LMATS_OLD=LMATS
    LREAL_OLD=LREAL
    LMATS=1 !trick to use build_gf for a specific input value of frequency
    LREAL=1
    !
    !
    if(allocated(wm))deallocate(wm)
    if(allocated(vm))deallocate(vm)
    if(allocated(wr))deallocate(wr)
    if(allocated(tau))deallocate(tau)
    !
    if(.not.allocated(wm))allocate(wm(Lmats))
    if(.not.allocated(vm))allocate(vm(0:Lmats))          !bosonic frequencies
    if(.not.allocated(wr))allocate(wr(Lreal))
    if(.not.allocated(tau))allocate(tau(0:Ltau))
    !
    wm(1)=omega
    wr(1)=omega
    !
    if(allocated(impGmats))deallocate(impGmats)
    if(allocated(impGreal))deallocate(impGreal)
    !
    allocate(impGmats(Nlat,Nlat,Nspin,Nspin,Norb,Norb,Lmats))
    allocate(impGreal(Nlat,Nlat,Nspin,Nspin,Norb,Norb,Lmats))
    allocate(t_k(Nlat,Nlat,Nspin,Nspin,Norb,Norb))  
    allocate(tmp_mat(Nlat*Nspin*Norb,Nlat*Nspin*Norb))
    !
    impGmats=zero    
    impGreal=zero
    t_k=zero
    tmp_mat=zero
    !    
    call TB_build_kgrid([10,10],kgrid) !!!!!!DIVIDI PER NUMERO SITI, RBZ
    !
    kgrid=kgrid/2.d0
    write(LOGfile,"(A)")"THIS IS A TEST:"
    !
    call build_gf_normal()
    !
    call generate_hcluster()
    !
    do ii=1,100
       t_k=t_prime
       t_k(1,4,1,1,1,1)=t_k(1,4,1,1,1,1)+t*exp(xi*kgrid(ii,2))
       t_k(4,1,1,1,1,1)=t_k(4,1,1,1,1,1)+t*exp(-xi*kgrid(ii,2))
       t_k(1,2,1,1,1,1)=t_k(1,2,1,1,1,1)+t*exp(xi*kgrid(ii,1))
       t_k(2,1,1,1,1,1)=t_k(2,1,1,1,1,1)+t*exp(-xi*kgrid(ii,1))
       !
       t_k(2,3,1,1,1,1)=t_k(2,3,1,1,1,1)+t*exp(-xi*kgrid(ii,2))
       t_k(3,2,1,1,1,1)=t_k(3,2,1,1,1,1)+t*exp(xi*kgrid(ii,2))
       t_k(3,4,1,1,1,1)=t_k(3,4,1,1,1,1)+t*exp(xi*kgrid(ii,1))
       t_k(4,3,1,1,1,1)=t_k(4,3,1,1,1,1)+t*exp(-xi*kgrid(ii,1))        
       !
       tmp_mat=eye(Nlat*Nspin*Norb)+matmul(vca_nnn2lso_reshape(t_k-t_prime,Nlat,Nspin,Norb),vca_nnn2lso_reshape(impGmats(:,:,:,:,:,:,1),Nlat,Nspin,Norb))  !GMATS AT 0 TEMP, ELSE I NEED GREAL
       ASD=ASD+log(abs(dble(det(tmp_mat))))
    enddo
    !
    call deallocate_grids
    LMATS=LMATS_OLD
    LREAL=LREAL_OLD
    !
  end function sum_kmesh


end MODULE VCA_OMEGA













