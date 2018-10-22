!#TODO;
!-debug!!


MODULE VCA_OMEGA
  !
  USE VCA_GF_SHARED
  USE VCA_GF_NORMAL
  USE SCIFOR
  !
  implicit none
  private 

  public                                            :: sum_kmesh
  public                                            :: test_ksum
  public                                            :: reconstruct_g
  public                                            :: frequency_integration
  public                                            :: frequency_integration_sample
  complex(8),allocatable,dimension(:,:)             :: tmp_mat,asd,lol
  complex(8),allocatable,dimension(:,:,:,:,:,:)     :: gfprime,wut ![Nlat][Nlat][Nspin][Nspin][Norb][Norb]
  complex(8),allocatable,dimension(:,:,:,:,:,:,:)   :: gftest ![Nlat][Nlat][Nspin][Nspin][Norb][Norb]



contains

 
  !+------------------------------------------------------------------+
  !PURPOSE  : Function that does the k-sum in the RBZ of the VCA variational function(al)
  !+------------------------------------------------------------------+
  function sum_kmesh(omega) result(out_1)
    integer                                                  :: ii,jj,kk
    real(8)                                                  :: omega
    real(8)                                                  :: out_1
    !
    out_1=0.d0
    !
    !
    if(allocated(tmp_mat))deallocate(tmp_mat)
    if(allocated(asd))deallocate(asd)
    if(allocated(lol))deallocate(lol)
    if(allocated(wut))deallocate(wut)
    if(allocated(gfprime))deallocate(gfprime)
    !
    allocate(tmp_mat(Nlat*Nspin*Norb,Nlat*Nspin*Norb))
    allocate(asd(Nlat*Nspin*Norb,Nlat*Nspin*Norb))
    allocate(lol(Nlat*Nspin*Norb,Nlat*Nspin*Norb))
    allocate(wut(Nlat,Nlat,Nspin,Nspin,Norb,Norb))
    allocate(gfprime(Nlat,Nlat,Nspin,Nspin,Norb,Norb))
    !
    tmp_mat=zero
    gfprime=zero
    !
    !    
    call vca_gf_cluster(xi*omega,gfprime)
    !
    do ii=1,size(impHk,7)
       asd=vca_nnn2lso_reshape(gfprime,Nlat,Nspin,Norb)
       wut=impHloc-impHk(:,:,:,:,:,:,ii)
       lol=vca_nnn2lso_reshape(wut,Nlat,Nspin,Norb) 
       tmp_mat=eye(Nlat*Nspin*Norb)+matmul(lol,asd) !!FIXME: OCIO!
       out_1=out_1+log(abs(det(tmp_mat)))
    enddo
    out_1=out_1/size(impHk,7) !*(pi**Ndim)
    !
    deallocate(tmp_mat)
    deallocate(gfprime)
    deallocate(asd)
    deallocate(lol)
    deallocate(wut)
    return
    !
  end function sum_kmesh

  !+------------------------------------------------------------------+
  !PURPOSE  : Do the frequency sum at zero T
  !+------------------------------------------------------------------+

  function frequency_integration() result(out_2)
    integer                         :: inf
    real(8)                         :: out_2
    !
    out_2=0.d0
    call quad(sum_kmesh,a=0.0d0,inf=1,verbose=.true.,result=out_2)!!FIXME
    out_2=2*out_2/pi  ! FIXME: PERCHÈ DIAVOLO?
    return
  end function frequency_integration

  function frequency_integration_sample() result(out_2)
      integer                          :: N,i
      real(8),dimension(:),allocatable :: x,func
      real(8)                          :: out_2,a,b
      !
      a=0.0001d0  !!FIXME
      b=1000.d0
      N=1000
      out_2=0.d0
      allocate(x(N),func(N))
      x = linspace(a,b,N)
      do i=1,N
         func(i) = sum_kmesh(x(i))
      enddo
      call quad(func,a,b,Ninterp=3,key=6,epsabs=0d0,epsrel=1d-4,verbose=.true.,result=out_2)
      out_2=out_2/pi
      deallocate(x,func)
      return
  end function frequency_integration_sample








 !+------------------------------------------------------------------+
  !PURPOSE  : DEBUG
  !+------------------------------------------------------------------+

  subroutine reconstruct_g
    character(len=64) :: suffix
    integer           :: ilat,jlat,iorb,ispin,ifreq
    allocate(gftest(Nlat,Nlat,Nspin,Nspin,Norb,Norb,Lmats))
    gftest=0
    !
    !
    if(.not.allocated(wm))allocate(wm(Lmats))
    wm     = pi/beta*real(2*arange(1,Lmats)-1,8)
    !
    do ifreq=1,size(wm)
        call vca_gf_cluster(xi*wm(ifreq),gftest(:,:,:,:,:,:,ifreq))
    enddo
    !
    do ilat=1,Nlat
       do jlat=1,Nlat
          do iorb=1,Norb
             do ispin=1,Nspin
                suffix="_Isite"//str(ilat,4)//"_Jsite"//str(jlat,4)//"_l"//str(iorb)//str(iorb)//"_s"//str(ispin)
                call splot("Gtest"//reg(suffix)//"_iw"//reg(file_suffix)//".vca"   ,wm,gftest(ilat,jlat,ispin,ispin,iorb,iorb,:))
             enddo
          enddo
       enddo
    enddo
    !
    if(allocated(wm))deallocate(wm)
    if(allocated(gftest))deallocate(gftest)
    !
  end subroutine reconstruct_g

  subroutine test_ksum
    integer                                                  :: ii,jj,kk
    real(8)                                                  :: omega
    real(8)                                                  :: out_1
    !
    out_1=0.d0
    omega=1.0d0
    !
    !
    if(allocated(tmp_mat))deallocate(tmp_mat)
    if(allocated(gfprime))deallocate(gfprime)
    !
    allocate(tmp_mat(Nlat*Nspin*Norb,Nlat*Nspin*Norb))
    allocate(gfprime(Nlat,Nlat,Nspin,Nspin,Norb,Norb))
    !
    tmp_mat=zero
    gfprime=zero
    !    
    call vca_gf_cluster(xi*omega,gfprime)
    do ii=1,size(impHk,7)
       print*,"ITERATION",ii
       !print*,"OMEGA",xi*omega
       !print*,"AAAAAAAAA",matmul(vca_nnn2lso_reshape(impHk(:,:,:,:,:,:,ii)-impHloc,Nlat,Nspin,Norb),vca_nnn2lso_reshape(gfprime,Nlat,Nspin,Norb))
       !print*,"BBBBBBBBB",vca_nnn2lso_reshape(impHk(:,:,:,:,:,:,ii)-impHloc,Nlat,Nspin,Norb)
       tmp_mat=eye(Nlat*Nspin*Norb)+matmul(vca_nnn2lso_reshape(impHk(:,:,:,:,:,:,ii)-impHloc,Nlat,Nspin,Norb),vca_nnn2lso_reshape(gfprime,Nlat,Nspin,Norb))
       out_1=out_1+log(abs(REAL(det(tmp_mat))))
       print*,"CCCCCCCCC",(vca_nnn2lso_reshape(impHk(:,:,:,:,:,:,ii)-impHloc,Nlat,Nspin,Norb))
    enddo
    out_1=out_1/size(impHk,7) !*(pi**Ndim)
    !
    deallocate(tmp_mat)
    deallocate(gfprime)
    return
    !
  end subroutine test_ksum



end MODULE VCA_OMEGA













