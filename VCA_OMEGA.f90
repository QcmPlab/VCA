!#TODO;
!-check get_Sft_potential
!-check Gmatrix type and functions
!-this has to give only omega (GLOBAL VARIABLE) no minimize loop inside main


MODULE VCA_OMEGA
  !
  USE VCA_GF_SHARED
  USE VCA_GF_NORMAL
  USE SCIFOR
  !
  implicit none
  private 

  public                                            :: sum_kmesh
  public                                            :: frequency_integration
  public                                            :: frequency_integration_sample
  complex(8),allocatable,dimension(:,:)             :: tmp_mat
  complex(8),allocatable,dimension(:,:,:,:,:,:)     :: gfprime ![Nlat][Nlat][Nspin][Nspin][Norb][Norb]



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
    if(allocated(gfprime))deallocate(gfprime)
    !
    allocate(tmp_mat(Nlat*Nspin*Norb,Nlat*Nspin*Norb))
    allocate(gfprime(Nlat,Nlat,Nspin,Nspin,Norb,Norb))
    !
    tmp_mat=zero
    gfprime=zero
    !    
    call vca_gf_cluster(xi*omega,gfprime)
    !
    do ii=1,size(impHk,7)
       tmp_mat=eye(Nlat*Nspin*Norb)+matmul(vca_nnn2lso_reshape(impHk(:,:,:,:,:,:,ii)-impHloc,Nlat,Nspin,Norb),vca_nnn2lso_reshape(gfprime,Nlat,Nspin,Norb))  !GMATS AT 0 TEMP, ELSE I NEED GREAL
       out_1=out_1+log(abs(REAL(det(tmp_mat))))
    enddo
    out_1=out_1*pi*pi/size(impHk,7) ! FIXME: careful, only works in 2d
    return
    !
  end function sum_kmesh

  !+------------------------------------------------------------------+
  !PURPOSE  : Do the frequency sum at zero T
  !+------------------------------------------------------------------+

  function frequency_integration() result(out_2)
    integer                         :: ii,ispin,iorb,inf
    real(8)                         :: out_2
    !
    out_2=0.d0
    call quad(sum_kmesh,a=0.d0,inf=1,verbose=.true.,result=out_2)
    out_2=out_2/pi
    return
  end function frequency_integration

  function frequency_integration_sample() result(out_2)
      integer                          :: N,i
      real(8),dimension(:),allocatable :: x,func
      real(8)                          :: out_2,a,b
      !
      a=0.d0
      b=1000.d0
      N=9999
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


end MODULE VCA_OMEGA













