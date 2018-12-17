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

  public                                            :: test_ksum
  public                                            :: reconstruct_g


  complex(8),allocatable,dimension(:,:,:,:,:,:)     :: gfprime ![Nlat][Nlat][Nspin][Nspin][Norb][Norb]
  
  complex(8),allocatable,dimension(:,:,:,:,:,:,:)   :: gftest ![Nlat][Nlat][Nspin][Nspin][Norb][Norb]



contains

 
  !+------------------------------------------------------------------+
  !PURPOSE  : Function that does the k-sum in the RBZ of the VCA variational function(al)
  !+------------------------------------------------------------------+
  function sum_kmesh(omega) result(out_1)
    integer                                                  :: ii,jj,kk
    real(8)                                                  :: omega
    real(8)                                                  :: out_1
    complex(8),allocatable,dimension(:,:)                    :: tmp_mat
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
    !    
    call vca_gf_cluster(xi*omega,gfprime)
    !
    do ii=1,size(impHk,7)
       tmp_mat=eye(Nlat*Nspin*Norb)+matmul(vca_nnn2lso_reshape(impHloc-impHk(:,:,:,:,:,:,ii),Nlat,Nspin,Norb),vca_nnn2lso_reshape(gfprime,Nlat,Nspin,Norb))
       out_1=out_1+log(abs(det(tmp_mat)))
    enddo
    out_1=out_1/size(impHk,7)
    !
    deallocate(tmp_mat)
    deallocate(gfprime)
    return
    !
  end function sum_kmesh

  function sum_kmesh_embedded(omega) result(out_1)
    integer                                                  :: ii,jj,kk
    real(8)                                                  :: omega
    real(8)                                                  :: out_1
    complex(8),allocatable,dimension(:,:)                    :: tmp_mat
    !
    out_1=0.d0
    !
    !
    if(allocated(tmp_mat))deallocate(tmp_mat)
    !
    allocate(tmp_mat(Nlat*Nspin*Norb*(Nbath+1),Nlat*Nspin*Norb*(Nbath+1)))
    !
    tmp_mat=zero
    !
    do ii=1,size(embeddedHk,3)
       tmp_mat=eye(Nlat*Nspin*Norb*(Nbath+1))+matmul(embeddedHloc-embeddedHk(:,:,ii),build_embedded_gf(omega))
       out_1=out_1+log(abs(det(tmp_mat)))
       !print*,out_1
    enddo
    out_1=out_1/size(embeddedHk,3)
    !
    deallocate(tmp_mat)
    return
    !
  end function sum_kmesh_embedded

  !+------------------------------------------------------------------+
  !PURPOSE  : Do the frequency sum at zero T
  !+------------------------------------------------------------------+

  function frequency_integration() result(out_2)
    integer                         :: inf
    real(8)                         :: out_2,spin_multiplicity
    !
    out_2=0.d0
    spin_multiplicity=3.d0-Nspin 
    !
    !
    if(Nbath > 0)then
      write(LOGfile,"(A)")"Calculating Omega with embedded matrices"
      call quad(sum_kmesh_embedded,a=0.0d0,inf=1,verbose=(verbose>=3),result=out_2)
    else
      write(LOGfile,"(A)")"Calculating Omega with original matrices"
      call quad(sum_kmesh,a=0.0d0,inf=1,verbose=(verbose>=3),result=out_2)
    endif
    out_2=spin_multiplicity*out_2/pi 
    return
  end function frequency_integration

  function frequency_integration_sample() result(out_2)
      integer                          :: N,i
      real(8),dimension(:),allocatable :: x,func
      real(8)                          :: out_2,a,b,spin_multiplicity
      !
      a=0.0001d0
      b=9999.d0
      N=1000
      out_2=0.d0
      spin_multiplicity=2.d0
      allocate(x(N),func(N))
      x = linspace(a,b,N)
      do i=1,N
         func(i) = sum_kmesh_embedded(x(i))
      enddo
      call quad(func,a,b,Ninterp=3,key=6,epsabs=0d0,epsrel=1d-4,verbose=.true.,result=out_2)
      out_2=spin_multiplicity*out_2/pi 
      deallocate(x,func)
      return
  end function frequency_integration_sample





 !+------------------------------------------------------------------+
  !PURPOSE  : DEBUG
  !+------------------------------------------------------------------+

  subroutine reconstruct_g
    complex(8),allocatable,dimension(:,:)                    :: tmp_mat
    character(len=64)                                        :: suffix
    integer                                                  :: ilat,jlat,iorb,ispin,ifreq
    !
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
    complex(8),allocatable,dimension(:,:)                    :: tmp_mat
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













