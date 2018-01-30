MODULE VCA_GREENS_FUNCTIONS
  !
  USE SF_CONSTANTS, only:one,xi,zero,pi
  USE SF_TIMER  
  USE SF_IOTOOLS, only: str,free_unit,reg,free_units,txtfy,splot
  USE SF_ARRAYS,  only: arange,linspace
  USE SF_LINALG,  only: inv,inv_sym,inv_her,eye
  USE SF_SP_LINALG, only: sp_lanc_tridiag
  !
  USE VCA_VARS_GLOBAL
  USE VCA_DIAG
  USE VCA_SETUP
  USE VCA_EIGENSPACE
  USE VCA_AUX_FUNX
  USE VCA_IO
  USE VCA_BATH_FUNCTIONS
  USE VCA_HAMILTONIAN_MATVEC
  !
  implicit none
  private 



  !Lanczos shared variables
  !=========================================================
  real(8),dimension(:),pointer                :: state_vec
  complex(8),dimension(:),pointer             :: state_cvec
  real(8)                                     :: state_e



  public :: build_gf_cluster



contains


  !+------------------------------------------------------------------+
  !PURPOSE  : Build the GF of the cluster 
  !+------------------------------------------------------------------+
  subroutine build_gf_cluster()
    impGmats=zero
    impGreal=zero
    !
    impSmats = zero
    impSreal = zero
    !
    impG0mats=zero
    impG0real=zero
    !
    !
    write(LOGfile,"(A)")"Get impurity Greens functions:"
    call build_gf_normal()
    call get_sigma_normal()
    !
    if(print_Sigma)call vca_print_impSigma()
    if(print_impG)call vca_print_impG()
    if(print_impG0)call vca_print_impG0()
    !
  end subroutine build_gf_cluster





  !+------------------------------------------------------------------+
  !PURPOSE  : Evaluate Green's functions
  !+------------------------------------------------------------------+
  subroutine build_gf_normal()
    integer :: iorb,ispin
    !
    if(allocated(impGmatrix))deallocate(impGmatrix)
    allocate(impGmatrix(Nlat,Nlat,Nspin,Nspin,Norb,Norb))
    !
    select case(vca_method)
    case ("full")
       !NORMAL: (default)
       call full_build_gf_normal()
    case ("lanc")
       !NORMAL: (default)    
       call lanc_build_gf_normal()
    case default
       stop "build_gf_normal error: vca_method != ['full','lanc']"
    end select
  end subroutine build_gf_normal


  !=========================================================
  !=========================================================
  include 'vca_full_gf_normal.f90'
  !
  include 'vca_lanc_gf_normal.f90'
  !=========================================================
  !=========================================================









  !+------------------------------------------------------------------+
  !                    SELF-ENERGY FUNCTIONS 
  !+------------------------------------------------------------------+
  subroutine get_sigma_normal
    integer                                                     :: i,ispin,iorb
    complex(8),dimension(Nlat,Nlat,Nspin,Nspin,Norb,Norb,Lmats) :: invG0mats,invGmats
    complex(8),dimension(Nlat,Nlat,Nspin,Nspin,Norb,Norb,Lreal) :: invG0real,invGreal
    complex(8),dimension(Nlat,Nlat)                             :: invGimp
    !
    call vca_allocate_time_freq_arrays()
    !
    invG0mats = zero
    invGmats  = zero
    invG0real = zero
    invGreal  = zero
    !
    !Get G0^-1
    invG0mats = invg0_bath_mats(dcmplx(0d0,wm(:)),vca_bath)
    invG0real = invg0_bath_real(dcmplx(wr(:),eps),vca_bath)
    !
    !
    !Get Gimp^-1
    invGmats=zero
    invGreal=zero
    do ispin=1,Nspin
       do iorb=1,Norb
          do i=1,Lmats
             invGimp = impGmats(:,:,ispin,ispin,iorb,iorb,i)
             call inv(invGimp)
             invGmats(:,:,ispin,ispin,iorb,iorb,i)=invGimp
          enddo
          !
          do i=1,Lreal
             invGimp = impGreal(:,:,ispin,ispin,iorb,iorb,i)
             call inv(invGimp)
             invGreal(:,:,ispin,ispin,iorb,iorb,i)=invGimp
          enddo
       enddo
    enddo
    !
    !Get Sigma functions: Sigma= G0^-1 - G^-1
    impSmats=zero
    impSreal=zero
    do ispin=1,Nspin
       do iorb=1,Norb
          impSmats(:,:,ispin,ispin,iorb,iorb,:) = invG0mats(:,:,ispin,ispin,iorb,iorb,:) - invGmats(:,:,ispin,ispin,iorb,iorb,:)
          impSreal(:,:,ispin,ispin,iorb,iorb,:) = invG0real(:,:,ispin,ispin,iorb,iorb,:) - invGreal(:,:,ispin,ispin,iorb,iorb,:)
       enddo
    enddo
    !
    !
    !Get G0and:
    impG0mats = g0and_bath_mats(dcmplx(0d0,wm(:)),vca_bath)
    impG0real = g0and_bath_real(dcmplx(wr(:),eps),vca_bath)
    !
    call vca_deallocate_time_freq_arrays()
    !
  end subroutine get_sigma_normal

















  !+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  !++++++++++++++++++COMPUTATIONAL ROUTINE: TQL2++++++++++++++++++++++++ 
  !+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  !---------------------------------------------------------------------
  ! PURPOSE computes all eigenvalues/vectors, real symmetric tridiagonal matrix.
  !    This subroutine finds the eigenvalues and eigenvectors of a symmetric
  !    tridiagonal matrix by the QL method.  The eigenvectors of a full
  !    symmetric matrix can also be found if TRED2 has been used to reduce this
  !    full matrix to tridiagonal form.
  !  Parameters:
  !    Input, integer ( kind = 4 ) N, the order of the matrix.
  !
  !    Input/output, real ( kind = 8 ) D(N).  On input, the diagonal elements of
  !    the matrix.  On output, the eigenvalues in ascending order.  If an error
  !    exit is made, the eigenvalues are correct but unordered for indices
  !    1,2,...,IERR-1.
  !
  !    Input/output, real ( kind = 8 ) E(N).  On input, E(2:N) contains the
  !    subdiagonal elements of the input matrix, and E(1) is arbitrary.
  !    On output, E has been destroyed.
  !
  !    Input, real ( kind = 8 ) Z(N,N).  On input, the transformation matrix
  !    produced in the reduction by TRED2, if performed.  If the eigenvectors of
  !    the tridiagonal matrix are desired, Z must contain the identity matrix.
  !    On output, Z contains the orthonormal eigenvectors of the symmetric
  !    tridiagonal (or full) matrix.  If an error exit is made, Z contains
  !    the eigenvectors associated with the stored eigenvalues.
  !
  !    Output, integer ( kind = 4 ) IERR, error flag.
  !    0, normal return,
  !    J, if the J-th eigenvalue has not been determined after
  !    30 iterations.
  !
  !---------------------------------------------------------------------
  subroutine tql2 ( n, d, e, z, ierr )
    integer :: n
    real(8) :: c
    real(8) :: c2
    real(8) :: c3
    real(8) :: d(n)
    real(8) :: dl1
    real(8) :: e(n)
    real(8) :: el1
    real(8) :: f
    real(8) :: g
    real(8) :: h
    integer ( kind = 4 ) i
    integer ( kind = 4 ) ierr
    integer ( kind = 4 ) ii
    integer ( kind = 4 ) j
    integer ( kind = 4 ) k
    integer ( kind = 4 ) l
    integer ( kind = 4 ) l1
    integer ( kind = 4 ) l2
    integer ( kind = 4 ) m
    integer ( kind = 4 ) mml
    real(8) :: p
    real(8) :: r
    real(8) :: s
    real(8) :: s2
    real(8) :: tst1
    real(8) :: tst2
    real(8) :: z(n,n)
    ierr = 0
    if ( n == 1 ) then
       return
    end if
    do i = 2, n
       e(i-1) = e(i)
    end do
    f = 0.0D+00
    tst1 = 0.0D+00
    e(n) = 0.0D+00
    do l = 1, n
       j = 0
       h = abs ( d(l) ) + abs ( e(l) )
       tst1 = max ( tst1, h )
       !
       !  Look for a small sub-diagonal element.
       !
       do m = l, n
          tst2 = tst1 + abs ( e(m) )
          if ( tst2 == tst1 ) then
             exit
          end if
       end do
       if ( m == l ) then
          go to 220
       end if
130    continue
       if ( 30 <= j ) then
          ierr = l
          return
       end if
       j = j + 1
       !
       !  Form shift.
       !
       l1 = l + 1
       l2 = l1 + 1
       g = d(l)
       p = ( d(l1) - g ) / ( 2.0D+00 * e(l) )
       r = pythag ( p, 1.0D+00 )
       d(l) = e(l) / ( p + sign ( r, p ) )
       d(l1) = e(l) * ( p + sign ( r, p ) )
       dl1 = d(l1)
       h = g - d(l)
       d(l2:n) = d(l2:n) - h
       f = f + h
       !
       !  QL transformation.
       !
       p = d(m)
       c = 1.0D+00
       c2 = c
       el1 = e(l1)
       s = 0.0D+00
       mml = m - l
       do ii = 1, mml
          c3 = c2
          c2 = c
          s2 = s
          i = m - ii
          g = c * e(i)
          h = c * p
          r = pythag ( p, e(i) )
          e(i+1) = s * r
          s = e(i) / r
          c = p / r
          p = c * d(i) - s * g
          d(i+1) = h + s * ( c * g + s * d(i) )
          !
          !  Form vector.
          !
          do k = 1, n
             h = z(k,i+1)
             z(k,i+1) = s * z(k,i) + c * h
             z(k,i) = c * z(k,i) - s * h
          end do
       end do
       p = - s * s2 * c3 * el1 * e(l) / dl1
       e(l) = s * p
       d(l) = c * p
       tst2 = tst1 + abs ( e(l) )
       if ( tst2 > tst1 ) then
          go to 130
       end if
220    continue
       d(l) = d(l) + f
    end do
    !
    !  Order eigenvalues and eigenvectors.
    !
    do ii = 2, n
       i = ii - 1
       k = i
       p = d(i)
       do j = ii, n
          if ( d(j) < p ) then
             k = j
             p = d(j)
          end if
       end do
       if ( k /= i ) then
          d(k) = d(i)
          d(i) = p
          do j = 1, n
             call r8_swap ( z(j,i), z(j,k) )
          end do
       end if
    end do
    return
  end subroutine tql2


  !---------------------------------------------------------------------
  ! PURPOSE: computes SQRT ( A * A + B * B ) carefully.
  !    The formula
  !    PYTHAG = sqrt ( A * A + B * B )
  !    is reasonably accurate, but can fail if, for example, A**2 is larger
  !    than the machine overflow.  The formula can lose most of its accuracy
  !    if the sum of the squares is very large or very small.
  !  Parameters:
  !    Input, real(8) :: A, B, the two legs of a right triangle.
  !    Output, real(8) :: PYTHAG, the length of the hypotenuse.
  !---------------------------------------------------------------------
  function pythag ( a, b )
    implicit none
    real(8) :: a
    real(8) :: b
    real(8) :: p
    real(8) :: pythag
    real(8) :: r
    real(8) :: s
    real(8) :: t
    real(8) :: u
    p = max ( abs ( a ), abs ( b ) )
    if ( p /= 0.0D+00 ) then
       r = ( min ( abs ( a ), abs ( b ) ) / p )**2
       do
          t = 4.0D+00 + r
          if ( t == 4.0D+00 ) then
             exit
          end if
          s = r / t
          u = 1.0D+00 + 2.0D+00 * s
          p = u * p
          r = ( s / u )**2 * r
       end do
    end if
    pythag = p
    return
  end function pythag

  !---------------------------------------------------------------------
  ! PURPOSE: swaps two R8's.
  !  Parameters:
  !    Input/output, real(8) :: X, Y.  On output, the values of X and
  !    Y have been interchanged.
  !---------------------------------------------------------------------
  subroutine r8_swap ( x, y )
    real(8) :: x
    real(8) :: y
    real(8) :: z
    z = x
    x = y
    y = z
    return
  end subroutine r8_swap





end MODULE VCA_GREENS_FUNCTIONS













