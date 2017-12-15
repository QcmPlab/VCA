MODULE VCA_SPARSE_MATRIX  !THIS VERSION CONTAINS ONLY COMPLX ELEMENT: (HERMITIAN MATRIX) 

  implicit none
  private

  complex(8),parameter :: zero=dcmplx(0d0,0d0)
  
  type sparse_element
     complex(8)                            :: cval !value of the entry: double complex
     integer                               :: col  !col connected to this compress value
     type(sparse_element),pointer          :: next !link to next entry in the row
  end type sparse_element
  public :: sparse_element

  type sparse_row
     integer                               :: size    !size of the list
     type(sparse_element),pointer          :: root    !head/root of the list\== list itself
  end type sparse_row
  public :: sparse_row

  type sparse_matrix
     integer                               :: Nrow
     logical                               :: status=.false.
     type(sparse_row),dimension(:),pointer :: row
  end type sparse_matrix
  public :: sparse_matrix


  !INIT SPARSE MATRICES (LL)
  interface sp_init_matrix
     module procedure sp_init_matrix_ll
  end interface sp_init_matrix
  public :: sp_init_matrix      !init the sparse matrix   !checked


  !DELETE SPARSE MATRIX (LL) OR ONE OF ITS ELEMENTS (LL)
  interface sp_delete_matrix
     module procedure sp_delete_matrix_ll
  end interface sp_delete_matrix
  public :: sp_delete_matrix    !delete the sparse matrix !checked
  public :: sp_delete_element   !delete n-th/last element !checked



  !GET NUMBER OF NON-ZERO ELEMENTS
  interface sp_get_nnz
     module procedure sp_get_nnz_ll
  end interface sp_get_nnz
  public :: sp_get_nnz


  !INSERT ELEMENTS (D,C) IN LL-SPARSE MATRIX
  interface sp_insert_element
     module procedure sp_insert_element_c
  end interface sp_insert_element
  public :: sp_insert_element   !insert an element        !checked



  !INSERT DIAGONAL ENTRY IN LL-SPARSE MATRIX
  interface sp_insert_diag
     module procedure sp_insert_diag_c
  end interface sp_insert_diag
  public :: sp_insert_diag      !insert a vector at diag  !checked



  !GET ELEMENTS ALONG THE DIAGONAL
  interface sp_get_diagonal
     module procedure sp_get_diagonal_c
  end interface sp_get_diagonal
  public :: sp_get_diagonal     !get diagonal elements    !checked



  !LOAD STANDARD MATRIX INTO SPARSE MATRICES
  interface sp_load_matrix
     module procedure sp_load_matrix_c
  end interface sp_load_matrix
  public :: sp_load_matrix      !create sparse from array !checked



  !DUMP SPARSE MATRIX INTO STANDARD MATRIX
  interface sp_dump_matrix
     module procedure sp_dump_matrix_c
  end interface sp_dump_matrix
  public :: sp_dump_matrix      !dump sparse into array   !checked



  !PRETTY PRINTING
  interface sp_print_matrix
     module procedure sp_print_matrix_ll
  end interface sp_print_matrix
  public :: sp_print_matrix     !print sparse             !checked


  !TEST
  public :: sp_test_symmetric



  !GET ELEMENT FROM SPARSE MATRIX
  public :: sp_get_element_c    !""                       !checked


  !INQUIRE IF ELEMENT EXISTS
  public :: sp_inquire_element  !inquire an element       !checked





contains       




  !+------------------------------------------------------------------+
  !PURPOSE:  initialize the sparse matrix list
  !+------------------------------------------------------------------+
  subroutine sp_init_matrix_ll(sparse,N)
    type(sparse_matrix),intent(inout) :: sparse
    integer                           :: i,N
    !put here a delete statement to avoid problems
    if(sparse%status)stop "sp_init_matrix: alreay allocate can not init"
    sparse%Nrow=N
    sparse%status=.true.
    allocate(sparse%row(N))
    do i=1,N
       allocate(sparse%row(i)%root)
       sparse%row(i)%root%next => null()
       sparse%row(i)%size=0
    end do
  end subroutine sp_init_matrix_ll








  !+------------------------------------------------------------------+
  !PURPOSE: delete an entire sparse matrix
  !+------------------------------------------------------------------+
  subroutine sp_delete_matrix_ll(sparse)    
    type(sparse_matrix),intent(inout) :: sparse
    integer                           :: i
    if(.not.sparse%status)stop "Warning SPARSE/sp_delete_matrix: sparse not allocated already."
    do i=1,sparse%Nrow
       call delete_row(sparse%row(i))
       deallocate(sparse%row(i)%root)
    end do
    deallocate(sparse%row)
    sparse%Nrow=0
    sparse%status=.false.
  end subroutine sp_delete_matrix_ll



  !+------------------------------------------------------------------+
  !PURPOSE: delete a single element at (i,j) from the sparse matrix
  !+------------------------------------------------------------------+
  subroutine sp_delete_element(matrix,i,j)
    type(sparse_matrix),intent(inout) :: matrix
    integer,intent(in)                :: i,j
    logical :: delete
    delete = delete_element_from_row(matrix%row(i),col=j)
    if(.not.delete)write(*,"(A,I3,I3)")"sp_delete_element: can not delete element in",i,j
  end subroutine sp_delete_element

  !+------------------------------------------------------------------+
  !PURPOSE: delete an entire row from the sparse matrix (private)
  !+------------------------------------------------------------------+
  subroutine delete_row(row)
    type(sparse_row),intent(inout) :: row
    type(sparse_element),pointer   :: p,c
    do
       p => row%root
       c => p%next
       if(.not.associated(c))exit  !empty list
       p%next => c%next !
       c%next=>null()
       deallocate(c)
    end do
  end subroutine delete_row




  !This shoud be better tested!
  !+------------------------------------------------------------------+
  !PURPOSE: delete a given element from a row of the sparse matrix (private)
  !+------------------------------------------------------------------+
  function delete_element_from_row(row,n,col) result(found)
    type(sparse_row),intent(inout)    :: row
    integer,optional                  :: n
    integer,optional                  :: col
    integer                           :: i,pos
    type(sparse_element),pointer      :: p,c
    logical                           :: found
    pos= row%size ; if(present(n))pos=n
    p => row%root
    c => p%next
    found = .false.
    if(present(col))then
       do 
          if(found .OR. .not.associated(c))return
          if(col == c%col)then
             found=.true.
             exit
          else
             p => c
             c => c%next
          endif
       end do
       if(found)then
          p%next => c%next !reallocate skipping the deleted link
          deallocate(c)           !free link
          row%size=row%size-1
       endif
    else
       do i=1,pos 
          if(.not.associated(c))return !empty list
          p => c
          c => c%next
       end do
       found=.true.
       p%next => c%next !reallocate skipping the deleted link
       deallocate(c)           !free link
       row%size=row%size-1
    endif
  end function delete_element_from_row










  !+------------------------------------------------------------------+
  !PURPOSE:  return total number of non-zero elements stored in sparse
  !+------------------------------------------------------------------+
  function sp_get_nnz_ll(sparse) result(Nnz)
    type(sparse_matrix) :: sparse
    integer             :: i
    integer             :: Nnz
    Nnz=0
    do i=1,sparse%Nrow
       Nnz=Nnz+sparse%row(i)%size
    enddo
  end function sp_get_nnz_ll











  !+------------------------------------------------------------------+
  !PURPOSE: insert an element value at position (i,j) in the sparse matrix
  !+------------------------------------------------------------------+
  subroutine sp_insert_element_c(sparse,value,i,j)
    type(sparse_matrix),intent(inout) :: sparse
    complex(8),intent(in)             :: value
    integer,intent(in)                :: i,j
    call insert_element_in_row_c(sparse%row(i),value,j)
  end subroutine sp_insert_element_c








  !+------------------------------------------------------------------+
  !PURPOSE: insert a vector of elements at the diagonal of the sparse matrix
  !+------------------------------------------------------------------+
  subroutine sp_insert_diag_c(sparse,diag)
    type(sparse_matrix),intent(inout)  :: sparse
    complex(8),intent(in),dimension(:) :: diag
    integer                            :: i
    if(size(diag)/=sparse%Nrow)stop "sp_insert_diag: error in dimensions"
    do i=1,size(diag)
       call insert_element_in_row_c(sparse%row(i),diag(i),i)
    enddo
  end subroutine sp_insert_diag_c








  !+------------------------------------------------------------------+
  !PURPOSE: insert an element in a given row (private) 
  !+------------------------------------------------------------------+
  subroutine insert_element_in_row_c(row,value,column)
    type(sparse_row),intent(inout)    :: row
    complex(8) ,intent(in)            :: value
    integer, intent(in)               :: column
    type(sparse_element),pointer      :: p,c
    logical :: iadd
    p => row%root
    c => p%next
    iadd = .false.                !check if column already exist
    do                            !traverse the list
       if(.not.associated(c))exit !empty list or end of the list
       if(c%col == column)then
          iadd=.true.
          exit
       endif
       !if(c%col > column)exit
       if(column <= c%col)exit
       p => c
       c => c%next
    end do
    if(iadd)then
       c%cval=c%cval + value
    else
       allocate(p%next)                !Create a new element in the list
       p%next%cval= value
       p%next%col = column
       row%size   = row%size+1
       if(.not.associated(c))then !end of the list special case (current=>current%next)
          p%next%next  => null()
       else
          p%next%next  => c      !the %next of the new node come to current
       end if
    endif
  end subroutine insert_element_in_row_c













  !+------------------------------------------------------------------+
  !PURPOSE: get the diagonal elements of the sparse matrix
  !+------------------------------------------------------------------+
  subroutine sp_get_diagonal_c(sparse,diag)
    type(sparse_matrix),intent(inout) :: sparse
    complex(8),dimension(:)           :: diag
    integer                           :: Ndim,i
    Ndim=size(diag);if(Ndim/=sparse%Nrow)stop "sp_get_diagonal: error in diag dimension." 
    do i=1,Ndim
       call get_element_from_row_c(sparse%row(i),diag(i),i)
    enddo
  end subroutine sp_get_diagonal_c











  !+------------------------------------------------------------------+
  !PURPOSE: get an element from position (i,j) of the sparse matrix
  !+------------------------------------------------------------------+
  function sp_get_element_c(sparse,i,j) result(value)
    type(sparse_matrix),intent(inout) :: sparse    
    integer,intent(in)                :: i,j
    complex(8)                        :: value
    call get_element_from_row_c(sparse%row(i),value,j)
  end function sp_get_element_c













  !+------------------------------------------------------------------+
  !PURPOSE: get an element from a given row of the matrix (private)
  !+------------------------------------------------------------------+
  subroutine get_element_from_row_c(row,value,column)
    type(sparse_row),intent(inout)    :: row
    complex(8)                        :: value
    integer, intent(in)               :: column
    type(sparse_element),pointer      :: c
    c => row%root%next
    value=cmplx(0.d0,0.d0,8)
    do                            !traverse the list
       if(.not.associated(c))return !empty list or end of the list
       if(c%col == column)exit
       c => c%next
    end do
    !
    value = c%cval
  end subroutine get_element_from_row_c












  !+------------------------------------------------------------------+
  !PURPOSE: check if a given element exists
  !+------------------------------------------------------------------+
  function sp_inquire_element(sparse,i,j) result(exist)
    type(sparse_matrix),intent(inout) :: sparse    
    integer,intent(in)                :: i,j
    logical                           :: exist
    exist = inquire_element_from_row(sparse%row(i),j)
  end function sp_inquire_element

  !+------------------------------------------------------------------+
  !PURPOSE: check if an element in a given row of the matrix exist (private)
  !+------------------------------------------------------------------+
  function inquire_element_from_row(row,column) result(exist)
    type(sparse_row),intent(inout)    :: row
    logical                           :: exist
    integer, intent(in)               :: column
    type(sparse_element),pointer      :: c
    c => row%root%next
    exist=.false.
    do                            !traverse the list
       if(.not.associated(c))return !empty list or end of the list
       if(c%col == column)exit
       c => c%next
    end do
    exist=.true.
  end function inquire_element_from_row




















  !+------------------------------------------------------------------+
  !PURPOSE: load a regular matrix (2dim array) into a sparse matrix
  !+------------------------------------------------------------------+
  subroutine sp_load_matrix_c(matrix,sparse)
    complex(8),dimension(:,:),intent(in)  :: matrix
    type(sparse_matrix),intent(inout)     :: sparse    
    integer                               :: i,j,Ndim1,Ndim2
    Ndim1=size(matrix,1)
    Ndim2=size(matrix,2)
    if(Ndim1/=Ndim2)print*,"Warning: SPARSE/load_matrix Ndim1.ne.Ndim2"
    if(sparse%Nrow /= Ndim1)stop "Warning SPARSE/load_matrix: dimensions error"
    do i=1,Ndim1
       do j=1,Ndim2
          if(matrix(i,j)/=cmplx(0.d0,0.d0,8))call sp_insert_element_c(sparse,matrix(i,j),i,j)
       enddo
    enddo
  end subroutine sp_load_matrix_c







  !+-----------------------------------------------------------------------------+!
  !PURPOSE:  test if a sparse matrix is symmetric 
  !+-----------------------------------------------------------------------------+
  subroutine sp_test_symmetric(sparse)
    type(sparse_matrix)                   :: sparse
    logical                               :: is_symmetric
    complex(8),dimension(:,:),allocatable :: cM
    integer                               :: Nrow,Ncol
    Nrow=sparse%Nrow
    Ncol=Nrow
    is_symmetric=.false.
    allocate(cM(Nrow,Ncol))
    call sp_dump_matrix(sparse,cM)
    if( maxval(abs(cM-conjg(transpose(cM))) ) < 1.d-12)is_symmetric=.true.
    if(is_symmetric)then
       write(*,"(A)")"Matrix IS Hermitian"
    else
       write(*,"(A)")"Matrix IS NOT Hermitian"
    endif
  end subroutine sp_test_symmetric






  !+------------------------------------------------------------------+
  !PURPOSE: dump a sparse matrix into a regular 2dim array
  !+------------------------------------------------------------------+
  subroutine sp_dump_matrix_c(sparse,matrix)
    type(sparse_matrix),intent(in)          :: sparse
    complex(8),dimension(:,:),intent(inout) :: matrix
    type(sparse_element),pointer            :: c
    integer                                 :: i,Ndim1,Ndim2
    Ndim1=size(matrix,1)
    Ndim2=size(matrix,2)
    !if(Ndim1/=Ndim2)print*,"Warning: SPARSE_MATRIX/sp_dump_matrix_d: Ndim1/=Ndim2"
    if(sparse%Nrow /= Ndim1)stop "Warning SPARSE/load_matrix: dimensions error"
    matrix=0.d0
    do i=1,Ndim1
       c => sparse%row(i)%root%next
       do while(associated(c))
          matrix(i,c%col) = c%cval
          c => c%next  !traverse list
       enddo
    enddo
  end subroutine sp_dump_matrix_c






















  !+------------------------------------------------------------------+
  !PURPOSE: pretty print a sparse matrix on a given unit using format fmt
  !+------------------------------------------------------------------+
  subroutine sp_print_matrix_ll(sparse,unit,fmt,full)
    type(sparse_matrix)            :: sparse
    integer,optional               :: unit
    integer                        :: i,j,unit_,Ns
    character(len=*),optional      :: fmt
    character(len=64)              :: fmt_
    logical,optional               :: full
    logical                        :: full_
    unit_=6;if(present(unit))unit_=unit
    fmt_='F6.2';if(present(fmt))fmt_=fmt
    full_=.false.;if(present(full))full_=full
    Ns=sparse%Nrow
    if(full_)then
       write(*,*)"Print sparse matrix (full mode < 100) ->",unit_
       do i=1,Ns
          write(unit_,"(100("//trim(fmt_)//",A1,"//trim(fmt_)//",2X))")(&
               real(sp_get_element_c(sparse,i,j)),",",imag(sp_get_element_c(sparse,i,j)),j=1,Ns)
       enddo
    else
       write(*,*)"Print sparse matrix (compact mode) ->",unit_
       do i=1,Ns
          call print_row_c(sparse%row(i),unit_,fmt_)
       enddo
    endif
    write(unit_,*)
  end subroutine sp_print_matrix_ll









  !+------------------------------------------------------------------+
  !PURPOSE: print an entire row of the sparse matrix (private)
  !+------------------------------------------------------------------+
  subroutine print_row_c(row,unit,fmt)
    type(sparse_row),intent(in)   :: row
    type(sparse_element),pointer  :: c
    integer                       :: count=0
    integer,optional :: unit
    integer          :: unit_
    character(len=*),optional :: fmt
    character(len=64)         :: fmt_
    unit_=6;if(present(unit))unit_=unit
    fmt_='F15.9';if(present(fmt))fmt_=fmt
    c => row%root%next   !assume is associated,ie list exists
    do
       if(.not.associated(c))exit
       count=count+1
       write(unit_,"(2"//trim(fmt_)//",A1,I3,3X)",advance='no')c%cval,',',c%col
       c => c%next  !traverse list
    end do
    write(unit_,*)
  end subroutine print_row_c













end module VCA_SPARSE_MATRIX
