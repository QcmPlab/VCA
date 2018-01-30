module VCA_DIAG
  USE VCA_INPUT_VARS
  USE VCA_VARS_GLOBAL
  USE VCA_EIGENSPACE
  USE VCA_SETUP
  USE VCA_HAMILTONIAN_MATVEC
  !
  USE SF_CONSTANTS
  USE SF_STAT
  USE SF_SP_LINALG
  USE SF_LINALG,  only: eigh
  USE SF_TIMER,   only: start_timer,stop_timer,eta
  USE SF_IOTOOLS, only: reg,free_unit
  USE SF_MISC,    only: assert_shape
  !
  implicit none
  private



  !>diagonalize cluster Hamiltonian
  public  :: diagonalize_cluster




contains

  !+-------------------------------------------------------------------+
  !>Setup the Hilbert space, create the Hamiltonian, get the
  ! GS, build the Green's functions calling all the necessary routines
  !+-------------------------------------------------------------------+
  subroutine diagonalize_cluster
    select case(vca_method)
    case ("full")
       call vca_diag_h
    case ("lanc")
       call vca_lanc_h
    case default
       stop "diagonalize_cluster error: vca_method != ['full','lanc']"
    end select
  end subroutine diagonalize_cluster








  !+-------------------------------------------------------------------+
  !PURPOSE  : diagonalize the Hamiltonian in each sector 
  !+------------------------------------------------------------------+
  subroutine vca_diag_h
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
    sector: do isector=1,Nsectors
       !
       Dim      = getdim(isector)
       !
       if(verbose==3)then
          nup  = getnup(isector)
          ndw  = getndw(isector)
          write(LOGfile,"(A,I4,A6,I2,A6,I2,A6,I15)")&
               "Solving sector:",isector,", nup:",nup,", ndw:",ndw,", dim=",getdim(isector)
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
    write(LOGfile,"(A)")"DIAG summary:"
    egs=minval(e0)
    open(free_unit(unit),file='egs'//reg(file_suffix)//".vca",position='append')
    do isector=1,Nsectors
       if(abs(e0(isector)-egs)>gs_threshold)cycle
       nup  = getnup(isector)
       ndw  = getndw(isector)
       write(LOGfile,"(A,F20.12,2I4)")'Egs =',e0(isector),nup,ndw
       write(unit,"(F20.12,2I4)")e0(isector),nup,ndw
    enddo
    close(unit)
    !
    !Get the partition function Z
    zeta_function=0.d0
    forall(isector=1:Nsectors)espace(isector)%e = espace(isector)%e - egs
    do isector=1,Nsectors
       dim=getdim(isector)
       do i=1,dim
          zeta_function=zeta_function+exp(-beta*espace(isector)%e(i))
       enddo
    enddo
    omega_potential = -1d0/beta*log(zeta_function)
    write(LOGfile,"(A,F20.12)")'Z     =',zeta_function
    write(LOGfile,"(A,F20.12)")'Omega =',omega_potential
    !
    return
  end subroutine vca_diag_h








  !+-------------------------------------------------------------------+
  !PURPOSE  : diagonalize the Hamiltonian in each sector and find the 
  ! spectrum DOUBLE COMPLEX
  !+------------------------------------------------------------------+
  subroutine vca_lanc_h
    integer                :: nup,ndw,isector,dim
    integer                :: isect,izero,sz,nt
    integer                :: i,j,iter,unit
    integer                :: Nitermax,Neigen,Nblock
    real(8)                :: oldzero,enemin,Ei,neigen_sector_error
    real(8),allocatable    :: eig_values(:)
    complex(8),allocatable :: eig_basis(:,:)
    integer                :: neigen_sector_prev(Nsectors)
    logical                :: lanc_solve,Tflag,converged_spectrum
    !
    if(state_list%status)call es_delete_espace(state_list)
    state_list=es_init_espace()
    oldzero=1000.d0
    write(LOGfile,"(A)")"Diagonalize impurity H:"
    call start_timer()
    !
    neigen_sector_prev=0
    converged_spectrum=.false.
    iter=0
    do while(.not.converged_spectrum.AND.iter<lanc_niter_spectrum)
       iter=iter+1
       call es_free_espace(state_list)
       print*,"Spectrum annealing iteration:",iter,"of ",lanc_niter_spectrum
       sector: do isector=1,Nsectors
          if(.not.twin_mask(isector))cycle sector !cycle loop if this sector should not be investigated
          Tflag    = twin_mask(isector).AND.diag_twin
          Tflag = Tflag.AND.(getnup(isector)/=getndw(isector))
          Dim      = getdim(isector)
          Neigen   = min(dim,neigen_sector(isector))
          Nitermax = min(dim,lanc_niter)
          Nblock   = min(dim,lanc_ncv_factor*Neigen + lanc_ncv_add)
          !
          lanc_solve  = .true.
          if(Neigen==dim)lanc_solve=.false.
          !
          if(verbose==3)then
             nup  = getnup(isector)
             ndw  = getndw(isector)
             write(LOGfile,"(A,I4,A6,I2,A6,I2,A6,I15,A12,3I6)")&
                  "Solving sector:",isector,", nup:",nup,", ndw:",ndw,", dim=",getdim(isector),", Lanc Info:",Neigen,Nitermax,Nblock
          elseif(verbose==1.OR.verbose==2)then
             call eta(iter,count(twin_mask),LOGfile)
          endif
          !
          call setup_Hv_sector(isector)
          if(lanc_solve)then
             if(allocated(eig_values))deallocate(eig_values)
             if(allocated(eig_basis))deallocate(eig_basis)
             allocate(eig_values(Neigen),eig_basis(Dim,Neigen))
             eig_values=0d0 ; eig_basis=zero
             call buildH_c()
             call sp_eigh(spHtimesV_cc,Dim,Neigen,Nblock,Nitermax,eig_values,eig_basis,tol=lanc_tolerance)
          else
             if(allocated(eig_values))deallocate(eig_values)
             if(allocated(eig_basis))deallocate(eig_basis)
             allocate(eig_values(Dim),eig_basis(Dim,dim))
             eig_values=0d0 ; eig_basis=zero
             call buildH_c(eig_basis)
             call eigh(eig_basis,eig_values,'V','U')
             if(dim==1)eig_basis(dim,dim)=one
          endif
          call delete_Hv_sector()
          !
          if(spH0%status)call sp_delete_matrix(spH0)
          !
          if(finiteT)then
             do i=1,Neigen
                call es_add_state(state_list,eig_values(i),eig_basis(1:dim,i),isector,twin=Tflag,size=lanc_nstates_total)
             enddo
          else
             do i=1,Neigen
                enemin = eig_values(i)
                if (enemin < oldzero-10.d0*gs_threshold)then
                   oldzero=enemin
                   call es_free_espace(state_list)
                   call es_add_state(state_list,enemin,eig_basis(1:dim,i),isector,twin=Tflag)
                elseif(abs(enemin-oldzero) <= gs_threshold)then
                   oldzero=min(oldzero,enemin)
                   call es_add_state(state_list,enemin,eig_basis(1:dim,i),isector,twin=Tflag)
                endif
             enddo
          endif
          unit=free_unit()
          open(unit,file="eigenvalues_list"//reg(file_suffix)//".vca",position='append',action='write')
          call print_eigenvalues_list(isector,eig_values(1:Neigen),unit)
          close(unit)
          !
          if(allocated(eig_values))deallocate(eig_values)
          if(allocated(eig_basis))deallocate(eig_basis)
          !
       enddo sector
       !
       call vca_lanc_analysis(.true.)
       if(finiteT)then
          neigen_sector_error = sum(abs(neigen_sector - neigen_sector_prev))/sum(abs(neigen_sector))
          converged_spectrum=(neigen_sector_error < lanc_spectrum_threshold)
          neigen_sector_prev=neigen_sector
       else
          !> at T=0 you only need to solve one: there is no annealing of the spectrum
          converged_spectrum=.true.
       endif
       !
    enddo
    !
    call stop_timer(LOGfile)
    !
  end subroutine vca_lanc_h



























  !###################################################################################################
  !
  !    POST-PROCESSING ROUTINES
  !
  !###################################################################################################
  !+-------------------------------------------------------------------+
  !PURPOSE  : analyse the spectrum and print some information after 
  !lanczos diagonalization. 
  !+------------------------------------------------------------------+
  subroutine vca_lanc_analysis(iprint)
    logical :: iprint
    integer             :: nup,ndw,sz,n,isector,dim
    integer             :: istate
    integer             :: i,unit
    integer             :: Nsize,NtoBremoved,nstates_below_cutoff
    integer             :: numgs
    real(8)             :: Egs,Ei,Ec,Etmp
    type(histogram)     :: h
    integer,allocatable :: list_sector(:),count_sector(:)
    !
    !POST PROCESSING:
    if(iprint)then
       unit=free_unit()
       open(unit,file="state_list"//reg(file_suffix)//".vca")
       call print_state_list(unit)
       close(unit)
       if(verbose>=2)call print_state_list(LOGfile)
    endif
    !
    zeta_function=0d0
    Egs = state_list%emin
    if(finiteT)then
       do i=1,state_list%size
          ei   = es_return_energy(state_list,i)
          zeta_function = zeta_function + exp(-beta*(Ei-Egs))
       enddo
    else
       zeta_function=real(state_list%size,8)
    end if
    !
    !
    numgs=es_return_gs_degeneracy(state_list,gs_threshold)
    if(numgs>Nsectors)stop "diag: too many gs"
    if(verbose>=2)then
       do istate=1,numgs
          isector = es_return_sector(state_list,istate)
          Egs     = es_return_energy(state_list,istate)
          nup  = getnup(isector)
          ndw  = getndw(isector)
          if(iprint)write(LOGfile,"(A,F20.12,2I4)")'Egs =',Egs,nup,ndw
       enddo
       if(iprint)write(LOGfile,"(A,F20.12)")'Z   =',zeta_function
    endif
    !
    !
    !    
    !get histogram distribution of the sector contributing to the evaluated spectrum:
    !go through states list and update the neigen_sector(isector) sector-by-sector
    if(finiteT)then
       !
       h = histogram_allocate(Nsectors)
       call lanc_build_histogram(h)
       if(iprint)then
          unit=free_unit()
          open(unit,file="histogram_states"//reg(file_suffix)//".vca",position='append')
          call histogram_print(h,unit)
          write(unit,*)""
          close(unit)
       endif
       call histogram_deallocate(h)
       !
       allocate(list_sector(state_list%size),count_sector(Nsectors))
       !get the list of actual sectors contributing to the list
       do i=1,state_list%size
          list_sector(i) = es_return_sector(state_list,i)
       enddo
       !count how many times a sector appears in the list
       do i=1,Nsectors
          count_sector(i) = count(list_sector==i)
       enddo
       !adapt the number of required Neig for each sector based on how many
       !appeared in the list.
       do i=1,Nsectors
          if(any(list_sector==i))then !if(count_sector(i)>1)then
             neigen_sector(i)=neigen_sector(i)+1
          else
             neigen_sector(i)=neigen_sector(i)-1
          endif
          !prevent Neig(i) from growing unbounded but 
          !try to put another state in the list from sector i
          if(neigen_sector(i) > count_sector(i))neigen_sector(i)=count_sector(i)+1
          if(neigen_sector(i) <= 0)neigen_sector(i)=1
       enddo
       !check if the number of states is enough to reach the required accuracy:
       !the condition to fullfill is:
       ! exp(-beta(Ec-Egs)) < \epsilon_c
       ! if this condition is violated then required number of states is increased
       ! if number of states is larger than those required to fullfill the cutoff: 
       ! trim the list and number of states.
       Egs  = state_list%emin
       Ec   = state_list%emax
       Nsize= state_list%size
       if(exp(-beta*(Ec-Egs)) > cutoff)then
          lanc_nstates_total=lanc_nstates_total + lanc_nstates_step
          if(iprint)write(LOGfile,"(A,I4)")"Increasing lanc_nstates_total:",lanc_nstates_total
       else
          ! !Find the energy level beyond which cutoff condition is verified & cut the list to that size
          if(iprint)write(LOGfile,*)
          isector = es_return_sector(state_list,state_list%size)
          Ei      = es_return_energy(state_list,state_list%size)
          do while ( exp(-beta*(Ei-Egs)) <= cutoff )
             if(verbose>=1.AND.iprint)write(LOGfile,"(A,I4,2x,2I4,2x,I5)")&
                  "Trimming state:",isector,getnup(isector),getndw(isector),state_list%size
             call es_pop_state(state_list)
             isector = es_return_sector(state_list,state_list%size)
             Ei      = es_return_energy(state_list,state_list%size)
          enddo
          if(verbose>=1.AND.iprint)then
             write(LOGfile,*)"Trimmed state list:"          
             call print_state_list(LOGfile)
          endif
          !
          lanc_nstates_total=max(state_list%size,lanc_nstates_step)+lanc_nstates_step
          write(LOGfile,"(A,I4)")"Adjusting lanc_nstates_total to:",lanc_nstates_total
          !
       endif
    endif
    write(LOGfile,*)""
    write(LOGfile,*)""
  end subroutine vca_lanc_analysis







  subroutine print_state_list(unit)
    integer :: nup,ndw,sz,n,isector
    integer :: istate
    integer :: unit
    real(8) :: Estate,Jz
    write(unit,"(A)")"# i       E_i           exp(-(E-E0)/T)       nup ndw  Sect     Dim"
    do istate=1,state_list%size
       Estate  = es_return_energy(state_list,istate)
       isector = es_return_sector(state_list,istate)
       nup   = getnup(isector)
       ndw   = getndw(isector)
       write(unit,"(i3,f18.12,2x,ES19.12,1x,2i3,3x,i3,i10)")&
            istate,Estate,exp(-beta*(Estate-state_list%emin)),nup,ndw,isector,getdim(isector)
    enddo
  end subroutine print_state_list


  subroutine print_eigenvalues_list(isector,eig_values,unit)
    integer              :: isector
    real(8),dimension(:) :: eig_values
    integer              :: unit,i
    write(unit,"(A7,A3,A3)")" # Sector","Nup","Ndw"
    write(unit,"(I4,2x,I3,I3)")isector,getnup(isector),getndw(isector)
    write(unit,*)""
  end subroutine print_eigenvalues_list





  subroutine lanc_build_histogram(hist)
    type(histogram) :: hist
    integer         :: i,isector
    call histogram_reset(hist)
    call histogram_set_range_uniform(hist,1d0,dble(Nsectors))
    do i=1,state_list%size
       isector = es_return_sector(state_list,i)
       call histogram_accumulate(hist,dble(isector),1d0)
    enddo
  end subroutine lanc_build_histogram


  ! !> this is unsafe: I am not adding check on allocation
  ! subroutine histogram_equality(h1,h2)
  !   type(histogram) :: h1
  !   type(histogram) :: h2
  !   h2%n = h1%n
  !   h2%range = h1%range
  !   h2%bin = h1%bin
  ! end subroutine histogram_equality


  ! function histogram_difference(h1,h2) result(error)
  !   type(histogram) :: h1
  !   type(histogram) :: h2
  !   real(8)         :: error,int1,int2
  !   integer         :: i
  !   real(8)         :: lower,upper,bin_value
  !   !
  !   int1=0d0
  !   do i=0,h1%n-1
  !      call histogram_get_range(h1,i,lower,upper)
  !      bin_value = histogram_get_value(h1,i)
  !      int1=int1+(upper-lower)*bin_value
  !   enddo
  !   !
  !   int2=0d0
  !   do i=0,h2%n-1
  !      call histogram_get_range(h2,i,lower,upper)
  !      bin_value = histogram_get_value(h2,i)
  !      int2=int2+(upper-lower)*bin_value
  !   enddo
  !   error = abs(int2-int1)/abs(int1+int2)
  ! end function histogram_difference



END MODULE VCA_DIAG













