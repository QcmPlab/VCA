  subroutine vca_gf_cluster_scalar(zeta,gf)
    complex(8)                                                          :: zeta
    complex(8),dimension(Nlat,Nlat,Nspin,Nspin,Norb,Norb),intent(inout) :: gf
    complex(8)                                                          :: green
    integer                                                             :: ispin
    integer                                                             :: ilat,jlat
    integer                                                             :: iorb,jorb
    integer                                                             :: iexc,Nexc
    integer                                                             :: ichan,Nchannel,istate,Nstates
    integer                                                             :: i,is,js
    real(8)                                                             :: weight,de
    !
    if(.not.allocated(impGmatrix))stop "vca_gf_cluster ERROR: impGmatrix not allocated!"
    !
    gf = zero
    !
    do ilat=1,Nlat
       do jlat=1,Nlat
          do iorb=1,Norb
             do ispin=1,Nspin
                !
                green = zero
                Nstates = size(impGmatrix(ilat,jlat,ispin,ispin,iorb,iorb)%state)
                do istate=1,Nstates
                  Nchannel = size(impGmatrix(ilat,jlat,ispin,ispin,iorb,iorb)%state(istate)%channel)
                  do ichan=1,Nchannel
                     Nexc  = size(impGmatrix(ilat,jlat,ispin,ispin,iorb,iorb)%state(istate)%channel(ichan)%poles)
                     if(Nexc .ne. 0)then
                       do iexc=1,Nexc
                          weight = impGmatrix(ilat,jlat,ispin,ispin,iorb,iorb)%state(istate)%channel(ichan)%weight(iexc)
                          de     = impGmatrix(ilat,jlat,ispin,ispin,iorb,iorb)%state(istate)%channel(ichan)%poles(iexc)
                          green = green + weight/(zeta-de)
                       enddo
                    endif
                  enddo
                enddo
                gf(ilat,jlat,ispin,ispin,iorb,iorb) = green
             enddo
          enddo
       enddo
    enddo
    do ispin=1,Nspin
       do iorb=1,Norb
          do ilat=1,Nlat
             do jlat=1,Nlat
                if(ilat==jlat)cycle
                gf(ilat,jlat,ispin,ispin,iorb,iorb) = 0.5d0*(gf(ilat,jlat,ispin,ispin,iorb,iorb) &
                     - gf(ilat,ilat,ispin,ispin,iorb,iorb) - gf(jlat,jlat,ispin,ispin,iorb,iorb))    
             enddo
          enddo
       enddo
    enddo
    !
  end subroutine vca_gf_cluster_scalar



  subroutine vca_gf_cluster_array(zeta,gf)
    complex(8),dimension(:)                                                        :: zeta
    complex(8),dimension(Nlat,Nlat,Nspin,Nspin,Norb,Norb,size(zeta)),intent(inout) :: gf
    complex(8),dimension(Nlat,Nlat,Nspin,Nspin,Norb,Norb)                          :: green
    integer                                                                        :: ispin
    integer                                                                        :: ilat,jlat
    integer                                                                        :: iorb,jorb
    integer                                                                        :: iexc,Nexc
    integer                                                                        :: ichan,Nchannel,istate,Nstates
    integer                                                                        :: i,is,js
    real(8)                                                                        :: weight,de
    !
    if(.not.allocated(impGmatrix))stop "vca_gf_cluster ERROR: impGmatrix not allocated!"
    !
    gf = zero
    do i=1,size(zeta)
       call vca_gf_cluster_scalar(zeta(i),green)
       gf(:,:,:,:,:,:,i) = green
    enddo
    !
  end subroutine vca_gf_cluster_array
