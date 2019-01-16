  subroutine vca_gf_cluster_scalar(zeta,gf)
    complex(8)                                                          :: zeta
    complex(8),dimension(Nlat,Nlat,Nspin,Nspin,Norb,Norb),intent(inout) :: gf
    complex(8)                                                          :: green
    integer                                                             :: ispin
    integer                                                             :: ilat,jlat
    integer                                                             :: iorb,jorb
    integer                                                             :: iexc,Nexc
    integer                                                             :: ichan,Nchannel
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
                Nchannel = size(impGmatrix(ilat,jlat,ispin,ispin,iorb,iorb)%channel)
                do ichan=1,Nchannel
                   Nexc  = size(impGmatrix(ilat,jlat,ispin,ispin,iorb,iorb)%channel(ichan)%poles)
                   do iexc=1,Nexc
                      weight = impGmatrix(ilat,jlat,ispin,ispin,iorb,iorb)%channel(ichan)%weight(iexc)
                      de     = impGmatrix(ilat,jlat,ispin,ispin,iorb,iorb)%channel(ichan)%poles(iexc)
                      green = green + weight/(zeta-de)
                   enddo
                enddo
                gf(ilat,jlat,ispin,ispin,iorb,iorb) = GS_MULT*green
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


  subroutine vca_gf_hyb_scalar(zeta,gf)
    complex(8)                                                                :: zeta
    complex(8),dimension(Nlat,Nlat*Nbath,Nspin,Nspin,Norb,Norb),intent(inout) :: gf
    complex(8)                                                                :: green
    integer                                                                   :: ispin
    integer                                                                   :: ilat,jlat
    integer                                                                   :: iorb,jorb
    integer                                                                   :: iexc,Nexc
    integer                                                                   :: ichan,Nchannel
    integer                                                                   :: i,is,js
    real(8)                                                                   :: weight,de
    !
    if(.not.allocated(hybGmatrix))stop "vca_gf_cluster ERROR: impGmatrix not allocated!"
    !
    gf = zero
    !
    do ilat=1,Nlat
       do jlat=1,Nlat*Nbath
          do iorb=1,Norb
             do ispin=1,Nspin
                !
                green = zero
                Nchannel = size(hybGmatrix(ilat,jlat,ispin,ispin,iorb,iorb)%channel)
                do ichan=1,Nchannel
                   Nexc  = size(hybGmatrix(ilat,jlat,ispin,ispin,iorb,iorb)%channel(ichan)%poles)
                   do iexc=1,Nexc
                      weight = hybGmatrix(ilat,jlat,ispin,ispin,iorb,iorb)%channel(ichan)%weight(iexc)
                      de     = hybGmatrix(ilat,jlat,ispin,ispin,iorb,iorb)%channel(ichan)%poles(iexc)
                      green = green + weight/(zeta-de)
                   enddo
                enddo
                gf(ilat,jlat,ispin,ispin,iorb,iorb) = green
             enddo
          enddo
       enddo
    enddo
    !
  end subroutine vca_gf_hyb_scalar

 subroutine vca_gf_bath_scalar(zeta,gf)
    complex(8)                                                                      :: zeta
    complex(8),dimension(Nlat*Nbath,Nlat*Nbath,Nspin,Nspin,Norb,Norb),intent(inout) :: gf
    complex(8)                                                                      :: green
    integer                                                                         :: ispin
    integer                                                                         :: ilat,jlat
    integer                                                                         :: iorb,jorb
    integer                                                                         :: iexc,Nexc
    integer                                                                         :: ichan,Nchannel
    integer                                                                         :: i,is,js
    real(8)                                                                         :: weight,de
    !
    if(.not.allocated(bathGmatrix))stop "vca_gf_cluster ERROR: impGmatrix not allocated!"
    !
    gf = zero
    !
    do ilat=1,Nlat*Nbath
       do jlat=1,Nlat*Nbath
          do iorb=1,Norb
             do ispin=1,Nspin
                !
                green = zero
                Nchannel = size(bathGmatrix(ilat,jlat,ispin,ispin,iorb,iorb)%channel)
                do ichan=1,Nchannel
                   Nexc  = size(bathGmatrix(ilat,jlat,ispin,ispin,iorb,iorb)%channel(ichan)%poles)
                   do iexc=1,Nexc
                      weight = bathGmatrix(ilat,jlat,ispin,ispin,iorb,iorb)%channel(ichan)%weight(iexc)
                      de     = bathGmatrix(ilat,jlat,ispin,ispin,iorb,iorb)%channel(ichan)%poles(iexc)
                      green = green + weight/(zeta-de)
                   enddo
                enddo
                gf(ilat,jlat,ispin,ispin,iorb,iorb) = green
             enddo
          enddo
       enddo
    enddo
    do ispin=1,Nspin
       do iorb=1,Norb
          do ilat=1,Nlat*Nbath
             do jlat=1,Nlat*Nbath
                if(ilat==jlat)cycle
                gf(ilat,jlat,ispin,ispin,iorb,iorb) = 0.5d0*(gf(ilat,jlat,ispin,ispin,iorb,iorb) &
                     - gf(ilat,ilat,ispin,ispin,iorb,iorb) - gf(jlat,jlat,ispin,ispin,iorb,iorb))    
             enddo
          enddo
       enddo
    enddo
    !
  end subroutine vca_gf_bath_scalar


  subroutine vca_gf_cluster_array(zeta,gf)
    complex(8),dimension(:)                                                        :: zeta
    complex(8),dimension(Nlat,Nlat,Nspin,Nspin,Norb,Norb,size(zeta)),intent(inout) :: gf
    complex(8),dimension(Nlat,Nlat,Nspin,Nspin,Norb,Norb)                          :: green
    integer                                                                        :: ispin
    integer                                                                        :: ilat,jlat
    integer                                                                        :: iorb,jorb
    integer                                                                        :: iexc,Nexc
    integer                                                                        :: ichan,Nchannel
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
