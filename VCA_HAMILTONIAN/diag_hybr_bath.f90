  
  if(allocated(diag_hybr))deallocate(diag_hybr)
  if(allocated(bath_diag))deallocate(bath_diag)

  allocate(diag_hybr(Nlat,Nlat_bath,Nspin,Norb,Norb_bath));diag_hybr=0d0
  allocate(bath_diag(Nlat_bath,Nspin,Norb_bath));bath_diag=0d0
  do ispin=1,Nspin
    do ilat_bath=1,Nlat_bath
      do iorb_bath=1,Norb_bath
        bath_diag(ilat_bath,ispin,iorb_bath)=vca_bath%h(ilat_bath,ilat_bath,ispin,ispin,iorb_bath,iorb_bath)
        do ilat=1,Nlat             
          do iorb=1,Norb
          diag_hybr(ilat,ilat_bath,ispin,iorb,iorb_bath)=vca_bath%v(ilat,ilat_bath,ispin,ispin,iorb,iorb_bath)
          enddo
        enddo
      enddo
    enddo
  enddo

