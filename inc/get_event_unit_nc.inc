C -*- fortran -*-
         nu_mode = zz(120)
         detect_nu = zz(118)
         zz(150) = 0             ! NC interactions
         iproc = 3
         include 'get_event.inc'
         do i = 1,nbins
            hErec_dat(i,iproc,detect_nu,nu_mode,iD) = event(i)
         enddo
