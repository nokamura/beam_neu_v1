C -*- fortran -*-
         nu_mode = zz(120)
         detect_nu = zz(118)

         zz(150) = 1             ! CC interactions
         if ((abs(detect_nu).eq.1).or.(abs(detect_nu).eq.2)) then
            iproc = abs(detect_nu)
c            call get_event(zz,iproc,rnevent_ren,evform,event,nevent)
            call get_event(zz,iproc,event,nevent)
            do i = 1,nbins
               hErec_dat(i,iproc,detect_nu,nu_mode,iD) = event(i)
            enddo
         endif

         zz(150) = 0             ! NC interactions
         iproc = 3
c         call get_event(zz,iproc,rnevent_ren,evform,event,nevent)
         call get_event(zz,iproc,event,nevent)
         do i = 1,nbins
            hErec_dat(i,iproc,detect_nu,nu_mode,iD) = event(i)
         enddo
