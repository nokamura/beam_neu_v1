C -*- fortran -*-
      call get_xsecfrac3(Ev,icc,1,int_mode,inutype,frac)
      rievent2 = rievent*frac
      z(1) = inutype
      z(2) = int_mode
      z(3) = Ev
      call MakeHisto1D(nout,fErec_dist,z,rievent2,nbins,x,evform
     &     ,serror,snmax,ihisto,event_mode,hevent_mode
     &     ,rnevent_mode,ierr)
      if (ierr.ne.0) then
         write(*,*) "get_nudist_mode: MakeHisto1D ierr = 1"
      endif
