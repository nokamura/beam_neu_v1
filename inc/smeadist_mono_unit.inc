C -*- fortran -*-
      do i = 1,nbins
         event(i) = 0d0
      enddo
      rtotevent = 0d0
      call get_xsecfrac(E,nutype,int_mode,frac)
      rnevent2 = rnevent*frac
      call get_nudist_mode(E,nutype,int_mode,evform,serror,snmax,x
     &     ,nbins,rnevent2,event_out,rnevent_out)
      do i = 1,nbins
         event(i) = event(i) +event_out(i)
      enddo
      rtotevent = rtotevent +rnevent_out