C -*- fortran -*-
      call get_xsecfrac3(Ev,icc,1,int_mode,inutype,frac)
      if ((inutype.gt.0).and.(int_mode.eq.3)) then   !  We ignore the neutrino CCRes-H interaction as an approximation. (To be considered later)
         frac = 0d0
      endif
      rievent2 = rievent*frac
      z(1) = inutype
      z(2) = int_mode
      z(3) = Ev
      if (ismear.eq.1) then  ! apply detector smearing effects.
         call MakeHisto1D(nout,fErec_dist,z,rievent2,nnbins,xx,evform
     &        ,serror,snmax,ihisto,event_mode,hevent_mode
     &        ,rnevent_mode,ierr)
         if (ierr.ne.0) then
            write(*,*) "get_nudist_mode: MakeHisto1D ierr = 1"
         endif
      elseif (ismear.eq.0) then
         do j = 1,nnbins
            event_mode(j) = 0d0
         enddo
         event_mode(i) = rievent2
      endif
