C -*- fortran -*-
      do i = 1,eventdim
         chisq_stat(i,iD) = 0d0
      enddo

      if ( iNnm.eq.1) then
         call ret_bins(2,evnbins(1,iD),evx(0,1,iD),nbins,x)
         chisq_stat(1,iD) = chi2_stat(total_event_dat(1,1,iD)
     &        ,total_event_th(1,1,iD),nbins,x,outnbins(1,iD)
     &        ,outx(0,1,iD))  
      endif
      if ( iNne.eq.1) then
         call ret_bins(1,evnbins(1,iD),evx(0,1,iD),nbins,x)
         chisq_stat(2,iD) = chi2_stat(total_event_dat(1,2,iD)
     &        ,total_event_th(1,2,iD),nbins,x,outnbins(2,iD)
     &        ,outx(0,2,iD))
      endif
      if ( iNam.eq.1) then
         call ret_bins(2,evnbins(1,iD),evx(0,1,iD),nbins,x)
         chisq_stat(3,iD) = chi2_stat(total_event_dat(1,3,iD)
     &        ,total_event_th(1,3,iD),nbins,x,outnbins(3,iD)
     &        ,outx(0,3,iD))
      endif
      if ( iNae.eq.1) then
         call ret_bins(1,evnbins(1,iD),evx(0,1,iD),nbins,x)
         chisq_stat(4,iD) = chi2_stat(total_event_dat(1,4,iD)
     &        ,total_event_th(1,4,iD),nbins,x,outnbins(4,iD)
     &        ,outx(0,4,iD))
      endif
