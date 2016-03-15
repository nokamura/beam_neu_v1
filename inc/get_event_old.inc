C -*- fortran -*-
            detect_nu = zz(118)
            call bining_basic(Emin,Emax,basic_binsize,nbins_basic,x_basic)
            call MakeHisto1D(nout,hfunc1D,zz,rnevent_ren,nbins_basic
     &           ,x_basic,evform,serror,snmax,ihisto,event1
     &           ,hevent1,nevent1,ierr) 
            if (zz(150).eq.1) then ! CC interactions
               if (iproc.eq.abs(detect_nu)) then
                  call get_nudist(detect_nu,evform,serror,snmax,ihisto
     &                 ,event1,x_basic,nbins_basic,event2,nevent2)
               endif
            elseif (zz(150).eq.0) then ! NC intractions
               if (iproc.eq.3) then
                  call get_1pi0dist(ipi0xsec,ipi0dist,detect_nu,evform
     &                 ,serror,snmax,ihisto,event1,nevent1,x_basic
     &                 ,nbins_basic,event2,nevent2)
               endif
            endif               
            ite = 0
            ibins = 1
            nevent = 0d0
            do i = 1,nbins
               event(i) = 0d0
            enddo
            do i = 1,nbins_basic
               ite = ite +1
               event(ibins) = event(ibins) +event2(i) 
               if (ite.eq.binsize_factor) then
                  ibins = ibins +1
                  ite = 0
               endif
            enddo
            nevent = nevent +nevent2
