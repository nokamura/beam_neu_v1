C -*- fortran -*-
         z_dat(15) = 0 ! N vs. sqrt(E)
         z(15) = z_dat(15)
         evform_dat = 2
         evform_th = 2

         include 'bining_x.inc'
c         include 'bining_sqrtx_Eres.inc'

         if (ifirst.eq.0) then
            call MakeHisto1D_stat(nout,hfunc1D,z_dat,nevent,nbins,x
c            call MakeHisto1D(nout,hfunc1D,z_dat,nevent,nbins,x
     &           ,evform_dat,serror,snmax,ihisto,event_dat,hevent_dat
     &           ,nevent_dat,ierr1) 
c            call smearing_nl(event_dat,x,nbins,rdx,Eres,Eres_nl,event2_dat)         
            call smearing_nl2(ismear,event_dat,x,nbins,Eres,Eres_nl,event2_dat)         
            do i = 1,nbins
               event2_dat(i) = event2_dat(i) +dsqrt(event2_dat(i))*gran()
               if (event2_dat(i).lt.0) event2_dat(i) = 0
            enddo
            ifirst = 1
         endif

         call MakeHisto1D(nout,hfunc1D,z,nevent,nbins,x
     &        ,evform_th,serror,snmax,ihisto,event_th,hevent_th
     &        ,nevent_th,ierr2)
c         call smearing_nl(event_th,x,nbins,rdx,Eres,Eres_nl,event2_th)         
         call smearing_nl2(ismear,event_dat,x,nbins,Eres,Eres_nl,event2_dat)         

         chisq = chi2_2(nout,event2_dat,event2_th,x,nbins,npar,z
     &        ,z_dat,error)  
