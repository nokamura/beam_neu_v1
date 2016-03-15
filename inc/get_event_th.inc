C -*- fortran -*-
            call MakeHisto1D(nout,hfunc1D,z,nevent,nbins,x
     &           ,evform_dat,serror,snmax,ihisto,event_tmp,hevent_tmp
     &           ,nevent_tmp,ierr1) 
c            call smearing_nl(event_tmp,x,nbins,rdx,Eres,Eres_nl,event2_tmp)         
            call smearing_nl2(ismear,event_tmp,x,nbins,Eres,Eres_nl,event_tmp2)         
            do i = 1,nbins
               event_th(i) = event_th(i) +event_tmp2(i)
            enddo
            nevent_th = nevent_th +nevent_tmp 