C -*- fortran -*-
            call MakeHisto1D(nout,hfunc1D,zz,nevent,nbins,x
     &           ,evform,serror,snmax,ihisto,event_tmp,hevent_tmp
     &           ,nevent_tmp,ierr) 
c            call smearing_nl(event_tmp,x,nbins,rdx,Eres,Eres_nl,event_tmp2)         
            call smearing_nl2(ismear,event_tmp,x,nbins,Eres,Eres_nl,event_tmp2)         
            do i = 1,nbins
               event(i) = event(i) +event_tmp2(i)
            enddo
            nevent = nevent +nevent_tmp 
