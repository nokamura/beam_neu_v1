C -*- fortran -*-
            write(*,*) "Events"
c            write(*,*)
            if (iNnm.eq.1) then
               nm = 1
               write(*,*) "nu_mu beam"
               write(*,*) "  mu-like:",nevent_dat(nm,iD)
               write(*,*) "    nbins =",outnbins(nm,iD)
               write(*,*) "  ",(outx(i,nm,iD),i=0,outnbins(nm,iD))
               write(*,*) "  ",(total_event_dat(i,nm,iD),i=1,evnbins(2,iD))
               write(*,*)
            endif
            if (iNne.eq.1) then
               nm = 2
               write(*,*) "  e-like:",nevent_dat(nm,iD)
               write(*,*) "    nbins =",outnbins(nm,iD)
               write(*,*) "  ",(outx(i,nm,iD),i=0,outnbins(nm,iD))
               write(*,*) "  ",(total_event_dat(i,nm,iD),i=1,evnbins(1,iD))
               write(*,*)
            endif
            if (iNam.eq.1) then
               nm = 3
               write(*,*) "anti nu_mu beam"
               write(*,*) "  mu-like:",nevent_dat(nm,iD)   
               write(*,*) "    nbins =",outnbins(nm,iD)
               write(*,*) "  ",(outx(i,nm,iD),i=0,outnbins(nm,iD))         
               write(*,*) "  ",(total_event_dat(i,nm,iD),i=1,evnbins(2,iD))      
               write(*,*)
            endif
            if (iNae.eq.1) then
               nm = 4
               write(*,*) "  e-like:",nevent_dat(nm,iD)
               write(*,*) "    nbins =",outnbins(nm,iD)
               write(*,*) "  ",(outx(i,nm,iD),i=0,outnbins(nm,iD))
               write(*,*) "  ",(total_event_dat(i,nm,iD),i=1,evnbins(1,iD))
               write(*,*)
            endif
c            write(*,*) "total:",nevent_tot(iD)
            write(*,*)
            

c$$$            write(*,*) "Event Histogram"
c$$$            if (iNnm.eq.1) then
c$$$               nm = 1
c$$$               write(*,*) "mu-like events ",
c$$$     &              "from nu_mu focusing beam (iNnm)"
c$$$               write(*,*) (total_event_dat(i,nm,iD),i=1,evnbins(2,iD))
c$$$               write(*,*)
c$$$               open(1,file="event_dat_Nnm.dat",status="replace")
c$$$               do i = 0,evnbins(2,iD)-1
c$$$                  write(1,*) evx(i,2,iD),total_event_dat(i+1,nm,iD)
c$$$               enddo
c$$$               close(1)
c$$$            endif               
c$$$            if (iNne.eq.1) then
c$$$               nm = 2
c$$$               write(*,*) "nu_e like events ",
c$$$     &              "from nu_mu focusing beam (iNne)"
c$$$               write(*,*) (total_event_dat(i,nm,iD),i=1,evnbins(1,iD))
c$$$               write(*,*)
c$$$               open(1,file="event_dat_Nne.dat",status="replace")
c$$$               do i = 0,evnbins(1,iD)-1
c$$$                  write(1,*) evx(i,1,iD),total_event_dat(i+1,nm,iD)
c$$$               enddo
c$$$               close(1)
c$$$            endif
c$$$            if (iNam.eq.1) then
c$$$               nm = 3
c$$$               write(*,*) "nu_mu-like events ",
c$$$     &              "from the anti-nu_mu focusing beam (iNam)"
c$$$               write(*,*) (total_event_dat(i,nm,iD),i=1,evnbins(2,iD))
c$$$               write(*,*)
c$$$               open(1,file="event_dat_Nam.dat",status="replace")
c$$$               do i = 0,evnbins(2,iD)-1
c$$$                  write(1,*) evx(i,2,iD),total_event_dat(i+1,nm,iD)
c$$$               enddo
c$$$               close(1)
c$$$            endif
c$$$            if (iNae.eq.1) then
c$$$               nm = 4
c$$$               write(*,*) "nu_e-like events ",
c$$$     &              "from the anti-nu_mu focusing beam (iNae)"
c$$$               write(*,*) (total_event_dat(i,nm,iD),i=1,evnbins(1,iD))
c$$$               write(*,*)
c$$$               open(1,file="event_dat_Nae.dat",status="replace")
c$$$               do i = 0,evnbins(1,iD)-1
c$$$                  write(1,*) evx(i,1,iD),total_event_dat(i+1,nm,iD)
c$$$               enddo
c$$$               close(1)
c$$$            endif            
c$$$
c$$$
c$$$            write(*,*) "Binning"
c$$$            if (iNnm.eq.1) then
c$$$               nm = 1
c$$$               write(*,*) "iNnm events"
c$$$               write(*,*) "nbins =",outnbins(nm,iD)
c$$$               write(*,*) (outx(i,nm,iD),i=0,outnbins(nm,iD))
c$$$               write(*,*)
c$$$            endif
c$$$            if (iNne.eq.1) then
c$$$               nm = 2
c$$$               write(*,*) "iNne events"
c$$$               write(*,*) "nbins =",outnbins(nm,iD)
c$$$               write(*,*) (outx(i,nm,iD),i=0,outnbins(nm,iD))
c$$$               write(*,*)
c$$$            endif
c$$$            if (iNam.eq.1) then
c$$$               nm = 3
c$$$               write(*,*) "iNam events"
c$$$               write(*,*) "nbins =",outnbins(nm,iD)
c$$$               write(*,*) (outx(i,nm,iD),i=0,outnbins(nm,iD))
c$$$               write(*,*)
c$$$            endif
c$$$            if (iNae.eq.1) then
c$$$               nm = 4
c$$$               write(*,*) "iNae events"
c$$$               write(*,*) "nbins =",outnbins(nm,iD)
c$$$               write(*,*) (outx(i,nm,iD),i=0,outnbins(nm,iD))
c$$$               write(*,*)
c$$$            endif

