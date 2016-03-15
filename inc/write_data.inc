C -*- fortran -*-
         open(1,file="event_total.dat",status="replace")
         do i = 0,evnbins(iev,iD)-1
            write(1,*) evx(i,iev,iD),hErec_total_dat(i+1)
         enddo
         close(1)

         open(1,file="event_pi0bg.dat",status="replace")
         do i = 0,evnbins(iev,iD)-1
            write(1,*) evx(i,iev,iD),hErec_pi0bg_dat(i+1,iD)
         enddo
         close(1)         
         open(1,file="event_pi0bg2.dat",status="replace")
         do i = 0,evnbins(iev,iD)-1
            write(1,*) evx(i,iev,iD),hErec_pi0bg2_dat(i+1,iD)
         enddo
         close(1)         

         open(1,file="event_nb.ne.ne.e.dat",status="replace")
         do i = 0,evnbins(iev,iD)-1
            write(1,*) evx(i,iev,iD),hErec_dat(i+1,1,1,1,iD,1)
         enddo
         close(1)
         open(1,file="event_nb.nm.nm.mu.dat",status="replace")
         do i = 0,evnbins(iev,iD)-1
            write(1,*) evx(i,iev,iD),hErec_dat(i+1,2,2,2,1)
         enddo
         close(1)
         open(1,file="event_nb.nm.nm.pi0.dat",status="replace")
         do i = 0,evnbins(iev,iD)-1
            write(1,*) evx(i,iev,iD),hErec_dat(i+1,3,2,2,iD,1)
         enddo
         close(1)
         open(1,file="event_nb.nm.nm.gam.dat",status="replace")
         do i = 0,evnbins(iev,iD)-1
            write(1,*) evx(i,iev,iD),hErec_dat(i+1,4,2,2,iD,1)
         enddo
         close(1)

         open(1,file="event_nb.nm.nm.e.dat",status="replace")
         do i = 0,evnbins(iev,iD)-1
            write(1,*) evx(i,iev,iD),hErec_dat(i+1,1,2,2,iD,1)
         enddo
         close(1)
         open(1,file="event_nb.nm.nm.mu.dat",status="replace")
         do i = 0,evnbins(iev,iD)-1
            write(1,*) evx(i,iev,iD),hErec_dat(i+1,2,2,2,1)
         enddo
         close(1)
         open(1,file="event_nb.nm.nm.pi0.dat",status="replace")
         do i = 0,evnbins(iev,iD)-1
            write(1,*) evx(i,iev,iD),hErec_dat(i+1,3,2,2,iD,1)
         enddo
         close(1)
         open(1,file="event_nb.nm.nm.gam.dat",status="replace")
         do i = 0,evnbins(iev,iD)-1
            write(1,*) evx(i,iev,iD),hErec_dat(i+1,4,2,2,iD,1)
         enddo
         close(1)

         open(1,file="event_nb.nm.ne.e.dat",status="replace")
         do i = 0,evnbins(iev,iD)-1
            write(1,*) evx(i,iev,iD),hErec_dat(i+1,1,1,2,iD,1)
         enddo
         close(1)
         open(1,file="event_nb.nm.ne.mu.dat",status="replace")
         do i = 0,evnbins(iev,iD)-1
            write(1,*) evx(i,iev,iD),hErec_dat(i+1,2,1,2,iD,1)
         enddo
         close(1)
         open(1,file="event_nb.nm.ne.pi0.dat",status="replace")
         do i = 0,evnbins(iev,iD)-1
            write(1,*) evx(i,iev,iD),hErec_dat(i+1,3,1,2,iD,1)
         enddo
         close(1)
         open(1,file="event_nb.nm.ne.gam.dat",status="replace")
         do i = 0,evnbins(iev,iD)-1
            write(1,*) evx(i,iev,iD),hErec_dat(i+1,3,3,2,iD,1)
         enddo
         close(1)

         open(1,file="event_ne.nm.e.dat",status="replace")
         do i = 0,evnbins(iev,iD)-1
            write(1,*) evx(i,iev,iD),hErec_dat(i+1,1,2,1,iD)
         enddo
         close(1)
         open(1,file="event_ne.nm.mu.dat",status="replace")
         do i = 0,evnbins(iev,iD)-1
            write(1,*) evx(i,iev,iD),hErec_dat(i+1,2,2,1,iD)
         enddo
         close(1)
         open(1,file="event_ne.nm.pi0.dat",status="replace")
         do i = 0,evnbins(iev,iD)-1
            write(1,*) evx(i,iev,iD),hErec_dat(i+1,3,2,1,iD)
         enddo
         close(1)
         open(1,file="event_ne.ne.e.dat",status="replace")
         do i = 0,evnbins(iev,iD)-1
            write(1,*) evx(i,iev,iD),hErec_dat(i+1,1,1,1,iD)
         enddo
         close(1)
         open(1,file="event_ne.ne.mu.dat",status="replace")
         do i = 0,evnbins(iev,iD)-1
            write(1,*) evx(i,iev,iD),hErec_dat(i+1,2,1,1,iD)
         enddo
         close(1)
         open(1,file="event_ne.ne.pi0.dat",status="replace")
         do i = 0,evnbins(iev,iD)-1
            write(1,*) evx(i,iev,iD),hErec_dat(i+1,3,1,1,iD)
         enddo
         close(1)
         open(1,file="event_ne.nt.pi0.dat",status="replace")
         do i = 0,evnbins(iev,iD)-1
            write(1,*) evx(i,iev,iD),hErec_dat(i+1,3,3,1,iD)
         enddo
         close(1)

         open(1,file="event_am.am.e.dat",status="replace")
         do i = 0,evnbins(iev,iD)-1
            write(1,*) evx(i,iev,iD),hErec_dat(i+1,1,-2,-2,iD)
         enddo
         close(1)
         open(1,file="event_am.am.mu.dat",status="replace")
         do i = 0,evnbins(iev,iD)-1
            write(1,*) evx(i,iev,iD),hErec_dat(i+1,2,-2,-2,iD)
         enddo
         close(1)
         open(1,file="event_am.am.pi0.dat",status="replace")
         do i = 0,evnbins(iev,iD)-1
            write(1,*) evx(i,iev,iD),hErec_dat(i+1,3,-2,-2,iD)
         enddo
         close(1)
         open(1,file="event_am.ae.e.dat",status="replace")
         do i = 0,evnbins(iev,iD)-1
            write(1,*) evx(i,iev,iD),hErec_dat(i+1,1,-1,-2,iD)
         enddo
         close(1)
         open(1,file="event_am.ae.mu.dat",status="replace")
         do i = 0,evnbins(iev,iD)-1
            write(1,*) evx(i,iev,iD),hErec_dat(i+1,2,-1,-2,iD)
         enddo
         close(1)
         open(1,file="event_am.ae.pi0.dat",status="replace")
         do i = 0,evnbins(iev,iD)-1
            write(1,*) evx(i,iev,iD),hErec_dat(i+1,3,-1,-2,iD)
         enddo
         close(1)
         open(1,file="event_am.at.pi0.dat",status="replace")
         do i = 0,evnbins(iev,iD)-1
            write(1,*) evx(i,iev,iD),hErec_dat(i+1,3,-3,-2,iD)
         enddo
         close(1)

         open(1,file="event_ae.am.e.dat",status="replace")
         do i = 0,evnbins(iev,iD)-1
            write(1,*) evx(i,iev,iD),hErec_dat(i+1,1,-2,-1,iD)
         enddo
         close(1)
         open(1,file="event_ae.am.mu.dat",status="replace")
         do i = 0,evnbins(iev,iD)-1
            write(1,*) evx(i,iev,iD),hErec_dat(i+1,2,-2,-1,iD)
         enddo
         close(1)
         open(1,file="event_ae.am.pi0.dat",status="replace")
         do i = 0,evnbins(iev,iD)-1
            write(1,*) evx(i,iev,iD),hErec_dat(i+1,3,-2,-1,iD)
         enddo
         close(1)
         open(1,file="event_ae.ae.e.dat",status="replace")
         do i = 0,evnbins(iev,iD)-1
            write(1,*) evx(i,iev,iD),hErec_dat(i+1,1,-1,-1,iD)
         enddo
         close(1)
         open(1,file="event_ae.ae.mu.dat",status="replace")
         do i = 0,evnbins(iev,iD)-1
            write(1,*) evx(i,iev,iD),hErec_dat(i+1,2,-1,-1,iD)
         enddo
         close(1)
         open(1,file="event_ae.ae.pi0.dat",status="replace")
         do i = 0,evnbins(iev,iD)-1
            write(1,*) evx(i,iev,iD),hErec_dat(i+1,3,-1,-1,iD)
         enddo
         close(1)
         open(1,file="event_ae.at.pi0.dat",status="replace")
         do i = 0,evnbins(iev,iD)-1
            write(1,*) evx(i,iev,iD),hErec_dat(i+1,3,-3,-1,iD)
         enddo
         close(1)
