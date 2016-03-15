      subroutine make_basic_data
C     ****************************************************
C     By Yoshitaro Takaesu @KIAS JAN 12 2014
C     ****************************************************
      implicit none
C     GLOBAL VARIABLES
      include 'inc/params.inc'
      include 'inc/main.inc'
      include 'inc/minfunc.inc'
C     CONSTANTS
C     ARGUMENTS 
C     LOCAL VARIABLES
      integer i
      integer binsize_factor_tmp,ihfunc_tmp,ihisto_tmp
      integer MH_tmp,nu_mode_tmp,detect_nu_tmp,ierr
      real*8 Emin_tmp,Emax_tmp,basic_binsize_tmp,L_tmp
      real*8 eventout(maxnbin),heventout(maxnbin),neventout
      real*8 event_tmp(maxnbin),hevent_tmp(maxnbin),nevent_tmp
      real*8 E,frac
      integer nbins_loc
      real*8 xmin,xmax,binsize_loc
      real*8 xl(0:maxnbin),yyl(0:maxnbin)
C     EXTERNAL FUNCTIONS
      real*8 hfunc1D
      external hfunc1D
C     ----------
C     BEGIN CODE
C     ----------
      basic_binsize_tmp = basic_binsize
      binsize_factor_tmp = binsize_factor
      MH_tmp = MHH
c      L_tmp = z_dat(116)
      L_tmp = L
      nu_mode_tmp = int(z_dat(120))
      detect_nu_tmp = int(z_dat(118))
c      ihfunc_tmp = int(z_dat(115))
      ihfunc_tmp = ihfunc
      ihisto_tmp = ihisto
      
      xmin = 0d0
      xmax = 5d0
      basic_binsize = 0.01
      binsize_factor = 5
      evform = evform_dat
c     
c     Flux	 
c     
c      z_dat(116) = 1d0 ! L
      L = 1d0 ! L
c     z_dat(116) = 295d0 ! L
c      z_dat(116) = 0.28d0
c      z_dat(115) = 1               ! ihfunc 
      ihfunc = 1               ! ihfunc 
      ihisto = 2
      z_dat(120) = 1  ! nu_e flux
      binsize_loc = basic_binsize*binsize_factor
      call bining_x(xmin,xmax,binsize_loc,nbins_loc,xl,yyl)
      call MakeHisto1D(nout,hfunc1D,z_dat,rnevent_ren,nbins_loc
     &     ,xl,evform,serror,snmax,ihisto,eventout
     &     ,heventout,neventout,ierr) 
      open(1,file="flux_ne.dat",status="replace")
      do i = 0,nbins_loc-1
         write(1,*) xl(i),eventout(i+1)
      enddo
      close(1)
      z_dat(120) = 2  ! nu_mode
      binsize_loc = basic_binsize*binsize_factor
      call bining_x(xmin,xmax,binsize_loc,nbins_loc,xl,yyl)
      call MakeHisto1D(nout,hfunc1D,z_dat,rnevent_ren,nbins_loc
     &     ,xl,evform,serror,snmax,ihisto,eventout
     &     ,heventout,neventout,ierr) 
      open(1,file="flux_nm.dat",status="replace")
      do i = 0,nbins_loc-1
         write(1,*) xl(i),eventout(i+1)
      enddo
      close(1)
      z_dat(120) = -1  ! nu_mode
      binsize_loc = basic_binsize*binsize_factor
      call bining_x(xmin,xmax,binsize_loc,nbins_loc,xl,yyl)
      call MakeHisto1D(nout,hfunc1D,z_dat,rnevent_ren,nbins_loc
     &     ,xl,evform,serror,snmax,ihisto,eventout
     &     ,heventout,neventout,ierr) 
      open(1,file="flux_ae.dat",status="replace")
      do i = 0,nbins_loc-1
         write(1,*) xl(i),eventout(i+1)
      enddo
      close(1)
      z_dat(120) = -2  ! nu_mode
      binsize_loc = basic_binsize*binsize_factor
      call bining_x(xmin,xmax,binsize_loc,nbins_loc,xl,yyl)
      call MakeHisto1D(nout,hfunc1D,z_dat,rnevent_ren,nbins_loc
     &     ,xl,evform,serror,snmax,ihisto,eventout
     &     ,heventout,neventout,ierr) 
      open(1,file="flux_am.dat",status="replace")
      do i = 0,nbins_loc-1
         write(1,*) xl(i),eventout(i+1)
      enddo
      close(1)
c     
c     Transition Probability
c     
      binsize_factor = 1
      binsize_loc = basic_binsize*binsize_factor
      call bining_x(xmin,xmax,binsize_loc,nbins_loc,xl,yyl)
c      z_dat(116) = L_tmp           ! L
      L = L_tmp           ! L
c      z_dat(115) = 2  ! ihfunc
      ihfunc = 2  ! ihfunc
      ihisto = 0
      z_dat(120) = 1 ! nu_mode
      z_dat(118) = 1 ! detect
      call MakeHisto1D(nout,hfunc1D,z_dat,rnevent_ren,nbins_loc
     &     ,xl,evform,serror,snmax,ihisto,eventout
     &     ,heventout,neventout,ierr) 
      open(1,file="prob_ne.ne.dat",status="replace")
      do i = 0,nbins_loc-1
         write(1,*) xl(i),eventout(i+1)
      enddo
      close(1)
      z_dat(118) = 2
      call MakeHisto1D(nout,hfunc1D,z_dat,rnevent_ren,nbins_loc
     &     ,xl,evform,serror,snmax,ihisto,eventout
     &     ,heventout,neventout,ierr) 
      open(1,file="prob_ne.nm.dat",status="replace")
      do i = 0,nbins_loc-1
         write(1,*) xl(i),eventout(i+1)
      enddo
      close(1)
      z_dat(118) = 3
      call MakeHisto1D(nout,hfunc1D,z_dat,rnevent_ren,nbins_loc
     &     ,xl,evform,serror,snmax,ihisto,eventout
     &     ,heventout,neventout,ierr) 
      open(1,file="prob_ne.nt.dat",status="replace")
      do i = 0,nbins_loc-1
         write(1,*) xl(i),eventout(i+1)
      enddo
      close(1)
      z_dat(120) = 2
      z_dat(118) = 1
      call MakeHisto1D(nout,hfunc1D,z_dat,rnevent_ren,nbins_loc
     &     ,xl,evform,serror,snmax,ihisto,eventout
     &     ,heventout,neventout,ierr) 
      open(1,file="prob_nm.ne.dat",status="replace")
      do i = 0,nbins_loc-1
         write(1,*) xl(i),eventout(i+1)
      enddo
      close(1)
      z_dat(118) = 2
      call MakeHisto1D(nout,hfunc1D,z_dat,rnevent_ren,nbins_loc
     &     ,xl,evform,serror,snmax,ihisto,eventout
     &     ,heventout,neventout,ierr) 
      open(1,file="prob_nm.nm.dat",status="replace")
      do i = 0,nbins_loc-1
         write(1,*) xl(i),eventout(i+1)
      enddo
      close(1)
      z_dat(118) = 3
      call MakeHisto1D(nout,hfunc1D,z_dat,rnevent_ren,nbins_loc
     &     ,xl,evform,serror,snmax,ihisto,eventout
     &     ,heventout,neventout,ierr) 
      open(1,file="prob_nm.nt.dat",status="replace")
      do i = 0,nbins_loc-1
         write(1,*) xl(i),eventout(i+1)
      enddo
      close(1)
      z_dat(120) = 3
      z_dat(118) = 1
      call MakeHisto1D(nout,hfunc1D,z_dat,rnevent_ren,nbins_loc
     &     ,xl,evform,serror,snmax,ihisto,eventout
     &     ,heventout,neventout,ierr) 
      open(1,file="prob_nt.ne.dat",status="replace")
      do i = 0,nbins_loc-1
         write(1,*) xl(i),eventout(i+1)
      enddo
      close(1)
      z_dat(118) = 2
      call MakeHisto1D(nout,hfunc1D,z_dat,rnevent_ren,nbins_loc
     &     ,xl,evform,serror,snmax,ihisto,eventout
     &     ,heventout,neventout,ierr) 
      open(1,file="prob_nt.nm.dat",status="replace")
      do i = 0,nbins_loc-1
         write(1,*) xl(i),eventout(i+1)
      enddo
      close(1)
      z_dat(118) = 3
      call MakeHisto1D(nout,hfunc1D,z_dat,rnevent_ren,nbins_loc
     &     ,xl,evform,serror,snmax,ihisto,eventout
     &     ,heventout,neventout,ierr) 
      open(1,file="prob_nt.nt.dat",status="replace")
      do i = 0,nbins_loc-1
         write(1,*) xl(i),eventout(i+1)
      enddo
      close(1)
      z_dat(120) = -1
      z_dat(118) = -1
      call MakeHisto1D(nout,hfunc1D,z_dat,rnevent_ren,nbins_loc
     &     ,xl,evform,serror,snmax,ihisto,eventout
     &     ,heventout,neventout,ierr) 
      open(1,file="prob_ae.ae.dat",status="replace")
      do i = 0,nbins_loc-1
         write(1,*) xl(i),eventout(i+1)
      enddo
      close(1)
      z_dat(118) = -2
      call MakeHisto1D(nout,hfunc1D,z_dat,rnevent_ren,nbins_loc
     &     ,xl,evform,serror,snmax,ihisto,eventout
     &     ,heventout,neventout,ierr) 
      open(1,file="prob_ae.am.dat",status="replace")
      do i = 0,nbins_loc-1
         write(1,*) xl(i),eventout(i+1)
      enddo
      close(1)
      z_dat(118) = -3
      call MakeHisto1D(nout,hfunc1D,z_dat,rnevent_ren,nbins_loc
     &     ,xl,evform,serror,snmax,ihisto,eventout
     &     ,heventout,neventout,ierr) 
      open(1,file="prob_ae.at.dat",status="replace")
      do i = 0,nbins_loc-1
         write(1,*) xl(i),eventout(i+1)
      enddo
      close(1)
      z_dat(120) = -2
      z_dat(118) = -1
      call MakeHisto1D(nout,hfunc1D,z_dat,rnevent_ren,nbins_loc
     &     ,xl,evform,serror,snmax,ihisto,eventout
     &     ,heventout,neventout,ierr) 
      open(1,file="prob_am.ae.dat",status="replace")
      do i = 0,nbins_loc-1
         write(1,*) xl(i),eventout(i+1)
      enddo
      close(1)
      z_dat(118) = -2
      call MakeHisto1D(nout,hfunc1D,z_dat,rnevent_ren,nbins_loc
     &     ,xl,evform,serror,snmax,ihisto,eventout
     &     ,heventout,neventout,ierr) 
      open(1,file="prob_am.am.dat",status="replace")
      do i = 0,nbins_loc-1
         write(1,*) xl(i),eventout(i+1)
      enddo
      close(1)
      z_dat(118) = -3
      call MakeHisto1D(nout,hfunc1D,z_dat,rnevent_ren,nbins_loc
     &     ,xl,evform,serror,snmax,ihisto,eventout
     &     ,heventout,neventout,ierr) 
      open(1,file="prob_am.at.dat",status="replace")
      do i = 0,nbins_loc-1
         write(1,*) xl(i),eventout(i+1)
      enddo
      close(1)
      z_dat(120) = -3
      z_dat(118) = -1
      call MakeHisto1D(nout,hfunc1D,z_dat,rnevent_ren,nbins_loc
     &     ,xl,evform,serror,snmax,ihisto,eventout
     &     ,heventout,neventout,ierr) 
      open(1,file="prob_at.ae.dat",status="replace")
      do i = 0,nbins_loc-1
         write(1,*) xl(i),eventout(i+1)
      enddo
      close(1)
      z_dat(118) = -2
      call MakeHisto1D(nout,hfunc1D,z_dat,rnevent_ren,nbins_loc
     &     ,xl,evform,serror,snmax,ihisto,eventout
     &     ,heventout,neventout,ierr) 
      open(1,file="prob_at.am.dat",status="replace")
      do i = 0,nbins_loc-1
         write(1,*) xl(i),eventout(i+1)
      enddo
      close(1)
      z_dat(118) = -3
      call MakeHisto1D(nout,hfunc1D,z_dat,rnevent_ren,nbins_loc
     &     ,xl,evform,serror,snmax,ihisto,eventout
     &     ,heventout,neventout,ierr) 
      open(1,file="prob_at.at.dat",status="replace")
      do i = 0,nbins_loc-1
         write(1,*) xl(i),eventout(i+1)
      enddo
      close(1)
c     
c     Xsec
c     
c      z_dat(115) = 3
      ihfunc = 3
      ihisto = 1
      z_dat(150) = 1               ! CC
      z_dat(118) = 1               ! detected neutrino
      binsize_loc = basic_binsize*binsize_factor
      call bining_x(xmin,xmax,binsize_loc,nbins_loc,xl,yyl)
      call MakeHisto1D(nout,hfunc1D,z_dat,rnevent_ren,nbins_loc
     &     ,xl,evform,serror,snmax,ihisto,eventout
     &     ,heventout,neventout,ierr) 
      open(1,file="xsec_cc_ne.dat",status="replace")
      do i = 0,nbins_loc-1
         write(1,*) xl(i),eventout(i+1)
      enddo
      close(1)
      z_dat(118) = 2               ! detected neutrino
      call MakeHisto1D(nout,hfunc1D,z_dat,rnevent_ren,nbins_loc
     &     ,xl,evform,serror,snmax,ihisto,eventout
     &     ,heventout,neventout,ierr) 
      open(1,file="xsec_cc_nm.dat",status="replace")
      do i = 0,nbins_loc-1
         write(1,*) xl(i),eventout(i+1)
      enddo
      close(1)
      z_dat(118) = -1              ! detected neutrino
      call MakeHisto1D(nout,hfunc1D,z_dat,rnevent_ren,nbins_loc
     &     ,xl,evform,serror,snmax,ihisto,eventout
     &     ,heventout,neventout,ierr) 
      open(1,file="xsec_cc_ae.dat",status="replace")
      do i = 0,nbins_loc-1
         write(1,*) xl(i),eventout(i+1)
      enddo
      close(1)
      z_dat(118) = -2              ! detected neutrino
      call MakeHisto1D(nout,hfunc1D,z_dat,rnevent_ren,nbins_loc
     &     ,xl,evform,serror,snmax,ihisto,eventout
     &     ,heventout,neventout,ierr) 
      open(1,file="xsec_cc_am.dat",status="replace")
      do i = 0,nbins_loc-1
         write(1,*) xl(i),eventout(i+1)
      enddo
      close(1)
      
      z_dat(150) = 0               ! NC
      z_dat(118) = 1               ! detected neutrino
      call MakeHisto1D(nout,hfunc1D,z_dat,rnevent_ren,nbins_loc
     &     ,xl,evform,serror,snmax,ihisto,eventout
     &     ,heventout,neventout,ierr) 
      open(1,file="xsec_nc_ne.dat",status="replace")
      do i = 0,nbins_loc-1
         write(1,*) xl(i),eventout(i+1)
      enddo
      close(1)
      z_dat(118) = 2               ! detected neutrino
      call MakeHisto1D(nout,hfunc1D,z_dat,rnevent_ren,nbins_loc
     &     ,xl,evform,serror,snmax,ihisto,eventout
     &     ,heventout,neventout,ierr) 
      open(1,file="xsec_nc_nm.dat",status="replace")
      do i = 0,nbins_loc-1
         write(1,*) xl(i),eventout(i+1)
      enddo
      close(1)
      z_dat(118) = -1              ! detected neutrino
      call MakeHisto1D(nout,hfunc1D,z_dat,rnevent_ren,nbins_loc
     &     ,xl,evform,serror,snmax,ihisto,eventout
     &     ,heventout,neventout,ierr) 
      open(1,file="xsec_nc_ae.dat",status="replace")
      do i = 0,nbins_loc-1
         write(1,*) xl(i),eventout(i+1)
      enddo
      close(1)
      z_dat(118) = -2              ! detected neutrino
      call MakeHisto1D(nout,hfunc1D,z_dat,rnevent_ren,nbins_loc
     &     ,xl,evform,serror,snmax,ihisto,eventout
     &     ,heventout,neventout,ierr) 
      open(1,file="xsec_nc_am.dat",status="replace")
      do i = 0,nbins_loc-1
         write(1,*) xl(i),eventout(i+1)
      enddo
      close(1)
      
      
c     
c     Flux*P*xsec
c     
c      z_dat(116) = 1000
      L = 1000
      z_dat(120) = 2
      z_dat(118) = 3
c      z_dat(111) = 1
      MHH = 1
      basic_binsize = 0.01
      binsize_factor = 5
c      z_dat(115) = 6               ! ihfunc 
      ihfunc = 4               ! ihfunc 
      ihisto = 2
      z_dat(150) = 0
      binsize_loc = basic_binsize*binsize_factor
      call bining_x(xmin,xmax,binsize_loc,nbins_loc,xl,yyl)
      call MakeHisto1D(nout,hfunc1D,z_dat,rnevent_ren,nbins_loc
     &     ,xl,evform,serror,snmax,ihisto,eventout
     &     ,heventout,neventout,ierr) 
      open(1,file="flux_P_xsec.dat",status="replace")
      do i = 0,nbins_loc-1
         write(1,*) xl(i),eventout(i+1)
      enddo
      close(1)
      
      
c     
c     Events after 1pi0 cut
c     
c      z_dat(116) = 1000
      L = 1000
      z_dat(120) = 2
      z_dat(118) = 3
      MHH = 1
      basic_binsize = 0.01
      binsize_factor = 5
c      z_dat(115) = 6               ! ihfunc 
      ihfunc = 6               ! ihfunc 
      ihisto = 2
      z_dat(150) = 0
      binsize_loc = basic_binsize*binsize_factor
      call bining_x(xmin,xmax,binsize_loc,nbins_loc,xl,yyl)
      call MakeHisto1D(nout,hfunc1D,z_dat,rnevent_ren,nbins_loc
     &     ,xl,evform,serror,snmax,ihisto,event_tmp
     &     ,hevent_tmp,nevent_tmp,ierr) 
      do i = 1,nbins_loc
         E = ( xl(i) +xl(i-1) )/2d0
         call get_xsecfrac3(E,int(z_dat(150)),0,5,int(z_dat(118)),frac)
         eventout(i) = frac*event_tmp(i)
      enddo
      open(1,file="pi0event.dat",status="replace")
      do i = 0,nbins_loc-1
         write(1,*) xl(i),eventout(i+1)
      enddo
      close(1)
      
c      z_dat(116) = L_tmp
      L = L_tmp
      z_dat(120) = nu_mode_tmp
      z_dat(118) = detect_nu_tmp
c      z_dat(111) = MH_tmp
      MHH = MH_tmp
      basic_binsize = basic_binsize_tmp
      binsize_factor = binsize_factor_tmp
c      z_dat(115) = ihfunc_tmp
      ihfunc = ihfunc_tmp
      ihisto = ihisto_tmp
      
      return
      end
