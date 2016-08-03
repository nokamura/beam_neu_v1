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
      real*8 Emin_tmp,Emax_tmp,basic_binsize_tmp,L_tmp,V_tmp,icc_tmp
      real*8 eventout(maxnbin),heventout(maxnbin),neventout,Y_tmp
      real*8 event_tmp(maxnbin),hevent_tmp(maxnbin),nevent_tmp
      real*8 E,frac,eventout2(maxnbin),iD_tmp,rho_tmp,oab_tmp
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
      L_tmp = L
      V_tmp = V
      Y_tmp = Y
      iD_tmp = iD
      icc_tmp = icc
      nu_mode_tmp = nu_mode
      detect_nu_tmp = detect
      ihfunc_tmp = ihfunc
      ihisto_tmp = ihisto
      rho_tmp = rho
      oab_tmp = oab      

      xmin = 0d0
      xmax = 5d0
      basic_binsize = 0.01
      binsize_factor = 5
      evform = evform_dat

c     Flux	 
c     
c      z_dat(116) = 1d0 ! L
      L = 1d0 ! L
c     z_dat(116) = 295d0 ! L
c      z_dat(116) = 0.28d0
c      z_dat(115) = 1               ! ihfunc 
      icc = 1 ! dummy
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
c     Xsec
c     
c      z_dat(115) = 3
      ihfunc = 3
      ihisto = 1
      ismear = 1

      icc = 1
c      z_dat(150) = 1               ! CC
      detect = 1
c      z_dat(118) = 1               ! detected neutrino
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
      detect = 2
c      z_dat(118) = 2               ! detected neutrino
      call MakeHisto1D(nout,hfunc1D,z_dat,rnevent_ren,nbins_loc
     &     ,xl,evform,serror,snmax,ihisto,eventout
     &     ,heventout,neventout,ierr) 
      open(1,file="xsec_cc_nm.dat",status="replace")
      do i = 0,nbins_loc-1
         write(1,*) xl(i),eventout(i+1)
      enddo
      close(1)
      detect = -1
c      z_dat(118) = -1              ! detected neutrino
      call MakeHisto1D(nout,hfunc1D,z_dat,rnevent_ren,nbins_loc
     &     ,xl,evform,serror,snmax,ihisto,eventout
     &     ,heventout,neventout,ierr) 
      open(1,file="xsec_cc_ae.dat",status="replace")
      do i = 0,nbins_loc-1
         write(1,*) xl(i),eventout(i+1)
      enddo
      close(1)
      detect = -2
c      z_dat(118) = -2              ! detected neutrino
      call MakeHisto1D(nout,hfunc1D,z_dat,rnevent_ren,nbins_loc
     &     ,xl,evform,serror,snmax,ihisto,eventout
     &     ,heventout,neventout,ierr) 
      open(1,file="xsec_cc_am.dat",status="replace")
      do i = 0,nbins_loc-1
         write(1,*) xl(i),eventout(i+1)
      enddo
      close(1)
      
      icc = 2
c      z_dat(150) = 0               ! NC
      detect = 1
c      z_dat(118) = 1               ! detected neutrino
      call MakeHisto1D(nout,hfunc1D,z_dat,rnevent_ren,nbins_loc
     &     ,xl,evform,serror,snmax,ihisto,eventout
     &     ,heventout,neventout,ierr) 
      open(1,file="xsec_nc_ne.dat",status="replace")
      do i = 0,nbins_loc-1
         write(1,*) xl(i),eventout(i+1)
      enddo
      close(1)
      detect = 2
c      z_dat(118) = 2               ! detected neutrino
      call MakeHisto1D(nout,hfunc1D,z_dat,rnevent_ren,nbins_loc
     &     ,xl,evform,serror,snmax,ihisto,eventout
     &     ,heventout,neventout,ierr) 
      open(1,file="xsec_nc_nm.dat",status="replace")
      do i = 0,nbins_loc-1
         write(1,*) xl(i),eventout(i+1)
      enddo
      close(1)
      detect = -1
c      z_dat(118) = -1              ! detected neutrino
      call MakeHisto1D(nout,hfunc1D,z_dat,rnevent_ren,nbins_loc
     &     ,xl,evform,serror,snmax,ihisto,eventout
     &     ,heventout,neventout,ierr) 
      open(1,file="xsec_nc_ae.dat",status="replace")
      do i = 0,nbins_loc-1
         write(1,*) xl(i),eventout(i+1)
      enddo
      close(1)
      detect = -2
c      z_dat(118) = -2              ! detected neutrino
      call MakeHisto1D(nout,hfunc1D,z_dat,rnevent_ren,nbins_loc
     &     ,xl,evform,serror,snmax,ihisto,eventout
     &     ,heventout,neventout,ierr) 
      open(1,file="xsec_nc_am.dat",status="replace")
      do i = 0,nbins_loc-1
         write(1,*) xl(i),eventout(i+1)
      enddo
      close(1)

c$$$c
c$$$c     Flux*P*xsec
c$$$c     
c$$$      L = 1000
c$$$      V = 1d0
c$$$      Y = 5d0
c$$$      MHH = 1
c$$$      basic_binsize = 0.01
c$$$      binsize_factor = 5
c$$$      ihfunc = 4              
c$$$      ihisto = 2
c$$$      icc = 1
c$$$      oab = 0.5
c$$$      rho = 3.0
c$$$      iD = 3
c$$$
c$$$      binsize_loc = basic_binsize*binsize_factor
c$$$      call bining_x(xmin,xmax,binsize_loc,nbins_loc,xl,yyl)
c$$$
c$$$      do i = 0,nbins_loc-1
c$$$         eventout2(i+1) = 0d0
c$$$      enddo      
c$$$c$$$      nu_mode = 1
c$$$c$$$      detect = 1
c$$$c$$$      call MakeHisto1D(nout,hfunc1D,z_dat,rnevent_ren,nbins_loc
c$$$c$$$     &     ,xl,evform,serror,snmax,ihisto,eventout
c$$$c$$$     &     ,heventout,neventout,ierr) 
c$$$c$$$      do i = 0,nbins_loc-1
c$$$c$$$         eventout2(i+1) = eventout2(i+1) +eventout(i+1)
c$$$c$$$      enddo
c$$$      nu_mode = 1
c$$$      detect = 2
c$$$      call MakeHisto1D(nout,hfunc1D,z_dat,rnevent_ren,nbins_loc
c$$$     &     ,xl,evform,serror,snmax,ihisto,eventout
c$$$     &     ,heventout,neventout,ierr) 
c$$$      do i = 0,nbins_loc-1
c$$$         eventout2(i+1) = eventout2(i+1) +eventout(i+1)
c$$$      enddo
c$$$      
c$$$c$$$      nu_mode = -1
c$$$c$$$      detect = -1
c$$$c$$$      call MakeHisto1D(nout,hfunc1D,z_dat,rnevent_ren,nbins_loc
c$$$c$$$     &     ,xl,evform,serror,snmax,ihisto,eventout
c$$$c$$$     &     ,heventout,neventout,ierr) 
c$$$c$$$      do i = 0,nbins_loc-1
c$$$c$$$         eventout2(i+1) = eventout2(i+1) +eventout(i+1)
c$$$c$$$      enddo
c$$$c$$$      nu_mode = -2
c$$$c$$$      detect = -1
c$$$c$$$      call MakeHisto1D(nout,hfunc1D,z_dat,rnevent_ren,nbins_loc
c$$$c$$$     &     ,xl,evform,serror,snmax,ihisto,eventout
c$$$c$$$     &     ,heventout,neventout,ierr) 
c$$$c$$$      do i = 0,nbins_loc-1
c$$$c$$$         eventout2(i+1) = eventout2(i+1) +eventout(i+1)
c$$$c$$$      enddo
c$$$      
c$$$      
c$$$      open(1,file="flux_P_xsec.dat",status="replace")
c$$$      do i = 0,nbins_loc-1
c$$$         write(1,*) xl(i),eventout2(i+1)
c$$$      enddo
c$$$      close(1)      

      basic_binsize = basic_binsize_tmp
      binsize_factor = binsize_factor_tmp
      MHH = MH_tmp
      L = L_tmp
      V = V_tmp
      Y = Y_tmp
      icc = icc_tmp
      nu_mode = nu_mode_tmp
      detect = detect_nu_tmp
      ihfunc = ihfunc_tmp
      ihisto = ihisto_tmp
      oab = oab_tmp
      rho = rho_tmp
      iD = iD_tmp
      
      return
      end
