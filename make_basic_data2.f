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
      real*8 eventout(maxnbin),heventout(maxnbin),neventout
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

c
c     Flux*P*xsec
c     
      L = 1000
      V = 1d0
      MHH = 1
      basic_binsize = 0.01
      binsize_factor = 5
      ihfunc = 4              
      ihisto = 2
      icc = 1
      oab = 0.5
      rho = 3.0
      iD = 3

      binsize_loc = basic_binsize*binsize_factor
      call bining_x(xmin,xmax,binsize_loc,nbins_loc,xl,yyl)

      do i = 0,nbins_loc-1
         eventout2(i+1) = 0d0
      enddo      
c$$$      nu_mode = 1
c$$$      detect = 1
c$$$      call MakeHisto1D(nout,hfunc1D,z_dat,rnevent_ren,nbins_loc
c$$$     &     ,xl,evform,serror,snmax,ihisto,eventout
c$$$     &     ,heventout,neventout,ierr) 
c$$$      do i = 0,nbins_loc-1
c$$$         eventout2(i+1) = eventout2(i+1) +eventout(i+1)
c$$$      enddo
      nu_mode = 1
      detect = 2
      call MakeHisto1D(nout,hfunc1D,z_dat,rnevent_ren,nbins_loc
     &     ,xl,evform,serror,snmax,ihisto,eventout
     &     ,heventout,neventout,ierr) 
      do i = 0,nbins_loc-1
         eventout2(i+1) = eventout2(i+1) +eventout(i+1)
      enddo
      
c$$$      nu_mode = -1
c$$$      detect = -1
c$$$      call MakeHisto1D(nout,hfunc1D,z_dat,rnevent_ren,nbins_loc
c$$$     &     ,xl,evform,serror,snmax,ihisto,eventout
c$$$     &     ,heventout,neventout,ierr) 
c$$$      do i = 0,nbins_loc-1
c$$$         eventout2(i+1) = eventout2(i+1) +eventout(i+1)
c$$$      enddo
c$$$      nu_mode = -2
c$$$      detect = -1
c$$$      call MakeHisto1D(nout,hfunc1D,z_dat,rnevent_ren,nbins_loc
c$$$     &     ,xl,evform,serror,snmax,ihisto,eventout
c$$$     &     ,heventout,neventout,ierr) 
c$$$      do i = 0,nbins_loc-1
c$$$         eventout2(i+1) = eventout2(i+1) +eventout(i+1)
c$$$      enddo
      
      
      open(1,file="flux_P_xsec.dat",status="replace")
      do i = 0,nbins_loc-1
         write(1,*) xl(i),eventout2(i+1)
      enddo
      close(1)      

      basic_binsize = basic_binsize_tmp
      binsize_factor = binsize_factor_tmp
      MHH = MH_tmp
      L = L_tmp
      V = V_tmp
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
