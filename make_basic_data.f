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
      integer MH_tmp,nu_mode_tmp,detect_nu_tmp,ierr,ismear_tmp
      integer nbins_loc
      real*8 Emin_tmp,Emax_tmp,basic_binsize_tmp,L_tmp,V_tmp,icc_tmp
      real*8 eventout(maxnbin),heventout(maxnbin),neventout,Y_tmp
      real*8 event_tmp(maxnbin),hevent_tmp(maxnbin),nevent_tmp
      real*8 E,frac,eventout2(maxnbin),iD_tmp,rho_tmp,oab_tmp
      real*8 eventout_nm(maxnbin),eventout_am(maxnbin)
      real*8 eventout_ne(maxnbin),eventout_ae(maxnbin)
      real*8 xmin,xmax,binsize_loc
      real*8 xl(0:maxnbin),yyl(0:maxnbin)
      character*2 exp
      character*12 file_name      
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
      ismear_tmp = ismear

      xmin = 0d0
      xmax = 5.9d0
      basic_binsize = 0.01
      binsize_factor = 5
      evform = evform_dat


CCC
CCC     Flux	 
CCC     
      icc = 1 ! dummy
      detect = 1 ! dummy
      ihfunc = 1               ! ihfunc 
      ihisto = 2
      binsize_loc = basic_binsize*binsize_factor
      call bining_x(xmin,xmax,binsize_loc,nbins_loc,xl,yyl)

      L = 1 ! km

      oab = 0d0
      exp = "00"
      call write_flux_det(exp,nbins_loc,xl,yyl)
      oab = 0.5d0
      exp = "05"
      call write_flux_det(exp,nbins_loc,xl,yyl)
      oab = 0.6d0
      exp = "06"
      call write_flux_det(exp,nbins_loc,xl,yyl)
      oab = 0.8d0
      exp = "08"
      call write_flux_det(exp,nbins_loc,xl,yyl)
      oab = 0.9d0
      exp = "09"
      call write_flux_det(exp,nbins_loc,xl,yyl)
      oab = 1.0d0
      exp = "10"
      call write_flux_det(exp,nbins_loc,xl,yyl)
      oab = 1.1d0
      exp = "11"
      call write_flux_det(exp,nbins_loc,xl,yyl)
      oab = 1.2d0
      exp = "12"
      call write_flux_det(exp,nbins_loc,xl,yyl)
      oab = 1.3d0
      exp = "13"
      call write_flux_det(exp,nbins_loc,xl,yyl)
      oab = 1.4d0
      exp = "14"
      call write_flux_det(exp,nbins_loc,xl,yyl)
      oab = 1.5d0
      exp = "15"
      call write_flux_det(exp,nbins_loc,xl,yyl)
      oab = 2.0d0
      exp = "20"
      call write_flux_det(exp,nbins_loc,xl,yyl)
      oab = 2.3d0
      exp = "23"
      call write_flux_det(exp,nbins_loc,xl,yyl)
      oab = 2.5d0
      exp = "25"
      call write_flux_det(exp,nbins_loc,xl,yyl)
      oab = 3.0d0
      exp = "30"
      call write_flux_det(exp,nbins_loc,xl,yyl)

      if (iSK.eq.1) then
         exp = "SK"  ! two charactors
         L = SL
         oab = SOAB
         call write_flux_det(exp,nbins_loc,xl,yyl)         
      endif

      if (iOki.eq.1) then
         exp = "Ok"  ! two charactors
         L = OL
         oab = OOAB
         call write_flux_det(exp,nbins_loc,xl,yyl)         
      endif

      if (iKr.eq.1) then
         exp = "Kr"  ! two charactors
         L = KL
         oab = KOAB         
         call write_flux_det(exp,nbins_loc,xl,yyl)
      endif


CCC     
CCC     Transition Probability
CCC     
      ihfunc = 2
      ihisto = 0
      ismear = 0
c      binsize_loc = basic_binsize*binsize_factor
c      binsize_loc = 0.005
      nbins_loc = 1001
c      call bining_x(xmin,xmax,binsize_loc,nbins_loc,xl,yyl)
      call bining_log10x(0.1d0,xmax,nbins_loc,xl,yyl)

      if (iSK.eq.1) then
         L = SL
         rho = Srho

         nu_mode = 1
         detect = 1
         call MakeHisto1D(nout,hfunc1D,z_dat,rnevent_ren,nbins_loc
     &        ,xl,evform,serror,snmax,ihisto,eventout
     &        ,heventout,neventout,ierr) 
         open(1,file="prob_ne.ne_SK.dat",status="replace")
         do i = 0,nbins_loc-1
            write(1,*) xl(i),eventout(i+1)
         enddo
         close(1)
         detect = 2
         call MakeHisto1D(nout,hfunc1D,z_dat,rnevent_ren,nbins_loc
     &        ,xl,evform,serror,snmax,ihisto,eventout
     &        ,heventout,neventout,ierr) 
         open(1,file="prob_ne.nm_SK.dat",status="replace")
         do i = 0,nbins_loc-1
            write(1,*) xl(i),eventout(i+1)
         enddo
         close(1)
         detect = 3
         call MakeHisto1D(nout,hfunc1D,z_dat,rnevent_ren,nbins_loc
     &        ,xl,evform,serror,snmax,ihisto,eventout
     &        ,heventout,neventout,ierr) 
         open(1,file="prob_ne.nt_SK.dat",status="replace")
         do i = 0,nbins_loc-1
            write(1,*) xl(i),eventout(i+1)
         enddo
         close(1)
         
         nu_mode = 2
         detect = 1
         call MakeHisto1D(nout,hfunc1D,z_dat,rnevent_ren,nbins_loc
     &        ,xl,evform,serror,snmax,ihisto,eventout
     &        ,heventout,neventout,ierr) 
         open(1,file="prob_nm.ne_SK.dat",status="replace")
         do i = 0,nbins_loc-1
            write(1,*) xl(i),eventout(i+1)
         enddo
         close(1)
         detect = 2
         call MakeHisto1D(nout,hfunc1D,z_dat,rnevent_ren,nbins_loc
     &        ,xl,evform,serror,snmax,ihisto,eventout
     &        ,heventout,neventout,ierr) 
         open(1,file="prob_nm.nm_SK.dat",status="replace")
         do i = 0,nbins_loc-1
            write(1,*) xl(i),eventout(i+1)
         enddo
         close(1)
         detect = 3
         call MakeHisto1D(nout,hfunc1D,z_dat,rnevent_ren,nbins_loc
     &        ,xl,evform,serror,snmax,ihisto,eventout
     &        ,heventout,neventout,ierr) 
         open(1,file="prob_nm.nt_SK.dat",status="replace")
         do i = 0,nbins_loc-1
            write(1,*) xl(i),eventout(i+1)
         enddo
         close(1)
         
         nu_mode = -1
         detect = -1
         call MakeHisto1D(nout,hfunc1D,z_dat,rnevent_ren,nbins_loc
     &        ,xl,evform,serror,snmax,ihisto,eventout
     &        ,heventout,neventout,ierr) 
         open(1,file="prob_ae.ae_SK.dat",status="replace")
         do i = 0,nbins_loc-1
            write(1,*) xl(i),eventout(i+1)
         enddo
         close(1)
         detect = -2
         call MakeHisto1D(nout,hfunc1D,z_dat,rnevent_ren,nbins_loc
     &        ,xl,evform,serror,snmax,ihisto,eventout
     &        ,heventout,neventout,ierr) 
         open(1,file="prob_ae.am_SK.dat",status="replace")
         do i = 0,nbins_loc-1
            write(1,*) xl(i),eventout(i+1)
         enddo
         close(1)
         detect = -3
         call MakeHisto1D(nout,hfunc1D,z_dat,rnevent_ren,nbins_loc
     &        ,xl,evform,serror,snmax,ihisto,eventout
     &        ,heventout,neventout,ierr) 
         open(1,file="prob_ae.at_SK.dat",status="replace")
         do i = 0,nbins_loc-1
            write(1,*) xl(i),eventout(i+1)
         enddo
         close(1)
         
         nu_mode = -2
         detect = -1
         call MakeHisto1D(nout,hfunc1D,z_dat,rnevent_ren,nbins_loc
     &        ,xl,evform,serror,snmax,ihisto,eventout
     &        ,heventout,neventout,ierr) 
         open(1,file="prob_am.ae_SK.dat",status="replace")
         do i = 0,nbins_loc-1
            write(1,*) xl(i),eventout(i+1)
         enddo
         close(1)
         detect = -2
         call MakeHisto1D(nout,hfunc1D,z_dat,rnevent_ren,nbins_loc
     &        ,xl,evform,serror,snmax,ihisto,eventout
     &        ,heventout,neventout,ierr) 
         open(1,file="prob_am.am_SK.dat",status="replace")
         do i = 0,nbins_loc-1
            write(1,*) xl(i),eventout(i+1)
         enddo
         close(1)
         detect = -3
         call MakeHisto1D(nout,hfunc1D,z_dat,rnevent_ren,nbins_loc
     &        ,xl,evform,serror,snmax,ihisto,eventout
     &        ,heventout,neventout,ierr) 
         open(1,file="prob_am.at_SK.dat",status="replace")
         do i = 0,nbins_loc-1
            write(1,*) xl(i),eventout(i+1)
         enddo
         close(1)
      endif

      if (iKr.eq.1) then
         L = KL
         rho = Krho

         nu_mode = 1
         detect = 1
         call MakeHisto1D(nout,hfunc1D,z_dat,rnevent_ren,nbins_loc
     &        ,xl,evform,serror,snmax,ihisto,eventout
     &        ,heventout,neventout,ierr) 
         open(1,file="prob_ne.ne_Kr.dat",status="replace")
         do i = 0,nbins_loc-1
            write(1,*) xl(i),eventout(i+1)
         enddo
         close(1)
         detect = 2
         call MakeHisto1D(nout,hfunc1D,z_dat,rnevent_ren,nbins_loc
     &        ,xl,evform,serror,snmax,ihisto,eventout
     &        ,heventout,neventout,ierr) 
         open(1,file="prob_ne.nm_Kr.dat",status="replace")
         do i = 0,nbins_loc-1
            write(1,*) xl(i),eventout(i+1)
         enddo
         close(1)
         detect = 3
         call MakeHisto1D(nout,hfunc1D,z_dat,rnevent_ren,nbins_loc
     &        ,xl,evform,serror,snmax,ihisto,eventout
     &        ,heventout,neventout,ierr) 
         open(1,file="prob_ne.nt_Kr.dat",status="replace")
         do i = 0,nbins_loc-1
            write(1,*) xl(i),eventout(i+1)
         enddo
         close(1)
         
         nu_mode = 2
         detect = 1
         call MakeHisto1D(nout,hfunc1D,z_dat,rnevent_ren,nbins_loc
     &        ,xl,evform,serror,snmax,ihisto,eventout
     &        ,heventout,neventout,ierr) 
         open(1,file="prob_nm.ne_Kr.dat",status="replace")
         do i = 0,nbins_loc-1
            write(1,*) xl(i),eventout(i+1)
         enddo
         close(1)
         detect = 2
         call MakeHisto1D(nout,hfunc1D,z_dat,rnevent_ren,nbins_loc
     &        ,xl,evform,serror,snmax,ihisto,eventout
     &        ,heventout,neventout,ierr) 
         open(1,file="prob_nm.nm_Kr.dat",status="replace")
         do i = 0,nbins_loc-1
            write(1,*) xl(i),eventout(i+1)
         enddo
         close(1)
         detect = 3
         call MakeHisto1D(nout,hfunc1D,z_dat,rnevent_ren,nbins_loc
     &        ,xl,evform,serror,snmax,ihisto,eventout
     &        ,heventout,neventout,ierr) 
         open(1,file="prob_nm.nt_Kr.dat",status="replace")
         do i = 0,nbins_loc-1
            write(1,*) xl(i),eventout(i+1)
         enddo
         close(1)
         
         nu_mode = -1
         detect = -1
         call MakeHisto1D(nout,hfunc1D,z_dat,rnevent_ren,nbins_loc
     &        ,xl,evform,serror,snmax,ihisto,eventout
     &        ,heventout,neventout,ierr) 
         open(1,file="prob_ae.ae_Kr.dat",status="replace")
         do i = 0,nbins_loc-1
            write(1,*) xl(i),eventout(i+1)
         enddo
         close(1)
         detect = -2
         call MakeHisto1D(nout,hfunc1D,z_dat,rnevent_ren,nbins_loc
     &        ,xl,evform,serror,snmax,ihisto,eventout
     &        ,heventout,neventout,ierr) 
         open(1,file="prob_ae.am_Kr.dat",status="replace")
         do i = 0,nbins_loc-1
            write(1,*) xl(i),eventout(i+1)
         enddo
         close(1)
         detect = -3
         call MakeHisto1D(nout,hfunc1D,z_dat,rnevent_ren,nbins_loc
     &        ,xl,evform,serror,snmax,ihisto,eventout
     &        ,heventout,neventout,ierr) 
         open(1,file="prob_ae.at_Kr.dat",status="replace")
         do i = 0,nbins_loc-1
            write(1,*) xl(i),eventout(i+1)
         enddo
         close(1)
         
         nu_mode = -2
         detect = -1
         call MakeHisto1D(nout,hfunc1D,z_dat,rnevent_ren,nbins_loc
     &        ,xl,evform,serror,snmax,ihisto,eventout
     &        ,heventout,neventout,ierr) 
         open(1,file="prob_am.ae_Kr.dat",status="replace")
         do i = 0,nbins_loc-1
            write(1,*) xl(i),eventout(i+1)
         enddo
         close(1)
         detect = -2
         call MakeHisto1D(nout,hfunc1D,z_dat,rnevent_ren,nbins_loc
     &        ,xl,evform,serror,snmax,ihisto,eventout
     &        ,heventout,neventout,ierr) 
         open(1,file="prob_am.am_Kr.dat",status="replace")
         do i = 0,nbins_loc-1
            write(1,*) xl(i),eventout(i+1)
         enddo
         close(1)
         detect = -3
         call MakeHisto1D(nout,hfunc1D,z_dat,rnevent_ren,nbins_loc
     &        ,xl,evform,serror,snmax,ihisto,eventout
     &        ,heventout,neventout,ierr) 
         open(1,file="prob_am.at_Kr.dat",status="replace")
         do i = 0,nbins_loc-1
            write(1,*) xl(i),eventout(i+1)
         enddo
         close(1)
      endif
         

CCC     
CCC     Xsec
CCC     
c      z_dat(115) = 3
      ihfunc = 3
      ihisto = 1
      ismear = 1
      binsize_loc = basic_binsize*binsize_factor
      call bining_x(xmin,xmax,binsize_loc,nbins_loc,xl,yyl)

      icc = 1
c      z_dat(150) = 1               ! CC
      detect = 1
c      z_dat(118) = 1               ! detected neutrino
      call MakeHisto1D(nout,hfunc1D,z_dat,rnevent_ren,nbins_loc
     &     ,xl,evform,serror,snmax,ihisto,eventout_ne
     &     ,heventout,neventout,ierr) 
      detect = 2
c      z_dat(118) = 2               ! detected neutrino
      call MakeHisto1D(nout,hfunc1D,z_dat,rnevent_ren,nbins_loc
     &     ,xl,evform,serror,snmax,ihisto,eventout_nm
     &     ,heventout,neventout,ierr) 
      detect = -1
c      z_dat(118) = -1              ! detected neutrino
      call MakeHisto1D(nout,hfunc1D,z_dat,rnevent_ren,nbins_loc
     &     ,xl,evform,serror,snmax,ihisto,eventout_ae
     &     ,heventout,neventout,ierr) 
      detect = -2
c      z_dat(118) = -2              ! detected neutrino
      call MakeHisto1D(nout,hfunc1D,z_dat,rnevent_ren,nbins_loc
     &     ,xl,evform,serror,snmax,ihisto,eventout_am
     &     ,heventout,neventout,ierr) 

      open(1,file="xsec_ccqe.dat",status="replace")
      write(1,*) "# Neutrino CCQE cross section data"
      write(1,*) "# Columns: Enu [GeV], xsecs for nu_e, nu_mu"
     &     ,", bar nu_e, bar nu_mu []"
      write(1,*) " "
      do i = 0,nbins_loc-1
         write(1,*) (xl(i) +xl(i+1))/2d0,eventout_ne(i+1)
     &        ,eventout_nm(i+1),eventout_ae(i+1),eventout_am(i+1)
      enddo
      close(1)
      

      icc = 2
c      z_dat(150) = 0               ! NC
      detect = 1
c      z_dat(118) = 1               ! detected neutrino
      call MakeHisto1D(nout,hfunc1D,z_dat,rnevent_ren,nbins_loc
     &     ,xl,evform,serror,snmax,ihisto,eventout_ne
     &     ,heventout,neventout,ierr) 
      detect = 2
c      z_dat(118) = 2               ! detected neutrino
      call MakeHisto1D(nout,hfunc1D,z_dat,rnevent_ren,nbins_loc
     &     ,xl,evform,serror,snmax,ihisto,eventout_nm
     &     ,heventout,neventout,ierr) 
      detect = -1
c      z_dat(118) = -1              ! detected neutrino
      call MakeHisto1D(nout,hfunc1D,z_dat,rnevent_ren,nbins_loc
     &     ,xl,evform,serror,snmax,ihisto,eventout_ae
     &     ,heventout,neventout,ierr) 
      detect = -2
c      z_dat(118) = -2              ! detected neutrino
      call MakeHisto1D(nout,hfunc1D,z_dat,rnevent_ren,nbins_loc
     &     ,xl,evform,serror,snmax,ihisto,eventout_am
     &     ,heventout,neventout,ierr) 
      open(1,file="xsec_nc.dat",status="replace")
      write(1,*) "# Neutrino NC cross section data"
      write(1,*) "# Columns: Enu [GeV], xsecs for nu_e, nu_mu"
     &     ,", bar nu_e, bar nu_mu []"
      write(1,*) " "
      do i = 0,nbins_loc-1
         write(1,*) (xl(i) +xl(i+1))/2d0,eventout_ne(i+1)
     &        ,eventout_nm(i+1),eventout_ae(i+1),eventout_am(i+1)
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
      ismear = ismear_tmp
      
      return
      end


      subroutine write_flux_det(exp,nbins_loc,xl,yyl)
C     ****************************************************
C     By Yoshitaro Takaesu @ U.Okayama  JAN 05 2016
C     ****************************************************
      implicit none
C     GLOBAL VARIABLES
      include 'inc/params.inc'
C     CONSTANTS
C     ARGUMENTS       
      integer nbins_loc
      real*8 xl(0:maxnbin),yyl(0:maxnbin)
      character*2 exp
C     LOCAL VARIABLES
      character*12 file_name
C     EXTERNAL FUNCTIONS
C     ----------
C     BEGIN CODE
C     ----------
      beam = 1
      file_name = "fluxn_"//exp//".dat"
      call write_flux(file_name,nbins_loc,xl,yyl)
      
      beam = -1
      file_name = "fluxa_"//exp//".dat"
      call write_flux(file_name,nbins_loc,xl,yyl)
      
      return
      end


      subroutine write_flux(file_name,nbins_loc,xl,yyl)
C     ****************************************************
C     By Yoshitaro Takaesu @ U.Okayama JAN 05 2016
C     ****************************************************
      implicit none
C     GLOBAL VARIABLES
      include 'inc/params.inc'
      include 'inc/main.inc'
      include 'inc/minfunc.inc'
C     CONSTANTS
C     ARGUMENTS       
      integer nbins_loc
      real*8 xl(0:maxnbin),yyl(0:maxnbin)
      character*12 file_name
C     LOCAL VARIABLES
      integer i,ierr
      real*8 heventout(maxnbin),neventout
      real*8 eventout_nm(maxnbin),eventout_am(maxnbin)
      real*8 eventout_ne(maxnbin),eventout_ae(maxnbin)
C     EXTERNAL FUNCTIONS
      real*8 hfunc1D
      external hfunc1D
C     ----------
C     BEGIN CODE
C     ----------
      nu_mode = 1               ! nu_e flux
      call MakeHisto1D(nout,hfunc1D,z_dat,rnevent_ren,nbins_loc
     &     ,xl,evform,serror,snmax,ihisto,eventout_ne
     &     ,heventout,neventout,ierr) 
      nu_mode = 2               ! nu_mu flux
      call MakeHisto1D(nout,hfunc1D,z_dat,rnevent_ren,nbins_loc
     &     ,xl,evform,serror,snmax,ihisto,eventout_nm
     &     ,heventout,neventout,ierr) 
      nu_mode = -1              ! anti-nu_e flux
      call MakeHisto1D(nout,hfunc1D,z_dat,rnevent_ren,nbins_loc
     &     ,xl,evform,serror,snmax,ihisto,eventout_ae
     &     ,heventout,neventout,ierr) 
      nu_mode = -2              ! anti-nu_mu flux
      call MakeHisto1D(nout,hfunc1D,z_dat,rnevent_ren,nbins_loc
     &     ,xl,evform,serror,snmax,ihisto,eventout_am
     &     ,heventout,neventout,ierr) 
      
      open(1,file=file_name,status="replace")
      write(1,*) "# neutrino flux data"
      write(1,*) "# L = ", L, "km, OAB = ", oab,"deg."
      write(1,*) "# Columns: Enu [GeV], fluxes for nu_e, nu_mu"
     &     ,", bar nu_e, bar nu_mu [1/cm^2/50MeV/10^{21}POT]"
      write(1,*) " "
      do i = 0,nbins_loc-1
         write(1,*) (xl(i) +xl(i+1))/2d0,eventout_ne(i+1)
     &        ,eventout_nm(i+1),eventout_ae(i+1),eventout_am(i+1)
      enddo
      close(1)
      
      return
      end
