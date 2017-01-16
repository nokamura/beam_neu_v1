      subroutine make_basic_data
C     ****************************************************
C     By Yoshitaro Takaesu @KIAS JAN 12 2014
C     ****************************************************
      implicit none
C     GLOBAL VARIABLES
      include 'inc/params.inc'
      include 'inc/main.inc'
      include 'inc/minfunc.inc'
C     LOCAL VARIABLES
      integer i,iout
      integer nbins_loc,ierr
      real*8 eventout(maxnbin),heventout(maxnbin),neventout
      real*8 E,frac,eventout2(maxnbin)
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
CCC   save global parameters
      call save_params_to_tmp

      evform = evform_dat
      call output_flux
      call output_prob
      call output_xsec
c      call output_flux_xsec_prob

      xmin = 0d0
      xmax = 5.9d0
      basic_binsize = 0.01
      binsize_factor = 5

c
c     Flux*P*xsec
c     
      L = 295
      V = 22.5d0
      Y = 5d0
c      MHH = 1
c      basic_binsize = 0.01
c      binsize_factor = 5
      ismear = 0
      ihfunc = 4
      ihisto = 2
      icc = 1
c      oab = 0.5
c      rho = 3.0
c      iD = 3

      binsize_loc = basic_binsize*binsize_factor
      call bining_x(xmin,xmax,binsize_loc,nbins_loc,xl,yyl)

CCC   nu_e events
      detect = 1
      do i = 0,nbins_loc-1
         eventout2(i+1) = 0d0
      enddo      

CCCCCC from nu_e flux
      nu_mode = 1
      call wrap_MakeHisto1D(nbins_loc,xl,eventout)
      do i = 0,nbins_loc-1
         eventout2(i+1) = eventout2(i+1) +eventout(i+1)         
      enddo
C$$$CCCCCC from nu_mu flux
C$$$      nu_mode = 2
C$$$      call wrap_MakeHisto1D(nbins_loc,xl,eventout)
C$$$      do i = 0,nbins_loc-1
C$$$         eventout2(i+1) = eventout2(i+1) +eventout(i+1)
C$$$      enddo
C$$$CCCCCC from anti-nu_e flux
C$$$      nu_mode = -1
C$$$      call wrap_MakeHisto1D(nbins_loc,xl,eventout)
C$$$      do i = 0,nbins_loc-1
C$$$         eventout2(i+1) = eventout2(i+1) +eventout(i+1)
C$$$      enddo
C$$$CCCCCC from anti-nu_mu flux
C$$$      nu_mode = -2
C$$$      call wrap_MakeHisto1D(nbins_loc,xl,eventout)
C$$$      do i = 0,nbins_loc-1
C$$$         eventout2(i+1) = eventout2(i+1) +eventout(i+1)
C$$$      enddo

CCCCCC output event numbers
      open(1,file="flux_P_xsec_ne.dat",status="replace")
      do i = 0,nbins_loc-1
c         write(1,*) xl(i),eventout2(i+1)
         write(1,*) xl(i),eventout(i+1)
      enddo
      close(1)      

CCC   reset global parameters to their initial values
      call reset_tmp_to_params
      
      return
      end


      subroutine save_params_to_tmp
C     ****************************************************
C     By Yoshitaro Takaesu @ U.Okayama JAN 16 2017
C     ****************************************************
      implicit none
C     GLOBAL VARIABLES
      include 'inc/params.inc'

      basic_binsize_tmp = basic_binsize
      binsize_factor_tmp = binsize_factor
      MH_tmp = MHH
      L_tmp = L
      V_tmp = V
      Y_tmp = Y
      icc_tmp = icc
      nu_mode_tmp = nu_mode
      detect_nu_tmp = detect
      ihfunc_tmp = ihfunc
      ihisto_tmp = ihisto
      rho_tmp = rho
      oab_tmp = oab      
      ismear_tmp = ismear
    
      return
      end


      subroutine reset_tmp_to_params
C     ****************************************************
C     By Yoshitaro Takaesu @ U.Okayama JAN 16 2017
C     ****************************************************
      implicit none
C     GLOBAL VARIABLES
      include 'inc/params.inc'

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
      rho = rho_tmp
      oab = oab_tmp
      ismear = ismear_tmp
    
      return
      end


      subroutine output_flux
C     ****************************************************
C     By Yoshitaro Takaesu @ U.Okayama JAN 16 2017
C     ****************************************************
      implicit none
C     GLOBAL VARIABLES
      include 'inc/params.inc'
C     LOCAL VARIABLES
      integer nbins_loc
      character*2 exp
      double precision xmin,xmax,binsize_loc
      double precision xl(0:maxnbin),yyl(0:maxnbin)

      icc = 1 ! dummy
      detect = 1 ! dummy
      ihfunc = 1               ! ihfunc 
      ihisto = 2

      xmin = 0d0
      xmax = 5.9d0
      basic_binsize = 0.01
      binsize_factor = 5
      binsize_loc = basic_binsize*binsize_factor
      call bining_x(xmin,xmax,binsize_loc,nbins_loc,xl,yyl)

      L = 1 ! km

      oab = 0d0
      exp = "00"
      call write_flux_det(exp,nbins_loc,xl)
      oab = 0.5d0
      exp = "05"
      call write_flux_det(exp,nbins_loc,xl)
      oab = 0.6d0
      exp = "06"
      call write_flux_det(exp,nbins_loc,xl)
      oab = 0.8d0
      exp = "08"
      call write_flux_det(exp,nbins_loc,xl)
      oab = 0.9d0
      exp = "09"
      call write_flux_det(exp,nbins_loc,xl)
      oab = 1.0d0
      exp = "10"
      call write_flux_det(exp,nbins_loc,xl)
      oab = 1.1d0
      exp = "11"
      call write_flux_det(exp,nbins_loc,xl)
      oab = 1.2d0
      exp = "12"
      call write_flux_det(exp,nbins_loc,xl)
      oab = 1.3d0
      exp = "13"
      call write_flux_det(exp,nbins_loc,xl)
      oab = 1.4d0
      exp = "14"
      call write_flux_det(exp,nbins_loc,xl)
      oab = 1.5d0
      exp = "15"
      call write_flux_det(exp,nbins_loc,xl)
      oab = 2.0d0
      exp = "20"
      call write_flux_det(exp,nbins_loc,xl)
      oab = 2.3d0
      exp = "23"
      call write_flux_det(exp,nbins_loc,xl)
      oab = 2.5d0
      exp = "25"
      call write_flux_det(exp,nbins_loc,xl)
      oab = 3.0d0
      exp = "30"
      call write_flux_det(exp,nbins_loc,xl)

      if (iSK.eq.1) then
         exp = "SK"  ! two charactors
         L = SL
         oab = SOAB
         call write_flux_det(exp,nbins_loc,xl)         
      endif

      if (iOki.eq.1) then
         exp = "Ok"  ! two charactors
         L = OL
         oab = OOAB
         call write_flux_det(exp,nbins_loc,xl)         
      endif

      if (iKr.eq.1) then
         exp = "Kr"  ! two charactors
         L = KL
         oab = KOAB         
         call write_flux_det(exp,nbins_loc,xl)
      endif

CCC   reset global parameters to their initial values
      call reset_tmp_to_params

      return
      end


      subroutine write_flux(iout,nbins_loc,xl)
C     ********************************************************
C     By Yoshitaro Takaesu @ U.Okayama JAN 05 2016
C     Modified: JAN 06 2016: Remove open close statements
C     ********************************************************
      implicit none
C     GLOBAL VARIABLES
      include 'inc/params.inc'
C     CONSTANTS
C     ARGUMENTS       
      integer nbins_loc
      real*8 xl(0:maxnbin)
C     LOCAL VARIABLES
      integer i,iout
      real*8 eventout_nm(maxnbin),eventout_am(maxnbin)
      real*8 eventout_ne(maxnbin),eventout_ae(maxnbin)
C     EXTERNAL FUNCTIONS
C     ----------
C     BEGIN CODE
C     ----------
      nu_mode = 1               ! nu_e flux
      call wrap_MakeHisto1D(nbins_loc,xl,eventout_ne)
      nu_mode = 2               ! nu_mu flux
      call wrap_MakeHisto1D(nbins_loc,xl,eventout_nm)
      nu_mode = -1              ! anti-nu_e flux
      call wrap_MakeHisto1D(nbins_loc,xl,eventout_ae)
      nu_mode = -2              ! anti-nu_mu flux
      call wrap_MakeHisto1D(nbins_loc,xl,eventout_am)
      
      write(iout,*) "# neutrino flux data"
      write(iout,*) "# L = ", L, "km, OAB = ", oab,"deg."
      write(iout,*) "# Columns: Enu [GeV], fluxes for nu_e, nu_mu"
     &     ,", bar nu_e, bar nu_mu [1/cm^2/50MeV/10^{21}POT]"
      write(iout,*) " "
      do i = 0,nbins_loc-1
         write(iout,*) (xl(i) +xl(i+1))/2d0,eventout_ne(i+1)
     &        ,eventout_nm(i+1),eventout_ae(i+1),eventout_am(i+1)
      enddo
      
      return
      end


      subroutine output_prob
C     ****************************************************
C     By Yoshitaro Takaesu @ U.Okayama JAN 16 2017
C     ****************************************************
      implicit none
C     GLOBAL VARIABLES
      include 'inc/params.inc'
C     LOCAL VARIABLES
      integer nbins_loc
      character*2 exp
      double precision xmin,xmax,binsize_loc
      double precision xl(0:maxnbin),yyl(0:maxnbin)

      ihfunc = 2
      ihisto = 0
      ismear = 0

      xmin = 0.1d0
      xmax = 5.9d0
c      binsize_loc = basic_binsize*binsize_factor
c      call bining_x(xmin,xmax,binsize_loc,nbins_loc,xl,yyl)
      nbins_loc = 1001
      call bining_log10x(xmin,xmax,nbins_loc,xl,yyl)

      if (iSK.eq.1) then
         exp = "SK"
         L = SL
         rho = Srho
         call write_prob_det(exp,nbins_loc,xl)        
      endif

      if (iKr.eq.1) then
         exp = "Kr"
         L = KL
         rho = Krho
         call write_prob_det(exp,nbins_loc,xl)        
      endif

CCC   reset global parameters to their initial values
      call reset_tmp_to_params

      return
      end


      subroutine write_flux_det(exp,nbins_loc,xl)
C     ****************************************************
C     By Yoshitaro Takaesu @ U.Okayama  JAN 05 2016
C     ****************************************************
      implicit none
C     GLOBAL VARIABLES
      include 'inc/params.inc'
C     CONSTANTS
C     ARGUMENTS       
      integer nbins_loc
      real*8 xl(0:maxnbin)
      character*2 exp
C     LOCAL VARIABLES
C     EXTERNAL FUNCTIONS
C     ----------
C     BEGIN CODE
C     ----------
      beam = 1
      open(1,file="fluxn_"//exp//".dat",status="replace")
      call write_flux(1,nbins_loc,xl)
      close(1)
      
      beam = -1
      open(1,file="fluxa_"//exp//".dat",status="replace")
      call write_flux(1,nbins_loc,xl)
      close(1)

      return
      end


      subroutine write_prob_det(exp,nbins_loc,xl)
C     ****************************************************
C     By Yoshitaro Takaesu @ U.Okayama  JAN 05 2016
C     ****************************************************
      implicit none
C     GLOBAL VARIABLES
      include 'inc/params.inc'
C     ARGUMENTS       
      integer nbins_loc
      character*2 exp
      double precision xl(0:maxnbin)
C     ----------
C     BEGIN CODE
C     ----------
      nu_mode = 1
      detect = 1
      call write_prob(exp,nbins_loc,xl)
      detect = 2
      call write_prob(exp,nbins_loc,xl)
      detect = 3
      call write_prob(exp,nbins_loc,xl)
      
      nu_mode = 2
      detect = 1
      call write_prob(exp,nbins_loc,xl)
      detect = 2
      call write_prob(exp,nbins_loc,xl)
      detect = 3
      call write_prob(exp,nbins_loc,xl)
      
      nu_mode = -1
      detect = -1
      call write_prob(exp,nbins_loc,xl)
      detect = -2
      call write_prob(exp,nbins_loc,xl)
      detect = -3
      call write_prob(exp,nbins_loc,xl)
      
      nu_mode = -2
      detect = -1
      call write_prob(exp,nbins_loc,xl)
      detect = -2
      call write_prob(exp,nbins_loc,xl)
      detect = -3
      call write_prob(exp,nbins_loc,xl)
      
      return
      end


      subroutine write_prob(exp,nbins_loc,xl)
C     ****************************************************
C     By Yoshitaro Takaesu @ U.Okayama  JAN 05 2016
C     ****************************************************
      implicit none
C     GLOBAL VARIABLES
      include 'inc/params.inc'
C     ARGUMENTS       
      integer nbins_loc
      character*2 exp
      double precision xl(0:maxnbin)
C     LOCAL VARIABLES
      integer i
      character*2 cnu_mode,cnu_detect
      double precision eventout(maxnbin)
C     ----------
C     BEGIN CODE
C     ----------
      if (nu_mode.eq.1) cnu_mode = "ne"
      if (nu_mode.eq.2) cnu_mode = "nm"
      if (nu_mode.eq.-1) cnu_mode = "ae"
      if (nu_mode.eq.-2) cnu_mode = "am"

      if (detect.eq.1) cnu_detect = "ne"
      if (detect.eq.2) cnu_detect = "nm"
      if (detect.eq.3) cnu_detect = "nt"
      if (detect.eq.-1) cnu_detect = "ae"
      if (detect.eq.-2) cnu_detect = "am"
      if (detect.eq.-3) cnu_detect = "at"

      call wrap_MakeHisto1D(nbins_loc,xl,eventout)
      open(1,file="prob_"//cnu_mode//"."//cnu_detect//"_"//exp//".dat"
     &     ,status="replace")
      do i = 0,nbins_loc-1
         write(1,*) xl(i),eventout(i+1)
      enddo
      close(1)

      return
      end


      subroutine output_xsec
C     ****************************************************
C     By Yoshitaro Takaesu @ U.Okayama JAN 16 2017
C     ****************************************************
      implicit none
C     GLOBAL VARIABLES
      include 'inc/params.inc'
C     LOCAL VARIABLES
      integer nbins_loc,iout
      character*2 exp
      double precision xmin,xmax,binsize_loc
      double precision xl(0:maxnbin),yyl(0:maxnbin)
C     ----------
C     BEGIN CODE
C     ----------
      ihfunc = 3
      ihisto = 1
      ismear = 1

      xmin = 0d0
      xmax = 5.9d0
      basic_binsize = 0.01
      binsize_factor = 5
      binsize_loc = basic_binsize*binsize_factor
      call bining_x(xmin,xmax,binsize_loc,nbins_loc,xl,yyl)

      iout = 1
CCC   CC cross section
      icc = 1
      xsec_mode = 0
      open(iout,file="xsec_ccqetot.dat",status="replace")
      write(iout,*) "# Neutrino total CCQE cross section of a water"
     &     ," target"
      call write_xsec(iout,nbins_loc,xl)
      close(iout)

      xsec_mode = 1
      open(iout,file="xsec_ccqe.dat",status="replace")
      write(iout,*) "# Neutrino CCQE cross section (after CCQE cut)"
     &     ," of a water target"
      call write_xsec(iout,nbins_loc,xl)
      close(iout)

      xsec_mode = 2
      open(iout,file="xsec_ccres.dat",status="replace")
      write(iout,*) "# Neutrino CC resonant cross section"
     &     ," (after CCQE cut) of a water target"
      call write_xsec(iout,nbins_loc,xl)
      close(iout)

CCC   NC cross section
      icc = 2
      iout = 2
      xsec_mode = 0
      open(iout,file="xsec_nctot.dat",status="replace")
      write(iout,*) "# Neutrino total NC cross section of a water"
     &     ," target"
      call write_xsec(iout,nbins_loc,xl)
      close(iout)

      xsec_mode = 1
      open(iout,file="xsec_ncqe.dat",status="replace")
      write(iout,*) "# Neutrino NCQE cross section of a water"
     &     ," target"
      call write_xsec(iout,nbins_loc,xl)
      close(iout)

      xsec_mode = 2
      open(iout,file="xsec_ncres.dat",status="replace")
      write(iout,*) "# Neutrino NC resonant cross section of a water"
     &     ," target"
      call write_xsec(iout,nbins_loc,xl)
      close(iout)

      xsec_mode = 3
      open(iout,file="xsec_ncdi.dat",status="replace")
      write(iout,*) "# Neutrino NCCI cross section of a water"
     &     ," target"
      call write_xsec(iout,nbins_loc,xl)
      close(iout)

      xsec_mode = 4
      open(iout,file="xsec_ncco.dat",status="replace")
      write(iout,*) "# Neutrino NCCoh + NCDiff cross section of a water"
     &     ," target"
      call write_xsec(iout,nbins_loc,xl)
      close(iout)

      return
      end


      subroutine write_xsec(iout,nbins_loc,xl)
C     ****************************************************
C     By Yoshitaro Takaesu @ U.Okayama JAN 06 2016
C     ****************************************************
      implicit none
C     GLOBAL VARIABLES
      include 'inc/params.inc'
C     CONSTANTS
C     ARGUMENTS       
      integer nbins_loc
      real*8 xl(0:maxnbin)
C     LOCAL VARIABLES
      integer i,iout
      real*8 eventout_nm(maxnbin),eventout_am(maxnbin)
      real*8 eventout_ne(maxnbin),eventout_ae(maxnbin)
C     EXTERNAL FUNCTIONS
C     ----------
C     BEGIN CODE
C     ----------
      detect = 1                ! nu_e
      call wrap_MakeHisto1D(nbins_loc,xl,eventout_ne)
      detect = 2                ! nu_mu
      call wrap_MakeHisto1D(nbins_loc,xl,eventout_nm)
      detect = -1               ! anti-nu_e
      call wrap_MakeHisto1D(nbins_loc,xl,eventout_ae)
      detect = -2               ! anti-nu_mu
      call wrap_MakeHisto1D(nbins_loc,xl,eventout_am)
      write(iout,*) "# Columns: Enu [GeV], xsecs for nu_e, nu_mu"
     &     ,", bar nu_e, bar nu_mu [cm^2/kton]"
      write(iout,*) " "
      do i = 0,nbins_loc-1
         write(iout,*) (xl(i) +xl(i+1))/2d0,eventout_ne(i+1)
     &        ,eventout_nm(i+1),eventout_ae(i+1),eventout_am(i+1)
      enddo
      
      return
      end


      subroutine wrap_MakeHisto1D(nbins_loc,xl,eventout)
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
      real*8 xl(0:maxnbin)
      real*8 eventout(maxnbin)
C     LOCAL VARIABLES
      integer ierr
      real*8 heventout(maxnbin),neventout
C     EXTERNAL FUNCTIONS
      real*8 hfunc1D
      external hfunc1D
C     ----------
C     BEGIN CODE
C     ----------
      call MakeHisto1D(nout,hfunc1D,z_dat,rnevent_ren,nbins_loc
     &     ,xl,evform,serror,snmax,ihisto,eventout
     &     ,heventout,neventout,ierr) 

      return
      end
