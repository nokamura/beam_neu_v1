      subroutine get_event_all(zz,hErec)
C     ****************************************************
C     By Yoshitaro Takaesu @KIAS JAN 23 2014
C     ****************************************************
      implicit none
C     GLOBAL VARIABLES
      include 'inc/params.inc'
      include 'inc/main.inc'
      include 'inc/minfunc.inc'
C     CONSTANTS
C     ARGUMENTS 
      real*8 zz(zdim),hErec(maxnbin,imaxint,nmode,-3:3,-3:3,ndetect
     &     ,-1:1)
C     LOCAL VARIABLES 
C     EXTERNAL FUNCTIONS
C     ----------
C     BEGIN CODE
C     ----------
      if (r_nu.gt.0) then
         beam = 1         ! nu_mu focusing beam
         if (iSK.eq.1) then
            iD = 1
            rho = Srho
            V = SV
            L = SL
            oab = SOAB
            call get_event_allmode(zz,hErec(1,1,1,-3,-3,iD,beam))
         endif
         if (iOki.eq.1) then
            iD = 2
            rho = Orho
            V = OV
            L = OL
            oab = OOAB
            call get_event_allmode(zz,hErec(1,1,1,-3,-3,iD,beam))
         endif
         if (iKr.eq.1) then
            iD = 3
            rho = Krho
            V = KV
            L = KL
            oab = KOAB
            call get_event_allmode(zz,hErec(1,1,1,-3,-3,iD,beam))
         endif
      endif
      if (r_anu.gt.0) then
         beam = -1        ! anti-nu_mu focusing beam
         if (iSK.eq.1) then
            iD = 1
            rho = Srho
            V = SV
            L = SL
            oab = SOAB
            call get_event_allmode(zz,hErec(1,1,1,-3,-3,iD,beam))
         endif
         if (iOki.eq.1) then
            iD = 2
            rho = Orho
            V = OV
            L = OL
            oab = OOAB
            call get_event_allmode(zz,hErec(1,1,1,-3,-3,iD,beam))
         endif
         if (iKr.eq.1) then
            iD = 3
            rho = Krho
            V = KV
            L = KL
            oab =KOAB
            call get_event_allmode(zz,hErec(1,1,1,-3,-3,iD,beam))
         endif
      endif
      
      return
      end
