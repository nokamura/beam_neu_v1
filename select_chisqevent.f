      subroutine select_chisqevent(hErec_in.hErec_out)
C     ****************************************************
C     By Yoshitaro Takaesu @KIAS JAN 7 2014
C     ****************************************************
      implicit none
C     GLOBAL VARIABLES
      include 'inc/params.inc'
      include 'inc/main.inc'
C     CONSTANTS
C     ARGUMENTS 
      real*8 hErec_in(maxnbin,nmode_rec,-3:3,-3:3,ndetect,-1:1)
      real*8 hErec_out(maxnbin,2,ndetect,-1:1)
C     LOCAL VARIABLES 
      integer i
      real*8 emin_SK,emax_SK,mmin_SK,mmax_SK
      real*8 emin_Oki,emax_Oki,mmin_Oki,mmax_Oki
      real*8 emin_Kr,emax_Kr,mmin_Kr,mmax_Kr
C     EXTERNAL FUNCTIONS
C     ----------
C     BEGIN CODE
C     ----------
      emin_SK = 0.4d0
      emax_SK = 1.2d0
      mmin_SK = 0.4d0
      mmax_SK = 5d0
      emin_Oki = 0.4d0
      emax_Oki = 2.8d0
      mmin_Oki = 0.4d0
      mmax_Oki = 5d0
      emin_Kr = 0.4d0
      emax_Kr = 2.8d0
      mmin_Kr = 0.4d0
      mmax_Kr = 5d0

      do i = 0,nbins
         if ((x(i).ge.emin_SK).and.(x(i).le.emax_SK) then
            hErec_out(ipos) = hErec_in(i+1)

      return
      end
