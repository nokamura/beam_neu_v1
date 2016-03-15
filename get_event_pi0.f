      subroutine get_event_pi0(z,hErec)
C     ****************************************************
C     By Yoshitaro Takaesu @KIAS JAN 7 2014
C     ****************************************************
      implicit none
C     GLOBAL VARIABLES
      include 'inc/params.inc'
C     CONSTANTS
C     ARGUMENTS 
      real*8 z(zdim),hErec(maxnbin,imaxint,nmode)
C     LOCAL VARIABLES 
      integer iproc
c      real*8 nevent_tmp
C     EXTERNAL FUNCTIONS
C     ----------
C     BEGIN CODE
C     ----------
      if (iNC.eq.1) then
         icc = 2                ! NC interactions
         iproc = 3              ! NC events for ismear = 0, NCpi0_bg for ismear = 1 
         call get_event(z,iproc,hErec(1,1,iproc))
      endif

      return
      end
