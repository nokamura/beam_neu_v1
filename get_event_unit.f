      subroutine get_event_unit(z,hErec)
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
      if ((iCCQE.eq.1).or.(iCCRes.eq.1)) then
         icc = 1                ! CC interactions
         if ((abs(detect).eq.1).or.(abs(detect).eq.2)) then
            iproc = abs(detect)
c            call get_event(z,iproc,hErec(1,1,iproc),nevent_tmp)
            call get_event(z,iproc,hErec(1,1,iproc))
         endif
      endif
      if (iNC.eq.1) then
         icc = 2                ! NC interactions
         iproc = 3              ! NC events for ismear = 0, NCpi0_bg for ismear = 1 
c         call get_event(z,iproc,hErec(1,1,iproc),nevent_tmp)
         call get_event(z,iproc,hErec(1,1,iproc))
      endif

      return
      end
