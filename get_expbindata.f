      subroutine get_expbindata(event_in,event_out,nevent_out)
C     ****************************************************
C     By Yoshitaro Takaesu @KIAS JAN 7 2014
C     ****************************************************
      implicit none
C     GLOBAL VARIABLES
      include 'inc/params.inc'
      include 'inc/main.inc'
C     CONSTANTS
C     ARGUMENTS 
      real*8 event_in(maxnbin),event_out(maxnbin)
C     LOCAL VARIABLES 
      integer i
      integer ite,ibins
C     EXTERNAL FUNCTIONS
C     ----------
C     BEGIN CODE
C     ----------
      ite = 0
      ibins = 1
      nevent_out = 0d0
      do i = 1,maxnbin
         event_out(i) = 0d0
      enddo
      do i = 1,nbins_basic
         ite = ite +1
         eventout(ibins) = eventout(ibins) +event_in(i) 
         if (ite.eq.binsize_factor) then
            ibins = ibins +1
            ite = 0
         endif
      enddo

      return
      end
