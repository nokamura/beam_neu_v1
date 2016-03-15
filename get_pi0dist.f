      subroutine get_pi0dist(event,x,nbins,event_pi0)
C     ****************************************************
C     By Yoshitaro Takaesu @KIAS AUG 15 2013
C
C     event should be prepared with 50 MeV bin. 
C     (This restriction may be removed in future)
C     ****************************************************
      implicit none
C     GLOBAL VARIABLES
      include 'pi0dist/pi0dist0025.inc'
      include 'pi0dist/pi0dist0075.inc'
      include 'pi0dist/pi0dist0125.inc'
      include 'pi0dist/pi0dist0175.inc'
C     CONSTANTS
C     ARGUMENTS 
      integer nbins
      real*8 event,event_pi0,x
C     LOCAL VARIABLES 
      integer pi0_event
      real*8 pi0_binwidth
      parameter (pi0_event=1000,pi0_binwidth=0.05d0)
C     EXTERNAL FUNCTIONS
C     ----------
C     BEGIN CODE
C     ----------
      do i = 1,nbins
         event_pi0(i) = 0d0
      enddo
      do i = 1,nbins
         bin_nevent = event(i)
         scale = bin_nevent/pi0_event
         event_pi0(i) = event_pi0(i) +pi0dist0025
         
         

      return
      end
