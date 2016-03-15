      subroutine ret_bins(itype,evnbins,evx,nbins,x)
C     ****************************************************
C     By Yoshitaro Takaesu @KIAS MAY 27  2013
C     ****************************************************
      implicitnone
C     CONSTANTS
      integer maxnbin
      parameter (maxnbin=20000)
C     GLOBAL VARIABLES
C     ARGUMENTS 
      integer itype,evnbins(2),nbins
      real*8 x(0:maxnbin),evx(0:maxnbin,2)
C     LOCAL VARIABLES 
      integer i
C     EXTERNAL FUNCTIONS
C     ----------
C     BEGIN CODE
C     ----------
      nbins = evnbins(itype)
      do i = 0,nbins
         x(i) = evx(i,itype)
      enddo

      return
      end
