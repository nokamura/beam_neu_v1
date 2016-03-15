      subroutine interpolate(x1,x2,x,y1,y2,y)
C     ****************************************************
C     By Yoshitaro Takaesu @KIAS May 21 2013
C     ****************************************************
      implicit none
C     GLOBAL VARIABLES
C     CONSTANTS
C     ARGUMENTS 
      real*8 x1,x2,x,y1,y2,y
C     LOCAL VARIABLES 
      real*8 a
C     EXTERNAL FUNCTIONS
C     ----------
C     BEGIN CODE
C     ----------
      a = (y2-y1)/(x2-x1)
      y = y1 +a*(x-x1)

      return
      end
