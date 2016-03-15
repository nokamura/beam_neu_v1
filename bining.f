      subroutine bining_basic(xmin,xmax,binsize,nbins,x)
C     ****************************************************
C     By Yoshitaro Takaesu @KIAS JAN 7 2014
C     ****************************************************
      implicit none
C     GLOBAL VARIABLES
C     CONSTANTS
C     ARGUMENTS 
      integer i
      integer nbins,maxnbins
      parameter (maxnbins=100000)
      real*8 xmin,xmax,binsize,x(0:maxnbins)
C     LOCAL VARIABLES 
C     EXTERNAL FUNCTIONS
C     ----------
C     BEGIN CODE
C     ----------
      nbins = int( ( xmax -xmin ) / binsize ) ! nbins should be less than 100000
      do i = 0,nbins
         x(i) = xmin +binsize*i
      enddo
      if (xmax -x(nbins).gt.1d-7) then
         nbins = nbins +1
         x(nbins) = xmax
      endif

      return
      end


      subroutine bining_x(xmin,xmax,binsize,nbins,x,yy)
C     ****************************************************
C     By Yoshitaro Takaesu @KIAS JAN 12 2014
C
C     Calculate borders of each bins from the lower and 
C     upper limits and binsize.
C     ****************************************************
      implicit none
C     GLOBAL VARIABLES
C     CONSTANTS
C     ARGUMENTS 
      integer i
      integer nbins,maxnbins
      parameter (maxnbins=100000)
      real*8 xmin,xmax,binsize,x(0:maxnbins),yy(0:maxnbins)
C     LOCAL VARIABLES 
C     EXTERNAL FUNCTIONS
C     ----------
C     BEGIN CODE
C     ----------
      nbins = int( ( xmax -xmin ) / binsize ) ! nbins should be less than 100000
      if (nbins.eq.0) then
         nbins = 1
         x(0) = xmin
         x(1) = xmax
         yy(0) = x(0)
         yy(1) = x(1)
      else
         do i = 0,nbins
            x(i) = xmin +binsize*i
            yy(i) = x(i)
         enddo
         if (xmax-x(nbins).gt.1d-7) then
            nbins = nbins+1
            x(nbins) = xmax
            yy(nbins) = x(nbins)
         endif
      endif
      
      return
      end
