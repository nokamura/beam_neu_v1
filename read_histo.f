      real*8 function read_histo(tbl,nraw,ncol,binw,icol,E)
C     ****************************************************
C     By Yoshitaro Takaesu @KIAS May.21 2013
C     ****************************************************
      implicitnone
C     ARGUMENTS 
      integer nraw,ncol,icol
      real*8 E,tbl(nraw,ncol),binw
C     LOCAL VARIABLES
      integer i
      real*8 value,xmin,xmax,ymin,ymax
C     ----------
C     BEGIN CODE
C     ----------
      do i = 1,nraw
         xmax = tbl(i,1) +0.025d0
         if (E.lt.xmax+1d-5) then
            value = tbl(i,icol)/binw
            exit
         endif
      enddo
      
      read_histo = value

      return
      end
