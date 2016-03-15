      subroutine classify_recevent(nbins,dat,like)
C     ****************************************************
C     By Yoshitaro Takaesu @KIAS JAN 31 2014
C     ****************************************************
      implicit none
C     GLOBAL VARIABLES
      include 'inc/params.inc'
C     CONSTANTS
C     ARGUMENTS 
      integer nbins
      real*8 dat(maxnbin,imaxint,nmode_rec,-3:3,-3:3,ndetect,-1:1)
      real*8 like(maxnbin,2,ndetect,-1:1)
C     LOCAL VARIABLES 
      integer i,j,k,l1,l2,l3,l4
      integer idebug
C     EXTERNAL FUNCTIONS
C     ----------
C     BEGIN CODE
C     ----------
c     check input data
      idebug = 1
      if (idebug.eq.1) then
         do l2 = -1,1,2
            do l1 = 1,ndetect
               do k = -3,3
                  do j = -3,3
                     do l4 = 1,5
                        do l3 = 1,imaxint
                           do i = 1,nbins      
                              if (dat(i,l3,l4,j,k,l1,l2)*0.ne.0) then
                                 write(*,*) "ERROR:classify_recevent.f"
     &                                ,":Input data is invalid"
                                 write(*,*) i,l3,l4,j,k,l1,l2
                                 stop
                              endif
                           enddo
                        enddo
                     enddo
                  enddo
               enddo
            enddo
         enddo
      endif

      do l1 = -1,1,2
         do k = 1,ndetect
            do j = 1,nmode_like
               do i = 1,nbins
                  like(i,j,k,l1) = 0d0
               enddo
            enddo
         enddo
      enddo   

      do l2 = -1,1,2
         do l1 = 1,ndetect
            do k = -3,3
               do j = -3,3
                  do l3 = 1,imaxint
                     do i = 1,nbins
C     e-like events
                        like(i,1,l1,l2) = like(i,1,l1,l2)
     &                       +dat(i,l3,1,j,k,l1,l2)
     &                       +dat(i,l3,3,j,k,l1,l2)
     &                       +dat(i,l3,4,j,k,l1,l2)
     &                       +dat(i,l3,5,j,k,l1,l2)
C     mu-like  events
                        like(i,2,l1,l2) = like(i,2,l1,l2)
     &                       +dat(i,l3,2,j,k,l1,l2)
     &                       +dat(i,l3,6,j,k,l1,l2)
                     enddo
                  enddo
               enddo
            enddo
         enddo
      enddo

      return
      end
