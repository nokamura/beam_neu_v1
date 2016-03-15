      subroutine get_recevent(z,nbins,hErec_raw,hErec_rec)
C     ****************************************************
C     By Yoshitaro Takaesu @KIAS JAN 7 2014
C     ****************************************************
      implicit none
C     GLOBAL VARIABLES
      include 'inc/params.inc'
C     CONSTANTS
C     ARGUMENTS 
      integer nbins
      real*8 z(zdim)
      real*8 hErec_raw(maxnbin,imaxint,nmode,-3:3,-3:3,ndetect,-1:1)
      real*8 hErec_rec(maxnbin,imaxint,nmode_rec,-3:3,-3:3,ndetect,-1:1)
C     LOCAL VARIABLES 
      integer i,j,k,l1,l2,l3
C     EXTERNAL FUNCTIONS
C     ----------
C     BEGIN CODE
C     ----------
      do l2 = -1,1,2
         do l1 = 1,ndetect
            do k = -3,3
               do j = -3,3
                  do l3 = 1,imaxint
                     do i = 1,nbins
C     e > e events
                       hErec_rec(i,l3,1,j,k,l1,l2) = 
     &                       (1d0 -Pe2m)*effe*hErec_raw(i,l3,1,j,k,l1
     &                       ,l2)
C     mu > mu  events
                        hErec_rec(i,l3,2,j,k,l1,l2) = 
     &                       (1d0 -Pm2e)*effmu*hErec_raw(i,l3,2,j,k,l1
     &                       ,l2)
C     Just for bokkkeeping of NC-pi0 and NC-gamma events
                        hErec_rec(i,l3,3,j,k,l1,l2) = hErec_raw(i,l3,3,j
     &                       ,k,l1,l2)
                        hErec_rec(i,l3,4,j,k,l1,l2) = hErec_raw(i,l3,4,j
     &                       ,k,l1,l2)
C     mu > e miss-ID events
                        hErec_rec(i,l3,5,j,k,l1,l2) = 
     &                       Pm2e*effmu*hErec_raw(i,l3,2,j,k,l1,l2)
C     e > mu miss-ID events
                        hErec_rec(i,l3,6,j,k,l1,l2) = 
     &                       Pe2m*effe*hErec_raw(i,l3,1,j,k,l1,l2)
                     enddo
                  enddo
               enddo
            enddo
         enddo
      enddo

      return
      end
