      subroutine prepare_events
C     ****************************************************
C     By Yoshitaro Takaesu @KIAS JAN 19 2014
C     ****************************************************
      implicit none
C     GLOBAL VARIABLES
      include 'inc/params.inc'
      include 'inc/main.inc'
C     CONSTANTS
C     ARGUMENTS 
C     LOCAL VARIABLES 
      integer i,j,k,l1,l2,l3,l4
C     EXTERNAL FUNCTIONS
C     ----------
C     BEGIN CODE
C     ----------
      do k = -1,1,2
         do j = 1,ndetect
            nevent_e(j,k) = 0d0
            nevent_mu(j,k) = 0d0
            nevent_pi0(j,k) = 0d0
            nevent_gam(j,k) = 0d0
            nevent_e2e(j,k) = 0d0
            nevent_e2mu(j,k) = 0d0
            nevent_mu2mu(j,k) = 0d0
            nevent_mu2e(j,k) = 0d0
            do i = 1,nbins
               hErec_e(i,j,k) = 0d0
               hErec_mu(i,j,k) = 0d0
               hErec_pi0(i,j,k) = 0d0
               hErec_gam(i,j,k) = 0d0
               hErec_e2e(i,j,k) = 0d0
               hErec_e2mu(i,j,k) = 0d0
               hErec_mu2mu(i,j,k) = 0d0
               hErec_mu2e(i,j,k) = 0d0
            enddo
         enddo
      enddo
      do l3 = -1,1,2
         do l2 = 1,ndetect
            do l1 = -3,3
               do k = -3,3
                  do l4 = 1,imaxint
                     do i = 1,nbins
                        hErec_e(i,l2,l3) = hErec_e(i,l2,l3) 
     &                    +hErec_raw_dat(i,l4,1,k,l1,l2,l3)
                        nevent_e(l2,l3) = nevent_e(l2,l3) 
     &                       +hErec_raw_dat(i,l4,1,k,l1,l2,l3)
                        hErec_mu(i,l2,l3) = hErec_mu(i,l2,l3) 
     &                       +hErec_raw_dat(i,l4,2,k,l1,l2,l3)
                        nevent_mu(l2,l3) = nevent_mu(l2,l3) 
     &                       +hErec_raw_dat(i,l4,2,k,l1,l2,l3)
                        hErec_pi0(i,l2,l3) = hErec_pi0(i,l2,l3) 
     &                       +hErec_raw_dat(i,l4,3,k,l1,l2,l3)
                        nevent_pi0(l2,l3) = nevent_pi0(l2,l3) 
     &                       +hErec_raw_dat(i,l4,3,k,l1,l2,l3)
                        hErec_gam(i,l2,l3) = hErec_gam(i,l2,l3) 
     &                       +hErec_raw_dat(i,l4,4,k,l1,l2,l3)
                        nevent_gam(l2,l3) = nevent_gam(l2,l3) 
     &                       +hErec_raw_dat(i,l4,4,k,l1,l2,l3)
                        hErec_e2e(i,l2,l3) = hErec_e2e(i,l2,l3) 
     &                       +hErec_dat(i,l4,1,k,l1,l2,l3)
                        nevent_e2e(l2,l3) = nevent_e2e(l2,l3) 
     &                       +hErec_dat(i,l4,1,k,l1,l2,l3)
                        hErec_e2mu(i,l2,l3) = hErec_e2mu(i,l2,l3) 
     &                       +hErec_dat(i,l4,6,k,l1,l2,l3)
                        nevent_e2mu(l2,l3) = nevent_e2mu(l2,l3) 
     &                       +hErec_dat(i,l4,6,k,l1,l2,l3)
                        hErec_mu2mu(i,l2,l3) = hErec_mu2mu(i,l2,l3) 
     &                       +hErec_dat(i,l4,2,k,l1,l2,l3)
                        nevent_mu2mu(l2,l3) = nevent_mu2mu(l2,l3) 
     &                       +hErec_dat(i,l4,2,k,l1,l2,l3)
                        hErec_mu2e(i,l2,l3) = hErec_mu2e(i,l2,l3) 
     &                       +hErec_dat(i,l4,5,k,l1,l2,l3)
                        nevent_mu2e(l2,l3) = nevent_mu2e(l2,l3) 
     &                       +hErec_dat(i,l4,5,k,l1,l2,l3)
                     enddo
                  enddo
               enddo
            enddo
         enddo
      enddo
      do l3 = -1,1,2
         do l2 = 1,ndetect
            do i = 1,nbins
               nevent_elike(l2,l3) = nevent_elike(l2,l3) 
     &              +hErec_like_dat(i,1,l2,l3)
               nevent_mulike(l2,l3) = nevent_mulike(l2,l3) 
     &              +hErec_like_dat(i,2,l2,l3)
            enddo
         enddo
      enddo

      return
      end
