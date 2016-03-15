C -*- fortran -*-
      do k = -1,1,2
         do j = 1,ndetect
            nevent_e(j,k) = 0d0
            nevent_mu(j,k) = 0d0
            nevent_pi0(j,k) = 0d0
            nevent_gam(j,k) = 0d0
            do i = 1,nbins
               hErec_e(i,j,k) = 0d0
               hErec_mu(i,j,k) = 0d0
               hErec_pi0(i,j,k) = 0d0
               hErec_gam(i,j,k) = 0d0
            enddo
         enddo
      enddo
      do l3 = -1,1,2
         do l2 = 1,ndetect
            do l1 = -3,3
               do k = -3,3
                  do i = 1,nbins
                     hErec_e(i,l2,l3) = hErec_e(i,l2,l3) 
     &                    +hErec_dat(i,1,k,l1,l2,l3)
                     nevent_e(l2,l3) = nevent_e(l2,l3) 
     &                    +hErec_dat(i,1,k,l1,l2,l3)
                     hErec_mu(i,l2,l3) = hErec_mu(i,l2,l3) 
     &                    +hErec_dat(i,2,k,l1,l2,l3)
                     nevent_mu(l2,l3) = nevent_mu(l2,l3) 
     &                    +hErec_dat(i,2,k,l1,l2,l3)
                     hErec_pi0(i,l2,l3) = hErec_pi0(i,l2,l3) 
     &                    +hErec_dat(i,3,k,l1,l2,l3)
                     nevent_pi0(l2,l3) = nevent_pi0(l2,l3) 
     &                    +hErec_dat(i,3,k,l1,l2,l3)
                     hErec_gam(i,l2,l3) = hErec_gam(i,l2,l3) 
     &                    +hErec_dat(i,4,k,l1,l2,l3)
                     nevent_gam(l2,l3) = nevent_gam(l2,l3) 
     &                    +hErec_dat(i,4,k,l1,l2,l3)
                  enddo
               enddo
            enddo
         enddo
      enddo
      do k = -1,1,2
         do j = 1,ndetect
            nevent_elike(j,k) = nevent_e(j,k) +nevent_pi0(j,k) 
     &           +nevent_gam(j,k)
            nevent_mulike(j,k) = nevent_mu(j,k)
            do i = 1,nbins
               hErec_elike(i,j,k) = hErec_e(i,j,k) +hErec_pi0(i,j,k) 
     &              +hErec_gam(i,j,k)
               hErec_mulike(i,j,k) = hErec_mu(i,j,k)
            enddo
         enddo
      enddo  

