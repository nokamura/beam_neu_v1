      subroutine apply_Ereccut(iD,hErec_in,hErec_out)
C     ****************************************************
C     By Yoshitaro Takaesu @KIAS JAN 7 2014
C     ****************************************************
      implicit none
C     GLOBAL VARIABLES
      include 'inc/params.inc'
      include 'inc/main.inc'
C     CONSTANTS
C     ARGUMENTS 
      integer iD
      real*8 hErec_in(maxnbin,nmode_like,ndetect,-1:1)
      real*8 hErec_out(maxnbin,nmode_like,ndetect,-1:1)
C     LOCAL VARIABLES 
      integer i,bfactor,ibeam,ipos,icount,lmode,idebug
      real*8 e_min,e_max,m_min,m_max,y_dum(0:maxnbin),binsize_chi2
C     EXTERNAL FUNCTIONS
C     ----------
C     BEGIN CODE
C     ----------
c     check input data
      idebug = 1
      if (idebug.eq.1) then
         do ibeam = -1,1,2
            do lmode = 1,nmode_like
               do i = 1,maxnbin
                  if ((hErec_in(i,lmode,iD,ibeam)*0).ne.0) then
                     write(*,*) "ERROR:apply_Ereccut.f"
     &                    ,":Input data is invalid"
                     stop
                  endif
               enddo
            enddo
         enddo
      endif

      if (iD.eq.1) then
         e_min = Semin
         e_max = Semax
         m_min = Smmin
         m_max = Smmax
         bfactor = Sbfactor
      elseif(iD.eq.2) then
         e_min = Oemin
         e_max = Oemax
         m_min = Ommin
         m_max = Ommax
         bfactor = Obfactor
      elseif(iD.eq.3) then
         e_min = Kemin
         e_max = Kemax
         m_min = Kmmin
         m_max = Kmmax
         bfactor = Kbfactor
      endif         

      binsize_chi2 = basic_binsize*bfactor
      lmode = 1
      call bining_x(e_min,e_max,binsize_chi2,nbins_chi2(lmode,iD)
     &     ,x_chi2(0,lmode,iD),y_dum)
      ibeam = 1
      call apply_Ereccut_unit(e_min,e_max,bfactor,nbins,x
     &     ,hErec_in(1,lmode,iD,ibeam),hErec_out(1,lmode,iD,ibeam))
      ibeam = -1
      call apply_Ereccut_unit(e_min,e_max,bfactor,nbins,x
     &     ,hErec_in(1,lmode,iD,ibeam),hErec_out(1,lmode,iD,ibeam))

      lmode = 2
      call bining_x(m_min,m_max,binsize_chi2,nbins_chi2(lmode,iD)
     &     ,x_chi2(0,lmode,iD),y_dum)
      ibeam = 1
      call apply_Ereccut_unit(m_min,m_max,bfactor,nbins,x
     &     ,hErec_in(1,lmode,iD,ibeam),hErec_out(1,lmode,iD,ibeam))
      ibeam = -1
      call apply_Ereccut_unit(m_min,m_max,bfactor,nbins,x
     &     ,hErec_in(1,lmode,iD,ibeam),hErec_out(1,lmode,iD,ibeam))

      return
      end


      subroutine apply_Ereccut_unit(min,max,bfactor,nbins,xx
     &     ,hErec_in,hErec_out)
C     ****************************************************
C     By Yoshitaro Takaesu @KIAS Feb 13 2014
C     ****************************************************
      implicit none
C     GLOBAL VARIABLES
      include 'inc/params.inc'
C     CONSTANTS
C     ARGUMENTS 
      integer nbins,bfactor
      real*8 min,max,xx(0:maxnbin)
      real*8 hErec_in(maxnbin),hErec_out(maxnbin)
C     LOCAL VARIABLES 
      integer i,ipos,icount
C     EXTERNAL FUNCTIONS
C     ----------
C     BEGIN CODE
C     ----------
      do i = 1,nbins
         hErec_out(i) = 0d0
      enddo
      ipos = 1
      icount = 0
      do i = 0,nbins-1
         if ((xx(i).ge.min).and.(xx(i+1).le.(max+1d-3))) then
            icount = icount +1
            hErec_out(ipos) = hErec_out(ipos) +hErec_in(i+1)
            if (icount.eq.bfactor) then
               icount = 0
               ipos = ipos +1
            endif
         endif
      enddo

      return
      end

