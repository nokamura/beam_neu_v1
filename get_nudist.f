      subroutine get_nudist(inutype,event_in,event_out)
c     &     ,rnevent_out)
C     ****************************************************
C     By Yoshitaro Takaesu @KIAS OCT 11 2013
C
C     modified: JAN 21 2014: iCCQE, iCCRes added
C                            params.inc, minfunc.inc added
C     Mar 27 2014: separated output for each int_mode 
C     ****************************************************
      implicit none
C     GLOBAL VARIABLES
      include 'inc/params.inc'
      include 'inc/minfunc.inc'
      include 'inc/get_event.inc'
C     CONSTANTS
C     ARGUMENTS 
      integer inutype
      real*8 event_in(maxnbin),event_out(maxnbin,imaxint)
      real*8 rnevent_out
C     LOCAL VARIABLES 
      integer i,j
      integer int_mode,ierr
      real*8 Ev,rievent,hevent_out(nbins_basic),event_mode(nbins_basic)
      real*8 frac,rievent2,z(zdim),hevent_mode(nbins_basic),rnevent_mode
C     EXTERNAL FUNCTIONS
      real*8 fErec_dist
      external fErec_dist
C     ----------
C     BEGIN CODE
C     ----------
      do j = 1,imaxint
         do i = 1,maxnbin
            event_out(i,j) = 0d0
         enddo
      enddo
      do i = 1,nbins_basic
         Ev = ( x_basic(i) +x_basic(i-1) )/2d0
         rievent = event_in(i)
         if (iCCQE.eq.1) then  ! include CCQE events
            do int_mode = 1,2            
               call get_nudist_mode2(inutype,int_mode,Ev,rievent,i
     &              ,event_mode)
               do j = 1,nbins_basic
                  event_out(j,int_mode) = event_out(j,int_mode) 
     &                 +event_mode(j)*fxsec_CCQE
               enddo
            enddo
         endif
         if (iCCRes.eq.1) then ! include Res events
            do int_mode = 3,4
               call get_nudist_mode2(inutype,int_mode,Ev,rievent,i
     &              ,event_mode)
               do j = 1,nbins_basic
                  event_out(j,int_mode) = event_out(j,int_mode) 
     &                 +event_mode(j)*fxsec_CCRes
               enddo
            enddo
         endif         
      enddo

      return
      end


      subroutine get_nudist_mode2(inutype,int_mode,Ev,rievent,ibin
     &     ,event_mode)
C     ****************************************************
C     By Yoshitaro Takaesu @KIAS OCT 11 2013
C
C     ****************************************************
      implicit none
C     GLOBAL VARIABLES
      include 'inc/params.inc'
      include 'inc/main.inc'
      include 'inc/minfunc.inc'
      include 'inc/get_event.inc'
C     CONSTANTS
C     ARGUMENTS 
      integer inutype,int_mode,ibin
      real*8 Ev,rievent,event_mode(nbins_basic)
C     LOCAL VARIABLES 
      integer i,j
      integer ierr
      real*8 hevent_mode(nbins_basic),rnevent_mode
      real*8 z(3),frac,rievent2
C     EXTERNAL FUNCTIONS
      real*8 fErec_dist
      external fErec_dist
C     ----------
C     BEGIN CODE
C     ----------
      call get_xsecfrac3(Ev,icc,1,int_mode,inutype,frac)
      if ((inutype.gt.0).and.(int_mode.eq.3)) then   !  We ignore the neutrino CCRes-H interaction as an approximation. (To be considered later)
         frac = 0d0
      endif
      if (iCCQE.eq.1) then
c         frac = 0.5d0           ! when using xsec_CCQE for CCQE event
         continue
      endif
      rievent2 = rievent*frac
      z(1) = inutype
      z(2) = int_mode
      z(3) = Ev
      call MakeHisto1D(nout,fErec_dist,z,rievent2,nbins_basic,x_basic
     &     ,evform,serror,snmax,ihisto,event_mode,hevent_mode
     &     ,rnevent_mode,ierr)
      if (ierr.ne.0) then
         write(*,*) "get_nudist_mode: MakeHisto1D ierr = 1"
      endif
         
      return
      end
