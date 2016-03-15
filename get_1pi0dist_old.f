      subroutine get_1pi0dist(ixsec,idist,inutype,evform,serror,snmax
     &     ,ihisto,dumm,event_in,nevent_in,x,nbins,event_out
     &     ,rnevent_out)
C     ****************************************************
C     By Yoshitaro Takaesu @KIAS OCT 26 2013
C     ****************************************************
      implicit none
C     GLOBAL VARIABLES
C     CONSTANTS
C     ARGUMENTS 
      integer nbins,inutype,evform,snmax,ihisto,idist,ixsec
      real*8 x(0:nbins),event_in(nbins),event_out(nbins),serror
      real*8 rnevent_out,dumm
C     LOCAL VARIABLES 
      integer i,j,k
      integer int_mode,max_mode,icc,nout,ierr
      parameter (max_mode=10)
      real*8 Ev,frac,event_mode(nbins),hevent_mode(nbins),rievent
      real*8 rievent2,rnevent_mode
      real*8 event_mode2(nbins),hevent_mode2(nbins),rnevent_mode2
      real*8 z(200),hevent_out(nbins),nevent_in
C     EXTERNAL FUNCTIONS
      real*8 fErec_1pi0dist,fpi0mom,fpi0mom_old,fpi0_Erec
      real*8 fpi0_Erec_old
      external fErec_1pi0dist,fpi0mom,fpi0mom_old,fpi0_Erec
      external fpi0_Erec_old
C     ----------
C     BEGIN CODE
C     ----------
      nout = 6
      icc = 2 ! NC interaction
      do i = 1,nbins
         event_out(i) = 0d0
      enddo
      do i = 1,nbins
         Ev = ( x(i) +x(i-1) )/2d0
         rievent = event_in(i)
c         do int_mode = 1,4
         do int_mode = 2,2 ! NCRS
c         do int_mode = 4,4 ! NCCO
            call get_1pi0dist_mode(inutype,ixsec,idist,int_mode,Ev
     &           ,rievent,i,event_mode)
            do j = 1,nbins
               event_out(j) = event_out(j) +event_mode(j)
            enddo
         enddo
      enddo
      rnevent_out = 0d0
      do i = 1,nbins
         rnevent_out = rnevent_out +event_out(i)
      enddo

      return
      end


      subroutine get_1pi0dist_mode(inutype,ixsec,idist,int_mode,Ev
     &     ,rievent,ibin,event_mode)
C     ****************************************************
C     By Yoshitaro Takaesu @KIAS Mar 12 2014
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
      integer inutype,int_mode,ibin,ixsec,idist
      real*8 Ev,rievent,event_mode(nbins_basic)
C     LOCAL VARIABLES 
      integer i,j,k
      integer ierr
      real*8 hevent_mode(nbins_basic),rnevent_mode
      real*8 z(3),frac,rievent2,fpi
      real*8 hevent_mode2(nbins_basic),rnevent_mode2
      real*8 event_mode2(nbins_basic)
C     EXTERNAL FUNCTIONS
      real*8 fpi0_Erec,fpi0mom,fErec_1pi0dist,fpi0_Erec_old,fpi0mom_old
      external fpi0_Erec,fpi0mom,fErec_1pi0dist,fpi0_Erec_old
      external fpi0mom_old
C     ----------
C     BEGIN CODE
C     ----------
c      int_mode = 5 ! for test
      call get_xsecfrac3(Ev,icc,ixsec,int_mode,inutype,frac)
      if (int_mode.eq.2) then
         fpi = fxsec_pirs
      elseif (int_mode.eq.4) then
         fpi = fxsec_pico
      else
         fpi = 1d0
      endif
      rievent2 = rievent*frac*fpi
      z(1) = inutype
      z(2) = Ev
      z(3) = int_mode
      if (ixsec.eq.0) then      ! NC 1pi0
         if (idist.eq.0) then   ! Erec dist
c     rievent2 = -1
c            rievent = 1
c            do j = 1,nbins
c               event_mode(j) = 0d0
c            enddo
c            rnevent_mode = 0d0
c     do j = 5,5
c            do j = 1,1
c               z(3) = j 
            call MakeHisto1D(nout,fpi0_Erec,z,rievent2,nbins,x
     &           ,evform,serror,snmax,ihisto,event_mode,hevent_mode
     &           ,rnevent_mode,ierr)
c               call MakeHisto1D(nout,fpi0_Erec,z,rievent2,nbins,x
c     &              ,evform,serror,snmax,ihisto,event_mode2
c     &              ,hevent_mode2,rnevent_mode2,ierr)
c               do k = 1,nbins
c                  event_mode(k) = event_mode(k) 
c     &                 +rievent*event_mode2(k)
c               enddo
c            enddo
c            do j = 1,nbins
c               rnevent_mode = rnevent_mode +event_mode(j)
c            enddo
         elseif (idist.eq.1) then ! momentum dist
            call MakeHisto1D(nout,fpi0mom,z,rievent2,nbins,x
     &           ,evform,serror,snmax,ihisto,event_mode,hevent_mode
     &           ,rnevent_mode,ierr)
         endif
      elseif (ixsec.eq.1) then  ! NC pi0-bg (polfit middID)
         if (idist.eq.0) then   ! Erec dist
            call MakeHisto1D(nout,fErec_1pi0dist,z,rievent2,nbins,x
     &           ,evform,serror,snmax,ihisto,event_mode,hevent_mode
     &           ,rnevent_mode,ierr)
         endif
      elseif (ixsec.eq.2) then  ! NC pi0-bg (old_func middID)
         if (idist.eq.0) then   ! Erec dist
c            z(3) = 1
            call MakeHisto1D(nout,fpi0_Erec_old,z,rievent2,nbins,x
     &           ,evform,serror,snmax,ihisto,event_mode,hevent_mode
     &           ,rnevent_mode,ierr)
         elseif (idist.eq.1) then ! momentum dist
            call MakeHisto1D(nout,fpi0mom_old,z,rievent2,nbins,x
     &           ,evform,serror,snmax,ihisto,event_mode,hevent_mode
     &           ,rnevent_mode,ierr)
         endif
      endif
      if (ierr.ne.0) then
         write(97,*) "get_1pi0dist: MakeHisto1D ierr = 1"
      endif         
      
      return
      end
      
