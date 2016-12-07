      subroutine get_1pi0dist(evform,fxsec_r,event_in,nevent_in,x,nbins
     &     ,event_out)
C     ****************************************************
C     By Yoshitaro Takaesu @KIAS OCT 26 2013
C     ****************************************************
      implicit none
C     GLOBAL VARIABLES
      include 'inc/params.inc' 
      real*8 ratioRS,ratioCO
      common /ratioRS/ ratioRS,ratioCO
C     CONSTANTS
C     ARGUMENTS 
      integer nbins,evform,fxsec_r
      real*8 x(0:nbins),event_in(maxnbin),event_out(maxnbin,imaxint)
C     LOCAL VARIABLES 
      integer i,j,k
      integer int_mode,max_mode,ierr
      parameter (max_mode=10)
      real*8 Ev,frac,event_mode(nbins),hevent_mode(nbins),rievent
      real*8 rievent2,rnevent_mode
      real*8 event_mode2(nbins),hevent_mode2(nbins),rnevent_mode2
      real*8 z(zdim),hevent_out(nbins),nevent_in
      real*8 Nncrs,Nncco,r0
C     EXTERNAL FUNCTIONS
      real*8 fErec_1pi0dist,fpi0mom,fpi0mom_old,fpi0_Erec
      real*8 fpi0_Erec_old
      external fErec_1pi0dist,fpi0mom,fpi0mom_old,fpi0_Erec
      external fpi0_Erec_old
C     ----------
C     BEGIN CODE
C     ----------
      do j = 1,imaxint
         do i = 1,nbins
            event_out(i,j) = 0d0
         enddo
      enddo
      do i = 1,nbins
         Ev = ( x(i) +x(i-1) )/2d0
         rievent = event_in(i)
c     Compute ratio between NCRS and NCCO 
         if (ipi0unc.eq.0) then
            ratioRS = 0d0 ! dummy parameter in this case. not be used
            ratioCO = 0d0 ! dummy parameter in this case. not be used
         elseif (ipi0unc.eq.1) then
            int_mode = 2
            call get_xsecfrac3(Ev,icc,ipi0xsec,int_mode,detect,frac)
            Nncrs = frac*rievent
            int_mode = 4
            call get_xsecfrac3(Ev,icc,ipi0xsec,int_mode,detect,frac)
            Nncco = frac*rievent
            if (Nncrs.le.0) then
               ratioRS = 0d0
               ratioCO = 1d0
            else
               r0 = Nncco/Nncrs         
               ratioRS = (1d0 +r0)/(1d0 +fxsec_r*r0) 
               ratioCO = fxsec_r*(1d0 +r0)/(1d0 +fxsec_r*r0) 
            endif
         elseif (ipi0unc.eq.2) then
            ratioRS = 0d0 ! dummy parameter in this case. not be used
            ratioCO = 0d0 ! dummy parameter in this case. not be used
         endif
c     Calculte distributions
cccc
cccc  NOTICE: Only the NC1pi0BG (polfit and old) mode supports the separate int_mode output.
cccc
         if ((ipi0xsec.eq.1).or.(ipi0xsec.eq.2)) then 
c         if (ipi0xsec.eq.1) then 
            do int_mode = 1,imaxint
               call get_1pi0dist_mode(int_mode,Ev,rievent,i,event_mode)
               do j = 1,nbins
                  event_out(j,int_mode) = event_out(j,int_mode) 
     &                 +event_mode(j)
               enddo
            enddo
         else
            int_mode = 1
            call get_1pi0dist_mode(int_mode,Ev,rievent,i,event_mode)
            do j = 1,nbins
               event_out(j,int_mode) = event_out(j,int_mode) 
     &              +event_mode(j)
            enddo
         endif
      enddo

      return
      end


      subroutine get_1pi0dist_mode(int_mode,Ev,rievent,ibin,event_mode)
C     ****************************************************
C     By Yoshitaro Takaesu @KIAS Mar 12 2014
C     Modified: Dec 7 2016: make dependent on the ismear switch
C
C     ****************************************************
      implicit none
C     GLOBAL VARIABLES
      include 'inc/params.inc'
      include 'inc/main.inc'
      include 'inc/minfunc.inc'
      include 'inc/get_event.inc'
      real*8 ratioRS,ratioCO
      common /ratioRS/ ratioRS,ratioCO
C     CONSTANTS
C     ARGUMENTS 
      integer int_mode,ibin
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
CCC   Initialization
      ierr = 0
      do i = 1,nbins_basic
         event_mode(i) = 0d0
      enddo

CCC   For resonant 1pi0BG
      if (int_mode.eq.2) then 
CCC      For independent uncertainties for NCRS 1pi0BG
         if (ipi0unc.eq.0) then 
            fpi = fxsec_pirs
CCC      For correlated uncertainties for NCRS 1pi0BG
         elseif (ipi0unc.eq.1) then 
            fpi = ratioRS
CCC      No uncertainties for NCRS 1pi0BG
         elseif (ipi0unc.eq.2) then 
            fpi = 1d0
         endif

CCC   For coherent 1pi0BG
      elseif (int_mode.eq.4) then 
CCC      For independent uncertainties for NCCO 1pi0BG
         if (ipi0unc.eq.0) then
            fpi = fxsec_pico
CCC      For correlated uncertainties for NCCO 1pi0BG
         elseif (ipi0unc.eq.1) then
            fpi = ratioCO
CCC      No uncertainties for NCCO 1pi0BG
         elseif (ipi0unc.eq.2) then 
            fpi = 1d0
         endif

CCC   For other possibilities
      else
         fpi = 1d0
      endif      

CCC   Calculate event number fraction of a mode
      call get_xsecfrac3(Ev,icc,ipi0xsec,int_mode,detect,frac)

      rievent2 = rievent*frac*fpi
      z(1) = detect
      z(2) = Ev
      z(3) = int_mode

CCC   NC 1pi0 distribution
      if (ipi0xsec.eq.0) then      

CC       Erec distribution
         if (ipi0dist.eq.0) then   
            write(97,*) "ERROR:get_1pi0dist.f:fpi0_Erec still",
     &           " does not support the separate int_mode output."
            stop
c            call MakeHisto1D(nout,fpi0_Erec,z,rievent2,nbins,x
c     &           ,evform,serror,snmax,ihisto,event_mode,hevent_mode
c     &           ,rnevent_mode,ierr)

CC       Momentum distribution
         elseif (ipi0dist.eq.1) then 
            write(97,*) "ERROR:get_1pi0dist.f:fpi0_mom still",
     &           " does not support the separate int_mode output."
            stop
c            call MakeHisto1D(nout,fpi0mom,z,rievent2,nbins,x
c     &           ,evform,serror,snmax,ihisto,event_mode,hevent_mode
c     &           ,rnevent_mode,ierr)
         endif

CCC   NC pi0-bg (polfit missID) distribution
      elseif (ipi0xsec.eq.1) then  

CC       Erec distribution
         if (ipi0dist.eq.0) then   
            if (ismear.eq.0) then
               event_mode(ibin) = rievent2
            elseif (ismear.eq.1) then
               call MakeHisto1D(nout,fErec_1pi0dist,z,rievent2,nbins,x
     &              ,evform,serror,snmax,ihisto,event_mode,hevent_mode
     &              ,rnevent_mode,ierr)
            endif
         endif

CC       Momentum distribution
         if (ipi0dist.eq.1) then
            write(97,*) "ERROR:get_1pi0dist.f:momentum distribution",
     &           " for NC pi0-bg (polfit) is not supported."
            stop
         endif

CCC   NC pi0-bg (old_func missID) distribution
      elseif (ipi0xsec.eq.2) then  

CC       Erec distribution
         if (ipi0dist.eq.0) then
            if (ismear.eq.0) then
               event_mode(ibin) = rievent2
            elseif (ismear.eq.1) then
               call MakeHisto1D(nout,fpi0_Erec_old,z,rievent2,nbins,x
     &              ,evform,serror,snmax,ihisto,event_mode,hevent_mode
     &              ,rnevent_mode,ierr)
            endif

CC       Momentum distribution
         elseif (ipi0dist.eq.1) then ! momentum dist
            write(97,*) "ERROR:get_1pi0dist.f:momentum distribution",
     &           " for NC pi0-bg (old) is not supported."
            stop
c            call MakeHisto1D(nout,fpi0mom_old,z,rievent2,nbins,x
c     &           ,evform,serror,snmax,ihisto,event_mode,hevent_mode
c     &           ,rnevent_mode,ierr)
         endif

      endif

CCC   Show ERROE message
      if (ierr.ne.0) then
         write(97,*) "get_1pi0dist: MakeHisto1D ierr = 1"
      endif         
      
      return
      end
      
