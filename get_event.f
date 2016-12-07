c      subroutine get_event(z,iproc,eventout,neventout)
      subroutine get_event(z,iproc,eventout)
C     **********************************************************
C     By Yoshitaro Takaesu @KIAS JAN 7 2014
C     Modified: Dec 7 2016: Remove the ismear switch. It is put
C                           in the get_nudist and get_1pi0dist
C                           subroutines directly.
C     **********************************************************
      implicit none
C     GLOBAL VARIABLES
      include 'inc/params.inc'
      include 'inc/main.inc'
      include 'inc/minfunc.inc'
      include 'inc/get_event.inc'
C     CONSTANTS
C     ARGUMENTS 
      integer iproc
      real*8 z(zdim)
      real*8 eventout(maxnbin,imaxint)
C     LOCAL VARIABLES 
      integer i,j
      integer ierr,ite,ibins
      real*8 event_tmp1(maxnbin),hevent_tmp1(maxnbin),nevent_tmp1
      real*8 event_tmp2(maxnbin,imaxint),nevent_tmp2
C     EXTERNAL FUNCTIONS
      real*8 hfunc1D
      external hfunc1D
C     ----------
C     BEGIN CODE
C     ----------
CCC   Set basic bins
      call bining_basic(Emin,Emax,basic_binsize,nbins_basic,x_basic)

CCC   Calculate unsmeared event distributions
      call MakeHisto1D(nout,hfunc1D,z,rnevent_ren,nbins_basic
     &     ,x_basic,evform,serror,snmax,ihisto,event_tmp1
     &     ,hevent_tmp1,nevent_tmp1,ierr) 

CCC   Initialize the "event_tmp2" array
      do j = 1,imaxint
         do i = 1,nbins_basic
            event_tmp2(i,j) = 0d0
         enddo
      enddo

CCC   Calculate smeared event distributions
c      if (ismear.eq.0) then
c         do i = 1,nbins_basic
c            event_tmp2(i,1) = event_tmp1(i)
c         enddo
c         nevent_tmp2 = nevent_tmp1
c      elseif (ismear.eq.1) then

CCC      CC interactions
         if (icc.eq.1) then  
            if (iproc.eq.abs(detect)) then
               if (detect.gt.0) then
                  fxsec_CCQE = z(7)
                  fxsec_CCRes = z(51)
               elseif (detect.lt.0) then
                  fxsec_CCQE = z(8)
                  fxsec_CCRes = z(52)
               endif
               call get_nudist(detect,event_tmp1,event_tmp2
     &              ,nevent_tmp2)
            endif

CCC      NC intractions
         elseif (icc.eq.2) then 
            if (iproc.eq.3) then
               if (ipi0unc.eq.0) then
                  fxsec_pirs = z(54)
                  fxsec_pico = z(55)
                  fxsec_r = 0d0 ! dymmy parameter in this case. not be used
               elseif (ipi0unc.eq.1) then
                  fxsec_r = z(56)
               endif
               call get_1pi0dist(evform,fxsec_r,event_tmp1
     &              ,nevent_tmp1,x_basic,nbins_basic,event_tmp2)
            endif
         endif      
c      endif

CCC   Make output "eventout" array by combining adjuscent bins to ensure >= 10 events in all the bins.
      do j = 1,imaxint
         ite = 0
         ibins = 1
         do i = 1,nbins
            eventout(i,j) = 0d0
         enddo
         do i = 1,nbins_basic
            ite = ite +1
            eventout(ibins,j) = eventout(ibins,j) +event_tmp2(i,j) 
            if (ite.eq.binsize_factor) then
               ibins = ibins +1
               ite = 0
            endif
         enddo
      enddo
      
      return
      end
