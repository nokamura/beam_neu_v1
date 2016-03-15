      subroutine MakeHisto1D(nout,f,z,nevent,nbins,x,evform,error,nmax
     &     ,mode,event,hevent,nevent_th,ierr)
C     ****************************************************
C     By Yoshitaro Takaesu @KIAS AUG 25 2012
C
C     mode: 0:no histogram 1:histo 2:left-point histo 3:center-point histo
C     ****************************************************
      implicit none

C     ARGUMENTS 
      integer nbins,nout,evform,nmax,mode,ierr
      real*8 f,x(0:nbins),error,z(200),nevent
      real*8 event(nbins),hevent(nbins),nevent_th
C     LOCAL VARIABLES 
      integer i,j
      integer nmerge
      real*8 sumy,xi,y(nbins),dy(nbins),rnevent,hy(nbins)
      real*8 dx(nbins),xxmin,xxmax,mx,cx
C     EXTERNAL FUNCTIONS
      external f
C     ----------
C     BEGIN CODE
C     ----------
      sumy = 0d0
      ierr = 0
      nevent_th = 0d0
      do i = 1,nbins
         event(i) = 0d0
         hevent(i) = 0d0
      enddo

      if (mode.eq.0) then
         do i = 1,nbins
            cx = ( x(i) +x(i-1) )/2d0
            event(i) = f(cx,z)
            hevent(i) = 0d0
            nevent_th = 0d0
         enddo
      elseif (mode.eq.1) then
         do i = 1,nbins
            dx(i) = x(i) -x(i-1)
            call hsimp1D(f,x(i-1),x(i),z,y(i),error,nmax,ierr)
            if (ierr.ne.0) then
               write(nout,*) "MakeHisto1D:",
     &              " Integration does not converge"
               return
            endif
            sumy = sumy +y(i)
         enddo
         
         if (nevent.lt.0d0) then
            sumy = 1d0
            rnevent = 1d0
         else
            rnevent = nevent
         endif

         do i = 1,nbins
            if (sumy.eq.0d0) then
               event(i) = 0d0
            else
               if (evform.eq.1) then
                  event(i) = idint( y(i)*rnevent/sumy )
               elseif (evform.eq.2) then
                  event(i) = y(i)*rnevent/sumy
               endif
            endif
            hevent(i) = event(i)/dx(i)
            nevent_th = nevent_th +event(i)
         enddo
      elseif (mode.eq.2) then
         do i = 1,nbins
            dx(i) = x(i) -x(i-1)         
            mx = x(i-1) +dx(i)/2d0
            y(i) = f(mx,z)*dx(i)
            sumy = sumy +y(i)
         enddo
         if (nevent.lt.0d0) then
            sumy = 1d0
            rnevent = 1d0
         else
            rnevent = nevent
         endif
         do i = 1,nbins
            if (sumy.eq.0d0) then
               event(i) = 0d0
            else
               if (evform.eq.1) then
                  event(i) = idint( y(i)*rnevent/sumy )
               elseif (evform.eq.2) then
                  event(i) = y(i)*rnevent/sumy
               endif
            endif
            hevent(i) = event(i)/dx(i)
            nevent_th = nevent_th +event(i) 
         enddo
      elseif (mode.eq.3) then
         do i = 1,nbins
            dx(i) = x(i) -x(i-1)
            xi = x(i-1)
            hevent(i) = f(xi,z)
            event(i) = hevent(i)*dx(i)
            nevent_th = nevent_th +event(i)
         enddo
      endif

      return
      end
