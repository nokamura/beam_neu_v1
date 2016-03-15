      subroutine write_results_mode1
C     ****************************************************
C     By Yoshitaro Takaesu @KIAS JAN 7 2014
C
C     Output chi^2 minimization results to files
C     ****************************************************
      implicit none
C     GLOBAL VARIABLES
      include 'inc/params.inc'
      real*8 dchisqmin,fedm,chisqmin_true,errdef,istat
      common /minuit/ dchisqmin,fedm,chisqmin_true,errdef,istat
C     CONSTANTS
C     ARGUMENTS 
C     LOCAL VARIABLES 
      integer i,iflag,ierr
      real*8 name(maxnparam),pval(maxnparam),perr(maxnparam)
      real*8 plo(maxnparam),phi(maxnparam)
C     EXTERNAL FUNCTIONS
      real*8 minfunc
      external minfunc
C     ----------
C     BEGIN CODE
C     ----------
c$$$  iD = 3
c$$$  nm = 2
c$$$  iev = 1
c$$$  
c$$$  write(*,*) (evx(i,iev,iD),i=0,evnbins(iev,iD))
c$$$  write(*,*) "data"
c$$$  write(*,*) (total_event_dat(i+1,nm,iD),i=0,evnbins(iev,iD)-1)
c$$$  write(*,*) "theory"
c$$$  write(*,*) (total_event_th(i+1,nm,iD),i=0,evnbins(iev,iD)-1)
c$$$  write(*,*)
c$$$  endif
c$$$  
c$$$  do i = 1,ndetect
c$$$  nevent_tot(i) = 0d0
c$$$  enddo
c$$$  do i = 1,ndetect
c$$$  do k = 1,eventdim
c$$$  nevent_tot(i) = nevent_tot(i) +nevent_dat(k,i)
c$$$  enddo
c$$$  enddo
c$$$  write(*,*) nevent_dat(2,1)
c$$$  if (iSK.eq.1) then
c$$$  write(*,*) "SK"
c$$$  write(*,*) 
c$$$  iD = 1
c$$$  include 'inc/output.inc'
c$$$  endif
c$$$  if (iOki.eq.1) then
c$$$  write(*,*) "Oki"
c$$$  write(*,*)
c$$$  iD = 2
c$$$  include 'inc/output.inc'
c$$$  endif
c$$$  if (iKr.eq.1) then
c$$$  write(*,*) "Korea"
c$$$  write(*,*)
c$$$  iD = 3
c$$$  include 'inc/output.inc'
c$$$  endif
c$$$  sum_nevent_tot = 0d0
c$$$  do i =1,ndetect
c$$$  sum_nevent_tot = sum_nevent_tot +nevent_tot(i)
c$$$  enddo
c$$$         write(*,*) "Total Events = ",sum_nevent_tot
c$$$  write(*,*) 
c$$$  write(*,*) "Delta Chi^2_min =",dchisqmin," +- ",fedm
c$$$  write(*,*)      
      return
      end
