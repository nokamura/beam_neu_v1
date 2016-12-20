      program main
      implicitnone
      include 'inc/params.inc'
      include 'inc/main.inc'

      character*10 name(10),iname
      integer i,j,k,l1,l2,l3
      integer ierr,lencm
      integer iflag,iD,nm,iev
      real*8 arg(10),pval(maxnparam),perr(maxnparam),plo(maxnparam)
      real*8 phi(maxnparam),gint
      real*8 chisqmin_wrong
      real*8 z(zdim),dchisq,grad,futil
      real*8 mean_nh,error_nh,mean_error_nh,error_error_nh
      real*8 mean_ih,error_ih,mean_error_ih,error_error_ih
      real*8 mean_dchi2min_nh,mean_dchi2min_ih
      real*8 sum_nevent_tot
      integer lench,time
      real*8 SigmaProb
      external minfunc,lench,SigmaProb

      open(97,file="my.log",status="replace")
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
CCCCCCCCCCCCCCCC   Initialization  CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      call hmffrd
      call bsffrd
c      call gran_init(time())
      call system_clock(count=time)
      call gran_init(time)
c      call gran_init(200)
      call write_input
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
CCCCCCCCCCCCCCCC   Main program  CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
c      ifirst = 0
CCCCCCCCCCCCCCCC  Event distributions  CCCCCCCCCCCCCCCCCCCCCCCCCCCC
      call init_minuit(npari,nparx)
      call set_params
      arg(1) = 0d0
      call mnexcm(minfunc,'SET PRINTOUT',arg,1,ierr,0)     
      idata = 1
      call minfunc(npari,grad,dchisq,z,iflag,0)            
      write(*,*) "first iteration"
CCCCCCCCCCCCCCCC  chi^2 minimization  CCCCCCCCCC
      if (mode.eq.1) then
         write(*,*) "[Delta-Chi2 analysis]"
         idata = -1
         call mnexcm(minfunc,'MIGRAD',arg,0,ierr,0)
         call mnstat(chisqmin_wrong,fedm,errdef,npari,nparx,istat)
         call write_params
         dchisqmin = chisqmin_wrong
         if (dchisqmin.gt.0d0) then
            sensitivity = SigmaProb(dsqrt(dchisqmin))
         else
            sensitivity = 0d0
         endif     
      endif
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
CCCCCCCCCCCCCCCC   Output Resulted Event Numbers  CCCCCCCCCCCCCCCCCCC
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
CCCCCCCCCCCCCCCC  Event distributions  CCCCCCCCCCCCCCCCCCCCCCCCCCCC
      open(1,file="dchi2.dat",status="replace")
      write(1,*) dchisqmin
      close(1)
      call prepare_events
      call display_results
      if (mode.eq.0) then
c         call write_data_oscsum(1)
         call write_data_oscsum_nuanu(1)
      elseif (mode.eq.1) then
c         call write_chi2event_dat
c         call write_chi2event_th
         call write_data_oscsum(1)
         call write_data_oscsum(2)
         call write_chi2_stat      
         call write_pulls
c         continue
      endif
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      close(97)      
      end
