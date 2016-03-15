      program test
      implicitnone
      include 'inc/const.inc'
      integer i,ii
      integer bfactor,nbins
      real*8 min,max,xmin,xmax,binsize,x(0:1000),yy(0:1000)
      real*8 datain(1000),dataout(1000)
      
      min = 1.4d0
      max = 8.4d0
      bfactor = 3
      xmin = 1.0d0
      xmax = 10d0
      binsize = 0.1d0

      call bining_x(xmin,xmax,binsize,nbins,x,yy)
      do i = 1,nbins
         datain(i) = i*0.1d0
      enddo
      call apply_Ereccut_unit(min,max,bfactor,nbins,x,datain,dataout)
      write(*,*) "datain"
      write(*,*) (datain(ii),ii=1,nbins)
      write(*,*) "dataout"
      write(*,*) (dataout(ii),ii=1,nbins)      

      return
      end
