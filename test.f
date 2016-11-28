      program test
      implicitnone
      include 'inc/const.inc'
      integer i,ii
      integer bfactor,nbins
      real*8 min,max,xmin,xmax,binsize,x(0:1000),yy(0:1000)
      real*8 datain(1000),dataout(1000),ans,prob
      external prob
      
      min = 1.4d0
      max = 8.4d0
      bfactor = 3
      xmin = 1.0d0
      xmax = 10d0
      binsize = 0.1d0

c$$$      call bining_x(xmin,xmax,binsize,nbins,x,yy)
c$$$      do i = 1,nbins
c$$$         datain(i) = i*0.1d0
c$$$      enddo
c$$$      call apply_Ereccut_unit(min,max,bfactor,nbins,x,datain,dataout)
c$$$      write(*,*) "datain"
c$$$      write(*,*) (datain(ii),ii=1,nbins)
c$$$      write(*,*) "dataout"
c$$$      write(*,*) (dataout(ii),ii=1,nbins)      

c 210  .1  1040.  0.85  1.  0.085  7.5E-05  0.002545  0.  0. -1
      ans = prob(2,1,0.1d0,1040d0,0.85d0,1.d0,0.085d0,7.5d-5,2.545d-3,0d0,0d0,-1)
      write(*,*) ans

      return
      end
