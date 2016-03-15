      integer ifirst
      common /first/ ifirst

      integer outnbins(nmode_like,ndetect)
      real*8 outx(0:maxnbin,nmode_like,ndetect)
      common /bins/ outx,outnbins

      integer evnbins(eventdim,ndetect)
      real*8 evx(0:maxnbin,eventdim,ndetect)
      common /ev/ evx,evnbins

      real*8 total_event_dat(maxnbin,eventdim,ndetect,-1:1)
      real*8 nevent_dat(eventdim,ndetect,-1:1)
      common /event_dat/ total_event_dat,nevent_dat

      real*8 total_event_th(maxnbin,eventdim,ndetect,-1:1)
      real*8 nevent_th(eventdim,ndetect,-1:1)
      common /event_th/ total_event_th,nevent_th

      integer nbins
      real*8 x(0:maxnbin),yy(0:maxnbin),binsize
      common /bining/ x,yy,binsize,nbins

      real*8 hErec_dat(maxnbin,imaxint,nmode_rec,-3:3,-3:3,ndetect,-1:1)
      real*8 hErec_th(maxnbin,imaxint,nmode_rec,-3:3,-3:3,ndetect,-1:1)
      real*8 hErec_raw_dat(maxnbin,imaxint,nmode,-3:3,-3:3,ndetect,-1:1)
      real*8 hErec_raw_th(maxnbin,imaxint,nmode,-3:3,-3:3,ndetect,-1:1)
      real*8 hErec_like_dat(maxnbin,2,ndetect,-1:1)
      real*8 hErec_like_th(maxnbin,2,ndetect,-1:1)
      real*8 hErec_chi2_dat(maxnbin,2,ndetect,-1:1)
      real*8 hErec_chi2_th(maxnbin,2,ndetect,-1:1)
      real*8 simple_dist(maxnbin,imaxint)
      common /hErec/ hErec_dat,hErec_th,hErec_raw_dat,hErec_raw_th
     &	     ,hErec_like_dat,hErec_like_th,hErec_chi2_dat
     &	     ,hErec_chi2_th,simple_dist

      integer idata
      real*8 ymin,ymax
      common /par/ ymin,ymax,idata

      real*8 hErec_e(maxnbin,ndetect,-1:1),nevent_e(ndetect,-1:1)
      real*8 hErec_mu(maxnbin,ndetect,-1:1),nevent_mu(ndetect,-1:1)
      real*8 hErec_pi0(maxnbin,ndetect,-1:1),nevent_pi0(ndetect,-1:1)
      real*8 hErec_gam(maxnbin,ndetect,-1:1),nevent_gam(ndetect,-1:1)
      real*8 hErec_e2e(maxnbin,ndetect,-1:1),nevent_e2e(ndetect,-1:1)
      real*8 hErec_e2mu(maxnbin,ndetect,-1:1),nevent_e2mu(ndetect,-1:1)
      real*8 hErec_mu2mu(maxnbin,ndetect,-1:1),nevent_mu2mu(ndetect,-1:1)
      real*8 hErec_mu2e(maxnbin,ndetect,-1:1),nevent_mu2e(ndetect,-1:1)
      real*8 nevent_elike(ndetect,-1:1),nevent_mulike(ndetect,-1:1)
      common /display/ hErec_e,hErec_mu,hErec_pi0,hErec_gam,hErec_e2e
     &	     ,hErec_e2mu,hErec_mu2mu,hErec_mu2e
     &	     ,nevent_e,nevent_mu,nevent_pi0,nevent_gam,nevent_e2e,nevent_e2mu
     &	     ,nevent_mu2mu,nevent_mu2e,nevent_elike,nevent_mulike 

      integer istat,npari,nparx
      real*8 fedm,errdef
      common /minuit/ fedm,errdef,npari,nparx,istat

      real*8 dchisqmin,sensitivity
      common /sensitivity/ dchisqmin,sensitivity

      integer nbins_chi2(nmode_like,ndetect)
      real*8 x_chi2(0:maxnbin,nmode_like,ndetect)
      common /cchi2/ nbins_chi2,x_chi2

      integer max_npar
      parameter (max_npar=100)	      
      integer parflag(max_npar)
      real*8 z_dat(zdim),error(max_npar)
      common /zdat/ z_dat,error,parflag

      real*8 chi2_stat_bin(maxnbin,nmode_like,ndetect,-1:1),pull_fact(max_npar)
      common /chi2_stat_bin/ chi2_stat_bin,pull_fact