      integer ndiv,mode,ifixL,ifluc,snmax,ismear,ihisto,ipi0unc,ihypo,iminevents
      common /par_sysi/ ndiv,mode,ifixL,ifluc,snmax,ismear,ihisto,ipi0unc,ihypo,iminevents

      integer iminuit,iinput,ithatm,idmatm2,iCCQE,iCCRes,iNC,ipi0dist,ipi0xsec
      common /par_sysi2/ iminuit,iinput,ithatm,idmatm2,iCCQE,iCCRes,iNC,ipi0dist,ipi0xsec

      integer ihfunc,iNnm,iNne,iNam,iNae,iSK,iOki,iKr,ioct,ieff,imissid
      common /par_sysi3/ ihfunc,iNnm,iNne,iNam,iNae,iSK,iOki,iKr,ioct,ieff,imissid

      real*8 Srho,SV,SL,SOAB,Orho,OV,OL,OOAB,Krho,KV,KL,KOAB
      common /par_sys4/ Srho,SV,SL,SOAB,Orho,OV,OL,OOAB,Krho,KV,KL,KOAB

      integer Sbfactor,Obfactor,Kbfactor
      real*8 Semin,Semax,Smmin,Smmax
      real*8 Oemin,Oemax,Ommin,Ommax
      real*8 Kemin,Kemax,Kmmin,Kmmax
      common /par_sys5/ Semin,Semax,Smmin,Smmax,Oemin,Oemax,Ommin,Ommax
     &	     ,Kemin,Kemax,Kmmin,Kmmax,Sbfactor,Obfactor,Kbfactor

      real*8 effe,effmu,Pe2m,Pm2e
      common /par_sys6/ effe,effmu,Pe2m,Pm2e

      integer binsize_factor
      real*8 L,serror,basic_binsize
      common /par_sysr/ L,serror,basic_binsize,binsize_factor

      integer nu_mode,beam,detect
      real*8 OAB,V,Y,r_nu,r_anu,Emin,Emax
      common /par_exp/ OAB,V,Y,r_nu,r_anu,Emin,Emax
      common /par_exp2/ nu_mode,beam,detect

      integer MH
      real*8 s2sol_2,thatm,s2rct_2,dm21_2,dmatm_2,dCP,rho
      common /par_phy/ s2sol_2,thatm,s2rct_2,dm21_2,dmatm_2,dCP,rho,MH

      real*8 fs2sol_2(4),fthatm(4),fs2rct_2(4),fdm21_2(4),fdmatm_2(4),fdCP(4)
      common /par_fit/ fs2sol_2,fthatm,fs2rct_2,fdm21_2,fdmatm_2,fdCP

      real*8 fCCQEn(3),fCCQEa(3),fCCResn(3),fCCResa(3),fpi0(3),fpirs(3),fpico(3)
      real*8 frpi0(3)
      common /err_cor/ fCCQEn,fCCQEa,fCCResn,fCCResa,fpi0,fpirs,fpico,frpi0

      real*8 fS(3),fSrho(3),fSnmn(3),fSnen(3),fSnma(3),fSnea(3)
      common /err_sk/ fS,fSrho,fSnmn,fSnen,fSnma,fSnea

      real*8 fSamn(3),fSaen(3),fSama(3),fSaea(3)
      common /err_sk2/ fSamn,fSaen,fSama,fSaea

      real*8 Seffe(3),Seffmu(3),SPe2m(3),SPm2e(3)  
      common /err_sk3/ Seffe,Seffmu,SPe2m,SPm2e  

      real*8 fO(3),fOrho(3),fOnmn(3),fOnen(3),fOnma(3),fOnea(3)
      common /err_ok/ fO,fOrho,fOnmn,fOnen,fOnma,fOnea

      real*8 fOamn(3),fOaen(3),fOama(3),fOaea(3)
      common /err_ok2/ fOamn,fOaen,fOama,fOaea

      real*8 Oeffe(3),Oeffmu(3),OPe2m(3),OPm2e(3)  
      common /err_ok3/ Oeffe,Oeffmu,OPe2m,OPm2e  

      real*8 fK(3),fKrho(3),fKnmn(3),fKnen(3),fKnma(3),fKnea(3)
      common /err_kre/ fK,fKrho,fKnmn,fKnen,fKnma,fKnea

      real*8 fKamn(3),fKaen(3),fKama(3),fKaea(3)
      common /err_kre/ fKamn,fKaen,fKama,fKaea

      real*8 Keffe(3),Keffmu(3),KPe2m(3),KPm2e(3)  
      common /err_kr3/ Keffe,Keffmu,KPe2m,KPm2e  

      integer icc,MHH,xsec_mode
      common /other/ icc,MHH,xsec_mode

      integer maxnbin,nout,zdim,eventdim,maxnparam
      integer ndetect,nmode,nmode_rec,nmode_like,imaxint
      parameter (maxnbin=1000,zdim=200,eventdim=4,nout=6)
      parameter (maxnparam=100,ndetect=3)
      parameter (nmode=4) ! 1:e 2:mu 3:pi0 4:gamma
      parameter (nmode_rec=6) ! 1:e>e 2:mu>mu 3:pi0 4:gamma 5:mu>e 6:e>mu
      parameter (nmode_like=2) ! 1:e-like 2:mu-like
      parameter (imaxint=4)