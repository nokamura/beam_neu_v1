      subroutine init_minuit(npari,nparx)
C     ****************************************************
C     By Yoshitaro Takaesu @KIAS DEC 26 2013
C
C     Initialize and set parameters for MINUIT
C     ****************************************************
      implicit none
C     GLOBAL VARIABLES
      include 'inc/params.inc'
      real*8 chisqmin,fedm,errdef,npari,nparx,istat
c      common /minuit/ dchisqmin,fedm,chisqmin_true,errdef,istat
C     CONSTANTS
C     ARGUMENTS 
C     LOCAL VARIABLES 
      integer iflag,ierr
C     EXTERNAL FUNCTIONS
      real*8 minfunc
      external minfunc
C     ----------
C     BEGIN CODE
C     ----------
      call mninit(5,6,7)
      call mnparm(1,'s2sol_2',fs2sol_2(1),fs2sol_2(2),0d0,0d0,ierr)
      if (ithatm.eq.1) then
         call mnparm(2,'s2atm_2',fthatm(1),fthatm(2),0d0,0d0,ierr)
      elseif (ithatm.eq.2) then
         call mnparm(2,'satm_2',fthatm(1),fthatm(2),0d0,0d0,ierr)
      else
         write(*,*) "ERROR:main.f:Invalid ithatm"
         write(*,*) "Stopping program..."
         stop
      endif
      call mnparm(3,'s2rct_2',fs2rct_2(1),fs2rct_2(2),0d0,0d0,ierr)
      call mnparm(4,'dm21_2',fdm21_2(1),fdm21_2(2),0d0,0d0,ierr)
      if (idmatm2.eq.1) then
         call mnparm(5,'dm31_2',fdmatm_2(1),fdmatm_2(2),0d0,0d0,ierr)
      elseif (idmatm2.eq.2) then
         call mnparm(5,'dm32_2',fdmatm_2(1),fdmatm_2(2),0d0,0d0,ierr)
      endif
      call mnparm(6,'dCP',fdCP(1),fdCP(2),0d0,0d0,ierr)
      call mnparm(7,'fCCQEn',fCCQEn(1),fCCQEn(2),0d0,0d0,ierr)
      call mnparm(8,'fCCQEa',fCCQEa(1),fCCQEa(2),0d0,0d0,ierr)
      call mnparm(9,'fS',fS(1),fS(2),0d0,0d0,ierr)
      call mnparm(10,'fSrho',fSrho(1),fSrho(2),0d0,0d0,ierr)
      call mnparm(11,'fSnmn',fSnmn(1),fSnmn(2),0d0,0d0,ierr)
      call mnparm(12,'fSnen',fSnen(1),fSnen(2),0d0,0d0,ierr)
      call mnparm(13,'fSnma',fSnma(1),fSnma(2),0d0,0d0,ierr)
      call mnparm(14,'fSnea',fSnea(1),fSnea(2),0d0,0d0,ierr)
      call mnparm(15,'fSamn',fSamn(1),fSamn(2),0d0,0d0,ierr)
      call mnparm(16,'fSaen',fSaen(1),fSaen(2),0d0,0d0,ierr)
      call mnparm(17,'fSama',fSama(1),fSama(2),0d0,0d0,ierr)
      call mnparm(18,'fSaea',fSaea(1),fSaea(2),0d0,0d0,ierr)
      if (ieff.eq.0) then
         call mnparm(19,'Seffe',1d0,Seffe(2),0d0,0d0,ierr)
         call mnparm(20,'Seffmu',1d0,Seffmu(2),0d0,0d0,ierr)
      elseif (ieff.eq.1) then
         call mnparm(19,'Seffe',Seffe(1),Seffe(2),0d0,0d0,ierr)
         call mnparm(20,'Seffmu',Seffmu(1),Seffmu(2),0d0,0d0,ierr)
      endif
      if (imissid.eq.0) then
         call mnparm(21,'SPe2m',0d0,SPe2m(2),0d0,0d0,ierr)
         call mnparm(22,'SPm2e',0d0,SPm2e(2),0d0,0d0,ierr)
      elseif (imissid.eq.1) then
         call mnparm(21,'SPe2m',SPe2m(1),SPe2m(2),0d0,0d0,ierr)
         call mnparm(22,'SPm2e',SPm2e(1),SPm2e(2),0d0,0d0,ierr)
      endif
      call mnparm(23,'fO',fO(1),fO(2),0d0,0d0,ierr)
      call mnparm(24,'fOrho',fOrho(1),fOrho(2),0d0,0d0,ierr)
      call mnparm(25,'fOnmn',fOnmn(1),fOnmn(2),0d0,0d0,ierr)
      call mnparm(26,'fOnen',fOnen(1),fOnen(2),0d0,0d0,ierr)
      call mnparm(27,'fOnma',fOnma(1),fOnma(2),0d0,0d0,ierr)
      call mnparm(28,'fOnea',fOnea(1),fOnea(2),0d0,0d0,ierr)
      call mnparm(29,'fOamn',fOamn(1),fOamn(2),0d0,0d0,ierr)
      call mnparm(30,'fOaen',fOaen(1),fOaen(2),0d0,0d0,ierr)
      call mnparm(31,'fOama',fOama(1),fOama(2),0d0,0d0,ierr)
      call mnparm(32,'fOaea',fOaea(1),fOaea(2),0d0,0d0,ierr)
      if (ieff.eq.0) then
         call mnparm(33,'Oeffe',1d0,Oeffe(2),0d0,0d0,ierr)
         call mnparm(34,'Oeffmu',1d0,Oeffmu(2),0d0,0d0,ierr)
      elseif (ieff.eq.1) then
         call mnparm(33,'Oeffe',Oeffe(1),Oeffe(2),0d0,0d0,ierr)
         call mnparm(34,'Oeffmu',Oeffmu(1),Oeffmu(2),0d0,0d0,ierr)
      endif
      if (imissid.eq.0) then
         call mnparm(35,'OPe2m',0d0,OPe2m(2),0d0,0d0,ierr)
         call mnparm(36,'OPm2e',0d0,OPm2e(2),0d0,0d0,ierr)
      elseif (imissid.eq.1) then
         call mnparm(35,'OPe2m',OPe2m(1),OPe2m(2),0d0,0d0,ierr)
         call mnparm(36,'OPm2e',OPm2e(1),OPm2e(2),0d0,0d0,ierr)
      endif
      call mnparm(37,'fK',fK(1),fK(2),0d0,0d0,ierr)
      call mnparm(38,'fKrho',fKrho(1),fKrho(2),0d0,0d0,ierr)
      call mnparm(39,'fKnmn',fKnmn(1),fKnmn(2),0d0,0d0,ierr)
      call mnparm(40,'fKnen',fKnen(1),fKnen(2),0d0,0d0,ierr)
      call mnparm(41,'fKnma',fKnma(1),fKnma(2),0d0,0d0,ierr)
      call mnparm(42,'fKnea',fKnea(1),fKnea(2),0d0,0d0,ierr)
      call mnparm(43,'fKamn',fKamn(1),fKamn(2),0d0,0d0,ierr)
      call mnparm(44,'fKaen',fKaen(1),fKaen(2),0d0,0d0,ierr)
      call mnparm(45,'fKama',fKama(1),fKama(2),0d0,0d0,ierr)
      call mnparm(46,'fKaea',fKaea(1),fKaea(2),0d0,0d0,ierr)
      if (ieff.eq.0) then
         call mnparm(47,'Keffe',1d0,Keffe(2),0d0,0d0,ierr)
         call mnparm(48,'Keffmu',1d0,Keffmu(2),0d0,0d0,ierr)
      elseif (ieff.eq.1) then
         call mnparm(47,'Keffe',Keffe(1),Keffe(2),0d0,0d0,ierr)
         call mnparm(48,'Keffmu',Keffmu(1),Keffmu(2),0d0,0d0,ierr)
      endif
      if (imissid.eq.0) then
         call mnparm(49,'KPe2m',0d0,KPe2m(2),0d0,0d0,ierr)
         call mnparm(50,'KPm2e',0d0,KPm2e(2),0d0,0d0,ierr)
      elseif (imissid.eq.1) then
         call mnparm(49,'KPe2m',KPe2m(1),KPe2m(2),0d0,0d0,ierr)
         call mnparm(50,'KPm2e',KPm2e(1),KPm2e(2),0d0,0d0,ierr)
      endif
      call mnparm(51,'fCCResn',fCCResn(1),fCCResn(2),0d0,0d0,ierr)
      call mnparm(52,'fCCResa',fCCResa(1),fCCResa(2),0d0,0d0,ierr)
      call mnparm(53,'fpi0',fpi0(1),fpi0(2),0d0,0d0,ierr)
      call mnparm(54,'fpirs',fpirs(1),fpirs(2),0d0,0d0,ierr)
      call mnparm(55,'fpico',fpico(1),fpico(2),0d0,0d0,ierr)
      call mnparm(56,'frpi0',frpi0(1),frpi0(2),0d0,0d0,ierr)
c     call mncomd(minfunc,'FIX 2',iflag,0)
      if (fs2sol_2(3).eq.0) then
         call mncomd(minfunc,'FIX 1',iflag,0)
      endif
      if (fthatm(3).eq.0) then
         call mncomd(minfunc,'FIX 2',iflag,0)
      endif
      if (fs2rct_2(3).eq.0) then
         call mncomd(minfunc,'FIX 3',iflag,0)
      endif
      if (fdm21_2(3).eq.0) then
         call mncomd(minfunc,'FIX 4',iflag,0)
      endif
      if (fdmatm_2(3).eq.0) then
         call mncomd(minfunc,'FIX 5',iflag,0)
      endif
      if (fdCP(3).eq.0) then
         call mncomd(minfunc,'FIX 6',iflag,0)
      endif
      if (iSK.eq.0) then
         call mncomd(minfunc,'FIX 9',iflag,0)
         call mncomd(minfunc,'FIX 10',iflag,0)
         call mncomd(minfunc,'FIX 11',iflag,0)
         call mncomd(minfunc,'FIX 12',iflag,0)
         call mncomd(minfunc,'FIX 13',iflag,0)
         call mncomd(minfunc,'FIX 14',iflag,0)
         call mncomd(minfunc,'FIX 15',iflag,0)
         call mncomd(minfunc,'FIX 16',iflag,0)
         call mncomd(minfunc,'FIX 17',iflag,0)
         call mncomd(minfunc,'FIX 18',iflag,0)
         call mncomd(minfunc,'FIX 19',iflag,0)
         call mncomd(minfunc,'FIX 20',iflag,0)
         call mncomd(minfunc,'FIX 21',iflag,0)
         call mncomd(minfunc,'FIX 22',iflag,0)
      endif
      if (iOki.eq.0) then
         call mncomd(minfunc,'FIX 23',iflag,0)
         call mncomd(minfunc,'FIX 24',iflag,0)
         call mncomd(minfunc,'FIX 25',iflag,0)
         call mncomd(minfunc,'FIX 26',iflag,0)
         call mncomd(minfunc,'FIX 27',iflag,0)
         call mncomd(minfunc,'FIX 28',iflag,0)
         call mncomd(minfunc,'FIX 29',iflag,0)
         call mncomd(minfunc,'FIX 30',iflag,0)
         call mncomd(minfunc,'FIX 31',iflag,0)
         call mncomd(minfunc,'FIX 32',iflag,0)
         call mncomd(minfunc,'FIX 33',iflag,0)
         call mncomd(minfunc,'FIX 34',iflag,0)
         call mncomd(minfunc,'FIX 35',iflag,0)
         call mncomd(minfunc,'FIX 36',iflag,0)
      endif
      if (iKr.eq.0) then
         call mncomd(minfunc,'FIX 37',iflag,0)
         call mncomd(minfunc,'FIX 38',iflag,0)
         call mncomd(minfunc,'FIX 39',iflag,0)
         call mncomd(minfunc,'FIX 40',iflag,0)
         call mncomd(minfunc,'FIX 41',iflag,0)
         call mncomd(minfunc,'FIX 42',iflag,0)
         call mncomd(minfunc,'FIX 43',iflag,0)
         call mncomd(minfunc,'FIX 44',iflag,0)
         call mncomd(minfunc,'FIX 45',iflag,0)
         call mncomd(minfunc,'FIX 46',iflag,0)
         call mncomd(minfunc,'FIX 47',iflag,0)
         call mncomd(minfunc,'FIX 48',iflag,0)
         call mncomd(minfunc,'FIX 49',iflag,0)
         call mncomd(minfunc,'FIX 50',iflag,0)
      endif
      if (ieff.eq.0) then
         call mncomd(minfunc,'FIX 19',iflag,0)
         call mncomd(minfunc,'FIX 20',iflag,0)
         call mncomd(minfunc,'FIX 33',iflag,0)
         call mncomd(minfunc,'FIX 34',iflag,0)
         call mncomd(minfunc,'FIX 47',iflag,0)
         call mncomd(minfunc,'FIX 48',iflag,0)
      endif
      if (imissid.eq.0) then
         call mncomd(minfunc,'FIX 21',iflag,0)
         call mncomd(minfunc,'FIX 22',iflag,0)
         call mncomd(minfunc,'FIX 35',iflag,0)
         call mncomd(minfunc,'FIX 36',iflag,0)
         call mncomd(minfunc,'FIX 49',iflag,0)
         call mncomd(minfunc,'FIX 50',iflag,0)
      endif
      if (ipi0unc.eq.0) then
         call mncomd(minfunc,'FIX 56',iflag,0)
      elseif (ipi0unc.eq.1) then
         call mncomd(minfunc,'FIX 54',iflag,0)
         call mncomd(minfunc,'FIX 55',iflag,0)
      elseif (ipi0unc.eq.2) then
         call mncomd(minfunc,'FIX 54',iflag,0)
         call mncomd(minfunc,'FIX 55',iflag,0)
         call mncomd(minfunc,'FIX 56',iflag,0)
      endif
      if (fpi0(2).eq.0) then
         call mncomd(minfunc,'FIX 53',iflag,0)
         fpi0(2) = 1
      endif
      if (fpirs(2).eq.0) then
         call mncomd(minfunc,'FIX 54',iflag,0)
         fpirs(2) = 1
      endif
      if (fpico(2).eq.0) then
         call mncomd(minfunc,'FIX 55',iflag,0)
         fpico(2) = 1
      endif

c      call mncomd(minfunc,'FIX 55',iflag,0)
c     For re-evaluation paper
      call mncomd(minfunc,'FIX 21',iflag,0)  ! SPe2m
      call mncomd(minfunc,'FIX 35',iflag,0)  ! OPe2m
      call mncomd(minfunc,'FIX 49',iflag,0)  ! KPe2m
c
c$$$      call mncomd(minfunc,'FIX 39',iflag,0) ! Kfnmn
c$$$      call mncomd(minfunc,'FIX 40',iflag,0) ! Kfnen
c$$$      call mncomd(minfunc,'FIX 41',iflag,0) ! Kfnma
c$$$      call mncomd(minfunc,'FIX 42',iflag,0) ! Kfnea
c$$$      call mncomd(minfunc,'FIX 43',iflag,0) ! Kfamn
c$$$      call mncomd(minfunc,'FIX 44',iflag,0) ! Kfaen
c$$$      call mncomd(minfunc,'FIX 45',iflag,0) ! Kfama
c$$$      call mncomd(minfunc,'FIX 46',iflag,0) ! Kfaea
c
c      call mncomd(minfunc,'FIX 33',iflag,0)  ! Oeffe
c      call mncomd(minfunc,'FIX 34',iflag,0)  ! Oeffm
c      call mncomd(minfunc,'FIX 35',iflag,0)  ! OPm2e
c      call mncomd(minfunc,'FIX 36',iflag,0)  ! OPe2m
c      call mncomd(minfunc,'FIX 47',iflag,0)  ! Keffe
c      call mncomd(minfunc,'FIX 48',iflag,0)  ! Keffm
c      call mncomd(minfunc,'FIX 49',iflag,0)  ! KPm2e
c      call mncomd(minfunc,'FIX 50',iflag,0)  ! KPe2m
c
c      call mncomd(minfunc,'FIX 20',iflag,0)  ! Seffm
c
c      call mncomd(minfunc,'FIX 19',iflag,0)  ! Seffe

      call mnstat(chisqmin,fedm,errdef,npari,nparx,istat)      

      return
      end
