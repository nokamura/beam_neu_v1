      subroutine write_minuit_results
C     ****************************************************
C     By Yoshitaro Takaesu @KIAS DEC 26 2013
C
C     Output MINUIT results to files or the std_output
C     ****************************************************
      implicit none
C     GLOBAL VARIABLES
      include 'inc/params.inc'
      include 'inc/main.inc'
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
      open(21,file='dchi2.dat',status='replace')
      open(22,file='params.dat',status='replace')
      do i = 1,nparx
         call mnpout(i,name(i),pval(i),perr(i),plo(i),phi(i)
     &        ,ierr)
      enddo
      write(21,*) dchisqmin
      write(22,'(3e10.3,152e13.5)') L,dchisqmin,fedm
     &     ,pval(1),perr(1),fs2sol_2(2),(pval(1)-fs2sol_2(1))/fs2sol_2(2)
     &     ,pval(2),perr(2),fthatm(2),(pval(2)-fthatm(1))/fthatm(2)
     &     ,pval(3),perr(3),fs2rct_2(2),(pval(3)-fs2rct_2(1))/fs2rct_2(2)
     &     ,pval(4),perr(4),fdm21_2(2),(pval(4)-fdm21_2(1))/fdm21_2(2)
     &     ,pval(5),perr(5),fdmatm_2(2),(pval(5)-fdmatm_2(1))/fdmatm_2(2)
     &     ,pval(6),perr(6),fdCP(2),(pval(6)-fdCP(1))/fdCP(2)
     &     ,pval(7),perr(7),fCCQEn(2),(pval(7)-fCCQEn(1))/fCCQEn(2)
     &     ,pval(8),perr(8),fCCQEa(2),(pval(8)-fCCQEa(1))/fCCQEa(2)
     &     ,pval(9),perr(9),fS(2),(pval(9)-fS(1))/fS(2)
     &     ,pval(10),perr(10),fSrho(2),(pval(10)-fSrho(1))/fSrho(2)
     &     ,pval(11),perr(11),fSnmn(2),(pval(11)-fSnmn(1))/fSnmn(2)
     &     ,pval(12),perr(12),fSnen(2),(pval(12)-fSnen(1))/fSnen(2)
     &     ,pval(13),perr(13),fSnma(2),(pval(13)-fSnma(1))/fSnma(2)
     &     ,pval(14),perr(14),fSnea(2),(pval(14)-fSnea(1))/fSnea(2)
     &     ,pval(15),perr(15),fSamn(2),(pval(15)-fSamn(1))/fSamn(2)
     &     ,pval(16),perr(16),fSaen(2),(pval(16)-fSaen(1))/fSaen(2)
     &     ,pval(17),perr(17),fSama(2),(pval(17)-fSama(1))/fSama(2)
     &     ,pval(18),perr(18),fSaea(2),(pval(18)-fSaea(1))/fSaea(2)
     &     ,pval(19),perr(19),fO(2),(pval(19)-fO(1))/fO(2)
     &     ,pval(20),perr(20),fOrho(2),(pval(20)-fOrho(1))/fOrho(2)
     &     ,pval(21),perr(21),fOnmn(2),(pval(21)-fOnmn(1))/fOnmn(2)
     &     ,pval(22),perr(22),fOnen(2),(pval(22)-fOnen(1))/fOnen(2)
     &     ,pval(23),perr(23),fOnma(2),(pval(23)-fOnma(1))/fOnma(2)
     &     ,pval(24),perr(24),fOnea(2),(pval(24)-fOnea(1))/fOnea(2)
     &     ,pval(25),perr(25),fOamn(2),(pval(25)-fOamn(1))/fOamn(2)
     &     ,pval(26),perr(26),fOaen(2),(pval(26)-fOaen(1))/fOaen(2)
     &     ,pval(27),perr(27),fOama(2),(pval(27)-fOama(1))/fOama(2)
     &     ,pval(28),perr(28),fOaea(2),(pval(28)-fOaea(1))/fOaea(2)
     &     ,pval(29),perr(29),fK(2),(pval(29)-fK(1))/fK(2)
     &     ,pval(30),perr(30),fKrho(2),(pval(30)-fKrho(1))/fKrho(2)
     &     ,pval(31),perr(31),fKnmn(2),(pval(31)-fKnmn(1))/fKnmn(2)
     &     ,pval(32),perr(32),fKnen(2),(pval(32)-fKnen(1))/fKnen(2)
     &     ,pval(33),perr(33),fKnma(2),(pval(33)-fKnma(1))/fKnma(2)
     &     ,pval(34),perr(34),fKnea(2),(pval(34)-fKnea(1))/fKnea(2)
     &     ,pval(35),perr(35),fKamn(2),(pval(35)-fKamn(1))/fKamn(2)
     &     ,pval(36),perr(36),fKaen(2),(pval(36)-fKaen(1))/fKaen(2)
     &     ,pval(37),perr(37),fKama(2),(pval(37)-fKama(1))/fKama(2)
     &     ,pval(38),perr(38),fKaea(2),(pval(38)-fKaea(1))/fKaea(2)
      
      write(*,'(4x,a14,e12.5,a4,e9.2)') "Delta-Chi2  = "
     &     ,dchisqmin," +- ",fedm
      write(*,'(4x,a14,e12.5,a4,e9.2)') "(sin2*sol)^2 = "
     &     ,pval(1)," +- ",perr(1)
      if (ithatm.eq.1) then
         write(*,'(4x,a14,e12.5,a4,e9.2)') "(sin2*atm)^2 = "
     &        ,pval(2)," +- ",perr(2)
      elseif (ithatm.eq.2) then
         write(*,'(4x,a14,e12.5,a4,e9.2)') "(sinatm)^2 = "
     &        ,pval(2)," +- ",perr(2)
      else
         write(*,*) "ERROR:main.f:Invalid ithatm"
         write(*,*) "Stopping program..."
         stop
      endif
      write(*,'(4x,a14,e12.5,a4,e9.2)') "(sin2*rct)^2 = "
     &     ,pval(3)," +- ",perr(3)
      write(*,'(4x,a14,e12.5,a4,e9.2)') "dm21^2      = "
     &     ,pval(4)," +- ",perr(4)
      if (idmatm2.eq.1) then
         write(*,'(4x,a14,e12.5,a4,e9.2)') "dm31^2      = "
     &        ,pval(5)," +- ",perr(5)
      elseif (idmatm2.eq.2) then
         write(*,'(4x,a14,e12.5,a4,e9.2)') "dm32^2      = "
     &        ,pval(5)," +- ",perr(5)
      endif
      write(*,'(4x,a14,e12.5,a4,e9.2)') "dCP      = "
     &     ,pval(6)," +- ",perr(6)
      write(*,'(4x,a14,e12.5,a4,e9.2)') "fCCQEn = "
     &     ,pval(7)," +- ",perr(7)
      write(*,'(4x,a14,e12.5,a4,e9.2)') "fCCQEa = "
     &     ,pval(8)," +- ",perr(8)
      if (iSK.eq.1) then
         write(*,'(4x,a14,e12.5,a4,e9.2)') "fS = "
     &        ,pval(9)," +- ",perr(9)
         write(*,'(4x,a14,e12.5,a4,e9.2)') "fSrho = "
     &        ,pval(10)," +- ",perr(10)
         write(*,'(4x,a14,e12.5,a4,e9.2)') "fSnmn = "
     &        ,pval(11)," +- ",perr(11)
         write(*,'(4x,a14,e12.5,a4,e9.2)') "fSnen = "
     &        ,pval(12)," +- ",perr(12)
         write(*,'(4x,a14,e12.5,a4,e9.2)') "fSnma = "
     &        ,pval(13)," +- ",perr(13)
         write(*,'(4x,a14,e12.5,a4,e9.2)') "fSnea = "
     &        ,pval(14)," +- ",perr(14)
         write(*,'(4x,a14,e12.5,a4,e9.2)') "fSamn = "
     &        ,pval(15)," +- ",perr(15)
         write(*,'(4x,a14,e12.5,a4,e9.2)') "fSaen = "
     &        ,pval(16)," +- ",perr(16)
         write(*,'(4x,a14,e12.5,a4,e9.2)') "fSama = "
     &        ,pval(17)," +- ",perr(17)
         write(*,'(4x,a14,e12.5,a4,e9.2)') "fSaea = "
     &        ,pval(18)," +- ",perr(18)
      endif
      if (iOki.eq.1) then
         write(*,'(4x,a14,e12.5,a4,e9.2)') "fO = "
     &        ,pval(19)," +- ",perr(19)
         write(*,'(4x,a14,e12.5,a4,e9.2)') "fOrho = "
     &        ,pval(20)," +- ",perr(20)
         write(*,'(4x,a14,e12.5,a4,e9.2)') "fOnmn = "
     &        ,pval(21)," +- ",perr(21)
         write(*,'(4x,a14,e12.5,a4,e9.2)') "fOnen = "
     &        ,pval(22)," +- ",perr(22)
         write(*,'(4x,a14,e12.5,a4,e9.2)') "fOnma = "
     &        ,pval(23)," +- ",perr(23)
         write(*,'(4x,a14,e12.5,a4,e9.2)') "fOnea = "
     &        ,pval(24)," +- ",perr(24)
         write(*,'(4x,a14,e12.5,a4,e9.2)') "fOamn = "
     &        ,pval(25)," +- ",perr(25)
         write(*,'(4x,a14,e12.5,a4,e9.2)') "fOaen = "
     &        ,pval(26)," +- ",perr(26)
         write(*,'(4x,a14,e12.5,a4,e9.2)') "fOama = "
     &        ,pval(27)," +- ",perr(27)
         write(*,'(4x,a14,e12.5,a4,e9.2)') "fOaea = "
     &        ,pval(28)," +- ",perr(28)
      endif
      if (iKr.eq.1) then
         write(*,'(4x,a14,e12.5,a4,e9.2)') "fK = "
     &        ,pval(29)," +- ",perr(29)
         write(*,'(4x,a14,e12.5,a4,e9.2)') "fKrho = "
     &        ,pval(30)," +- ",perr(30)
         write(*,'(4x,a14,e12.5,a4,e9.2)') "fKnmn = "
     &        ,pval(31)," +- ",perr(31)
         write(*,'(4x,a14,e12.5,a4,e9.2)') "fKnen = "
     &        ,pval(32)," +- ",perr(32)
         write(*,'(4x,a14,e12.5,a4,e9.2)') "fKnma = "
     &        ,pval(33)," +- ",perr(33)
         write(*,'(4x,a14,e12.5,a4,e9.2)') "fKnea = "
     &        ,pval(34)," +- ",perr(34)
         write(*,'(4x,a14,e12.5,a4,e9.2)') "fKamn = "
     &        ,pval(35)," +- ",perr(35)
         write(*,'(4x,a14,e12.5,a4,e9.2)') "fKaen = "
     &        ,pval(36)," +- ",perr(36)
         write(*,'(4x,a14,e12.5,a4,e9.2)') "fKama = "
     &        ,pval(37)," +- ",perr(37)
         write(*,'(4x,a14,e12.5,a4,e9.2)') "fKaea = "
     &        ,pval(38)," +- ",perr(38)
      endif
      
      call mncomd(minfunc,'SET OUTPUTFILE 6',iflag,0)
      call mncomd(minfunc,'SHOW COVARIANCE',iflag,0)
      
      close(21)
      close(22)
      
      return
      end
