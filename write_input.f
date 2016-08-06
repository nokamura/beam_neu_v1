      subroutine write_input
C     ****************************************************
C     By Yoshitaro Takaesu @KIAS DEC 26 2013
C
C     Output input parameters in params.card
C     ****************************************************
      implicit none
C     GLOBAL VARIABLES
      include 'inc/params.inc'
C     CONSTANTS
C     ARGUMENTS 
C     LOCAL VARIABLES 
C     EXTERNAL FUNCTIONS
C     ----------
C     BEGIN CODE
C     ----------
      write(*,*) "********** Input Parameters **********"
      write(*,*)
      write(*,*) "[ System Parameters ]"
      write(*,*) "mode =",mode
      write(*,*) "iSK =",iSK
      write(*,*) "  Srho =",Srho
      write(*,*) "  SV =",SV
      write(*,*) "  SL =",SL
      write(*,*) "  SOAB =",SOAB
      write(*,*) "iOki =",iOki
      write(*,*) "  Orho =",Orho
      write(*,*) "  OV =",OV
      write(*,*) "  OL =",OL
      write(*,*) "  OOAB =",OOAB
      write(*,*) "iKr =",iKr
      write(*,*) "  Krho =",Krho
      write(*,*) "  KV =",KV
      write(*,*) "  KL =",KL
      write(*,*) "  KOAB =",KOAB
      write(*,*) "iNnm =",iNnm      
      write(*,*) "iNne =",iNne
      write(*,*) "iNam =",iNam
      write(*,*) "iNae =",iNae      
      write(*,*) "serror =",serror      
      write(*,*) "snmax =",snmax      
      write(*,*) "ihisto =",ihisto
      write(*,*) "ismear =",ismear
      write(*,*) "iCCQE =",iCCQE
      write(*,*) "iCCRes =",iCCRes
      write(*,*) "ipi0dist =",ipi0dist
      write(*,*) "ipi0xsec =",ipi0xsec
      write(*,*) "ihfunc =",ihfunc
      write(*,*) "ieff =",ieff
      write(*,*) "imissid =",imissid
      write(*,*) "ihypo =",ihypo
      write(*,*) "iminevents =",iminevents
      write(*,*) "basic_binsize =",basic_binsize
      write(*,*)
      write(*,*) "[ Experimental Parameters ]"
      write(*,*) "beam =",beam
      write(*,*) "Y =",Y," [years]"
      write(*,*) "r_nu =",r_nu
      write(*,*) "r_anu =",r_anu
      write(*,*) "Emin =",Emin," [GeV]"
      write(*,*) "Emax =",Emax," [GeV]"
      write(*,*) "binsize_factor =",binsize_factor
      write(*,*) "binsize =",basic_binsize*binsize_factor
      write(*,*)
      write(*,*) "[ Physics Parameters ]"
      write(*,*) "iinput =",iinput
      write(*,*) "ithatm =",ithatm
      write(*,*) "ioct =",ioct
      write(*,*) "idmatm2 =",idmatm2
      write(*,*) "MH =",MH
      write(*,'(a11,e12.5)') "sin(2*th_sol)^2 = ",s2sol_2
      if (fs2sol_2(3).eq.0) then
         write(*,'(a24)') "sin(2*th_sol)^2 is fixed"
      endif
      if (fs2sol_2(4).eq.0) then
         write(*,'(a31,a9)') "sin(2*th_sol)^2 is not included",
     &        " in chi^2"
      endif
      if (ithatm.eq.1) then
         write(*,'(a11,e12.5)') "sin(2*th_atm)^2 = ",thatm
         if (fthatm(3).eq.0) then
            write(*,'(a24)') "sin(2*th_atm)^2 is fixed"
         endif
         if (fthatm(4).eq.0) then
            write(*,'(a31,a9)') "sin(2*th_atm)^2 is not included",
     &           " in chi^2"
         endif
      elseif (ithatm.eq.2) then
         write(*,'(a11,e12.5)') "sin(th_atm)^2 = ",thatm
         if (fthatm(3).eq.0) then
            write(*,'(a22)') "sin(th_atm)^2 is fixed"
         endif
         if (fthatm(4).eq.0) then
            write(*,'(a29,a9)') "sin(th_atm)^2 is not included",
     &           " in chi^2"
         endif
      else
         write(*,*) "ERROR:main.f:Invalid ithatm"
         write(*,*) "Stopping program..."
         stop
      endif
      write(*,'(a11,e12.5)') "sin(2*th_rct)^2 = ",s2rct_2
         if (fs2rct_2(3).eq.0) then
            write(*,'(a24)') "sin(2*th_rct)^2 is fixed"
         endif
         if (fs2rct_2(4).eq.0) then
            write(*,'(a31,a9)') "sin(2*th_rct)^2 is not included",
     &           " in chi^2"
         endif
      write(*,'(a11,e12.5,a7)') "dm21^2 = ",dm21_2," [eV^2]"
         if (fdm21_2(3).eq.0) then
            write(*,'(a15)') "dm21^2 is fixed"
         endif
         if (fdm21_2(4).eq.0) then
            write(*,'(a22,a9)') "dm21^2 is not included",
     &           " in chi^2"
         endif
      if (idmatm2.eq.1) then
         write(*,'(a11,e12.5,a7)') "dm31^2 = ",dmatm_2," [eV^2]"
         if (fdmatm_2(3).eq.0) then
            write(*,'(a15)') "dm31^2 is fixed"
         endif
         if (fdmatm_2(4).eq.0) then
            write(*,'(a22,a9)') "dm31^2 is not included",
     &           " in chi^2"
         endif
      elseif (idmatm2.eq.2) then
         write(*,'(a11,e12.5,a7)') "dm32^2 = ",dmatm_2," [eV^2]"
         if (fdmatm_2(3).eq.0) then
            write(*,'(a15)') "dm32^2 is fixed"
         endif
         if (fdmatm_2(4).eq.0) then
            write(*,'(a22,a9)') "dm32^2 is not included",
     &           " in chi^2"
         endif
      endif
      write(*,'(a11,e12.5,a6)') "delta_CP = ",dCP," [rad]"
      if (fdCP(3).eq.0) then
         write(*,'(a14)') "fixed delta_CP"
      endif
      if (fdCP(4).eq.0) then
         write(*,'(a24,a9)') "delta_CP is not included",
     &        " in chi^2"
      endif
      write(*,*) 
      write(*,*) "*************************************"

c$$$      write(*,*) "[ Systematic Errors ]"
c$$$      write(*,'(a11,e12.5)') "fCCQEn = ",fCCQEn(1)
c$$$      write(*,'(a11,e12.5)') "fCCQEa = ",fCCQEa(1)
c$$$      if (iSK.eq.1) then
c$$$         write(*,'(a11,e12.5)') "fS = ",fS(1)
c$$$         write(*,'(a11,e12.5)') "fSrho = ",fSrho(1)
c$$$         write(*,'(a11,e12.5)') "fSnmn = ",fSnmn(1)
c$$$         write(*,'(a11,e12.5)') "fSnen = ",fSnen(1)
c$$$         write(*,'(a11,e12.5)') "fSnma = ",fSnma(1)
c$$$         write(*,'(a11,e12.5)') "fSnea = ",fSnea(1)
c$$$         write(*,'(a11,e12.5)') "fSamn = ",fSamn(1)
c$$$         write(*,'(a11,e12.5)') "fSaen = ",fSaen(1)
c$$$         write(*,'(a11,e12.5)') "fSama = ",fSama(1)
c$$$         write(*,'(a11,e12.5)') "fSaea = ",fSaea(1)
c$$$      endif
c$$$      if(iOki.eq.1) then
c$$$         write(*,'(a11,e12.5)') "fO = ",fO(1)
c$$$         write(*,'(a11,e12.5)') "fOrho = ",fOrho(1)
c$$$         write(*,'(a11,e12.5)') "fOnmn = ",fOnmn(1)
c$$$         write(*,'(a11,e12.5)') "fOnen = ",fOnen(1)
c$$$         write(*,'(a11,e12.5)') "fOnma = ",fOnma(1)
c$$$         write(*,'(a11,e12.5)') "fOnea = ",fOnea(1)
c$$$         write(*,'(a11,e12.5)') "fO = ",fOamn(1)
c$$$         write(*,'(a11,e12.5)') "fO = ",fOaen(1)
c$$$         write(*,'(a11,e12.5)') "fOama = ",fOama(1)
c$$$         write(*,'(a11,e12.5)') "fOaea = ",fOaea(1)
c$$$      endif
c$$$      if(iKr.eq.1) then
c$$$         write(*,'(a11,e12.5)') "fK = ",fK(1)
c$$$         write(*,'(a11,e12.5)') "fKrho = ",fKrho(1)
c$$$         write(*,'(a11,e12.5)') "fKnmn = ",fKnmn(1)
c$$$         write(*,'(a11,e12.5)') "fKnen = ",fKnen(1)
c$$$         write(*,'(a11,e12.5)') "fKnma = ",fKnma(1)
c$$$         write(*,'(a11,e12.5)') "fKnea = ",fKnea(1)
c$$$         write(*,'(a11,e12.5)') "fK = ",fKamn(1)
c$$$         write(*,'(a11,e12.5)') "fK = ",fKaen(1)
c$$$         write(*,'(a11,e12.5)') "fKama = ",fKama(1)
c$$$         write(*,'(a11,e12.5)') "fKaea = ",fKaea(1)
c$$$      endif

      write(*,*) ""
      write(*,*) "( Simpson Integration Parameters )"
      write(*,*) "Accuracy:",serror
      write(*,*) "Max division = 2^",snmax+1
      write(*,*) ""

      if (mode.eq.0) then
         if (beam.eq.1) then
            write(*,*) "beam mode: nu_mu focusing beam"      
         elseif (beam.eq.-1) then
            write(*,*) "beam mode: anti-nu_mu focusing beam"      
         endif
         if (detect.eq.1) then
            write(*,*) "detecting neutrino: nu_e"      
         elseif (detect.eq.-1) then
            write(*,*) "detecting neutrino: anti-nu_e"      
         elseif (detect.eq.2) then
            write(*,*) "detecting neutrino: nu_mu"      
         elseif (detect.eq.-2) then
            write(*,*) "detecting neutrino: anti-nu_mu"      
         elseif (detect.eq.3) then
            write(*,*) "detecting neutrino: nu_tau"      
         elseif (detect.eq.-3) then
            write(*,*) "detecting neutrino: anti-nu_tau"      
         endif
      endif
      write(*,*) "Ev range:",Emin," -",Emax," [GeV]"
      if (mode.eq.1) then
         write(*,*) ""
         write(*,*) "Run in minimization mode"
         write(*,*) ""
      endif
      
      return
      end
