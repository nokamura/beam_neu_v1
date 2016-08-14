      subroutine set_params
C     ****************************************************
C     By Yoshitaro Takaesu @KIAS JAN 12 2014
C
C     Initialize variables in minfunc
C     ****************************************************
      implicit none
C     GLOBAL VARIABLES
      include 'inc/params.inc'
      include 'inc/main.inc'
C     CONSTANTS
C     ARGUMENTS 
      real*8 z(zdim)
C     LOCAL VARIABLES 
      integer i,j,k,l1,l2
C     EXTERNAL FUNCTIONS
C     ----------
C     BEGIN CODE
C     ----------
      if (nparx.gt.100) write(*,*) "ERROR: nparx is larger than 100"
      if (ieff.eq.0) then
         Seffe(1) = 1d0
         Seffmu(1) = 1d0
         Oeffe(1) = 1d0
         Oeffmu(1) = 1d0
         Keffe(1) = 1d0
         Keffmu(1) = 1d0
      endif
      if (imissid.eq.0) then
         SPe2m(1) = 0d0
         SPm2e(1) = 0d0
         OPe2m(1) = 0d0
         OPm2e(1) = 0d0
         KPe2m(1) = 0d0
         KPm2e(1) = 0d0
      endif

      do i = 1,zdim
         z_dat(i) = 0d0
      enddo
      z_dat(1) = s2sol_2
      z_dat(2) = thatm
      z_dat(3) = s2rct_2
      z_dat(4) = dm21_2
      z_dat(5) = dmatm_2
      z_dat(6) = dCP
      z_dat(7) = fCCQEn(1)
      z_dat(8) = fCCQEa(1)
      z_dat(9) = fS(1)
      z_dat(10) = fSrho(1)
      z_dat(11) = fSnmn(1)
      z_dat(12) = fSnen(1)
      z_dat(13) = fSnma(1)
      z_dat(14) = fSnea(1)
      z_dat(15) = fSamn(1)
      z_dat(16) = fSaen(1)
      z_dat(17) = fSama(1)
      z_dat(18) = fSaea(1)
      z_dat(19) = Seffe(1)
      z_dat(20) = Seffmu(1)
      z_dat(21) = SPe2m(1)
      z_dat(22) = SPm2e(1)
      z_dat(23) = fO(1)
      z_dat(24) = fOrho(1)
      z_dat(25) = fOnmn(1)
      z_dat(26) = fOnen(1)
      z_dat(27) = fOnma(1)
      z_dat(28) = fOnea(1)
      z_dat(29) = fOamn(1)
      z_dat(30) = fOaen(1)
      z_dat(31) = fOama(1)
      z_dat(32) = fOaea(1)
      z_dat(33) = Oeffe(1)
      z_dat(34) = Oeffmu(1)
      z_dat(35) = OPe2m(1)
      z_dat(36) = OPm2e(1)
      z_dat(37) = fK(1)
      z_dat(38) = fKrho(1)
      z_dat(39) = fKnmn(1)
      z_dat(40) = fKnen(1)
      z_dat(41) = fKnma(1)
      z_dat(42) = fKnea(1)
      z_dat(43) = fKamn(1)
      z_dat(44) = fKaen(1)
      z_dat(45) = fKama(1)
      z_dat(46) = fKaea(1)
      z_dat(47) = Keffe(1)
      z_dat(48) = Keffmu(1)
      z_dat(49) = KPe2m(1)
      z_dat(50) = KPm2e(1)
      z_dat(51) = fCCResn(1)
      z_dat(52) = fCCResa(1)
      z_dat(53) = fpi0(1)
      z_dat(54) = fpirs(1)
      z_dat(55) = fpico(1)
      z_dat(56) = frpi0(1)

      error(1) = fs2sol_2(2)
      error(2) = fthatm(2)
      error(3) = fs2rct_2(2)
      error(4) = fdm21_2(2)
      error(5) = fdmatm_2(2)
      error(6) = fdCP(2)
      error(7) = fCCQEn(2)
      error(8) = fCCQEa(2)
      error(9) = fS(2)
      error(10) = fSrho(2)
      error(11) = fSnmn(2)
      error(12) = fSnen(2)
      error(13) = fSnma(2)
      error(14) = fSnea(2)
      error(15) = fSamn(2)
      error(16) = fSaen(2)
      error(17) = fSama(2)
      error(18) = fSaea(2)
      error(19) = Seffe(2)
      error(20) = Seffmu(2)
      error(21) = SPe2m(2)
      error(22) = SPm2e(2)
      error(23) = fO(2)
      error(24) = fOrho(2)
      error(25) = fOnmn(2)
      error(26) = fOnen(2)
      error(27) = fOnma(2)
      error(28) = fOnea(2)
      error(29) = fOamn(2)
      error(30) = fOaen(2)
      error(31) = fOama(2)
      error(32) = fOaea(2)
      error(33) = Oeffe(2)
      error(34) = Oeffmu(2)
      error(35) = OPe2m(2)
      error(36) = OPm2e(2)
      error(37) = fK(2)
      error(38) = fKrho(2)
      error(39) = fKnmn(2)
      error(40) = fKnen(2)
      error(41) = fKnma(2)
      error(42) = fKnea(2)
      error(43) = fKamn(2)
      error(44) = fKaen(2)
      error(45) = fKama(2)
      error(46) = fKaea(2)
      error(47) = Keffe(2)
      error(48) = Keffmu(2)
      error(49) = KPe2m(2)
      error(50) = KPm2e(2)
      error(51) = fCCResn(2)
      error(52) = fCCResa(2)
      error(53) = fpi0(2)
      error(54) = fpirs(2)
      error(55) = fpico(2)
      error(56) = frpi0(2)

c      do i = 1,nparx
      do i = 1,max_npar
         parflag(i) = 1
c         parflag(i) = 0
      enddo
c     Specify which parameters are NOT included in the chi^2 function
      if (fs2sol_2(4).eq.0) then
         parflag(1) = 0        
      endif
      if (fthatm(4).eq.0) then
         parflag(2) = 0        
      endif
      if (fs2rct_2(4).eq.0) then
         parflag(3) = 0        
      endif
      if (fdm21_2(4).eq.0) then
         parflag(4) = 0        
      endif
      if (fdmatm_2(4).eq.0) then
         parflag(5) = 0        
      endif
      if (fdCP(4).eq.0) then
         parflag(6) = 0         ! CP
      endif
      if (ipi0unc.eq.0) then
         parflag(56) = 0
      elseif (ipi0unc.eq.1) then
         parflag(54) = 0
         parflag(55) = 0
      endif

      return
      end
