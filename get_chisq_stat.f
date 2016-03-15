      subroutine get_chisq_stat
C     ****************************************************
C     By Yoshitaro Takaesu @KIAS Feb 17 2014
C     ****************************************************
      implicit none
C     GLOBAL VARIABLES
      include 'inc/params.inc'
      include 'inc/minfunc.inc'
C     CONSTANTS
C     ARGUMENTS 
      integer i,j,k
      integer ilike
C     LOCAL VARIABLES 
C     EXTERNAL FUNCTIONS
C     ----------
C     BEGIN CODE
C     ----------
      do k = -1,1,2
         do j = 1,ndetect
            do i = 1,nmode_like
               chisq_stat(i,j,k) = 0d0
            enddo
         enddo
      enddo
      if (iSK.eq.1) then
         iD = 1
         if (iNne.eq.1) then
            beam = 1
            ilike = 1
            call get_chisq_stat_unit(ilike)
         endif
         if (iNnm.eq.1) then
            beam = 1
            ilike = 2
            call get_chisq_stat_unit(ilike)
         endif
         if (iNae.eq.1) then
            beam = -1
            ilike = 1
            call get_chisq_stat_unit(ilike)
         endif
         if (iNam.eq.1) then
            beam = -1
            ilike = 2
            call get_chisq_stat_unit(ilike)
         endif
      endif
      if (iOki.eq.1) then
         iD = 2
         if (iNne.eq.1) then
            beam = 1
            ilike = 1
            call get_chisq_stat_unit(ilike)
         endif
         if (iNnm.eq.1) then
            beam = 1
            ilike = 2
            call get_chisq_stat_unit(ilike)
         endif
         if (iNae.eq.1) then
            beam = -1
            ilike = 1
            call get_chisq_stat_unit(ilike)
         endif
         if (iNam.eq.1) then
            beam = -1
            ilike = 2
            call get_chisq_stat_unit(ilike)
         endif
      endif
      if (iKr.eq.1) then
         iD = 3
         if (iNne.eq.1) then
            beam = 1
            ilike = 1
            call get_chisq_stat_unit(ilike)
         endif
         if (iNnm.eq.1) then
            beam = 1
            ilike = 2
            call get_chisq_stat_unit(ilike)
         endif
         if (iNae.eq.1) then
            beam = -1
            ilike = 1
            call get_chisq_stat_unit(ilike)
         endif
         if (iNam.eq.1) then
            beam = -1
            ilike = 2
            call get_chisq_stat_unit(ilike)
         endif
      endif      

      return
      end


      subroutine get_chisq_stat_unit(ilike)
C     ****************************************************
C     By Yoshitaro Takaesu @KIAS JAN 13 2014
C     ****************************************************
      implicit none
C     GLOBAL VARIABLES
      include 'inc/params.inc'
      include 'inc/main.inc'
      include 'inc/minfunc.inc'
C     CONSTANTS
C     ARGUMENTS 
      integer ilike
C     LOCAL VARIABLES 
C     EXTERNAL FUNCTIONS
c      real*8 chi2_stat
c      external chi2_stat
C     ----------
C     BEGIN CODE
C     ----------
c$$$      chisq_stat(ilike,iD,beam) = 
c$$$     &     chi2_stat(hErec_chi2_dat(1,ilike,iD,beam)
c$$$     &     ,hErec_chi2_th(1,ilike,iD,beam),nbins_chi2(ilike,iD)
c$$$     &     ,x_chi2(0,ilike,iD),outnbins(ilike,iD),outx(0,ilike,iD))

      call get_chi2_stat_bin(hErec_chi2_dat(1,ilike,iD,beam)
     &     ,hErec_chi2_th(1,ilike,iD,beam),nbins_chi2(ilike,iD)
     &     ,x_chi2(0,ilike,iD),iminevents,chi2_stat_bin(1,ilike,iD,beam)
     &     ,chisq_stat(ilike,iD,beam),outnbins(ilike,iD)
     &     ,outx(0,ilike,iD))

      return
      end
