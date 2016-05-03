      subroutine minfunc(npar,grad,chisq,z,iflag,futil)
C     ****************************************************
C     By Yoshitaro Takaesu @KIAS AUG 27  2012
C     ****************************************************
      implicitnone
C     GLOBAL VARIABLES
      include 'inc/params.inc'
      include 'inc/main.inc'
      include 'inc/minfunc.inc'
C     ARGUMENTS 
      integer npar,iflag
      real*8 z(zdim),grad,chisq
C     LOCAL VARIABLES 
      integer i,j,k,l1,l2,l3,ilike
C     EXTERNAL FUNCTIONS
      real*8 futil,chi2_para
      external futil,chi2_para
C     ----------
C     BEGIN CODE
C     ----------
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
CCCCCCCCCCCCCCCC   Initialization  CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      call init_minfunc(z)
c      call init_minfunc
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
CCCCCCCCCCCCCCCC   Data Taking  CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      if (idata.eq.1) then
         call make_basic_data
c         call make_prob_data
         call get_event_all(z_dat,hErec_raw_dat)
         call get_recevent(z_dat,nbins,hErec_raw_dat,hErec_dat)
         call classify_recevent(nbins,hErec_dat,hErec_like_dat)
         if (iSK.eq.1) then
            iD = 1
            call apply_Ereccut(iD,hErec_like_dat,hErec_chi2_dat)
         endif
         if (iOki.eq.1) then            
            iD = 2
            call apply_Ereccut(iD,hErec_like_dat,hErec_chi2_dat)
         endif
         if (iKr.eq.1) then            
            iD = 3
            call apply_Ereccut(iD,hErec_like_dat,hErec_chi2_dat)
         endif
      endif
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
CCCCCCCCCCCCCCCC  Calculation of the Chi^2 function  CCCCCCCCCCCCCCCC
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      if (idata.eq.-1) then
         call get_event_all(z,hErec_raw_th)
         call get_recevent(z,nbins,hErec_raw_th,hErec_th)
         call classify_recevent(nbins,hErec_th,hErec_like_th)
         if (iSK.eq.1) then
            iD = 1
            call apply_Ereccut(iD,hErec_like_th,hErec_chi2_th)
         endif
         if (iOki.eq.1) then            
            iD = 2
            call apply_Ereccut(iD,hErec_like_th,hErec_chi2_th)
         endif
         if (iKr.eq.1) then            
            iD = 3
            call apply_Ereccut(iD,hErec_like_th,hErec_chi2_th)
         endif
         call get_chisq_stat
         chisq = 0d0
         do k = -1,1,2
            do i = 1,ndetect
               do j = 1,nmode_like
                  chisq = chisq +chisq_stat(j,i,k)
               enddo
            enddo
         enddo
         call get_pull(nparx,z,z_dat,pull_fact,error,parflag)
         chisq_param_tot = 0d0
         do i = 1,nparx
            chisq_param_tot = chisq_param_tot +pull_fact(i)**2
         enddo
c         chisq_param_tot = chi2_para(nparx,z,z_dat,error,parflag)
         chisq = chisq +chisq_param_tot
      endif
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      return
      end
      
