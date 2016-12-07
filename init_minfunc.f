      subroutine init_minfunc(z)
C     ****************************************************
C     By Yoshitaro Takaesu @KIAS JAN 12 2014
C
C     Initialize variables in minfunc
C     ****************************************************
      implicit none
C     GLOBAL VARIABLES
      include 'inc/params.inc'
      include 'inc/main.inc'
      include 'inc/minfunc.inc'
C     CONSTANTS
C     ARGUMENTS 
      real*8 z(zdim)
C     LOCAL VARIABLES 
      integer i,j,k,l1,l2
C     EXTERNAL FUNCTIONS
C     ----------
C     BEGIN CODE
C     ----------
      if (idata.eq.1) then 
         MHH = MH
      elseif (idata.eq.-1) then 
         MHH = ihypo*MH
      endif
      rnevent_ren = -1  ! -1:No Normalize events other:Normalize events
      evform_dat = 2    ! 1:integer event number 2:real event number
      evform_th = 2      
      evform = evform_dat 

      binsize = basic_binsize*binsize_factor
      call bining_x(Emin,Emax,binsize,nbins,x,yy)
      
      if (mode.eq.1) then !for minimization
         rnevent_ren = -1  ! -1:No Normalize events other:Normalize events
      endif

      nu_mode = 1 ! dummy nu_mode for make_basic_data routine
    
      return
      end
