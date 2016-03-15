      subroutine convert_th(s2sol_2,s2atm_2,s2rct_2,oct_atm,s212_2
     &     ,s223_2,s213_2,oct_23,ierror)
C     ****************************************************
C     By Yoshitaro Takaesu @KIAS MAY 30 2013
C
C     This subroutine converts theta_sol, theta_rct and theta_atm
C     to theta_12, theta_13 and theta_23.
C     ****************************************************
      implicit none
C     GLOBAL VARIABLES
C     CONSTANTS
C     ARGUMENTS 
      real*8 s2sol_2,s2atm_2,s2rct_2,s212_2,s223_2,s213_2
      integer oct_atm,oct_23,ierror
C     LOCAL VARIABLES 
      real*8 s2sol,satm_2,satm,srct_2,srct,s13,c13_2,s212,s23
C     EXTERNAL FUNCTIONS
C     ----------
C     BEGIN CODE
C     ----------
      s2sol = dsqrt(s2sol_2)
      satm_2 = 0.5*( 1d0 +oct_atm*dsqrt( 1d0 -s2atm_2 ) )
      satm = dsqrt(satm_2)
      srct_2 = 0.5*( 1d0 -dsqrt( 1d0 -s2rct_2 ) )
      srct = dsqrt(srct_2)

      s13 = srct
      c13_2 = 1d0 -s13**2 
      s212 = s2sol/c13_2
      s23 = satm/dsqrt(c13_2)

      if(s23**2.gt.0.5) then
         oct_23 = 1
      elseif(s23**2.le.0.5) then
         oct_23 = -1
      endif

      s212_2 = s212**2
      s223_2 = 4*s23**2*(1d0 -s23**2)
      s213_2 = s2rct_2      

      ierror = 0
      if (s212_2.gt.1) then
         ierror = 1
      elseif (s212_2.lt.0) then
         ierror = 1
      endif
      if (s223_2.gt.1) then
         ierror = 1
      elseif (s223_2.lt.0) then
         ierror = 1
      endif
      if (s213_2.gt.1) then
         ierror = 1
      elseif (s213_2.lt.0) then
         ierror = 1
      endif


      return
      end


      subroutine convert_thinv(s212_2,s223_2,s213_2,s2sol_2,s2atm_2
     &     ,s2rct_2)
C     ****************************************************
C     By Yoshitaro Takaesu @KIAS MAY 30 2013
C
C     This subroutine converts theta_sol, theta_rct and theta_atm
C     to theta_12, theta_13 and theta_23.
C     ****************************************************
      implicit none
C     GLOBAL VARIABLES
C     CONSTANTS
C     ARGUMENTS 
      real*8 s2sol_2,s2atm_2,s2rct_2,s212_2,s223_2,s213_2
C     LOCAL VARIABLES 
      real*8 s2sol,satm_2,satm,srct_2,srct,s13,c13_2,s212,s23
      real*8 s23_2,s13_2
C     EXTERNAL FUNCTIONS
C     ----------
C     BEGIN CODE
C     ----------
      s212 = dsqrt(s212_2)
      s23_2 = 0.5*( 1d0 -dsqrt( 1d0 -s223_2 ) )
      s23 = dsqrt(s23_2)
      s13_2 = 0.5*( 1d0 -dsqrt( 1d0 -s213_2 ) )
      s13 = dsqrt(s13_2)

      s2sol = s212*(1d0-s13_2) 
      satm = s23*dsqrt(1d0-s13_2) 

      s2sol_2 = s2sol**2
      s2atm_2 = 4*satm**2*(1d0 -satm**2)
      s2rct_2 = s213_2      

      return
      end
