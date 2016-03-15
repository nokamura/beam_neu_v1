      real*8 function smear_CCQE_mu(Enu,Erec)
C     ****************************************************
C     By Yoshitaro Takaesu @KIAS AUG 25 2012
C
C     Enu [GeV], Erec: reconstructed muon energy [GeV]
C     ****************************************************
      implicit none
C     GLOBAL VARIABLES
C     CONSTANTS
      include 'inc/const.inc'
C     ARGUMENTS 
      real*8 Enu,Erec
C     LOCAL VARIABLES 
      integer i
      real*8 x,xi,r(3),sigma(3),dE(3),sum_rsig,A,sum_gauss
C     EXTERNAL FUNCTIONS
C     ----------
C     BEGIN CODE
C     ----------
      x = Enu -1d0
      xi = dsqrt(Enu) -1d0

      r(1) = 1d0
      r(2) = 1.1d0 -0.96*x +0.44*x**2 -0.076*x**3 +0.0047*x**4
      r(3) = 0.365d0 -1.97*xi*( 1d0 -0.634*xi +0.464*xi**2 +0.293*xi**3
     &     -0.342*xi**4 )/( 1d0 +1.23*x )

      sigma(1) = ( 39.7d0 +68.5*xi )*1d-3
      sigma(2) = ( 82.7d0 -50.2*x +259*xi )*1d-3
      sigma(3) = ( 197d0 +486*x -606*xi +203*x*xi )*1d-3

      dE(1) = ( 35d0 -2.5*x*( 1d0 -1.2*x )/( 1d0 +x ) )*1d-3
      dE(2) = dE(1)
      dE(3) = ( sigma(3)*( 0.053d0 +0.033*x 
     &     *( ( 1d0 +5.27*x -8.67*x**2 +1.83*x**3 )
     &     /( 1d0 +3.65*x +4.35*x**2 ) ) ) )*1d-3

      sum_rsig = 0d0
      do i = 1,3
         sum_rsig = sum_rsig +r(i)*sigma(i)
      enddo

      A = dsqrt(2*pi)*sum_rsig

      sum_gauss = 0d0
      do i = 1,3
         sum_gauss = sum_gauss +r(i)
     &        *dexp( -( Erec -Enu +dE(i) )**2/( 2*sigma(i)**2 ) )
      enddo

      smear_CCQE_mu = 1d0/A*sum_gauss

      return
      end
