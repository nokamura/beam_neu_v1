      subroutine get_event_allmode(z,hErec)
C     ****************************************************
C     By Yoshitaro Takaesu @KIAS JAN 12 2014
C
C     Calculate events for all the neutrino beam and 
C     detected neutrino types.
C     ****************************************************
      implicit none
C     GLOBAL VARIABLES
      include 'inc/params.inc'
C     CONSTANTS
C     ARGUMENTS 
      real*8 z(zdim),hErec(maxnbin,imaxint,nmode,-3:3,-3:3)
C     LOCAL VARIABLES
C     EXTERNAL FUNCTIONS
C     ----------
C     BEGIN CODE
C     ----------
c     nu_e > nu_e
      nu_mode = 1 ! neutrino in the beam: 1:nu_e 2:nu_mu 3:anu_e 4:anu_mu
      detect = 1  ! detected neutrino type
      call get_event_nu(z,hErec(1,1,1,detect,nu_mode))
c     nu_e > nu_mu
      nu_mode = 1 ! neutrino in the beam: 1:nu_e 2:nu_mu 3:anu_e 4:anu_mu
      detect = 2  ! detected neutrino type
      call get_event_nu(z,hErec(1,1,1,detect,nu_mode))
c     nu_e > nu_tau
      nu_mode = 1 ! neutrino in the beam: 1:nu_e 2:nu_mu 3:anu_e 4:anu_mu
      detect = 3  ! detected neutrino type
      call get_event_nu(z,hErec(1,1,1,detect,nu_mode))
c     nu_mu > nu_e
      nu_mode = 2 ! neutrino in the beam: 1:nu_e 2:nu_mu 3:anu_e 4:anu_mu
      detect = 1  ! detected neutrino type
      call get_event_nu(z,hErec(1,1,1,detect,nu_mode))
c     nu_mu > nu_mu
      nu_mode = 2 ! neutrino in the beam: 1:nu_e 2:nu_mu 3:anu_e 4:anu_mu
      detect = 2  ! detected neutrino type
      call get_event_nu(z,hErec(1,1,1,detect,nu_mode))
c     nu_mu > nu_tau
      nu_mode = 2 ! neutrino in the beam: 1:nu_e 2:nu_mu 3:anu_e 4:anu_mu
      detect = 3  ! detected neutrino type
      call get_event_nu(z,hErec(1,1,1,detect,nu_mode))
c     anti-nu_e > anti-nu_e
      nu_mode = -1 ! neutrino in the beam: 1:nu_e 2:nu_mu 3:anu_e 4:anu_mu
      detect = -1  ! detected neutrino type
      call get_event_nu(z,hErec(1,1,1,detect,nu_mode))
c     anti-nu_e > anti-nu_mu
      nu_mode = -1 ! neutrino in the beam: 1:nu_e 2:nu_mu 3:anu_e 4:anu_mu
      detect = -2  ! detected neutrino type
      call get_event_nu(z,hErec(1,1,1,detect,nu_mode))
c     anti-nu_e > anti-nu_tau
      nu_mode = -1 ! neutrino in the beam: 1:nu_e 2:nu_mu 3:anu_e 4:anu_mu
      detect = -3  ! detected neutrino type
      call get_event_nu(z,hErec(1,1,1,detect,nu_mode))
c     anti-nu_mu > anti-nu_e
      nu_mode = -2 ! neutrino in the beam: 1:nu_e 2:nu_mu 3:anu_e 4:anu_mu
      detect = -1  ! detected neutrino type
      call get_event_nu(z,hErec(1,1,1,detect,nu_mode))
c     anti-nu_mu > anti-nu_mu
      nu_mode = -2 ! neutrino in the beam: 1:nu_e 2:nu_mu 3:anu_e 4:anu_mu
      detect = -2  ! detected neutrino type
      call get_event_nu(z,hErec(1,1,1,detect,nu_mode))
c     anti-nu_mu > anti-nu_tau
      nu_mode = -2 ! neutrino in the beam: 1:nu_e 2:nu_mu 3:anu_e 4:anu_mu
      detect = -3  ! detected neutrino type
c      call get_event_nu(z,hErec(1,1,1,detect,nu_mode))
      call get_event_nu(z,hErec(1,1,1,detect,nu_mode))

c     NC pi0 BG
      nu_mode = 1  ! dummy parameter
      detect = 1   ! dummy parameter
      call get_event_pi0(z,hErec(1,1,1,detect,nu_mode))

      nu_mode = -1 ! dummy parameter
      detect = -1  ! dummy parameter
      call get_event_pi0(z,hErec(1,1,1,detect,nu_mode))

      return
      end
