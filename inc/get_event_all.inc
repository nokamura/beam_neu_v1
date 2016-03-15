C -*- fortran -*-
c     nu_e > nu_e
         zz(120) = 1 ! neutrino in the beam
         zz(118) = 1 ! detected neutrino type: 1:nu_e 2:nu_mu 3:anu_e 4:anu_mu
         call get_event_unit(zz,iD)
c     nu_e > nu_mu
         zz(120) = 1 ! neutrino in the beam
         zz(118) = 2 ! detected neutrino type: 1:nu_e 2:nu_mu 3:anu_e 4:anu_mu
         call get_event_unit(zz,iD)
c     nu_e > nu_tau
         zz(120) = 1 ! neutrino in the beam
         zz(118) = 3 ! detected neutrino type: 1:nu_e 2:nu_mu 3:anu_e 4:anu_mu
         call get_event_unit(zz,iD)
c     nu_mu > nu_e
         zz(120) = 2 ! neutrino in the beam
         zz(118) = 1 ! detected neutrino type: 1:nu_e 2:nu_mu 3:anu_e 4:anu_mu
         call get_event_unit(zz,iD)
c     nu_mu > nu_mu
         zz(120) = 2 ! neutrino in the beam
         zz(118) = 2 ! detected neutrino type: 1:nu_e 2:nu_mu 3:anu_e 4:anu_mu
         call get_event_unit(zz,iD)
c     nu_mu > nu_tau
         zz(120) = 2 ! neutrino in the beam
         zz(118) = 3 ! detected neutrino type: 1:nu_e 2:nu_mu 3:anu_e 4:anu_mu
         call get_event_unit(zz,iD)
c     anti-nu_e > anti-nu_e
         zz(120) = -1 ! neutrino in the beam
         zz(118) = -1 ! detected neutrino type: 1:nu_e 2:nu_mu 3:anu_e 4:anu_mu
         call get_event_unit(zz,iD)
c     anti-nu_e > anti-nu_mu
         zz(120) = -1 ! neutrino in the beam
         zz(118) = -2 ! detected neutrino type: 1:nu_e 2:nu_mu 3:anu_e 4:anu_mu
         call get_event_unit(zz,iD)
c     anti-nu_e > anti-nu_tau
         zz(120) = -1 ! neutrino in the beam
         zz(118) = -3 ! detected neutrino type: 1:nu_e 2:nu_mu 3:anu_e 4:anu_mu
         call get_event_unit(zz,iD)
c     anti-nu_mu > anti-nu_e
         zz(120) = -2 ! neutrino in the beam
         zz(118) = -1 ! detected neutrino type: 1:nu_e 2:nu_mu 3:anu_e 4:anu_mu
         call get_event_unit(zz,iD)
c     anti-nu_mu > anti-nu_mu
         zz(120) = -2 ! neutrino in the beam
         zz(118) = -2 ! detected neutrino type: 1:nu_e 2:nu_mu 3:anu_e 4:anu_mu
         call get_event_unit(zz,iD)
c     anti-nu_mu > anti-nu_tau
         zz(120) = -2 ! neutrino in the beam
         zz(118) = -3 ! detected neutrino type: 1:nu_e 2:nu_mu 3:anu_e 4:anu_mu
         call get_event_unit(zz,iD)
