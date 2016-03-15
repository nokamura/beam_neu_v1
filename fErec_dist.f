      real*8 function fErec_dist(x,z)
c     E: incoming neutrino energy in GeV
      implicit none
      include 'inc/const.inc'
      real*8 x,z(200)
      integer int_mode,nutype
      real*8 E,EvG
      real*8 Erecdist_ae_ccqeh,Erecdist_ae_ccqeo,Erecdist_ae_ccrsh
      real*8 Erecdist_ae_ccrso
      real*8 Erecdist_am_ccqeh,Erecdist_am_ccqeo,Erecdist_am_ccrsh
      real*8 Erecdist_am_ccrso
      real*8 Erecdist_nm_ccqeo,Erecdist_nm_ccrso,Erecdist_ne_ccqeo
      real*8 Erecdist_ne_ccrso
      external Erecdist_am_ccqeh,Erecdist_am_ccqeo,Erecdist_am_ccrsh
      external Erecdist_am_ccrso
      external Erecdist_ae_ccqeh,Erecdist_ae_ccqeo,Erecdist_ae_ccrsh
      external Erecdist_ae_ccrso
      external Erecdist_nm_ccqeo,Erecdist_nm_ccrso,Erecdist_ne_ccqeo
      external Erecdist_ne_ccrso

      nutype = int(z(1)) ! incoming neutrino type: 1:nu_e 2:nu_mu 3:nu_tau -1:anti-nu_e -2: anti-nu_mu -3:anti-nu_tau 
      int_mode = int(z(2)) ! interaction mode: 1:CCQE-H 2:CCQE-O 3:RES-H 4:RES-O
      EvG = z(3) ! incoming nu energy in GeV
      E = x*1d3 ! x in MeV
      
      if (nutype.eq.1) then
         if (int_mode.eq.1) then ! CCQE H
            fErec_dist = 0d0    !   nutrino deos not interact with a proton via CCQE.
         elseif (int_mode.eq.2) then ! CCQE O
            fErec_dist = Erecdist_ne_ccqeo(EvG,E)
         elseif (int_mode.eq.3) then ! CCRS H
            fErec_dist = 0d0   !  we ignore this interaction (Under construction).
         elseif (int_mode.eq.4) then ! CCRS O
            fErec_dist = Erecdist_ne_ccrso(EvG,E)
         endif
      elseif (nutype.eq.2) then
         if (int_mode.eq.1) then ! CCQE H
            fErec_dist = 0d0   !   nutrino deos not interact with a proton via CCQE.
         elseif (int_mode.eq.2) then ! CCQE O
            fErec_dist = Erecdist_nm_ccqeo(EvG,E)
         elseif (int_mode.eq.3) then ! CCRS H
            fErec_dist = 0d0    !  we ignore this interaction (Under construction).
         elseif (int_mode.eq.4) then ! CCRS O
            fErec_dist = Erecdist_nm_ccrso(EvG,E)
         endif
      elseif (nutype.eq.3) then
         fErec_dist = 0d0
      elseif (nutype.eq.-1) then ! anti-nu_e
         if (int_mode.eq.1) then ! CCQE H
            fErec_dist = Erecdist_ae_ccqeh(EvG,E)
         elseif (int_mode.eq.2) then ! CCQE O
            fErec_dist = Erecdist_ae_ccqeo(EvG,E)
         elseif (int_mode.eq.3) then ! CCRS H
            fErec_dist = Erecdist_ae_ccrsh(EvG,E)
         elseif (int_mode.eq.4) then ! CCRS O
            fErec_dist = Erecdist_ae_ccrso(EvG,E)
         endif
      elseif (nutype.eq.-2) then
         if (int_mode.eq.1) then ! CCQE H
            fErec_dist = Erecdist_am_ccqeh(EvG,E)
         elseif (int_mode.eq.2) then ! CCQE O
            fErec_dist = Erecdist_am_ccqeo(EvG,E)
         elseif (int_mode.eq.3) then ! CCRS H
            fErec_dist = Erecdist_am_ccrsh(EvG,E)
         elseif (int_mode.eq.4) then ! CCRS O
            fErec_dist = Erecdist_am_ccrso(EvG,E)
         endif
      elseif (nutype.eq.-3) then
         fErec_dist = 0d0
      endif

      return
      end


      real*8 function Erecdist_nm_ccqeh(EvG,Erec)
c     Erec: Reconstructed neutrino energy [MeV]
c     EvG: incoming neutrino energy [GeV]
      implicit none
      real*8 EvG,Erec
      real*8 Ev,r1,r2,r3,sig1,sig2,sig3,dE1,dE2,dE3
      real*8 func
      Erecdist_nm_ccqeh = 0d0
      return
      end

      real*8 function Erecdist_nm_ccqeo(EvG,Erec)
c     Erec: Reconstructed neutrino energy [MeV]
c     EvG: incoming neutrino energy [GeV]
      implicit none
      real*8 EvG,Erec
      real*8 x,xi
      real*8 Ev,r1,r2,r3,sig1,sig2,sig3,dE1,dE2,dE3
      real*8 func
      Ev = EvG*1d3
      x = EvG -1d0
      xi = dsqrt(EvG) -1d0
      r1 = 1d0
      sig1 = 39.7d0 +68.5*xi
      dE1 = 35d0 -2.5*x*(1d0 -1.2*x)/(1d0 +x)
      r2 = 1.1d0 -0.96*x +0.44*x**2 -0.076*x**3 +0.0047*x**4
      sig2 = 82.7d0 -50.2*x +259*xi
      dE2 = dE1
      r3 = 0.365d0 -1.97*xi
     &     *( 1d0 -0.634*xi +0.464*xi**2 +0.293*xi**3 -0.342*xi**4 )
     &     /(1d0 +1.23*x)
      sig3 = 197d0 +486*x -606*xi +203*x*xi
      dE3 = sig3*( 0.053d0 +0.033*x
     &     *( (1d0 +5.27*x -8.67*x**2 +1.83*x**3)
     &     /(1d0 +3.65*x +4.35*x**2) ) ) 

      func = r1*dexp( -(Erec -Ev +dE1 )**2/(2*sig1**2) )
     &        +r2*dexp( -(Erec -Ev +dE2 )**2/(2*sig2**2) )
     &        +r3*dexp( -(Erec -Ev +dE3 )**2/(2*sig3**2) )      
      Erecdist_nm_ccqeo = func
      return
      end


      real*8 function Erecdist_nm_ccrsh(EvG,Erec)
c     Erec: Reconstructed neutrino energy [MeV]
c     EvG: incoming neutrino energy [GeV]
      implicit none
      real*8 EvG,Erec
      real*8 Ev,r1,r2,r3,r4,sig1,sig2,sig3,sig4,dE1,dE2,dE3,dE4
      real*8 func
      Erecdist_nm_ccrsh = 0d0
      return
      end

      real*8 function Erecdist_nm_ccrso(EvG,Erec)
c     Erec: Reconstructed neutrino energy [MeV]
c     EvG: incoming neutrino energy [GeV]
      implicit none
      real*8 EvG,Erec
      real*8 x,xi
      real*8 Ev,r1,r2,r3,r4,sig1,sig2,sig3,sig4,dE1,dE2,dE3,dE4
      real*8 func
      Ev = EvG*1d3
      if (EvG.le.1.2d0) then
         x = EvG -1d0
         xi = dsqrt(EvG) -1d0
         r1 = 1d0
         sig1 = 97.8d0 -3670*x +7410*xi +1540*x*xi
         dE1 = 382d0 -8170*x +16200*xi +3670*x*xi
         r2 = 0.42d0 +0.16*x -1.9*x**2
         sig2 = 98d0 +390*x +500*x**2
         dE2 = 579d0 -22600*x +45800*xi +10800*x*xi
         r3 = 0.1d0 -0.202*x/(1d0 +2.03*x)
         sig3 = 27d0 -22*x -33*x**2 -480*x**3 -1000*x**4
         dE3 = 210d0 -12*x -7.8*x**2 -170*x**3         
         func = r1*dexp( -(Erec -Ev +dE1 )**2/(2*sig1**2) )
     &        +r2*dexp( -(Erec -Ev +dE2 )**2/(2*sig2**2) )
     &        +r3*dexp( -(Erec -Ev +dE3 )**2/(2*sig3**2) )      
      elseif (EvG.gt.1.2d0) then
         x = EvG -2d0
         xi = dsqrt(EvG/2d0) -1d0
         r1 = 1d0
         sig1 = 110d0 +12*x +x**2
         dE1 = 350d0 -1.2*x +1.4*x**2
         r2 = 0.31d0 -0.18*x +0.22*x**2 -0.096*x**3 +0.015*x**4
         sig2 = 160d0 +36*x +2.4*x**2
         dE2 = 530d0 +20*x -4.9*x**2 +0.58*x**3
         r3 = 0.194d0 +0.316*x -1.37*xi -0.122*x*xi
         sig3 = 320d0 -1240*x +5640*xi +328*x*xi
         dE3 = 823d0 -( 3550*x -14100*xi -1840*x*xi )/(1d0 +0.921*x)
         r4 = 0.0369d0 -0.0728*x +0.254*xi +0.0356*x*xi
         sig4 = 177d0 +2580*x -9190*xi -747*x*xi
         dE4 = 1500d0 +990*x -1500*xi -760*x*xi
         func = r1*dexp( -(Erec -Ev +dE1 )**2/(2*sig1**2) )
     &        +r2*dexp( -(Erec -Ev +dE2 )**2/(2*sig2**2) )
     &        +r3*dexp( -(Erec -Ev +dE3 )**2/(2*sig3**2) )      
     &        +r4*dexp( -(Erec -Ev +dE4 )**2/(2*sig4**2) )      
      else
         func = 0d0
      endif
      Erecdist_nm_ccrso = func
      return
      end

      real*8 function Erecdist_ne_ccqeh(EvG,Erec)
c     Erec: Reconstructed neutrino energy [MeV]
c     EvG: incoming neutrino energy [GeV]
      implicit none
      real*8 EvG,Erec
      real*8 Ev,r1,r2,r3,sig1,sig2,sig3,dE1,dE2,dE3
      real*8 func
      Erecdist_ne_ccqeh = 0d0
      return
      end

      real*8 function Erecdist_ne_ccqeo(EvG,Erec)
c     Erec: Reconstructed neutrino energy [MeV]
c     EvG: incoming neutrino energy [GeV]
      implicit none
      real*8 EvG,Erec
      real*8 x,xi,ErecG
      real*8 Ev,r1,r2,r3,sig1,sig2,sig3,dE1,dE2,dE3
      real*8 func
      ErecG = Erec/1d3
      Ev = EvG*1d3
      x = EvG -1d0
      xi = dsqrt(EvG) -1d0
      r1 = 1d0
      sig1 = 55.5d0 -19.6*x +98.9*xi
      dE1 = 40d0 -0.99*x +3.3*x**2 -0.71*x**3 -2.2*x/(1d0 +x)
      r2 = 0.67d0 -0.58*x +0.58*x**2 -0.16*x**3 +0.019*x**4
      sig2 = 125d0 -51.3*x +201*xi
      dE2 = dE1
      r3 = 0.094d0 -0.040*x +0.031*x**2 -0.016*x**3 +0.0059*x**4
      sig3 = 273d0 -102*x +1560*xi +111*x*xi
      dE3 = sig3*( -0.16d0 +0.68*x -2.6*xi +1.1*x/(1d0 +x) )

      func = r1*dexp( -(Erec -Ev +dE1 )**2/(2*sig1**2) )
     &        +r2*dexp( -(Erec -Ev +dE2 )**2/(2*sig2**2) )
     &        +r3*dexp( -(Erec -Ev +dE3 )**2/(2*sig3**2) )      
      Erecdist_ne_ccqeo = func
      return
      end


      real*8 function Erecdist_ne_ccrsh(EvG,Erec)
c     Erec: Reconstructed neutrino energy [MeV]
c     EvG: incoming neutrino energy [GeV]
      implicit none
      real*8 EvG,Erec
      real*8 Ev,r1,r2,r3,r4,sig1,sig2,sig3,sig4,dE1,dE2,dE3,dE4
      real*8 func
      Erecdist_ne_ccrsh = 0d0
      return
      end


      real*8 function Erecdist_ne_ccrso(EvG,Erec)
c     Erec: Reconstructed neutrino energy [MeV]
c     EvG: incoming neutrino energy [GeV]
      implicit none
      real*8 EvG,Erec
      real*8 x,xi
      real*8 Ev,r1,r2,r3,r4,sig1,sig2,sig3,sig4,dE1,dE2,dE3,dE4
      real*8 func
      Ev = EvG*1d3
      if (EvG.le.1.2d0) then
         x = EvG -1d0
         xi = dsqrt(EvG) -1d0
         r1 = 1d0
         sig1 = 102d0 -2960*x +5940*xi +1300*x*xi
         dE1 = 372d0 -9530*x +18900*xi +4630*x*xi
         r2 = 0.37d0 -5.8*x +13*xi
         sig2 = 114d0 +370.0*x +244*x**2
         dE2 = 580.0d0 -15600*x +31400*xi +6070*x*xi
         r3 = 0.067d0 -0.068*x -0.37*x**2 -3.4*x**3
         sig3 = 27.1d0 -717*x +1390*xi +425*x*xi
         dE3 = sig3*( 7.6d0 +5.3*x -9.4*x**2 -24*x**3 )
         func = r1*dexp( -(Erec -Ev +dE1 )**2/(2*sig1**2) )
     &        +r2*dexp( -(Erec -Ev +dE2 )**2/(2*sig2**2) )
     &        +r3*dexp( -(Erec -Ev +dE3 )**2/(2*sig3**2) )      
      elseif (EvG.gt.1.2d0) then
         x = EvG -2d0
         xi = dsqrt(EvG/2d0) -1d0
         r1 = 1d0
         sig1 = 110d0 +14*x -0.50*x**2 
         dE1 = 360d0 -x -1.6*x**2 +0.95*x**3 -0.026*x**4
         r2 = 0.28d0 -0.092*x +0.23*x**2 -0.11*x**3 +0.017*x**4
         sig2 = 181d0 +2150*x +207*x**2 -8440*xi -1890*x*xi 
         dE2 = 510d0 -61.1*x +261*xi
         r3 = 0.196d0 -0.487*x +2.02*xi +0.148*x*xi
         sig3 = 334d0 -1750*x +7620*xi +556*x*xi
         dE3 = 796d0 -( 3260*x +115*x**2 -12900*xi -2070*x*xi )
     &        /(1d0 +0.94*x)
         r4 = 0.0298d0 +0.05*xi -0.213*xi**2 +0.23*xi**3
         sig4 = 222d0 +10300*x -4030*xi -3570*x*xi
         dE4 = 1500d0 -690*x +4600*xi 
         func = r1*dexp( -(Erec -Ev +dE1 )**2/(2*sig1**2) )
     &        +r2*dexp( -(Erec -Ev +dE2 )**2/(2*sig2**2) )
     &        +r3*dexp( -(Erec -Ev +dE3 )**2/(2*sig3**2) )      
     &        +r4*dexp( -(Erec -Ev +dE4 )**2/(2*sig4**2) )      
      else
         func = 0d0
      endif
      Erecdist_ne_ccrso = func
      return
      end


      real*8 function Erecdist_am_ccqeh(EvG,Erec)
c     Erec: Reconstructed neutrino energy [MeV]
c     EvG: incoming neutrino energy [GeV]
      implicit none
      real*8 EvG,Erec
      real*8 Ev,r1,r2,r3,sig1,sig2,sig3,dE1,dE2,dE3
      real*8 func

      r1 = 1d0
      sig1 = 3.20d0 +25.5*EvG
      dE1 = 0d0
      r2 = 0.106d0
      sig2 = -3.73d0 +52.7*EvG
      dE2 = 11.6d0 -17.0*EvG

      Ev = EvG*1d3              ! incoming nu energy in MeV
      if (Ev.le.700) then
         func = r1*dexp( -(Erec -Ev +dE1 )**2/(2*sig1**2) )
      elseif (Ev.gt.700) then
         func = r1*dexp( -(Erec -Ev +dE1 )**2/(2*sig1**2) )
     &        +r2*dexp( -(Erec -Ev +dE2 )**2/(2*sig2**2) )
      endif
      Erecdist_am_ccqeh = func

      return
      end


      real*8 function Erecdist_am_ccqeo(EvG,Erec)
c     Erec: Reconstructed neutrino energy [MeV]
c     EvG: incoming neutrino energy [GeV]
      implicit none
      real*8 EvG,Erec
      real*8 Ev,r1,r2,r3,sig1,sig2,sig3,dE1,dE2,dE3
      real*8 func

      r1 = 3.20d5 -2.16d5*EvG +5.62d4*EvG**2 -5.04d3*EvG**3
      sig1 = 11.3d0 +27.4*EvG -1.01*EvG**2
      dE1 = 43.5d0 -9.64*EvG +2.86*EvG**2 -0.229*EvG**3
      r2 = 2.05d5 -1.52d5*EvG +4.06d4*EvG**2 -3.60d3*EvG**3
      sig2 = 36.0d0 +49.0*EvG -3.31*EvG**2
      dE2 = 28.1d0 -2.89*EvG -0.482*EvG**2
      r3 = 1.10d4 -8.28d3*EvG +2.24d3*EvG**2 -198*EvG**3
      sig3 = -18.0d0 +241*EvG -56.8*EvG**2 +4.85*EvG**3
      dE3 = -4.84d0 -5.34*EvG -1.59*EvG**2

      Ev = EvG*1d3              ! incoming nu energy in MeV
      func = r1*dexp( -(Erec -Ev +dE1 )**2/(2*sig1**2) )
     &     +r2*dexp( -(Erec -Ev +dE2 )**2/(2*sig2**2) )
     &     +r3*dexp( -(Erec -Ev +dE3 )**2/(2*sig3**2) )
      Erecdist_am_ccqeo = func

      return
      end


      real*8 function Erecdist_am_ccrsh(EvG,Erec)
c     Erec: Reconstructed neutrino energy [MeV]
c     EvG: incoming neutrino energy [GeV]
      implicit none
      real*8 EvG,Erec
      real*8 Ev,r1,r2,r3,r4,sig1,sig2,sig3,sig4,dE1,dE2,dE3,dE4
      real*8 func

      r1 = 1.00d4*exp(-0.442*EvG) +7.84d4*exp(-2.90*EvG)
      sig1 = 18.0 +22.8*EvG
      dE1 = 212 +3.81*EvG -0.536*EvG**2
      r2 = 1.70d4*exp( -0.912*EvG )
      sig2 = 8.61 +13.5*EvG
      dE2 = 202 -1.63*EvG
      r3 = 8.93d3*Exp( -1.53*EvG )
      sig3 = 31.5 +46.6*EvG -6.27*EvG**2
      dE3 = 284 -6.23*EvG +16.1*EvG**2
      r4 = 624*( exp(-0.352*EvG) -exp(-0.434*EvG +0.0365))
      sig4 = -498 +763*EvG -203*EvG**2 +18.7*EvG**3
      dE4 = 655 +35*EvG

      Ev = EvG*1d3              ! incoming nu energy in MeV
      if (Ev.le.1000) then
         func = r1*dexp( -(Erec -Ev +dE1 )**2/(2*sig1**2) )
     &        +r2*dexp( -(Erec -Ev +dE2 )**2/(2*sig2**2) )
     &        +r3*dexp( -(Erec -Ev +dE3 )**2/(2*sig3**2) )
      elseif (Ev.gt.1000) then
         func = r1*dexp( -(Erec -Ev +dE1 )**2/(2*sig1**2) )
     &        +r2*dexp( -(Erec -Ev +dE2 )**2/(2*sig2**2) )
     &        +r3*dexp( -(Erec -Ev +dE3 )**2/(2*sig3**2) )
     &        +r4*dexp( -(Erec -Ev +dE4 )**2/(2*sig4**2) )
      endif
      Erecdist_am_ccrsh = func

      return
      end


      real*8 function Erecdist_am_ccrso(EvG,Erec)
c     Erec: Reconstructed neutrino energy [MeV]
c     EvG: incoming neutrino energy [GeV]
      implicit none
      real*8 EvG,Erec
      real*8 Ev,r1,r2,r3,r4,sig1,sig2,sig3,sig4,dE1,dE2,dE3,dE4
      real*8 func

      r1 = 1d0
      sig1 = 70.3d0 +13.9*EvG
      dE1 = 355d0
      r2 = 0.325d0 -0.0652*EvG +0.0142*EvG**2
      sig2 = 114d0 +40.1*EvG
      dE2 = 539d0/( 1d0 +1.53*exp(-2.76*EvG) )
      r3 = 0.0804d0 -0.00176*EvG +0.00194*EvG**2
      sig3 = 452*( exp(0.0105*EvG) -exp(1.98d0 -1.67*EvG) )
      dE3 = 850d0 +32*EvG
      r4 = -0.326d0 +0.229*EvG
      sig4 = 111*( exp(-1.67*EvG) -exp(-323*EvG) )
      dE4 = -214d0 +998*EvG

      Ev = EvG*1d3              ! incoming nu energy in MeV
      if (Ev.le.700) then
         func = r1*dexp( -(Erec -Ev +dE1 )**2/(2*sig1**2) )
     &        +r4*dexp( -(Erec -Ev +dE4 )**2/(2*sig4**2) )
      elseif ((Ev.gt.700).and.(Ev.le.1000)) then
         func = r1*dexp( -(Erec -Ev +dE1 )**2/(2*sig1**2) )
     &        +r2*dexp( -(Erec -Ev +dE2 )**2/(2*sig2**2) )
     &        +r4*dexp( -(Erec -Ev +dE4 )**2/(2*sig4**2) )
      elseif ((Ev.gt.1000).and.(Ev.le.1500)) then
         func = r1*dexp( -(Erec -Ev +dE1 )**2/(2*sig1**2) )
     &        +r2*dexp( -(Erec -Ev +dE2 )**2/(2*sig2**2) )
     &        +r3*dexp( -(Erec -Ev +dE3 )**2/(2*sig3**2) )
     &        +r4*dexp( -(Erec -Ev +dE4 )**2/(2*sig4**2) )
      elseif (Ev.gt.1500) then
         func = r1*dexp( -(Erec -Ev +dE1 )**2/(2*sig1**2) )
     &        +r2*dexp( -(Erec -Ev +dE2 )**2/(2*sig2**2) )
     &        +r3*dexp( -(Erec -Ev +dE3 )**2/(2*sig3**2) )
      endif
      Erecdist_am_ccrso = func

      return
      end


      real*8 function Erecdist_ae_ccqeh(EvG,Erec)
c     Erec: Reconstructed neutrino energy [MeV]
c     EvG: incoming neutrino energy [GeV]
      implicit none
      real*8 EvG,Erec
      real*8 Ev,r1,r2,r3,sig1,sig2,sig3,dE1,dE2,dE3
      real*8 func

      r1 = 1d0
      sig1 = 268*( dexp(-0.0225*EvG) -dexp(-0.135*EvG -0.0470d0) )
      dE1 = 18.3d0/( 1d0 +34.0*dexp(-1.02*EvG) )
      r2 = 4.47*( dexp(0.114*EvG -2.59d0) 
     &     -dexp(-0.698*EvG -2.36d0) )
      sig2 = 297*( 1d0 -dexp(-0.235*EvG) )
      dE2 = -7.26d0 -5.66*EvG
      r3 = 0.00532*dexp(0.620*EvG)
      sig3 = 4.41*( dexp(0.0584*EvG +4.27d0) 
     &     -dexp(-0.489*EvG +4.47d0) )
      dE3 = 665*( dexp(-0.261*EvG +0.0811d0) -dexp(-0.107*EvG) )
      
      Ev = EvG*1d3              ! incoming nu energy in MeV
      if (Ev.le.1000) then
         func = r1*dexp( -(Erec -Ev +dE1 )**2/(2*sig1**2) )
     &        +r2*dexp( -(Erec -Ev +dE2 )**2/(2*sig2**2) )
      elseif (Ev.gt.1000) then
         func = r1*dexp( -(Erec -Ev +dE1 )**2/(2*sig1**2) )
     &        +r2*dexp( -(Erec -Ev +dE2 )**2/(2*sig2**2) )
     &        +r3*dexp( -(Erec -Ev +dE3 )**2/(2*sig3**2) )
      endif
      Erecdist_ae_ccqeh = func

      return
      end

      real*8 function Erecdist_ae_ccqeo(EvG,Erec)
c     Erec: Reconstructed neutrino energy [MeV]
c     EvG: incoming neutrino energy [GeV]
      implicit none
      real*8 EvG,Erec
      real*8 Ev,r1,r2,r3,sig1,sig2,sig3,dE1,dE2,dE3
      real*8 func

      r1 = 1d0
      sig1 = 23.3d0 +25.8*EvG -1.54*EvG**2
      dE1 = 42.9d0 -11.0*EvG +4.42*EvG**2 -0.382*EvG**3
      r2 = 0.586d0 -0.258*EvG +0.122*EvG**2 -0.0137*EvG**3
      sig2 = 53.8d0 +46.5*EvG -2.79*EvG**2
      dE2 = 27.3d0 -6.43*EvG
      r3 = 0.0140d0 -0.00131*EvG +0.00502*EvG**2
      sig3 = 12.5d0 +250*EvG -66.1*EvG**2 +6.43*EvG**3
      dE3 = -176d0 +178*EvG -86.5*EvG**2 +15.5*EvG**3 -1.04*EvG**4
      
      Ev = EvG*1d3              ! incoming nu energy in MeV
      func = r1*dexp( -(Erec -Ev +dE1 )**2/(2*sig1**2) )
     &          +r2*dexp( -(Erec -Ev +dE2 )**2/(2*sig2**2) )
     &          +r3*dexp( -(Erec -Ev +dE3 )**2/(2*sig3**2) )
      Erecdist_ae_ccqeo = func

      return
      end


      real*8 function Erecdist_ae_ccrsh(EvG,Erec)
c     Erec: Reconstructed neutrino energy [MeV]
c     EvG: incoming neutrino energy [GeV]
      implicit none
      real*8 EvG,Erec
      real*8 Ev,r1,r2,r3,r4,sig1,sig2,sig3,sig4,dE1,dE2,dE3,dE4
      real*8 func

      r1 = 7.99d3*exp(-0.382*EvG) +5.27d4*exp(-1.89*EvG)
      sig1 = 17.0d0 +23.9*EvG -1.16*EvG**2
      dE1 = 206d0
      r2 = 1.81d4*exp(-3.14*EvG) +4.05d3*exp(-0.629*EvG)
      sig2 = 29.5d0 +46.2*EvG -1.70*EvG**2
      dE2 = 272d0 -3.13*exp(0.682*EvG)
      r3 = 93.7*exp(-0.186*EvG) +3.00d3*exp(-4.26*EvG)
      sig3 = -162d0 +345*EvG -42.6*EvG**2
      dE3 = 210d0 +387*EvG -60*EvG**2
      r4 = 97.8*exp(-0.379*EvG)
      sig4 = 62.6*exp(0.651*EvG)
      dE4 = -1.43d3 +1.66d3*EvG -118*EvG**2
      
      Ev = EvG*1d3              ! incoming nu energy in MeV
      if (Ev.le.1500) then
         func = r1*dexp( -(Erec -Ev +dE1 )**2/(2*sig1**2) )
     &        +r2*dexp( -(Erec -Ev +dE2 )**2/(2*sig2**2) )
     &        +r3*dexp( -(Erec -Ev +dE3 )**2/(2*sig3**2) )
      elseif (Ev.gt.1500) then
         func = r1*dexp( -(Erec -Ev +dE1 )**2/(2*sig1**2) )
     &        +r2*dexp( -(Erec -Ev +dE2 )**2/(2*sig2**2) )
     &        +r3*dexp( -(Erec -Ev +dE3 )**2/(2*sig3**2) )
     &        +r4*dexp( -(Erec -Ev +dE4 )**2/(2*sig4**2) )
      endif
      Erecdist_ae_ccrsh = func

      return
      end


      real*8 function Erecdist_ae_ccrso(EvG,Erec)
c     Erec: Reconstructed neutrino energy [MeV]
c     EvG: incoming neutrino energy [GeV]
      implicit none
      real*8 EvG,Erec
      real*8 Ev,r1,r2,r3,r4,sig1,sig2,sig3,sig4,dE1,dE2,dE3,dE4
      real*8 func

      r1 = 1d0
      sig1 = 72.9d0 +12.8*EvG
      dE1 = 348d0 +0.352*EvG
      r2 = 0.143d0 +0.0605*EvG
      sig2 = 105d0 +40.2*EvG
      dE2 = 443d0/( 1d0 -0.434*dexp(-0.902*EvG) )
      r3 = -0.0459d0 +0.105*EvG -0.0134*EvG**2
      sig3 = 147d0 +92.7*EvG
      dE3 = 753d0 +17.6*EvG

      Ev = EvG*1d3              ! incoming nu energy in MeV
      if (Ev.le.700) then
         func = r1*dexp( -(Erec -Ev +dE1 )**2/(2*sig1**2) )
      elseif ((Ev.gt.700).and.(Ev.le.1000)) then
         func = r1*dexp( -(Erec -Ev +dE1 )**2/(2*sig1**2) )
     &        +r2*dexp( -(Erec -Ev +dE2 )**2/(2*sig2**2) )
      elseif (Ev.gt.1000) then
         func = r1*dexp( -(Erec -Ev +dE1 )**2/(2*sig1**2) )
     &        +r2*dexp( -(Erec -Ev +dE2 )**2/(2*sig2**2) )
     &        +r3*dexp( -(Erec -Ev +dE3 )**2/(2*sig3**2) )
      endif
      Erecdist_ae_ccrso = func

      return
      end
