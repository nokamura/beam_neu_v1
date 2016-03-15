      real*8 function flux_react(E,P)
C     ****************************************************
C     By Yoshitaro Takaesu @KIAS Aug.31 2012
C     
C     Anti-electron-neutrino flux of energy E[MeV] 
C     from the reactors with the total thermal 
C     power of P[GW] in [ 1/s/MeV ] unit
C     ****************************************************
      implicitnone
C     
C     ARGUMENTS 
C     
      real*8 E,P
      real*8 a,EU235,EU238,EPu239,EPu241,fU235,fU238,fPu239,fPu241
      real*8 preflux,pi
C     ----------
C     BEGIN CODE
C     ----------
      pi = dacos(-1d0)

      fU235 = 0.58d0
      fU238 = 0.07d0
      fPu239 = 0.30d0
      fPu241 = 0.05d0
      EU235 = 201.7d0
      EU238 = 205.0d0
      EPu239 = 210.0d0
      EPu241 = 212.4d0

      a = 1d0 / ( fU235*EU235 +fU238*EU238 +fPu239*EPu239 
     &     +fPu241*EPu241 )

      preflux = fU235*dexp( 0.870d0 -0.160*E -0.091*E**2 )
     &      +fPu239*dexp( 0.896d0 -0.239*E -0.0981*E**2 )
     &      +fU238*dexp( 0.976d0 -0.162*E -0.0790*E**2 )
     &      +fPu241*dexp( 0.793d0 -0.080*E -0.1085*E**2 )

      flux_react = P*a*preflux*6.24*1d21
c      flux = a*preflux*1d12/1.602176487d0
c      flux = a*preflux*1d22/1.602176487d0


      return
      end


      real*8 function flux_beam(v1,oab,nu_mode,E,L)
C     ****************************************************
C     By Yoshitaro Takaesu @KIAS May.21 2013
C     
C     Neutrino beam flux in [ 1/GeV/1kton/10^21POT ] unit 
C     of energy E[GeV] at a L [km] point from J-PARK with
C     OAB angle = oab
C     ****************************************************
      implicitnone
C     GLOBAL VARIABLES
      integer oi
      include 'beam/oab00n.inc'
      include 'beam/oab05n.inc'
      include 'beam/oab06n.inc'
      include 'beam/oab09n.inc'
      include 'beam/oab10n.inc'
      include 'beam/oab14n.inc'
      include 'beam/oab15n.inc'
      include 'beam/oab20n.inc'
      include 'beam/oab25n.inc'
      include 'beam/oab30n.inc'
      include 'beam/oab00a.inc'
      include 'beam/oab05a.inc'
      include 'beam/oab06a.inc'
      include 'beam/oab09a.inc'
      include 'beam/oab10a.inc'
      include 'beam/oab14a.inc'
      include 'beam/oab15a.inc'
      include 'beam/oab20a.inc'
      include 'beam/oab25a.inc'
      include 'beam/oab30a.inc'
C     ARGUMENTS 
      integer v1,nu_mode
      real*8 E,L,oab
C     LOCAL VARIABLES
      integer i,icol
      real*8 flux
C     EXTERNAL FUNCS
      real*8 read_tbl,read_histo,binw
      external read_tbl,read_histo
C     ----------
C     BEGIN CODE
C     ----------
      if (abs(v1).ne.1) then
         write(*,*) "ERROR: flux_beam in funcs.f: Invalid beam mode"
         stop
      endif
      if (abs(nu_mode).ne.1) then
         if (abs(nu_mode).ne.2) then
            write(*,*) "ERROR: flux_beam in funcs.f:"
     &           ," Invalid neutrino mode"
            stop
         endif
      endif
      if (E.lt.0) then
         write(*,*) "ERROR: flux_beam in funcs.f:"
     &        ," E should be positive"
         stop
      endif

c      binw = 0.05d0
      binw = 1d0

      if (nu_mode.eq.1) icol = 4
      if (nu_mode.eq.-1) icol = 5
      if (nu_mode.eq.2) icol = 2
      if (nu_mode.eq.-2) icol = 3

c      if (abs(v1).eq.1) then
c         stop
c      elseif (v1.eq.1) then
      if (v1.eq.1) then
         if (oab.eq.0d0) flux = read_histo(oab00n,300,5,binw,icol,E)
         if (oab.eq.0.5d0) flux = read_histo(oab05n,300,5,binw,icol,E)
         if (oab.eq.0.6d0) flux = read_histo(oab06n,300,5,binw,icol,E)
         if (oab.eq.0.9d0) flux = read_histo(oab09n,300,5,binw,icol,E)
         if (oab.eq.1.0d0) flux = read_histo(oab10n,300,5,binw,icol,E)
         if (oab.eq.1.4d0) flux = read_histo(oab14n,300,5,binw,icol,E)
         if (oab.eq.1.5d0) flux = read_histo(oab15n,300,5,binw,icol,E)
         if (oab.eq.2.0d0) flux = read_histo(oab20n,300,5,binw,icol,E)
         if (oab.eq.2.5d0) flux = read_histo(oab25n,300,5,binw,icol,E)
         if (oab.eq.3.0d0) flux = read_histo(oab30n,300,5,binw,icol,E)
      elseif (v1.eq.-1) then
         if (oab.eq.0d0) flux = read_histo(oab00a,300,5,binw,icol,E)
         if (oab.eq.0.5d0) flux = read_histo(oab05a,300,5,binw,icol,E)
         if (oab.eq.0.6d0) flux = read_histo(oab06a,300,5,binw,icol,E)
         if (oab.eq.0.9d0) flux = read_histo(oab09a,300,5,binw,icol,E)
         if (oab.eq.1.0d0) flux = read_histo(oab10a,300,5,binw,icol,E)
         if (oab.eq.1.4d0) flux = read_histo(oab14a,300,5,binw,icol,E)
         if (oab.eq.1.5d0) flux = read_histo(oab15a,300,5,binw,icol,E)
         if (oab.eq.2.0d0) flux = read_histo(oab20a,300,5,binw,icol,E)
         if (oab.eq.2.5d0) flux = read_histo(oab25a,300,5,binw,icol,E)
         if (oab.eq.3.0d0) flux = read_histo(oab30a,300,5,binw,icol,E)
c         if (mod(oab,3).lt.0.1) flux = read_histo(oab30a,300,5,binw,icol,E)
c      elseif (abs(v1).eq.3) then
c         stop
      endif
      flux = flux/0.05d0  ! convert event number in 50 MeV bin to flux [1/GeV/cm^2/10^21POT]

      flux_beam = flux*(1d3/L)**2

      return
      end


      real*8 function xsec_IBD(E)
C     ****************************************************
C     By Yoshitaro Takaesu @KIAS Aug.17 2012
C     
C     Cross section between electron neutrino and detector in cm^2 unit.
C     Ref: hep-ph/9903554 eq.11
C
C     input
C           E: energy of the positron in MeV           
C     ****************************************************
      implicitnone
C
C     CONSTANTS
C
      real*8 me,mp,mn
C     
C     ARGUMENTS 
C     
      real*8 E,Ee
C     ----------
C     BEGIN CODE
C     ----------
      me = 0.510998910d0
      mp = 938.272d0
      mn = 939.565d0

      Ee = (E +mp)/2d0*(1d0 -(mn**2 -me**2)/(2*E*mp +mp**2))
c      Ee = E -(mn-mp)
      xsec_IBD = 0.0952d0*Ee*dsqrt(Ee**2 -me**2)*1d-42

      return
      end


      real*8 function xsec_CCQE(v2,E)
C     ****************************************************
C     By Yoshitaro Takaesu @KIAS MAY.23 2013
C     
C     The Effective Charged-Current QE Neutrino Cross Section 
C     with a water-filled 1kton detector
C     ****************************************************
      implicitnone
      integer oi
      include 'xsec/xsecCCqeN.inc'
C     CONSTANTS
C     ARGUMENTS 
      integer v2
      real*8 E
C     LOCALVARIABLES
      integer icol
      real*8 read_tbl,read_histo,binw
      external read_tbl,read_histo
C     ----------
C     BEGIN CODE
C     ----------
      if ((abs(v2).gt.3).or.(abs(v2).eq.0)) then      
         write(*,*) "ERROR: xsec_CCQE in funcs.f: Invalid detect mode"
      endif
      binw = 1d0

      if (v2.eq.1) icol = 2
      if (v2.eq.-1) icol = 3
      if (v2.eq.2) icol = 4
      if (v2.eq.-2) icol = 5
      if (v2.eq.3) icol = 6
      if (v2.eq.-3) icol = 7

c      xsec_CCQE = read_tbl(xsecCCqeN,300,7,icol,E)
      xsec_CCQE = read_histo(xsecCCqeN,300,7,binw,icol,E)/100d0

      return
      end


      real*8 function xsec_CC(v2,E)
C     ****************************************************
C     By Yoshitaro Takaesu @KIAS AUG.12 2013
C     
C     The Effective Charged-Current Neutrino Cross Section 
C     with a water-filled 1kton detector
C     ****************************************************
      implicitnone
      integer oi
      include 'xsec/xsecCCtotN.inc'
C     CONSTANTS
C     ARGUMENTS 
      integer v2
      real*8 E
C     LOCALVARIABLES
      integer icol
      real*8 read_tbl,read_histo,binw
      external read_tbl,read_histo
C     ----------
C     BEGIN CODE
C     ----------
      if ((abs(v2).gt.3).or.(abs(v2).eq.0)) then      
         write(*,*) "ERROR: xsec_CC in funcs.f: Invalid detect mode"
      endif
      binw = 1d0

      if (v2.eq.1) icol = 2
      if (v2.eq.-1) icol = 3
      if (v2.eq.2) icol = 4
      if (v2.eq.-2) icol = 5
      if (v2.eq.3) icol = 6
      if (v2.eq.-3) icol = 7

      xsec_CC = read_histo(xsecCCtotN,300,7,binw,icol,E)/100d0

      return
      end


      real*8 function xsec_NC(v2,E)
C     ****************************************************
C     By Yoshitaro Takaesu @KIAS MAY.23 2013
C
C     The Effective NC Neutrino Cross Section 
C     with a water-filled 1kton detector
C     ****************************************************
      implicitnone
      integer oi
      include 'xsec/xsecNCtotN.inc'
C     CONSTANTS
C     ARGUMENTS 
      integer v2
      real*8 E
C     LOCALVARIABLES
      integer icol
      real*8 read_tbl,read_histo,binw
      external read_tbl,read_histo
C     ----------
C     BEGIN CODE
C     ----------
      if ((abs(v2).gt.3).or.(abs(v2).eq.0)) then      
         write(*,*) "ERROR: xsec_NC in funcs.f: Invalid detect mode"
      endif
      binw = 1d0

      if (v2.eq.1) icol = 2
      if (v2.eq.-1) icol = 3
      if (v2.eq.2) icol = 4
      if (v2.eq.-2) icol = 5
      if (v2.eq.3) icol = 4
      if (v2.eq.-3) icol = 5

      xsec_NC = read_histo(xsecNCtotN,300,7,binw,icol,E)/100d0

      return
      end


      real*8 function read_tbl(tbl,nraw,ncol,icol,E)
C     ****************************************************
C     By Yoshitaro Takaesu @KIAS May.21 2013
C     ****************************************************
      implicitnone
C     ARGUMENTS 
      integer nraw,ncol,icol
      real*8 E,tbl(nraw,ncol)
C     LOCAL VARIABLES
      integer i
      real*8 value,xmin,xmax,ymin,ymax
C     ----------
C     BEGIN CODE
C     ----------
      do i = 1,nraw
         if (i.eq.1) then
            xmin = 0d0
            ymin = 0d0
         else
            xmin = tbl(i-1,1)
            ymin = tbl(i-1,icol)
         endif
         xmax = tbl(i,1)
         ymax = tbl(i,icol)
         if (E.lt.xmax) then
            call interpolate(xmin,xmax,E
     &           ,ymin,ymax,value)
            exit
         endif
      enddo
      
      read_tbl = value

      return
      end


      real*8 function prob_ee(a,param,error,sign,mode,unc_mode)
C     ****************************************************
C     By Yoshitaro Takaesu @KIAS Aug.17 2012
C     
C     P_ee
C
C     input
C           a: L/E in km/MeV
C     ****************************************************
      implicitnone
C
C     CONSTANTS
C
      real*8 s2sun_2,s23_2,s213_2,s12,c12,s13,c13
      real*8 s12_2,s13_2,unc_s2sun_2,unc_s213_2
      real*8 s2sun_2_eff,s213_2_eff
      real*8 dm13_2,dm12_2,dm23_2,unc_dm12_2,unc_dm13_2
      real*8 dm12_2_eff,dm13_2_eff,dm23_2_eff
      real*8 ue1,ue2,ue3,ue1ue2,ue1ue3,ue2ue3
      real*8 dim_fact
C     
C     ARGUMENTS 
C     
      integer sign,mode,unc_mode,icheck
      real*8 a,aa,param(4),error(4)
C     ----------
C     BEGIN CODE
C     ----------
      s2sun_2 = param(1)
      s213_2 = param(2)
      dm12_2 = param(3)
      dm13_2 = param(4)
      unc_s2sun_2 = error(1)
      unc_s213_2 = error(2)
      unc_dm12_2 = error(3)
      unc_dm13_2 = error(4)
      dm23_2 = sign*dm13_2 -dm12_2

      s2sun_2_eff = s2sun_2 +unc_mode*unc_s2sun_2
      s213_2_eff = s213_2 +unc_mode*unc_s213_2
      dm12_2_eff = dm12_2 
      dm13_2_eff = dm13_2 
      dm23_2_eff = dm23_2

      s13_2 = 0.5*s213_2_eff/( 1d0 +dsqrt(1d0 -s213_2_eff) )
      s13 = dsqrt(s13_2)
      c13 = dsqrt(1d0 -s13**2)
c      s12_2 = ( 1d0 -dsqrt(1d0 -s2sun_2_eff/c13**4) )/2d0
      s12_2 = ( 1d0 -dsqrt(1d0 -s2sun_2_eff) )/2d0
      s12 = dsqrt(s12_2)
      c12 = dsqrt(1d0 -s12**2)
      ue1 = c12*c13
      ue2 = s12*c13
      ue3 = s13

      dim_fact = 1d6/197.3269631d0
      aa = a*dim_fact

      if (mode.eq.0) then
         prob_ee = 1d0 -4*ue1**2*ue2**2*dsin(dm12_2_eff*aa/4d0)**2
     &        -4*ue1**2*ue3**2*dsin(dm13_2_eff*aa/4d0)**2
     &        -4*ue2**2*ue3**2*dsin(dm23_2_eff*aa/4d0)**2
c         Prob_ee = 1d0 -4*c13**4*s12**2*c12**2*dsin(dm12_2*a*1270d0)**2
c     &        -c12**2*s213_2*dsin(dm13_2*a*1270d0)**2
c     &        -s12**2*s213_2*dsin(dm23_2*a*1270d0)**2
      elseif (mode.eq.21) then
         prob_ee = 1d0 -4*ue1**2*ue2**2*dsin(dm12_2_eff*aa/4d0)**2
      elseif (mode.eq.31) then
         prob_ee = -4*ue1**2*ue3**2*dsin(dm13_2_eff*aa/4d0)**2
      elseif (mode.eq.32) then
         prob_ee = -4*ue2**2*ue3**2*dsin(dm23_2_eff*aa/4d0)**2
      endif

      return
      end


      real*8 function prob_approx1(v1,v2,E,L,s212_2,s223_2,s213_2,dm21_2
     &     ,dm31_2,dCP,rho)
C     ****************************************************
C     Approximated Oscillation Probability
C     eq.22 in 1209.2763

C     By Yoshitaro Takaesu @KIAS May.17 2013
C     ****************************************************
      implicitnone
C     CONSTANTS
C     ARGUMENTS
      integer v1,v2,beam_sign,detect_sign
      real*8 E,L,s212_2,s223_2,s213_2,dm21_2,dm31_2,dCP,rho
C     LOCAL VARIABLES
      real*8 rho2,s2sol_2,s2atm_2,s2rct_2,ssol_2,satm_2,srct_2
      real*8 catm_2,a,Del12,Del13,Amu,Bmu,csol_2,tatm_2
      real*8 Ae,Be,Ce,ssol,satm,srct,crct_2,csol,catm,crct,tatm      
      real*8 c2sol_2,c2atm_2,c2rct_2,c2sol,c2atm,c2rct
      real*8 s2sol,s2atm,s2rct,dCP2

      prob_approx1 = 0d0

      beam_sign = isign(1,v1)
      detect_sign = isign(1,v2) 
      if (beam_sign.ne.detect_sign) then
         return
      endif

      s2sol_2 = s212_2
      s2atm_2 = s223_2
      s2rct_2 = s213_2
      rho2 = beam_sign*rho ! [g/cm^3]
      dCP2 = beam_sign*dCP

      s2sol = dsqrt(s2sol_2)
      s2atm = dsqrt(s2atm_2)
      s2rct = dsqrt(s2rct_2)
      ssol_2 = 0.5*(1d0 -dsqrt( 1d0 -s2sol_2 ))
      satm_2 = 0.5*(1d0 -dsqrt( 1d0 -s2atm_2 ))
      srct_2 = 0.5*(1d0 -dsqrt( 1d0 -s2rct_2 ))
      ssol = dsqrt(ssol_2)
      satm = dsqrt(satm_2)
      srct = dsqrt(srct_2)
      
      c2sol_2 = 1d0 -s2sol_2
      c2atm_2 = 1d0 -s2atm_2
      c2rct_2 = 1d0 -s2rct_2
      c2sol = dsqrt(c2sol_2)
      c2atm = dsqrt(c2atm_2)
      c2rct = dsqrt(c2rct_2)
      csol_2 = 1d0 -ssol_2
      catm_2 = 1d0 -satm_2
      crct_2 = 1d0 -srct_2
      csol = dsqrt(csol_2)
      catm = dsqrt(catm_2)
      crct = dsqrt(crct_2)

      tatm_2 = satm_2/catm_2
      tatm = satm/catm
c      a = 2*dsqrt(2d0)*GF*E*ne
      a = 7.56d-5*rho2*E 
      Del12 = 2.534*dm21_2/E*L
      Del13 = 2.534*dm31_2/E*L

      Amu = -a*L/(Del13*E)*(1d0 -2*satm_2)/catm_2*srct_2

      Bmu = a*L/(4*E)*(1d0 -2*satm_2)/catm_2*srct_2 
     &     -Del12/2d0*( csol_2 +tatm_2*ssol_2*srct_2
     &     -dsqrt(tatm_2)*dsqrt(s2sol_2)*dsqrt(srct_2)*dcos(dCP2) )

      Ae = a*L/(Del13*E)*c2rct
     &     -Del12/2d0*s2sol/(tatm*srct)*dsin(dCP2)*( 1d0 +a*L/(2*Del13*E) )
     &     +Del12/4d0*( Del12 +a*L/(2*E) )
     &      *( s2sol/(tatm*srct)*dcos(dCP2) -2*ssol_2 )
     &     -1d0/2d0*(a*L/(2*E))**2 +3d0/4d0*(a*L/(Del13*E))**2

      Be = -a*L/(4*E) +Del12/4d0*( s2sol/(tatm*srct)*dcos(dCP2) -2*ssol_2 )
     &      *( 1d0 +a*L/(2*Del13*E) )
     &     +Del12/8d0*( Del12 +a*L/(2*E) )*s2sol/(tatm*srct)*dsin(dCP2)
     &     -1d0/Del13*(a*L/(2*E))**2

      Ce = (Del12)**2/4d0*s2sol_2*catm_2 
     &     -Del12/2d0*a*L/(2*E)*s2sol*s2atm*srct*dcos(dCP2)
     &     +(a*L/(2*E))**2*srct_2*satm_2
      
      if ( (abs(v1).eq.2).and.(abs(v2).eq.1) ) then
         prob_approx1 = 4*satm_2*crct_2*srct_2*( ( 1d0 +Ae )
     &        *dsin(Del13/2d0)**2 +Be*dsin(Del13) ) +Ce
      elseif ( (abs(v1).eq.2).and.(abs(v2).eq.2) ) then
         prob_approx1 = 1d0 -s2atm_2*(1d0 +Amu)*dsin( Del13/2d0 +Bmu )**2
      endif

      return
      end


      real*8 function prob_approx2(v1,v2,E,L,s212_2,s223_2,s213_2,dm21_2
     &     ,dm31_2,dCP,rho)
C     ****************************************************
C     Approximated Oscillation Probability
C     eq.29 in 1107.5857

C     By Yoshitaro Takaesu @KIAS May.17 2013
C     ****************************************************
      implicitnone
C     CONSTANTS
      include 'inc/const.inc'
C     ARGUMENTS
      integer v1,v2,beam_sign,detect_sign
      real*8 E,L,s212_2,s223_2,s213_2,dm21_2,dm31_2,dCP,rho
C     LOCAL VARIABLES
      real*8 rho2,s2sol_2,s2atm_2,s2rct_2,ssol_2,satm_2,srct_2
      real*8 catm_2,a,Del12,Del13,Amu,Bmu,csol_2,tatm_2
      real*8 Ae,Be,Ce,ssol,satm,srct,crct_2,csol,catm,crct,tatm      
      real*8 c2sol_2,c2atm_2,c2rct_2,c2sol,c2atm,c2rct
      real*8 s2sol,s2atm,s2rct,dCP2,y

      prob_approx2 = 0d0

      beam_sign = isign(1,v1)
      detect_sign = isign(1,v2) 
      if (beam_sign.ne.detect_sign) then
         return
      endif

      if ((abs(v1).ne.2).or.(abs(v2).ne.1)) then
         write(*,*) "ERROR:prob_approx2:",
     &        "Only P(nu_mu -> nu_e) is supported now."
         return
      endif

      s2sol_2 = s212_2
      s2atm_2 = s223_2
      s2rct_2 = s213_2
      rho2 = beam_sign*rho ! [g/cm^3]
      dCP2 = beam_sign*dCP

      s2sol = dsqrt(s2sol_2)
      s2atm = dsqrt(s2atm_2)
      s2rct = dsqrt(s2rct_2)
      ssol_2 = 0.5*(1d0 -dsqrt( 1d0 -s2sol_2 ))
      satm_2 = 0.5*(1d0 -dsqrt( 1d0 -s2atm_2 ))
      srct_2 = 0.5*(1d0 -dsqrt( 1d0 -s2rct_2 ))
      ssol = dsqrt(ssol_2)
      satm = dsqrt(satm_2)
      srct = dsqrt(srct_2)
      
      c2sol_2 = 1d0 -s2sol_2
      c2atm_2 = 1d0 -s2atm_2
      c2rct_2 = 1d0 -s2rct_2
      c2sol = dsqrt(c2sol_2)
      c2atm = dsqrt(c2atm_2)
      c2rct = dsqrt(c2rct_2)
      csol_2 = 1d0 -ssol_2
      catm_2 = 1d0 -satm_2
      crct_2 = 1d0 -srct_2
      csol = dsqrt(csol_2)
      catm = dsqrt(catm_2)
      crct = dsqrt(crct_2)

      tatm_2 = satm_2/catm_2
      tatm = satm/catm
c      a = 2*dsqrt(2d0)*GF*E*ne
      a = 7.56d-5*rho2*E 
      Del12 = 2.534*dm21_2/E*L
      Del13 = 2.534*dm31_2/E*L
      y = srct_2/0.04d0

c      Amu = -a*L/(Del13*E)*(1d0 -2*satm_2)/catm_2*srct_2

c      Bmu = a*L/(4*E)*(1d0 -2*satm_2)/catm_2*srct_2 
c     &     -Del12/2d0*( csol_2 +tatm_2*ssol_2*srct_2
c     &     -dsqrt(tatm_2)*dsqrt(s2sol_2)*dsqrt(srct_2)*dcos(dCP2) )

      Ae = 0.11d0*pi/Del13*L/295d0*rho2/3d0*(1d0 -0.08*y) 
     &     -0.014*(L/295d0*rho2/3d0)**2*(1d0 -0.61*(pi/Del13)**2)
     &     -0.29*dsin(dCP2)/dsqrt(y)*abs(Del13)/pi
     &     *(1d0 +0.054*pi/Del13*L/295d0*rho2/3d0)
     &     +0.014*(dcos(dCP2)/dsqrt(y) -0.11d0)*(Del13/pi)**2
     &     *(1d0 +1.7*pi/abs(Del13)*L/295d0*rho/3d0)

      Be = -0.084*L/295d0*rho2/3d0
     &     *(1d0 -0.08*y +0.11*pi/Del13*L/295d0*rho2/3d0)
     &     +0.14*(dcos(dCP2)/dsqrt(y) -0.12d0)*abs(Del13)/pi
     &     *(1d0 +0.054*pi/Del13*L/295d0*rho2/3d0)
     &     +0.0072*dsin(dCP2)/dsqrt(y)*(Del13/pi)**2
     &     *(1d0 +1.7*pi/abs(Del13)*L/295d0*rho2/3d0)

      Ce = 0.0011*(Del13/pi)**2 
     &     -0.0013*dsqrt(y)*dcos(dCP2)*abs(Del13)/pi*L/295d0
     &     +0.00036*y*(L/295d0)**2
      
      if ( (abs(v1).eq.2).and.(abs(v2).eq.1) ) then
         prob_approx2 = 4*satm_2*crct_2*srct_2*( ( 1d0 +Ae )
     &        *dsin(Del13/2d0)**2 +Be*dsin(Del13) ) +Ce
c      elseif ( (abs(v1).eq.2).and.(abs(v2).eq.2) ) then
c         prob_approx2 = 1d0 -s2atm_2*(1d0 +Amu)*dsin( Del13/2d0 +Bmu )**2
      endif

      return
      end


      real*8 function prob(v1,v2,E,L,s212_2,s223_2,s213_2,dm21_2,dm31_2
     &     ,dCP,rho,oct)
C     ****************************************************
C     By Yoshitaro Takaesu @KIAS May.17 2013
C     ****************************************************
      implicitnone
      include 'inc/const.inc'
C     CONSTANTS
C     ARGUMENTS
      integer v1,v2,beam_sign,detect_sign,nn,sort,oct
      real*8 E,L,s212_2,s223_2,s213_2,dm21_2,dm31_2,dCP,rho
C     LOCAL VARIABLES
      integer i,j,k
      real*8 rho2,a,lam(3),lam2(3)
      real*8 s12_2,s23_2,s13_2,s212,s223,s213,s12,s23,s13
      real*8 c212_2,c223_2,c213_2,c12_2,c23_2,c13_2,c212,c223,c213
      real*8 c12,c23,c13,t23_2,t23,dCP2,D21,D31,D32
      complex*16 U(3,3),Uinv(3,3),H(3,3),Ut(3,3),Utinv(3,3)
C     EXTERNAL FUNCTIONS
      real*8 sumre,sumim,Delt
      external sumre,sumim,Delt

      prob = 0d0

      if (E.le.0) then
         return
      endif
      if (v1.gt.0) then
         beam_sign = 1
      elseif (v1.lt.0) then
         beam_sign = -1
      endif
      if (v2.gt.0) then
         detect_sign = 1
      elseif (v2.lt.0) then
         detect_sign = -1
      endif
      if (beam_sign.ne.detect_sign) then
         return
      endif

      rho2 = beam_sign*rho ! [g/cm^3]
      dCP2 = beam_sign*dCP

      s212 = dsqrt(s212_2)
      s223 = dsqrt(s223_2)
      s213 = dsqrt(s213_2)
      s12_2 = 0.5*(1d0 -dsqrt( 1d0 -s212_2 ))
      if (oct.eq.1) then
         s23_2 = 0.5*(1d0 +dsqrt( 1d0 -s223_2 ))
      elseif (oct.eq.-1) then
         s23_2 = 0.5*(1d0 -dsqrt( 1d0 -s223_2 ))
      endif
      s13_2 = 0.5*(1d0 -dsqrt( 1d0 -s213_2 ))
      s12 = dsqrt(s12_2)
      s23 = dsqrt(s23_2)
      s13 = dsqrt(s13_2)
      
      c212_2 = 1d0 -s212_2
      c223_2 = 1d0 -s223_2
      c213_2 = 1d0 -s213_2
      c212 = dsqrt(c212_2)
      c223 = dsqrt(c223_2)
      c213 = dsqrt(c213_2)
      c12_2 = 1d0 -s12_2
      c23_2 = 1d0 -s23_2
      c13_2 = 1d0 -s13_2
      c12 = dsqrt(c12_2)
      c23 = dsqrt(c23_2)
      c13 = dsqrt(c13_2)

      t23_2 = s23_2/c23_2
      t23 = s23/c23
c      a = 2*dsqrt(2d0)*GF*E*ne
      a = 7.56d-5*rho2*E 
      
      U(1,1) = c12*c13
      U(1,2) = s12*c13
      U(1,3) = s13*cdexp(-ci*dCP2)
      U(2,1) = -s12*c23 -c12*s23*s13*cdexp(ci*dCP2)
      U(2,2) = c12*c23 -s12*s23*s13*cdexp(ci*dCP2)
      U(2,3) = s23*c13
      U(3,1) = s12*s23 -c12*c23*s13*cdexp(ci*dCP2)
      U(3,2) = -c12*s23 -s12*c23*s13*cdexp(ci*dCP2)
      U(3,3) = c23*c13
      
      do i = 1,3
         do j = 1,3
            Uinv(i,j) = (0d0,0d0)
            H(i,j) = (0d0,0d0)
            Utinv(i,j) = (0d0,0d0)
            Ut(i,j) = (0d0,0d0)
         enddo
      enddo   

      Uinv(1,1) = dconjg(U(1,1))
      Uinv(1,2) = dconjg(U(2,1))
      Uinv(1,3) = dconjg(U(3,1))
      Uinv(2,1) = dconjg(U(1,2))
      Uinv(2,2) = dconjg(U(2,2))
      Uinv(2,3) = dconjg(U(3,2))
      Uinv(3,1) = dconjg(U(1,3))
      Uinv(3,2) = dconjg(U(2,3))
      Uinv(3,3) = dconjg(U(3,3))

      H(1,1) = dm21_2*Uinv(2,1)*U(1,2)+dm31_2*Uinv(3,1)*U(1,3) +a
      H(1,2) = dm21_2*Uinv(2,2)*U(1,2) 
     &                    +dm31_2*Uinv(3,2)*U(1,3)
      H(1,3) = dm21_2*Uinv(2,3)*U(1,2) 
     &                    +dm31_2*Uinv(3,3)*U(1,3)
      H(2,1) = dm21_2*Uinv(2,1)*U(2,2) 
     &                    +dm31_2*Uinv(3,1)*U(2,3)
      H(2,2) = dm21_2*Uinv(2,2)*U(2,2) 
     &                    +dm31_2*Uinv(3,2)*U(2,3)
      H(2,3) = dm21_2*Uinv(2,3)*U(2,2) 
     &                    +dm31_2*Uinv(3,3)*U(2,3)
      H(3,1) = dm21_2*Uinv(2,1)*U(3,2) 
     &                    +dm31_2*Uinv(3,1)*U(3,3)
      H(3,2) = dm21_2*Uinv(2,2)*U(3,2) 
     &                    +dm31_2*Uinv(3,2)*U(3,3)
      H(3,3) = dm21_2*Uinv(2,3)*U(3,2) 
     &                    +dm31_2*Uinv(3,3)*U(3,3)

      nn = 3
      sort = 0
      call HEigensystem(nn,H,nn,lam,Utinv,nn,sort)

      Ut(1,1) = dconjg(Utinv(1,1))
      Ut(1,2) = dconjg(Utinv(2,1))
      Ut(1,3) = dconjg(Utinv(3,1))
      Ut(2,1) = dconjg(Utinv(1,2))
      Ut(2,2) = dconjg(Utinv(2,2))
      Ut(2,3) = dconjg(Utinv(3,2))
      Ut(3,1) = dconjg(Utinv(1,3))
      Ut(3,2) = dconjg(Utinv(2,3))
      Ut(3,3) = dconjg(Utinv(3,3))

      if (v1.eq.v2) then
         prob = 1d0 -4*sumre(abs(v1),abs(v2),nn,lam,Ut,E,L)       
     &        -2*sumim(abs(v1),abs(v2),nn,lam,Ut,E,L)
      else
         prob = -4*sumre(abs(v1),abs(v2),nn,lam,Ut,E,L)       
     &        -2*sumim(abs(v1),abs(v2),nn,lam,Ut,E,L)
      endif

      return
      end


      real*8 function sumre(a,b,n,lam,U,E,L)
C     ****************************************************
C     By Yoshitaro Takaesu @KIAS May.25 2013
C     ****************************************************
      implicitnone
      include 'inc/const.inc'
C     CONSTANTS
C     ARGUMENTS
      integer a,b,n
      real*8 lam(n),E,L
      complex*16 U(n,n)
C     LOCAL VARIABLES
      integer i,j
C     EXTERNAL FUNCTIONS
      real*8 Delt
      external Delt

      sumre = 0d0
      do j = 1,n-1
         do i = j+1,n
           sumre = sumre 
     &           +dble(dconjg(U(a,i))*U(b,i)*U(a,j)*dconjg(U(b,j)))
     &           *dsin( Delt(i,j,n,lam,E,L)/2d0 )**2
         enddo
      enddo

      return
      end


      real*8 function sumim(a,b,n,lam,U,E,L)
C     ****************************************************
C     By Yoshitaro Takaesu @KIAS May.25 2013
C     ****************************************************
      implicitnone
      include 'inc/const.inc'
C     CONSTANTS
C     ARGUMENTS
      integer a,b,n
      real*8 lam(n),E,L
      complex*16 U(n,n)
C     LOCAL VARIABLES
      integer i,j
C     EXTERNAL FUNCTIONS
      real*8 Delt
      external Delt

      sumim = 0d0
      do j = 1,n-1
         do i = j+1,n
           sumim = sumim
     &           +dimag(dconjg(U(a,i))*U(b,i)*U(a,j)*dconjg(U(b,j)))
     &           *dsin( Delt(i,j,n,lam,E,L) )
         enddo
      enddo

      return
      end


      real*8 function Delt(ii,jj,n,lam,E,L)
C     ****************************************************
C     By Yoshitaro Takaesu @KIAS May.25 2013
C
C     E [GeV], L [km]
C     ****************************************************
      implicitnone
      include 'inc/const.inc'
C     CONSTANTS
C     ARGUMENTS
      integer ii,jj,n
      real*8 lam(n),E,L
C     LOCAL VARIABLES
      
      Delt = (lam(jj) -lam(ii))/(2*E)*L*1d0/hbarc

      return
      end
