      real*8 function hfunc1D(x,z)
      implicit none
      include 'inc/const.inc'
      include 'inc/params.inc'
      include 'inc/minfunc.inc'
      integer i
      integer Ediv,ierr,oct_atm,oct_23,ierror,ierror2
      real*8 x,z(zdim),E,YY,fxsec_CCQE,fxsec_CCRes,fxsec_NCpi0
      real*8 s212_2,s223_2,s213_2,hdm21_2,hdm31_2,hdCP,hhdCP
      real*8 xmin,xmax,s2atm2,Enumin,Enumax,Enu,nuevent,ff
      real*8 hrho,xsec,fxsec,hfCCQEn,hfCCQEa,hfCCResn,hfCCResa,hfpi0
      real*8 fD,frho,fnmn,fnen,fnma,fnea,famn,faen,fama,faea
      real*8 frac,frac1,frac2,frac3,frac4,corr_fact
      real*8 flux_beam,xsec_NC,xsec_CC,xsec_CCQE,prob
      real*8 prob_approx1,prob_approx2,prob_approx3,smear_CCQE_mu
      real*8 flux
      external flux_beam,xsec_NC,xsec_CC,xsec_CCQE,prob
      external prob_approx1,prob_approx2,prob_approx3,smear_CCQE_mu      

c     oscillation parameter check
      if (z(1).gt.1) then ! sin(2*th_12)^2
         hfunc1D = 0d0
         return
      elseif (z(1).lt.0) then
         hfunc1D = 0d0
         return
      endif
      if (z(2).gt.1) then ! sin(2*th_23)^2 or sin(th_23)^2
         hfunc1D = 0d0
         return
      elseif (z(2).lt.0) then
         hfunc1D = 0d0
         return
      endif
      if (z(3).gt.1) then ! sin(2*th_13)^2
         hfunc1D = 0d0
         return
      elseif (z(3).lt.0) then
         hfunc1D = 0d0
         return
      endif
      If (iinput.eq.1) then
         if (ithatm.eq.1) then
            oct_atm = ioct 
            call convert_th(z(1),z(2),z(3),oct_atm,s212_2,s223_2,s213_2
     &           ,oct_23,ierror2)
            if (ierror2.eq.1) then
               hfunc1D = 0d0
               return
            endif
         elseif (ithatm.eq.2) then
            s2atm2 = 4*z(2)*(1d0 -z(2)) 
            if (z(2).le.0.5d0) then
               oct_atm = -1 ! This is an octant for sin(th_atm)^2
            elseif (z(2).gt.0.5d0) then
               oct_atm = 1 ! This is an octant for sin(th_atm)^2
            endif
            call convert_th(z(1),s2atm2,z(3),oct_atm,s212_2,s223_2
     &           ,s213_2,oct_23,ierror2)
            if (ierror2.eq.1) then
               hfunc1D = 0d0
               return
            endif
         else
            write(*,*) "ERROR:hfunc1D:Invalid ithatm (1/2)"
            write(*,*) "Stopping program..."
            stop
         endif
      elseif (iinput.eq.2) then
         s212_2 = z(1)
         if (ithatm.eq.1) then
            s223_2 = z(2)
            oct_23 = ioct
         elseif (ithatm.eq.2) then
            s223_2 = 4*z(2)*(1d0 -z(2))
            if (z(2).le.0.5d0) then
               oct_23 = -1
            elseif (z(2).gt.0.5d0) then
               oct_23 = 1
            endif
         else
            write(*,*) "ERROR:hfunc1D:Invalid ithatm (1/2)"
            write(*,*) "Stopping program..."
            stop
         endif
         s213_2 = z(3)
      else
         write(*,*) "ERROR:hfunc1D:Invalid iinput (1/2)"
         write(*,*) "Stopping program..."
         stop
      endif

      if (z(4).ge.0) then
         hdm21_2 = z(4)
      elseif (z(4).lt.0) then
         hfunc1D = 0d0
         return
      endif
      if (z(5).lt.0) then ! |delta_m_31^2| or |delta_m_32^2|
         hfunc1D = 0d0
         return
      endif
      if (idmatm2.eq.1) then
         hdm31_2 = MHH*z(5) 
      elseif (idmatm2.eq.2) then
         hdm31_2 = MHH*z(5) +hdm21_2 
      endif

      hdCP = z(6)/360d0*2*pi ! convert degree to radian
      if (z(7).ge.0) then
         hfCCQEn = z(7)
      elseif (z(7).lt.0) then
         hfunc1D = 0d0
         return
      endif
      if (z(8).ge.0) then
         hfCCQEa = z(8)
      elseif (z(8).lt.0) then
         hfunc1D = 0d0
         return
      endif
      if (z(51).ge.0) then
         hfCCResn = z(51)
      elseif (z(51).lt.0) then
         hfunc1D = 0d0
         return
      endif
      if (z(52).ge.0) then
         hfCCResa = z(52)
      elseif (z(52).lt.0) then
         hfunc1D = 0d0
         return
      endif
      if (z(53).ge.0) then
         hfpi0 = z(53)
      elseif (z(53).lt.0) then
         hfunc1D = 0d0
         return
      endif
      if (iD.eq.1) then
         fD = z(9)
         frho = z(10)
         fnmn = z(11)
         fnen = z(12)
         fnma = z(13)
         fnea = z(14)
         famn = z(15)
         faen = z(16)
         fama = z(17)
         faea = z(18)
         effe = z(19)
         effmu = z(20)
         Pe2m = abs(z(21))
         Pm2e = abs(z(22))
      elseif (iD.eq.2) then
         fD = z(23)
         frho = z(24)
         fnmn = z(25)
         fnen = z(26)
         fnma = z(27)
         fnea = z(28)
         famn = z(29)
         faen = z(30)
         fama = z(31)
         faea = z(32)
         effe = z(33)
         effmu = z(34)
         Pe2m = abs(z(35))
         Pm2e = abs(z(36))
      elseif (iD.eq.3) then
         fD = z(37)
         frho = z(38)
         fnmn = z(39)
         fnen = z(40)
         fnma = z(41)
         fnea = z(42)
         famn = z(43)
         faen = z(44)
         fama = z(45)
         faea = z(46)
         effe = z(47)
         effmu = z(48)
         Pe2m = abs(z(49))
         Pm2e = abs(z(50))
c$$$         effe = z(19) ! for re-evaluation paper
c$$$         effmu = z(20) ! for re-evaluation paper
c$$$         Pe2m = abs(z(21)) ! for re-evaluation paper
c$$$         Pm2e = abs(z(22)) ! for re-evaluation paper
      endif      
      ierror = 0
      if (fD.lt.0d0) ierror = 1
      if (frho.lt.0d0) ierror = 1
      if (fnmn.lt.0d0) ierror = 1
      if (fnen.lt.0d0) ierror = 1
      if (fnma.lt.0d0) ierror = 1
      if (fnea.lt.0d0) ierror = 1
      if (famn.lt.0d0) ierror = 1
      if (faen.lt.0d0) ierror = 1
      if (fama.lt.0d0) ierror = 1
      if (faea.lt.0d0) ierror = 1
      if ((effe.gt.2d0).or.(effe.lt.0d0)) then
         ierror = 1
      elseif (effe.gt.1d0) then
         effe = 2d0 -effe
      endif
      if ((effmu.gt.2d0).or.(effe.lt.0d0)) then
         ierror = 1
      elseif (effmu.gt.1d0) then
         effmu = 2d0 -effmu
      endif
      if ((Pe2m.gt.1d0).or.(Pe2m.lt.0d0)) ierror = 1
      if ((Pm2e.gt.1d0).or.(Pm2e.lt.0d0)) ierror = 1
      if (z(53).lt.0d0) ierror = 1
      if (z(54).lt.0d0) ierror = 1
      if (z(55).lt.0d0) ierror = 1
      if (z(56).lt.0d0) ierror = 1
      if (ierror.eq.1) then
         hfunc1D = 0d0
         return
      endif

      if (detect.gt.0) then
         fxsec_CCQE = hfCCQEn
         fxsec_CCRes = hfCCResn
      elseif (detect.lt.0) then
         fxsec_CCQE = hfCCQEa
         fxsec_CCRes = hfCCResa
      endif
      fxsec_NCpi0 = hfpi0
      if (beam.gt.0) then
         if (nu_mode.eq.1) ff = fnen
         if (nu_mode.eq.-1) ff = fnea
         if (nu_mode.eq.2) ff = fnmn
         if (nu_mode.eq.-2) ff = fnma
      elseif (beam.lt.0) then
         if (nu_mode.eq.1) ff = faen
         if (nu_mode.eq.-1) ff = faea
         if (nu_mode.eq.2) ff = famn
         if (nu_mode.eq.-2) ff = fama
      endif
      
      if ((r_nu.ge.0).and.(r_anu.ge.0)) then
         if ((r_nu.gt.0).or.(r_anu.gt.0)) then
            if (beam.eq.1) then
               YY = Y*r_nu/(r_nu + r_anu)
            elseif (beam.eq.-1) then
               YY = Y*r_anu/(r_nu + r_anu)
            endif
         else
            write(*,*) "ERROR: r_nu, r_anu = 0"
         endif
      else
         write(*,*) "ERROR: r_nu or r_anu are negative"
      endif

      E = x
      if (icc.eq.1) then
CCC     Fluxes
         flux = flux_beam(beam,oab,nu_mode,E,L)
CCC     Cross sections
         if (ismear.eq.0) then
            frac = 0d0
            if (iCCQE.eq.1) then
               call get_xsecfrac3(E,icc,1,1,detect,frac1)
               call get_xsecfrac3(E,icc,1,2,detect,frac2)
               frac = frac +fxsec_CCQE*(frac1 +frac2)
            endif
            if (iCCRes.eq.1) then
               call get_xsecfrac3(E,icc,1,3,detect,frac3)
               call get_xsecfrac3(E,icc,1,4,detect,frac4)
               frac = frac +fxsec_CCRes*(frac3 +frac4)
            endif
c            xsec = xsec_CC(detect,E)*frac  ! for re-evaluation paper
c            xsec = fxsec_CCQE*xsec_CCQE(detect,E)
c            xsec = xsec_CCQE(detect,E)
            xsec = xsec_CCQE(detect,E)*frac
         elseif (ismear.eq.1) then
            xsec = xsec_CC(detect,E) ! for re-evaluation paper
c            xsec = xsec_CCQE(detect,E)
         endif
      elseif (icc.eq.2) then
CCC     Fluxes ! Total flux is used for NC events   
         if (nu_mode.eq.1) then ! nu_flux
            if (beam.eq.1) then
               flux = fnen*flux_beam(beam,oab,1,E,L)
     &               +fnmn*flux_beam(beam,oab,2,E,L)
            elseif (beam.eq.-1) then
               flux = faen*flux_beam(beam,oab,1,E,L)
     &               +famn*flux_beam(beam,oab,2,E,L)
            endif
         elseif (nu_mode.eq.-1) then ! anti-nu flux
            if (beam.eq.1) then
               flux = fnea*flux_beam(beam,oab,-1,E,L)
     &              +fnma*flux_beam(beam,oab,-2,E,L)
            elseif (beam.eq.-1) then
               flux = faea*flux_beam(beam,oab,-1,E,L)
     &              +fama*flux_beam(beam,oab,-2,E,L)
            endif
         endif
c         flux = flux_beam(beam,oab,nu_mode,E,L)
CCC     Cross sections
         xsec = fxsec_NCpi0*xsec_NC(detect,E)
      endif
c      xsec_CCQE(detect,E)/(100*1d9*avog)
c      xsec = xsec_CCQE(detect,E)
c      xsec = xsec_CC(detect,E)/(33.6d30*8) ! nu_e-H2O xsec/neutron
c      xsec = xsec_CC(detect,E)/(33.6d30*10) ! bar{nu}_e-H2O xsec/proton
c      xsec = xsec_CC(detect,E)/(3.34d31*8)*4.89d33/18d0 ! nu_e-H2O xsec/neutron scaled to 1kton LS
c      xsec = xsec_CC(detect,E)/(3.34d31*10)*6.24d33/18d0 ! bar{nu}_e-H2O xsec/proton scaled to 1kton LS
      hrho = frho*rho
      if (ihfunc.eq.0) then
         if (icc.eq.1) then
            hfunc1D = ff*flux*xsec*fD*V*YY
     &           *prob(nu_mode,detect,E,L,s212_2,s223_2,s213_2,hdm21_2
     &           ,hdm31_2,hdCP,hrho,oct_23)
c     &           *prob_approx3(nu_mode,detect,E,L,s212_2,s223_2,s213_2
c     &           ,hdm21_2,hdm31_2,hdCP,hrho,oct_23)
         elseif (icc.eq.2) then
            hfunc1D = flux*xsec*fD*V*YY
c            hfunc1D = ff*flux*xsec*fD*V*YY
c     &           *prob(nu_mode,detect,E,L,s212_2,s223_2,s213_2,hdm21_2
c     &           ,hdm31_2,hdCP,hrho,oct_23)
         endif
      elseif (ihfunc.eq.1) then
         hfunc1D = flux
      elseif (ihfunc.eq.2) then
         hfunc1D = prob(nu_mode,detect,E,L,s212_2,s223_2,s213_2,hdm21_2
     &        ,hdm31_2,hdCP,hrho,oct_23)
      elseif (ihfunc.eq.3) then
         hfunc1D = xsec
      elseif (ihfunc.eq.4) then
         hfunc1D = flux*xsec*V*YY
     &        *prob(nu_mode,detect,E,L,s212_2,s223_2,s213_2,hdm21_2
     &        ,hdm31_2,hdCP,hrho,oct_23)
      endif

      return
      end
