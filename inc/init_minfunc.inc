C -*- fortran -*-
      do i = 1,zdim
         z_dat(i) = 0d0
      enddo

      z_dat(1) = s2sol_2
      z_dat(2) = thatm
      z_dat(3) = s2rct_2
      z_dat(4) = dm21_2
      z_dat(5) = dm31_2
      z_dat(6) = dCP
      z_dat(7) = fCCQEn(1)
      z_dat(8) = fCCQEa(1)
      z_dat(9) = fS(1)
      z_dat(10) = fSrho(1)
      z_dat(11) = fSnmn(1)
      z_dat(12) = fSnen(1)
      z_dat(13) = fSnma(1)
      z_dat(14) = fSnea(1)
      z_dat(15) = fSamn(1)
      z_dat(16) = fSaen(1)
      z_dat(17) = fSama(1)
      z_dat(18) = fSaea(1)
      z_dat(19) = fO(1)
      z_dat(20) = fOrho(1)
      z_dat(21) = fOnmn(1)
      z_dat(22) = fOnen(1)
      z_dat(23) = fOnma(1)
      z_dat(24) = fOnea(1)
      z_dat(25) = fOamn(1)
      z_dat(26) = fOaen(1)
      z_dat(27) = fOama(1)
      z_dat(28) = fOaea(1)
      z_dat(29) = fK(1)
      z_dat(30) = fKrho(1)
      z_dat(31) = fKnmn(1)
      z_dat(32) = fKnen(1)
      z_dat(33) = fKnma(1)
      z_dat(34) = fKnea(1)
      z_dat(35) = fKamn(1)
      z_dat(36) = fKaen(1)
      z_dat(37) = fKama(1)
      z_dat(38) = fKaea(1)
      error(1) = fs2sol_2(2)
      error(2) = fthatm(2)
      error(3) = fs2rct_2(2)
      error(4) = fdm21_2(2)
      error(5) = fdm31_2(2)
      error(6) = fdCP(2)
      error(7) = fCCQEn(2)
      error(8) = fCCQEa(2)
      error(9) = fS(2)
      error(10) = fSrho(2)
      error(11) = fSnmn(2)
      error(12) = fSnen(2)
      error(13) = fSnma(2)
      error(14) = fSnea(2)
      error(15) = fSamn(2)
      error(16) = fSaen(2)
      error(17) = fSama(2)
      error(18) = fSaea(2)
      error(19) = fO(2)
      error(20) = fOrho(2)
      error(21) = fOnmn(2)
      error(22) = fOnen(2)
      error(23) = fOnma(2)
      error(24) = fOnea(2)
      error(25) = fOamn(2)
      error(26) = fOaen(2)
      error(27) = fOama(2)
      error(28) = fOaea(2)
      error(29) = fK(2)
      error(30) = fKrho(2)
      error(31) = fKnmn(2)
      error(32) = fKnen(2)
      error(33) = fKnma(2)
      error(34) = fKnea(2)
      error(35) = fKamn(2)
      error(36) = fKaen(2)
      error(37) = fKama(2)
      error(38) = fKaea(2)

      z_dat(111) = MH                 ! NH/IH
      z_dat(112) = V                  ! Volume [kton]
      z_dat(113) = P                  ! Power [GW]
      z_dat(114) = Y                  ! Exposure time [year]
      z_dat(115) = ihfunc  ! hfunc1D mode
      z_dat(116) = L                  ! L [km]
      z_dat(117) = beam               ! beam mode
      z_dat(118) = detect             ! detected neutrino type
      z_dat(119) = OAB                ! Off Axis Beam Angle @SK
      z_dat(120) = nu_mode             ! nutrino in the beam
      z_dat(121) = rho                ! matter density [g/cm^3]
      z_dat(122) = 0  
      z_dat(123) = iinput 
      z_dat(124) = ithatm

      z(111) = hypo*z_dat(111)              
      z(112) = z_dat(112)
      z(113) = z_dat(113)
      z(114) = z_dat(114)
      z(115) = z_dat(115)
      z(116) = z_dat(116)
      z(117) = z_dat(117)
      z(118) = z_dat(118)
      z(119) = z_dat(119)
      z(120) = z_dat(120)
      z(121) = z_dat(121)
      z(122) = z_dat(122)
      z(123) = z_dat(123)
      z(124) = z_dat(124)
      rnevent_ren = -1  ! -1:No Normalize events other:Normalize events
      nnbins = 1000
      evform_dat = 2
      evform_th = 2      
      ipi0dist = 0 ! 0:NC pi0bg-Erec dist 1:NC pi0-mom dist

      do i = 1,ndetect
         do j = 1,eventdim
            do k = 1,maxnbin
               total_event_th(k,j,i) = 0d0
            enddo
         enddo
      enddo
      do i = 1,ndetect
         do j = 1,eventdim
            nevent_th(j,i) = 0d0
         enddo
      enddo
      do i = 1,ndetect
         do j = -3,3
            do k = -3,3
               do l1 = 1,nmode
                  do l2 = 1,maxnbin
                     hErec_dat(l2,l1,k,j,i) = 0d0
                     hErec_th(l2,l1,k,j,i) = 0d0
                  enddo
               enddo
            enddo
         enddo
      enddo

      iD = 1
      nm = 1
      iev = 1
      ymin = Emin
      ymax = Emax
      include 'bining_x.inc'
c     include 'inc/bining_sqrtx.inc'
c     include 'inc/bining_sqrtx_Eres.inc'
      evnbins(iev,iD) = nbins
      do i = 0,nbins
         evx(i,iev,iD) = x(i)
      enddo
      evform = evform_dat
      
      if (mode.eq.1) then !for minimization
         rnevent_ren = -1  ! -1:No Normalize events other:Normalize events
         do i = 1,nparx
            parflag(i) = 1
         enddo
         parflag(2) = 0
         parflag(5) = 0
         parflag(6) = 0
      endif
