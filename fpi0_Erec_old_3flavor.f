      real*8 function fpi0_Erec_old(x,z)
c     E: incoming neutrino energy in GeV
      implicit none
      include 'inc/const.inc'
      real*8 x,z(200)
      integer int_mode,nutype,icol
      real*8 E,Ev
      real*8 pi0_Erec_old_ae,pi0_Erec_old_am,pi0_Erec_old_ne
      real*8 pi0_Erec_old_nm
      external pi0_Erec_old_ae,pi0_Erec_old_am,pi0_Erec_old_ne
      external pi0_Erec_old_nm

      fpi0_Erec_old = 0d0
      nutype = int(z(1)) ! incoming neutrino type: 1:nu_e 2:nu_mu 3:nu_tau -1:anti-nu_e -2: anti-nu_mu -3:anti-nu_tau 
      Ev = z(2)*1d3 ! incoming nu energy in MeV
      icol = int(z(3)) +1
      E = x ! reconstructed nu energy in GeV
      
      if (nutype.eq.1) then
         fpi0_Erec_old = pi0_Erec_old_nm(icol,Ev,E)
      elseif (nutype.eq.2) then
         fpi0_Erec_old = pi0_Erec_old_nm(icol,Ev,E)
      elseif (nutype.eq.3) then
         fpi0_Erec_old = pi0_Erec_old_nm(icol,Ev,E)
      elseif (nutype.eq.-1) then
         fpi0_Erec_old = pi0_Erec_old_am(icol,Ev,E)
      elseif (nutype.eq.-2) then
         fpi0_Erec_old = pi0_Erec_old_am(icol,Ev,E)
      elseif (nutype.eq.-3) then
         fpi0_Erec_old = pi0_Erec_old_am(icol,Ev,E)
      endif

      return
      end


      real*8 function pi0_Erec_old_ne(icol,Ev,Erec)
c     Erec: Reconstructed neutrino energy [GeV]
c     Ev: incoming neutrino energy [MeV]
      implicit none
      integer oi
      include 'include_pi0_Erec_old_ne1100.inc'
      real*8 Ev,Erec
      character*5 cEv4
      integer i
      integer nrow,ncol,icol,inum,iEv
      real*8 binw,Evmin,dEv,iEvmax,f,res
      real*8 f1,f2
      real*8 read_histo
      external read_histo
      nrow = 120
      ncol = 2
c      binw = 0.05d0
      binw = 1d0 ! table is for dist probability, not for events/bin_width
      include 'inc/read_pi0_Erec_old_ne1100.inc'
      pi0_Erec_old_ne = f
      return
      end

      real*8 function pi0_Erec_old_nm(icol,Ev,Erec)
c     Erec: Reconstructed neutrino energy [MeV]
c     Ev: incoming neutrino energy [MeV]
      implicit none
      integer oi
      include 'include_pi0_Erec_old_nm1100.inc'
      real*8 Ev,Erec
      character*5 cEv4
      integer i
      integer nrow,ncol,icol,inum,iEv
      real*8 binw,Evmin,dEv,iEvmax,f,res
      real*8 f1,f2
      real*8 read_histo
      external read_histo
      nrow = 120
      ncol = 2
c      binw = 0.05d0
      binw = 1d0 ! table is for dist probability, not for events/bin_width
      include 'inc/read_pi0_Erec_old_nm1100.inc'
      pi0_Erec_old_nm = f
      return
      end

      real*8 function pi0_Erec_old_ae(icol,Ev,Erec)
c     Erec: Reconstructed neutrino energy [MeV]
c     Ev: incoming neutrino energy [MeV]
      implicit none
      integer oi
      include 'include_pi0_Erec_old_ae1100.inc'
      real*8 Ev,Erec
      character*5 cEv4
      integer i
      integer nrow,ncol,icol,inum,iEv
      real*8 binw,Evmin,dEv,iEvmax,f,res
      real*8 f1,f2
      real*8 read_histo
      external read_histo
      nrow = 120
      ncol = 2
c      binw = 0.05d0
      binw = 1d0 ! table is for dist probability, not for events/bin_width
      include 'inc/read_pi0_Erec_old_ae1100.inc'
      pi0_Erec_old_ae = f
      return
      end

      real*8 function pi0_Erec_old_am(icol,Ev,Erec)
c     Erec: Reconstructed neutrino energy [MeV]
c     Ev: incoming neutrino energy [MeV]
      implicit none
      integer oi
      include 'include_pi0_Erec_old_am1100.inc'
      real*8 Ev,Erec
      character*5 cEv4
      integer i
      integer nrow,ncol,icol,inum,iEv
      real*8 binw,Evmin,dEv,iEvmax,f,res
      real*8 f1,f2
      real*8 read_histo
      external read_histo      
      nrow = 120
      ncol = 2
c      binw = 0.05d0
      binw = 1d0 ! table is for dist probability, not for events/bin_width
      include 'inc/read_pi0_Erec_old_am1100.inc'
      pi0_Erec_old_am = f
      return
      end
