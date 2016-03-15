      real*8 function fErec_1pi0dist(x,z)
c     E: incoming neutrino energy in GeV
      implicit none
      include 'inc/const.inc'
      real*8 x,z(200)
      integer int_mode,nutype,icol
      real*8 E,Ev
      real*8 Erec_1pi0dist_ae,Erec_1pi0dist_am,Erec_1pi0dist_ne
      real*8 Erec_1pi0dist_nm
      external Erec_1pi0dist_ae,Erec_1pi0dist_am,Erec_1pi0dist_ne
      external Erec_1pi0dist_nm

      nutype = int(z(1)) ! incoming neutrino type: 1:nu_e 2:nu_mu 3:nu_tau -1:anti-nu_e -2: anti-nu_mu -3:anti-nu_tau 
      Ev = z(2)*1d3 ! incoming nu energy in MeV
      icol = int(z(3)) +1 ! interaction mode
      E = x ! reconstructed nu energy in GeV
      
      if (nutype.eq.1) then
c         fErec_1pi0dist = Erec_1pi0dist_ne(Ev,E)
         fErec_1pi0dist = Erec_1pi0dist_nm(icol,Ev,E)
      elseif (nutype.eq.2) then
         fErec_1pi0dist = Erec_1pi0dist_nm(icol,Ev,E)
      elseif (nutype.eq.3) then
         fErec_1pi0dist = Erec_1pi0dist_nm(icol,Ev,E)
      elseif (nutype.eq.-1) then
c         fErec_1pi0dist = Erec_1pi0dist_ae(Ev,E)
         fErec_1pi0dist = Erec_1pi0dist_am(icol,Ev,E)
      elseif (nutype.eq.-2) then
         fErec_1pi0dist = Erec_1pi0dist_am(icol,Ev,E)
      elseif (nutype.eq.-3) then
         fErec_1pi0dist = Erec_1pi0dist_am(icol,Ev,E)
      endif

      return
      end


      real*8 function Erec_1pi0dist_ne(icol,Ev,Erec)
c     Erec: Reconstructed neutrino energy [MeV]
c     Ev: incoming neutrino energy [MeV]
      implicit none
      integer oi
      include 'include_pi0bgincs_ne1100.inc'
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
c      icol = 2
      binw = 1d0 ! table is for dist probability, not for events/bin_width
      include 'inc/read_pi0bginc_ne1100.inc'
      Erec_1pi0dist_ne = f
      return
      end

      real*8 function Erec_1pi0dist_nm(icol,Ev,Erec)
c     Erec: Reconstructed neutrino energy [MeV]
c     Ev: incoming neutrino energy [MeV]
      implicit none
      integer oi
      include 'include_pi0bgincs_nm1100.inc'
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
c      icol = 2
      binw = 1d0 ! table is for dist probability, not for events/bin_width
      include 'inc/read_pi0bginc_nm1100.inc'
      Erec_1pi0dist_nm = f
      return
      end

      real*8 function Erec_1pi0dist_ae(icol,Ev,Erec)
c     Erec: Reconstructed neutrino energy [MeV]
c     Ev: incoming neutrino energy [MeV]
      implicit none
      integer oi
      include 'include_pi0bgincs_ae1100.inc'
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
c      icol = 2
      binw = 1d0 ! table is for dist probability, not for events/bin_width
      include 'inc/read_pi0bginc_ae1100.inc'
      Erec_1pi0dist_ae = f
      return
      end

      real*8 function Erec_1pi0dist_am(icol,Ev,Erec)
c     Erec: Reconstructed neutrino energy [MeV]
c     Ev: incoming neutrino energy [MeV]
      implicit none
      integer oi
      include 'include_pi0bgincs_am1100.inc'
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
c      icol = 2
      binw = 1d0 ! table is for dist probability, not for events/bin_width
      include 'inc/read_pi0bginc_am1100.inc'
      Erec_1pi0dist_am = f
      return
      end
