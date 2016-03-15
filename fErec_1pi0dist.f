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

      nutype = int(z(1)) ! incoming neutrino type
      Ev = z(2)*1d3 ! incoming nu energy in MeV
      icol = int(z(3)) +1 ! interaction mode 1:NCQE, 2:NCRS, 3:NCDI, 4:NCCO+NCDF
      E = x ! reconstructed nu energy in GeV
      
      if (nutype.gt.0) then
         fErec_1pi0dist = Erec_1pi0dist_nm(icol,Ev,E)
      elseif (nutype.lt.0) then
         fErec_1pi0dist = Erec_1pi0dist_am(icol,Ev,E)
      endif

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
c      ncol = 2
      ncol = 6
c      icol = 2
      binw = 1d0 ! table is for dist probability, not for events/bin_width
      include 'inc/read_pi0bginc_nm1100.inc'
      Erec_1pi0dist_nm = f
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
c      ncol = 2
      ncol = 6
c      icol = 2
      binw = 1d0 ! table is for dist probability, not for events/bin_width
      include 'inc/read_pi0bginc_am1100.inc'
      Erec_1pi0dist_am = f
      return
      end
