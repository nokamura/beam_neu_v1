      subroutine get_xsecfrac3(E,icc,ixsec,iproc,inutype,frac)
C     *******************************************************************
C     By Yoshitaro Takaesu @KIAS Mar 14 2014
C     For Separated CC and NC xsec folder
C     Arguments
C        E = incoming neutrino energy
C        icc = 1:CC 2:NC
C        ixsec = 0:NC 1pi0 1:NC pi0-bg (polfit) 2:NC pi0-bg (old)
C        iproc (icc=1) = 1:CCQE-H 2:CCQE-O 3:Res-H 4:Res-O
C        iproc (icc=2) = 1:NCQE 2:NCRes 3:NCDI 4:NCCoh+NCDiff 5:NC Total
C        inutype = incoming neutrino flavor
C     Output
C        frac
C     *******************************************************************
      implicit none
C     GLOBAL VARIABLES
      integer oi
      include 'xsecCC/nm_ratioCC.inc'
      include 'xsecCC/ne_ratioCC.inc'
      include 'xsecCC/am_ratioCC.inc'
      include 'xsecCC/ae_ratioCC.inc'
      include 'xsecNC/nm_ratioNC.inc'
      include 'xsecNC/am_ratioNC.inc'
c      include 'xsecNC/nm_ratioNC_1pi0.inc'
c      include 'xsecNC/am_ratioNC_1pi0.inc'
      include 'xsecNC/nm_ratioNC_old.inc'
      include 'xsecNC/am_ratioNC_old.inc'
C     CONSTANTS
C     ARGUMENTS 
      integer icc,iproc,inutype,ixsec
      real*8 E,frac
C     LOCAL VARIABLES 
      integer i,j
      integer icol,ncol
      real*8 binw,xsec_mode,xsec_total
C     EXTERNAL FUNCTIONS
      real*8 read_histo
      external read_histo
C     ----------
C     BEGIN CODE
C     ----------
      binw = 1d0
      icol = iproc +1
      frac = 0d0

      if (icc.eq.2) then ! NC interactions
         if (ixsec.eq.0) then ! NC 1pi0cut
            if (inutype.gt.0) then ! neutrino
c               frac = read_histo(nm_ratioNC_1pi0,120,10,binw,icol,E) 
               continue
            elseif (inutype.lt.0) then ! anti-nu
c               frac = read_histo(am_ratioNC_1pi0,120,10,binw,icol,E)
               continue
            endif
         elseif (ixsec.eq.1) then ! NC pi0-bg (polfit)
            if (inutype.gt.0) then ! nu
               frac = read_histo(nm_ratioNC,120,6,binw,icol,E)
            elseif (inutype.lt.0) then ! anti-nu
               frac = read_histo(am_ratioNC,120,6,binw,icol,E)
            endif
         elseif (ixsec.eq.2) then ! NC pi0-bg (old)
c            icol = 6
            if (inutype.gt.0) then ! nu
c               frac = read_histo(nm_ratioNC_old,120,10,binw,icol,E)
               frac = read_histo(nm_ratioNC_old,120,6,binw,icol,E)
            elseif (inutype.lt.0) then ! anti-nu
c               frac = read_histo(am_ratioNC_old,120,10,binw,icol,E)
               frac = read_histo(am_ratioNC_old,120,6,binw,icol,E)
            endif
         endif
      elseif (icc.eq.1) then ! CC interactions
            if (inutype.eq.1) then ! nu_e
               frac = read_histo(ne_ratioCC,120,5,binw,icol,E)
            elseif (inutype.eq.2) then ! nu_mu
               frac = read_histo(nm_ratioCC,120,5,binw,icol,E)
            elseif (inutype.eq.-1) then
               frac = read_histo(ae_ratioCC,120,5,binw,icol,E)
            elseif (inutype.eq.-2) then
               frac = read_histo(am_ratioCC,120,5,binw,icol,E)
            endif
      endif

      return
      end
