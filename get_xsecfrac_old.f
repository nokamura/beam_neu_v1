      subroutine get_xsecfrac(E,nutype,mode,frac)
C     ****************************************************
C     By Yoshitaro Takaesu @KIAS OCT 11 2013
C     ****************************************************
      implicit none
C     GLOBAL VARIABLES
      integer oi
      include 'xsec/ne_ratioCC.inc'
      include 'xsec/nm_ratioCC.inc'
      include 'xsec/ae_ratioCC.inc'
      include 'xsec/am_ratioCC.inc'
C     CONSTANTS
C     ARGUMENTS 
      integer nutype,mode 
      real*8 E,frac
C     LOCAL VARIABLES 
      integer i,j
      integer icol
      real*8 binw,xsec_mode,xsec_total
C     EXTERNAL FUNCTIONS
      real*8 read_histo
      external read_histo
C     ----------
C     BEGIN CODE
C     ----------
      binw = 1d0
      icol = mode +1

      if (nutype.eq.1) then
         frac = read_histo(ne_ratioCC,120,5,binw,icol,E)
      elseif (nutype.eq.2) then
         frac = read_histo(nm_ratioCC,120,5,binw,icol,E)
      elseif (nutype.eq.-1) then
         frac = read_histo(ae_ratioCC,120,5,binw,icol,E)
      elseif (nutype.eq.-2) then
         frac = read_histo(am_ratioCC,120,5,binw,icol,E)
      endif

      return
      end


      subroutine get_xsecfrac2(E,icc,iproc,inutype,frac)
C     ****************************************************
C     By Yoshitaro Takaesu @KIAS OCT 11 2013
C     ****************************************************
      implicit none
C     GLOBAL VARIABLES
      integer oi
      include 'xsec/ae_ratioCC.inc'
      include 'xsec/am_ratioCC.inc'
      include 'xsec/nm_ratioNC.inc'
      include 'xsec/am_ratioNC.inc'
      include 'xsec/ne_ratioNC.inc'
      include 'xsec/ae_ratioNC.inc'
C     CONSTANTS
C     ARGUMENTS 
      integer icc,iproc,inutype
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

      if (icc.eq.2) then ! NC interactions
         if (inutype.eq.1) then ! nu_e
            frac = read_histo(ne_ratioNC,120,10,binw,icol,E)
         elseif (inutype.eq.2) then ! nu_mu
            frac = read_histo(nm_ratioNC,120,10,binw,icol,E)
         elseif (inutype.eq.-1) then ! anti-nu_
            frac = read_histo(ae_ratioNC,120,10,binw,icol,E)
         elseif (inutype.eq.-2) then ! anti-nu_mu
            frac = read_histo(am_ratioNC,120,10,binw,icol,E)
         endif
      elseif (icc.eq.1) then ! CC interactions

         if (inutype.eq.1) then ! nu_e
            frac = 0d0
         elseif (inutype.eq.2) then ! nu_mu
            frac = 0d0
         elseif (inutype.eq.-1) then
            frac = read_histo(ae_ratioCC,120,5,binw,icol,E)
         elseif (inutype.eq.-2) then
            frac = read_histo(am_ratioCC,120,5,binw,icol,E)
         endif
      endif

      return
      end


      subroutine get_xsecfrac3(E,icc,ixsec,iproc,inutype,frac)
C     ****************************************************
C     By Yoshitaro Takaesu @KIAS OCT 11 2013
C     ****************************************************
      implicit none
C     GLOBAL VARIABLES
      integer oi
      include 'xsec/nm_ratioCC.inc'
      include 'xsec/ne_ratioCC.inc'
      include 'xsec/am_ratioCC.inc'
      include 'xsec/ae_ratioCC.inc'
      include 'xsec/nm_ratioNC.inc'
      include 'xsec/ne_ratioNC.inc'
      include 'xsec/am_ratioNC.inc'
      include 'xsec/ae_ratioNC.inc'
      include 'xsec/nm_ratioNC_1pi0.inc'
      include 'xsec/ne_ratioNC_1pi0.inc'
      include 'xsec/am_ratioNC_1pi0.inc'
      include 'xsec/ae_ratioNC_1pi0.inc'
      include 'xsec/nm_ratioNC_old.inc'
      include 'xsec/ne_ratioNC_old.inc'
      include 'xsec/am_ratioNC_old.inc'
      include 'xsec/ae_ratioNC_old.inc'
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
               frac = read_histo(nm_ratioNC_1pi0,120,10,binw,icol,E) 
            elseif (inutype.lt.0) then ! anti-nu
               frac = read_histo(am_ratioNC_1pi0,120,10,binw,icol,E)
            endif
         elseif (ixsec.eq.1) then ! NC pi0-bg (polfit)
            if (inutype.gt.0) then ! nu
               frac = read_histo(nm_ratioNC,120,10,binw,icol,E)
            elseif (inutype.lt.0) then ! anti-nu
               frac = read_histo(am_ratioNC,120,10,binw,icol,E)
            endif
         elseif (ixsec.eq.2) then ! NC pi0-bg (old)
            if (inutype.gt.0) then ! nu
               frac = read_histo(nm_ratioNC_old,120,10,binw,icol,E)
            elseif (inutype.lt.0) then ! anti-nu
               frac = read_histo(am_ratioNC_old,120,10,binw,icol,E)
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


      subroutine get_xsecfrac4(E,icc,ixsec,iproc,inutype,frac)
C     ****************************************************
C     By Yoshitaro Takaesu @KIAS Mar 14 2014
C     For Separated CC and NC xsec folder
C     ****************************************************
      implicit none
C     GLOBAL VARIABLES
      integer oi
      include 'xsecCC/nm_ratioCC.inc'
      include 'xsecCC/ne_ratioCC.inc'
      include 'xsecCC/am_ratioCC.inc'
      include 'xsecCC/ae_ratioCC.inc'
      include 'xsecNC/nm_ratioNC.inc'
      include 'xsecNC/ne_ratioNC.inc'
      include 'xsecNC/am_ratioNC.inc'
      include 'xsecNC/ae_ratioNC.inc'
      include 'xsecNC/nm_ratioNC_1pi0.inc'
      include 'xsecNC/ne_ratioNC_1pi0.inc'
      include 'xsecNC/am_ratioNC_1pi0.inc'
      include 'xsecNC/ae_ratioNC_1pi0.inc'
      include 'xsecNC/nm_ratioNC_old.inc'
      include 'xsecNC/ne_ratioNC_old.inc'
      include 'xsecNC/am_ratioNC_old.inc'
      include 'xsecNC/ae_ratioNC_old.inc'
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
               frac = read_histo(nm_ratioNC_1pi0,120,10,binw,icol,E) 
            elseif (inutype.lt.0) then ! anti-nu
               frac = read_histo(am_ratioNC_1pi0,120,10,binw,icol,E)
            endif
         elseif (ixsec.eq.1) then ! NC pi0-bg (polfit)
            if (inutype.gt.0) then ! nu
               frac = read_histo(nm_ratioNC,120,10,binw,icol,E)
            elseif (inutype.lt.0) then ! anti-nu
               frac = read_histo(am_ratioNC,120,10,binw,icol,E)
            endif
         elseif (ixsec.eq.2) then ! NC pi0-bg (old)
            if (inutype.gt.0) then ! nu
               frac = read_histo(nm_ratioNC_old,120,10,binw,icol,E)
            elseif (inutype.lt.0) then ! anti-nu
               frac = read_histo(am_ratioNC_old,120,10,binw,icol,E)
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
