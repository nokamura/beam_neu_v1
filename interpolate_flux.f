      subroutine interpolate_flux(neu,oab)
C     ****************************************************
C     By Yoshitaro Takaesu @Okayama U July.29 2016
C     
C     Make the flux table for a OAB by interpolating two 
C     flux tables for the larger and smaller OABs.
C     ****************************************************
      implicitnone    
C     GLOBAL VARIABLES
      integer oi
      include 'beam/oab00n.inc'
      include 'beam/oab05n.inc'
      include 'beam/oab10n.inc'
      include 'beam/oab15n.inc'
      include 'beam/oab20n.inc'
      include 'beam/oab25n.inc'
      include 'beam/oab30n.inc'
      include 'beam/oab00a.inc'
      include 'beam/oab05a.inc'
      include 'beam/oab10a.inc'
      include 'beam/oab15a.inc'
      include 'beam/oab20a.inc'
      include 'beam/oab25a.inc'
      include 'beam/oab30a.inc'
C     ARGUMENTS
      integer neu
      real*8 oab
C     LOCAL VARIABLES
      integer rows,cols,i,j
      parameter (rows = 300, cols=5)
      real*8 oabs(rows,cols),oabl(rows,cols),oabout(rows,cols)
      real*8 factor
C     ----------
C     BEGIN CODE
C     ----------
      if (neu.eq.1) then
         if ((oab.ge.0d0).and.(oab.lt.0.5d0)) then
            oabs = oab00n
            oabl = oab05n
            factor = oab/0.5d0
         elseif ((oab.ge.0.5d0).and.(oab.lt.1.0d0)) then
            oabs = oab05n
            oabl = oab10n
            factor = (oab -0.5d0)/0.5d0
         elseif ((oab.ge.1.0d0).and.(oab.lt.1.5d0)) then
            oabs = oab10n
            oabl = oab15n
            factor = (oab -1.0d0)/0.5d0
         elseif ((oab.ge.1.5d0).and.(oab.lt.2.0d0)) then
            oabs = oab15n
            oabl = oab20n
            factor = (oab -1.5d0)/0.5d0
         elseif ((oab.ge.2.0d0).and.(oab.lt.2.5d0)) then
            oabs = oab20n
            oabl = oab25n
            factor = (oab -2.0d0)/0.5d0
         elseif ((oab.ge.2.5d0).and.(oab.lt.3.0d0)) then
            oabs = oab25n
            oabl = oab30n
            factor = (oab -2.5d0)/0.5d0
         else
            write(*,*) "ERROR in interpolate_flux:"
     &           ,"invalid arugument oab: 0 <= oab <= 3.0"
         endif
      elseif (neu.eq.-1) then
         if ((oab.ge.0d0).and.(oab.lt.0.5d0)) then
            oabs = oab00a
            oabl = oab05a
            factor = (oab -0.0d0)/0.5d0
         elseif ((oab.ge.0.5d0).and.(oab.lt.1.0d0)) then
            oabs = oab05a
            oabl = oab10a
            factor = (oab -0.5d0)/0.5d0
         elseif ((oab.ge.1.0d0).and.(oab.lt.1.5d0)) then
            oabs = oab10a
            oabl = oab15a
            factor = (oab -1.0d0)/0.5d0
         elseif ((oab.ge.1.5d0).and.(oab.lt.2.0d0)) then
            oabs = oab15a
            oabl = oab20a
            factor = (oab -1.5d0)/0.5d0
         elseif ((oab.ge.2.0d0).and.(oab.lt.2.5d0)) then
            oabs = oab20a
            oabl = oab25a
            factor = (oab -2.0d0)/0.5d0
         elseif ((oab.ge.2.5d0).and.(oab.lt.3.0d0)) then
            oabs = oab25a
            oabl = oab30a
            factor = (oab -2.5d0)/0.5d0
         else
            write(*,*) "ERROR in interpolate_flux:"
     &           ,"invalid arugument oab: 0 <= oab <= 3.0"
         endif
      else
         write(*,*) "ERROR: interplate_flux: invalid argument neu."
     &        ," neu = 1 or -1"
      endif

      do i = 1,rows
         oabout(i,1) = oabs(i,1)
         do j = 2,5
            oabout(i,j) = factor*(oabl(i,j) -oabs(i,j)) +oabs(i,j)
         enddo
      enddo

      open(1,file="oabout.inc",status="replace") 
      write(1,'(a,a)') "C/*--------------------------------------------"
     &     ,"---------"
      write(1,'(a)') "C// Neutrino flux table at 1000 km"
      if (neu.eq.1) then
         write(1,'(a,a)') "C// Off-axis angle is OOAB deg."
     &        ," for nu_mu focusing beam"
      elseif (neu.eq.-1) then
         write(1,'(a,a)') "C// Off-axis angle is OOAB deg."
     &        ," for anti-nu_mu focusing beam"
      endif
      write(1,'(a)') "C// Enu(GeV), nu_mu, bar{nu_mu}, nu_e, bar{nu_e}"
      write(1,'(a,a)') "C----------------------------------------------"
     &     ,"---------*"
      write(1,'(6x,a,i0,a,i0,a)') "real*8 oabout(",rows,",",cols,")"
      do i = 1,rows
         write(1,'(6x,a,i0,a,4(e13.7,a),e13.7,a)') "data (oabout(",i
     &        ,",oi),oi=1,5) / ",oabout(i,1),", ",oabout(i,2),", "
     &        ,oabout(i,3),", ",oabout(i,4),", ",oabout(i,5)," /"
      enddo
      close(1)

      return
      end
