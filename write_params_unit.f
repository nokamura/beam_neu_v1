      subroutine write_params_unit(iout,ivar,chi2)
C     ****************************************************
C     By Yoshitaro Takaesu @KIAS Mar 7 2014
C
C     write down paramter fit info from the parameter number
C     ivar in the file specified by iout
C     ****************************************************
      implicit none
C     GLOBAL VARIABLES
      include 'inc/params.inc'
      include 'inc/main.inc'
C     CONSTANTS
C     ARGUMENTS 
      integer ivar,iout
C     LOCAL VARIABLES
      character*10 chnam
      integer ivarbl
      real*8 val,er,bnd1,bnd2,pull,chi2
C     EXTERNAL FUNCTIONS
C     ----------
C     BEGIN CODE
C     ----------
      call mnpout(ivar,chnam,val,er,bnd1,bnd2,ivarbl)
      if (ivarbl.gt.0) then         
         pull = (val -z_dat(ivar))/error(ivar)
         chi2 = pull**2
         write(iout,*) chnam,ivar,z_dat(ivar),error(ivar),val,er,pull
     &        ,chi2
      endif

      return
      end
