      subroutine write_params
C     ****************************************************
C     By Yoshitaro Takaesu @KIAS Mar 6 2014
C
C     write down paramter fit info
C     ****************************************************
      implicit none
C     GLOBAL VARIABLES
      include 'inc/params.inc'
      include 'inc/main.inc'
C     CONSTANTS
C     ARGUMENTS 
C     LOCAL VARIABLES
      integer i
      integer iout
      real*8 pull
      character*10 chnam
      integer ivarbl
      real*8 val,er,bnd1,bnd2,chi2_tot,chi2
C     EXTERNAL FUNCTIONS
C     ----------
C     BEGIN CODE
C     ----------
      iout = 1      
      open(iout,file="fit_params.dat",status="replace")
      chi2_tot = 0d0
      do i = 1,nparx
         call write_params_unit(iout,i,chi2)
         if (i.ne.6) then
            chi2_tot = chi2_tot +chi2
         endif
      enddo
      write(iout,*) "chi2_tot =",chi2_tot
      close(iout)

      open(iout,file="fit_fpico.dat",status="replace")
      call mnpout(55,chnam,val,er,bnd1,bnd2,ivarbl)
c      pull = (val -z_dat(55))/error(55)
      pull = (val -z_dat(53))/error(53)
      write(iout,*) pull
      close(iout)

      return
      end
