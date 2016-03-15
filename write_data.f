      subroutine write_data(itype)
C     ****************************************************
C     By Yoshitaro Takaesu @KIAS Mar 23 2014
C     ****************************************************
      implicit none
C     GLOBAL VARIABLES
      include 'inc/params.inc'
      include 'inc/main.inc'
C     CONSTANTS
C     ARGUMENTS
      integer itype
C     LOCAL VARIABLES 
      integer i
      integer ib,iD,inu,idet,iev 
C     EXTERNAL FUNCTIONS
C     ----------
C     BEGIN CODE
C     ----------
      do ib = -1,1,2
         do iD = 1,3
            do inu = -3,3
               do idet = -3,3
                  do iev = 1,2
                     if (inu*idet.gt.0) then
                        call write_data_unit(itype,ib,iD,inu,idet,iev)
                     endif
                  enddo
               enddo
            enddo
         enddo
      enddo

CCC  For Oscillation Probability output
c$$$      do ib = 1,1,2
c$$$         do iD = 3,3
c$$$            do inu = -2,2
c$$$               do idet = -2,2
c$$$                  do iev = 1,2
c$$$                     if ((inu*idet.gt.0).and.(iev.eq.abs(idet))) then
c$$$                        call write_data_unit(itype,ib,iD,inu,idet,iev)
c$$$                     endif
c$$$                  enddo
c$$$               enddo
c$$$            enddo
c$$$         enddo
c$$$      enddo
         
      return
      end

      subroutine write_data_unit(itype,ib,iD,inu,idet,iev)
C     ****************************************************
C     By Yoshitaro Takaesu @KIAS Mar 23 2014
C     ****************************************************
      implicit none
C     GLOBAL VARIABLES
c      include 'inc/params.inc'
c      include 'inc/main.inc'
C     CONSTANTS
C     ARGUMENTS
      integer itype,ib,iD,inu,idet,iev 
C     LOCAL VARIABLES 
      character*2 cb,cD,cnu,cdet
      character*11 file_mid
      character*15 file_app
c      integer i
C     EXTERNAL FUNCTIONS
C     ----------
C     BEGIN CODE
C     ----------
      if (inu*idet.eq.0) goto 100

      if (ib.eq.1) cb = "nb"
      if (ib.eq.-1) cb = "ab"
      if (iD.eq.1) cD = "sk"                  
      if (iD.eq.2) cD = "ok"
      if (iD.eq.3) cD = "kr"
      if (inu.eq.1) cnu = "ne"
      if (inu.eq.2) cnu = "nm"
      if (inu.eq.3) cnu = "nt"
      if (inu.eq.-1) cnu = "ae"
      if (inu.eq.-2) cnu = "am"
      if (inu.eq.-3) cnu = "at"
      if (idet.eq.1) cdet = "ne"
      if (idet.eq.2) cdet = "nm"
      if (idet.eq.3) cdet = "nt"
      if (idet.eq.-1) cdet = "ae"
      if (idet.eq.-2) cdet = "am"
      if (idet.eq.-3) cdet = "at"
      file_mid = cb//"."//cD//"."//cnu//"."//cdet

      if (itype.eq.1) then
         file_app = "dat_"//file_mid
      elseif (itype.eq.2) then
         file_app = "the_"//file_mid
      endif
      if (iev.eq.1) then
         open(1,file=file_app//".e2e.dat",status="replace")
      elseif (iev.eq.2) then
         open(1,file=file_app//".m2m.dat",status="replace")
      elseif (iev.eq.3) then
         open(1,file=file_app//".pi0.dat",status="replace")
      elseif (iev.eq.4) then
         open(1,file=file_app//".gam.dat",status="replace")
      elseif (iev.eq.5) then
         open(1,file=file_app//".m2e.dat",status="replace")
      elseif (iev.eq.6) then
         open(1,file=file_app//".e2m.dat",status="replace")
      endif
      call write_data_atm(1,itype,ib,iD,inu,idet,iev)
      close(1)        

 100  return
      end


      subroutine write_data_atm(ilun,itype,ib,iD,inu,idet,iev)
C     ****************************************************
C     By Yoshitaro Takaesu @KIAS Mar 23 2014
C     ****************************************************
      implicit none
C     GLOBAL VARIABLES
      include 'inc/params.inc'
      include 'inc/main.inc'
C     CONSTANTS
C     ARGUMENTS
      integer ilun,itype,ib,iD,inu,idet,iev
C     LOCAL VARIABLES 
      integer i,j
      real*8 ans(imaxint),tot
C     EXTERNAL FUNCTIONS
C     ----------
C     BEGIN CODE
C     ----------
      write(ilun,*) "# E_nu, total, CCQE-H, CCQE-O,"
     &     ,"CCRes-H, CCRes-O for CC events" 
      write(ilun,*) "# E_nu, total, NCQE, NCRS,"
     &     ," NCDI, NCCO+NCDF for NC events"
      do i = 0,nbins-1
         tot = 0d0
         do j = 1,imaxint
            ans(j) = 0d0
         enddo
         if (itype.eq.1) then
            do j = 1,imaxint
               ans(j) = ans(j) +hErec_dat(i+1,j,iev,idet,inu,iD,ib)
               tot = tot +ans(j)
            enddo
         elseif (itype.eq.2) then
            do j = 1,imaxint
               ans(j) = ans(j) +hErec_th(i+1,j,iev,idet,inu,iD,ib)
               tot = tot +ans(j)
            enddo
         endif         
         write(ilun,'(f5.3,f10.5,f10.5,f10.5,f10.5,f10.5)') 
     &        x(i),tot,ans(1),ans(2),ans(3),ans(4)
      enddo

      return
      end
