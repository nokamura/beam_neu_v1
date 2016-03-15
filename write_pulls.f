      subroutine write_pulls
C     ****************************************************
C     By Yoshitaro Takaesu @U.Tokyo Jan 27 2015
C
C     Output MINUIT results to files or the std_output
C     ****************************************************
      implicit none
C     GLOBAL VARIABLES
      include 'inc/params.inc'
      include 'inc/main.inc'
C     CONSTANTS
C     ARGUMENTS 
C     LOCAL VARIABLES 
      integer i,iflag,ierr
      real*8 name(maxnparam),pval(maxnparam),perr(maxnparam)
      real*8 plo(maxnparam),phi(maxnparam)
C     EXTERNAL FUNCTIONS
      real*8 minfunc
      external minfunc
C     ----------
C     BEGIN CODE
C     ----------
      open(1,file='pulls.dat',status='replace')
      write(1,*) pull_fact(1),pull_fact(2),pull_fact(3),pull_fact(4)
     &     ,pull_fact(5)      
      close(1)
      open(1,file='pulls_all.dat',status='replace')
      write(1,*) pull_fact(1),pull_fact(2),pull_fact(3),pull_fact(4)
     &     ,pull_fact(5),pull_fact(6),pull_fact(7),pull_fact(8)
     &     ,pull_fact(9),pull_fact(10),pull_fact(11),pull_fact(12)
     &     ,pull_fact(13),pull_fact(14),pull_fact(15),pull_fact(16)
     &     ,pull_fact(17),pull_fact(18),pull_fact(19),pull_fact(20)
     &     ,pull_fact(21),pull_fact(22),pull_fact(23),pull_fact(24)
     &     ,pull_fact(25),pull_fact(26),pull_fact(27),pull_fact(28)
     &     ,pull_fact(29),pull_fact(30),pull_fact(31),pull_fact(32)
     &     ,pull_fact(33),pull_fact(34),pull_fact(35),pull_fact(36)
     &     ,pull_fact(37),pull_fact(38),pull_fact(39),pull_fact(40)
     &     ,pull_fact(41),pull_fact(42),pull_fact(43),pull_fact(44)
     &     ,pull_fact(45),pull_fact(46),pull_fact(47),pull_fact(48)
     &     ,pull_fact(49),pull_fact(50),pull_fact(51),pull_fact(52)
     &     ,pull_fact(53),pull_fact(54),pull_fact(55),pull_fact(56)
      close(1)
      
      return
      end
