C -*- fortran -*-
               iD = 3
               yymin(1,iD) = 0d0 
               yymax(1,iD) = 2.8d0 
c               yymax(1,iD) = 5d0 
               yymin(2,iD) = 0d0
               yymax(2,iD) = 5d0
               binsize = 0.2

c               z_dat(116) = 1250d0 !L
c               z_dat(112) = 10d0 !V
               z_dat(116) = 1000d0 !L
               z_dat(112) = 100d0 !V
               if (OAB.eq.2) then
                  z_dat(119) = 1.5d0 !OAB
               elseif (OAB.eq.2.5) then
                  z_dat(119) = 1.0d0 !OAB
               elseif (OAB.eq.3) then
                  z_dat(119) = 0.5d0 !OAB
c                  z_dat(119) = 3d0 !OAB
               else 
                  write(*,*) "ERROR:setup_Oki.inc: invalid OAB"
                  write(*,*) "OAB: 2/2.5/3"
                  write(*,*) "Exiting Program..."
                  stop
               endif
               z_dat(121) = 2.90d0 !rho    
               z_dat(122) = iD !detector flag    

               z(116) = z_dat(116) !L
               z(112) = z_dat(112) !V
               z(119) = z_dat(119) !OAB
               z(121) = z_dat(121) !rho    
               z(122) = z_dat(122) !detector flag
