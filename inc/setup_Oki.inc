C -*- fortran -*-
               iD = 2
               yymin(1,iD) = 0.4d0 
               yymax(1,iD) = 2.4d0 
               yymin(2,iD) = 0.4d0
               yymax(2,iD) = 5d0
               binsize = 0.2

               z_dat(116) = 653d0 !L
               z_dat(112) = 100d0 !V
               if (OAB.eq.2) then
                  z_dat(119) = 0.6d0 !OAB
               elseif (OAB.eq.2.5) then
                  z_dat(119) = 0.9d0 !OAB
               elseif (OAB.eq.3) then
                  z_dat(119) = 1.4d0 !OAB
               else 
                  write(*,*) "ERROR:setup_Oki.inc: invalid OAB"
                  write(*,*) "OAB: 2/2.5/3"
                  write(*,*) "Exiting Program..."
                  stop
               endif
               z_dat(121) = 2.75d0 !rho    
               z_dat(122) = iD !detector flag    

               z(116) = z_dat(116) !L
               z(112) = z_dat(112) !V
               z(119) = z_dat(119) !OAB
               z(121) = z_dat(121) !rho    
               z(122) = z_dat(122) !rho    
