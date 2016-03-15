C -*- fortran -*-
      iD = 1
      yymin(1,iD) = 0.4d0 
      yymax(1,iD) = 1.2d0 
      yymin(2,iD) = 0.4d0
      yymax(2,iD) = 5d0
      binsize = 0.2
      
      z_dat(116) = 295d0        !L
      if ((iOki.eq.0).and.(iKorea.eq.0)) then
         z_dat(112) = 122.5d0   !V
      else
         z_dat(112) = 22.5d0    !V
      endif
      z_dat(119) = OAB          !OAB
      z_dat(121) = 2.60d0       !rho    
      z_dat(122) = iD
      
      z(116) = z_dat(116)       !L
      z(112) = z_dat(112)       !V
      z(119) = z_dat(119)       !OAB
      z(121) = z_dat(121)       !rho    
      z(122) = z_dat(122) 
