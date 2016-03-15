C -*- fortran -*-
            if (ifluc.eq.0) then
               do i = 1,nbins
                  sErec(i) = event(i)
               enddo
            elseif (ifluc.eq.1) then
               do i = 1,nbins
                  sErec(i) = event(i) +dsqrt(event(i))*gran()
                  if (event(i).lt.0) sErec(i) = 0
               enddo
            endif
            snevent = nevent
