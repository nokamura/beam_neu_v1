C -*- fortran -*-
      if (Ev.lt.200) then
         f = 0d0
      elseif (Ev.lt.300) then
         f1 = read_histo(pi0mom_old_nm200_1100,nrow,ncol,binw,icol,Erec)
         f2 = read_histo(pi0mom_old_nm300_1100,nrow,ncol,binw,icol,Erec)
         f = f1 +(f2 -f1)/1d2*(Ev -200d0)
      elseif (Ev.lt.400) then
         f1 = read_histo(pi0mom_old_nm300_1100,nrow,ncol,binw,icol,Erec)
         f2 = read_histo(pi0mom_old_nm400_1100,nrow,ncol,binw,icol,Erec)
         f = f1 +(f2 -f1)/1d2*(Ev -300d0)
      elseif (Ev.lt.500) then
         f1 = read_histo(pi0mom_old_nm400_1100,nrow,ncol,binw,icol,Erec)
         f2 = read_histo(pi0mom_old_nm500_1100,nrow,ncol,binw,icol,Erec)
         f = f1 +(f2 -f1)/1d2*(Ev -400d0)
      elseif (Ev.lt.600) then
         f1 = read_histo(pi0mom_old_nm500_1100,nrow,ncol,binw,icol,Erec)
         f2 = read_histo(pi0mom_old_nm600_1100,nrow,ncol,binw,icol,Erec)
         f = f1 +(f2 -f1)/1d2*(Ev -500d0)
      elseif (Ev.lt.700) then
         f1 = read_histo(pi0mom_old_nm600_1100,nrow,ncol,binw,icol,Erec)
         f2 = read_histo(pi0mom_old_nm700_1100,nrow,ncol,binw,icol,Erec)
         f = f1 +(f2 -f1)/1d2*(Ev -600d0)
      elseif (Ev.lt.800) then
         f1 = read_histo(pi0mom_old_nm700_1100,nrow,ncol,binw,icol,Erec)
         f2 = read_histo(pi0mom_old_nm800_1100,nrow,ncol,binw,icol,Erec)
         f = f1 +(f2 -f1)/1d2*(Ev -700d0)
      elseif (Ev.lt.900) then
         f1 = read_histo(pi0mom_old_nm800_1100,nrow,ncol,binw,icol,Erec)
         f2 = read_histo(pi0mom_old_nm900_1100,nrow,ncol,binw,icol,Erec)
         f = f1 +(f2 -f1)/1d2*(Ev -800d0)
      elseif (Ev.lt.1000) then
         f1 = read_histo(pi0mom_old_nm900_1100,nrow,ncol,binw,icol,Erec)
         f2 = read_histo(pi0mom_old_nm1000_1100,nrow,ncol,binw,icol,Erec)
         f = f1 +(f2 -f1)/1d2*(Ev -900d0)
      elseif (Ev.lt.1100) then
         f1 = read_histo(pi0mom_old_nm1000_1100,nrow,ncol,binw,icol,Erec)
         f2 = read_histo(pi0mom_old_nm1100_1100,nrow,ncol,binw,icol,Erec)
         f = f1 +(f2 -f1)/1d2*(Ev -1000d0)
      elseif (Ev.lt.1200) then
         f1 = read_histo(pi0mom_old_nm1100_1100,nrow,ncol,binw,icol,Erec)
         f2 = read_histo(pi0mom_old_nm1200_1100,nrow,ncol,binw,icol,Erec)
         f = f1 +(f2 -f1)/1d2*(Ev -1100d0)
      elseif (Ev.lt.1300) then
         f1 = read_histo(pi0mom_old_nm1200_1100,nrow,ncol,binw,icol,Erec)
         f2 = read_histo(pi0mom_old_nm1300_1100,nrow,ncol,binw,icol,Erec)
         f = f1 +(f2 -f1)/1d2*(Ev -1200d0)
      elseif (Ev.lt.1400) then
         f1 = read_histo(pi0mom_old_nm1300_1100,nrow,ncol,binw,icol,Erec)
         f2 = read_histo(pi0mom_old_nm1400_1100,nrow,ncol,binw,icol,Erec)
         f = f1 +(f2 -f1)/1d2*(Ev -1300d0)
      elseif (Ev.lt.1500) then
         f1 = read_histo(pi0mom_old_nm1400_1100,nrow,ncol,binw,icol,Erec)
         f2 = read_histo(pi0mom_old_nm1500_1100,nrow,ncol,binw,icol,Erec)
         f = f1 +(f2 -f1)/1d2*(Ev -1400d0)
      elseif (Ev.lt.1600) then
         f1 = read_histo(pi0mom_old_nm1500_1100,nrow,ncol,binw,icol,Erec)
         f2 = read_histo(pi0mom_old_nm1600_1100,nrow,ncol,binw,icol,Erec)
         f = f1 +(f2 -f1)/1d2*(Ev -1500d0)
      elseif (Ev.lt.1700) then
         f1 = read_histo(pi0mom_old_nm1600_1100,nrow,ncol,binw,icol,Erec)
         f2 = read_histo(pi0mom_old_nm1700_1100,nrow,ncol,binw,icol,Erec)
         f = f1 +(f2 -f1)/1d2*(Ev -1600d0)
      elseif (Ev.lt.1800) then
         f1 = read_histo(pi0mom_old_nm1700_1100,nrow,ncol,binw,icol,Erec)
         f2 = read_histo(pi0mom_old_nm1800_1100,nrow,ncol,binw,icol,Erec)
         f = f1 +(f2 -f1)/1d2*(Ev -1700d0)
      elseif (Ev.lt.1900) then
         f1 = read_histo(pi0mom_old_nm1800_1100,nrow,ncol,binw,icol,Erec)
         f2 = read_histo(pi0mom_old_nm1900_1100,nrow,ncol,binw,icol,Erec)
         f = f1 +(f2 -f1)/1d2*(Ev -1800d0)
      elseif (Ev.lt.2000) then
         f1 = read_histo(pi0mom_old_nm1900_1100,nrow,ncol,binw,icol,Erec)
         f2 = read_histo(pi0mom_old_nm2000_1100,nrow,ncol,binw,icol,Erec)
         f = f1 +(f2 -f1)/1d2*(Ev -1900d0)
      elseif (Ev.lt.2100) then
         f1 = read_histo(pi0mom_old_nm2000_1100,nrow,ncol,binw,icol,Erec)
         f2 = read_histo(pi0mom_old_nm2100_1100,nrow,ncol,binw,icol,Erec)
         f = f1 +(f2 -f1)/1d2*(Ev -2000d0)
      elseif (Ev.lt.2200) then
         f1 = read_histo(pi0mom_old_nm2100_1100,nrow,ncol,binw,icol,Erec)
         f2 = read_histo(pi0mom_old_nm2200_1100,nrow,ncol,binw,icol,Erec)
         f = f1 +(f2 -f1)/1d2*(Ev -2100d0)
      elseif (Ev.lt.2300) then
         f1 = read_histo(pi0mom_old_nm2200_1100,nrow,ncol,binw,icol,Erec)
         f2 = read_histo(pi0mom_old_nm2300_1100,nrow,ncol,binw,icol,Erec)
         f = f1 +(f2 -f1)/1d2*(Ev -2200d0)
      elseif (Ev.lt.2400) then
         f1 = read_histo(pi0mom_old_nm2300_1100,nrow,ncol,binw,icol,Erec)
         f2 = read_histo(pi0mom_old_nm2400_1100,nrow,ncol,binw,icol,Erec)
         f = f1 +(f2 -f1)/1d2*(Ev -2300d0)
      elseif (Ev.lt.2500) then
         f1 = read_histo(pi0mom_old_nm2400_1100,nrow,ncol,binw,icol,Erec)
         f2 = read_histo(pi0mom_old_nm2500_1100,nrow,ncol,binw,icol,Erec)
         f = f1 +(f2 -f1)/1d2*(Ev -2400d0)
      elseif (Ev.lt.2600) then
         f1 = read_histo(pi0mom_old_nm2500_1100,nrow,ncol,binw,icol,Erec)
         f2 = read_histo(pi0mom_old_nm2600_1100,nrow,ncol,binw,icol,Erec)
         f = f1 +(f2 -f1)/1d2*(Ev -2500d0)
      elseif (Ev.lt.2700) then
         f1 = read_histo(pi0mom_old_nm2600_1100,nrow,ncol,binw,icol,Erec)
         f2 = read_histo(pi0mom_old_nm2700_1100,nrow,ncol,binw,icol,Erec)
         f = f1 +(f2 -f1)/1d2*(Ev -2600d0)
      elseif (Ev.lt.2800) then
         f1 = read_histo(pi0mom_old_nm2700_1100,nrow,ncol,binw,icol,Erec)
         f2 = read_histo(pi0mom_old_nm2800_1100,nrow,ncol,binw,icol,Erec)
         f = f1 +(f2 -f1)/1d2*(Ev -2700d0)
      elseif (Ev.lt.2900) then
         f1 = read_histo(pi0mom_old_nm2800_1100,nrow,ncol,binw,icol,Erec)
         f2 = read_histo(pi0mom_old_nm2900_1100,nrow,ncol,binw,icol,Erec)
         f = f1 +(f2 -f1)/1d2*(Ev -2800d0)
      elseif (Ev.lt.3000) then
         f1 = read_histo(pi0mom_old_nm2900_1100,nrow,ncol,binw,icol,Erec)
         f2 = read_histo(pi0mom_old_nm3000_1100,nrow,ncol,binw,icol,Erec)
         f = f1 +(f2 -f1)/1d2*(Ev -2900d0)
      elseif (Ev.lt.3100) then
         f1 = read_histo(pi0mom_old_nm3000_1100,nrow,ncol,binw,icol,Erec)
         f2 = read_histo(pi0mom_old_nm3100_1100,nrow,ncol,binw,icol,Erec)
         f = f1 +(f2 -f1)/1d2*(Ev -3000d0)
      elseif (Ev.lt.3200) then
         f1 = read_histo(pi0mom_old_nm3100_1100,nrow,ncol,binw,icol,Erec)
         f2 = read_histo(pi0mom_old_nm3200_1100,nrow,ncol,binw,icol,Erec)
         f = f1 +(f2 -f1)/1d2*(Ev -3100d0)
      elseif (Ev.lt.3300) then
         f1 = read_histo(pi0mom_old_nm3200_1100,nrow,ncol,binw,icol,Erec)
         f2 = read_histo(pi0mom_old_nm3300_1100,nrow,ncol,binw,icol,Erec)
         f = f1 +(f2 -f1)/1d2*(Ev -3200d0)
      elseif (Ev.lt.3400) then
         f1 = read_histo(pi0mom_old_nm3300_1100,nrow,ncol,binw,icol,Erec)
         f2 = read_histo(pi0mom_old_nm3400_1100,nrow,ncol,binw,icol,Erec)
         f = f1 +(f2 -f1)/1d2*(Ev -3300d0)
      elseif (Ev.lt.3500) then
         f1 = read_histo(pi0mom_old_nm3400_1100,nrow,ncol,binw,icol,Erec)
         f2 = read_histo(pi0mom_old_nm3500_1100,nrow,ncol,binw,icol,Erec)
         f = f1 +(f2 -f1)/1d2*(Ev -3400d0)
      elseif (Ev.lt.3600) then
         f1 = read_histo(pi0mom_old_nm3500_1100,nrow,ncol,binw,icol,Erec)
         f2 = read_histo(pi0mom_old_nm3600_1100,nrow,ncol,binw,icol,Erec)
         f = f1 +(f2 -f1)/1d2*(Ev -3500d0)
      elseif (Ev.lt.3700) then
         f1 = read_histo(pi0mom_old_nm3600_1100,nrow,ncol,binw,icol,Erec)
         f2 = read_histo(pi0mom_old_nm3700_1100,nrow,ncol,binw,icol,Erec)
         f = f1 +(f2 -f1)/1d2*(Ev -3600d0)
      elseif (Ev.lt.3800) then
         f1 = read_histo(pi0mom_old_nm3700_1100,nrow,ncol,binw,icol,Erec)
         f2 = read_histo(pi0mom_old_nm3800_1100,nrow,ncol,binw,icol,Erec)
         f = f1 +(f2 -f1)/1d2*(Ev -3700d0)
      elseif (Ev.lt.3900) then
         f1 = read_histo(pi0mom_old_nm3800_1100,nrow,ncol,binw,icol,Erec)
         f2 = read_histo(pi0mom_old_nm3900_1100,nrow,ncol,binw,icol,Erec)
         f = f1 +(f2 -f1)/1d2*(Ev -3800d0)
      elseif (Ev.lt.4000) then
         f1 = read_histo(pi0mom_old_nm3900_1100,nrow,ncol,binw,icol,Erec)
         f2 = read_histo(pi0mom_old_nm4000_1100,nrow,ncol,binw,icol,Erec)
         f = f1 +(f2 -f1)/1d2*(Ev -3900d0)
      elseif (Ev.lt.4100) then
         f1 = read_histo(pi0mom_old_nm4000_1100,nrow,ncol,binw,icol,Erec)
         f2 = read_histo(pi0mom_old_nm4100_1100,nrow,ncol,binw,icol,Erec)
         f = f1 +(f2 -f1)/1d2*(Ev -4000d0)
      elseif (Ev.lt.4200) then
         f1 = read_histo(pi0mom_old_nm4100_1100,nrow,ncol,binw,icol,Erec)
         f2 = read_histo(pi0mom_old_nm4200_1100,nrow,ncol,binw,icol,Erec)
         f = f1 +(f2 -f1)/1d2*(Ev -4100d0)
      elseif (Ev.lt.4300) then
         f1 = read_histo(pi0mom_old_nm4200_1100,nrow,ncol,binw,icol,Erec)
         f2 = read_histo(pi0mom_old_nm4300_1100,nrow,ncol,binw,icol,Erec)
         f = f1 +(f2 -f1)/1d2*(Ev -4200d0)
      elseif (Ev.lt.4400) then
         f1 = read_histo(pi0mom_old_nm4300_1100,nrow,ncol,binw,icol,Erec)
         f2 = read_histo(pi0mom_old_nm4400_1100,nrow,ncol,binw,icol,Erec)
         f = f1 +(f2 -f1)/1d2*(Ev -4300d0)
      elseif (Ev.lt.4500) then
         f1 = read_histo(pi0mom_old_nm4400_1100,nrow,ncol,binw,icol,Erec)
         f2 = read_histo(pi0mom_old_nm4500_1100,nrow,ncol,binw,icol,Erec)
         f = f1 +(f2 -f1)/1d2*(Ev -4400d0)
      elseif (Ev.lt.4600) then
         f1 = read_histo(pi0mom_old_nm4500_1100,nrow,ncol,binw,icol,Erec)
         f2 = read_histo(pi0mom_old_nm4600_1100,nrow,ncol,binw,icol,Erec)
         f = f1 +(f2 -f1)/1d2*(Ev -4500d0)
      elseif (Ev.lt.4700) then
         f1 = read_histo(pi0mom_old_nm4600_1100,nrow,ncol,binw,icol,Erec)
         f2 = read_histo(pi0mom_old_nm4700_1100,nrow,ncol,binw,icol,Erec)
         f = f1 +(f2 -f1)/1d2*(Ev -4600d0)
      elseif (Ev.lt.4800) then
         f1 = read_histo(pi0mom_old_nm4700_1100,nrow,ncol,binw,icol,Erec)
         f2 = read_histo(pi0mom_old_nm4800_1100,nrow,ncol,binw,icol,Erec)
         f = f1 +(f2 -f1)/1d2*(Ev -4700d0)
      elseif (Ev.lt.4900) then
         f1 = read_histo(pi0mom_old_nm4800_1100,nrow,ncol,binw,icol,Erec)
         f2 = read_histo(pi0mom_old_nm4900_1100,nrow,ncol,binw,icol,Erec)
         f = f1 +(f2 -f1)/1d2*(Ev -4800d0)
      elseif (Ev.lt.5000) then
         f1 = read_histo(pi0mom_old_nm4900_1100,nrow,ncol,binw,icol,Erec)
         f2 = read_histo(pi0mom_old_nm5000_1100,nrow,ncol,binw,icol,Erec)
         f = f1 +(f2 -f1)/1d2*(Ev -4900d0)
      elseif (Ev.lt.5100) then
         f1 = read_histo(pi0mom_old_nm5000_1100,nrow,ncol,binw,icol,Erec)
         f2 = read_histo(pi0mom_old_nm5100_1100,nrow,ncol,binw,icol,Erec)
         f = f1 +(f2 -f1)/1d2*(Ev -5000d0)
      elseif (Ev.lt.5200) then
         f1 = read_histo(pi0mom_old_nm5100_1100,nrow,ncol,binw,icol,Erec)
         f2 = read_histo(pi0mom_old_nm5200_1100,nrow,ncol,binw,icol,Erec)
         f = f1 +(f2 -f1)/1d2*(Ev -5100d0)
      elseif (Ev.lt.5300) then
         f1 = read_histo(pi0mom_old_nm5200_1100,nrow,ncol,binw,icol,Erec)
         f2 = read_histo(pi0mom_old_nm5300_1100,nrow,ncol,binw,icol,Erec)
         f = f1 +(f2 -f1)/1d2*(Ev -5200d0)
      elseif (Ev.lt.5400) then
         f1 = read_histo(pi0mom_old_nm5300_1100,nrow,ncol,binw,icol,Erec)
         f2 = read_histo(pi0mom_old_nm5400_1100,nrow,ncol,binw,icol,Erec)
         f = f1 +(f2 -f1)/1d2*(Ev -5300d0)
      elseif (Ev.lt.5500) then
         f1 = read_histo(pi0mom_old_nm5400_1100,nrow,ncol,binw,icol,Erec)
         f2 = read_histo(pi0mom_old_nm5500_1100,nrow,ncol,binw,icol,Erec)
         f = f1 +(f2 -f1)/1d2*(Ev -5400d0)
      elseif (Ev.lt.5600) then
         f1 = read_histo(pi0mom_old_nm5500_1100,nrow,ncol,binw,icol,Erec)
         f2 = read_histo(pi0mom_old_nm5600_1100,nrow,ncol,binw,icol,Erec)
         f = f1 +(f2 -f1)/1d2*(Ev -5500d0)
      elseif (Ev.lt.5700) then
         f1 = read_histo(pi0mom_old_nm5600_1100,nrow,ncol,binw,icol,Erec)
         f2 = read_histo(pi0mom_old_nm5700_1100,nrow,ncol,binw,icol,Erec)
         f = f1 +(f2 -f1)/1d2*(Ev -5600d0)
      elseif (Ev.lt.5800) then
         f1 = read_histo(pi0mom_old_nm5700_1100,nrow,ncol,binw,icol,Erec)
         f2 = read_histo(pi0mom_old_nm5800_1100,nrow,ncol,binw,icol,Erec)
         f = f1 +(f2 -f1)/1d2*(Ev -5700d0)
      elseif (Ev.lt.5900) then
         f1 = read_histo(pi0mom_old_nm5800_1100,nrow,ncol,binw,icol,Erec)
         f2 = read_histo(pi0mom_old_nm5900_1100,nrow,ncol,binw,icol,Erec)
         f = f1 +(f2 -f1)/1d2*(Ev -5800d0)
      else
         f = 0d0
      endif
