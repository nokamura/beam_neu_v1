C -*- fortran -*-
      if (Ev.lt.200) then
         f = 0d0
      elseif (Ev.lt.300) then
         f1 = read_histo(pi0_erec_am200_1100,nrow,ncol,binw,icol,Erec)
         f2 = read_histo(pi0_erec_am300_1100,nrow,ncol,binw,icol,Erec)
         f = f1 +(f2 -f1)/1d2*(Ev -200d0)
      elseif (Ev.lt.400) then
         f1 = read_histo(pi0_erec_am300_1100,nrow,ncol,binw,icol,Erec)
         f2 = read_histo(pi0_erec_am400_1100,nrow,ncol,binw,icol,Erec)
         f = f1 +(f2 -f1)/1d2*(Ev -300d0)
      elseif (Ev.lt.500) then
         f1 = read_histo(pi0_erec_am400_1100,nrow,ncol,binw,icol,Erec)
         f2 = read_histo(pi0_erec_am500_1100,nrow,ncol,binw,icol,Erec)
         f = f1 +(f2 -f1)/1d2*(Ev -400d0)
      elseif (Ev.lt.600) then
         f1 = read_histo(pi0_erec_am500_1100,nrow,ncol,binw,icol,Erec)
         f2 = read_histo(pi0_erec_am600_1100,nrow,ncol,binw,icol,Erec)
         f = f1 +(f2 -f1)/1d2*(Ev -500d0)
      elseif (Ev.lt.700) then
         f1 = read_histo(pi0_erec_am600_1100,nrow,ncol,binw,icol,Erec)
         f2 = read_histo(pi0_erec_am700_1100,nrow,ncol,binw,icol,Erec)
         f = f1 +(f2 -f1)/1d2*(Ev -600d0)
      elseif (Ev.lt.800) then
         f1 = read_histo(pi0_erec_am700_1100,nrow,ncol,binw,icol,Erec)
         f2 = read_histo(pi0_erec_am800_1100,nrow,ncol,binw,icol,Erec)
         f = f1 +(f2 -f1)/1d2*(Ev -700d0)
      elseif (Ev.lt.900) then
         f1 = read_histo(pi0_erec_am800_1100,nrow,ncol,binw,icol,Erec)
         f2 = read_histo(pi0_erec_am900_1100,nrow,ncol,binw,icol,Erec)
         f = f1 +(f2 -f1)/1d2*(Ev -800d0)
      elseif (Ev.lt.1000) then
         f1 = read_histo(pi0_erec_am900_1100,nrow,ncol,binw,icol,Erec)
         f2 = read_histo(pi0_erec_am1000_1100,nrow,ncol,binw,icol,Erec)
         f = f1 +(f2 -f1)/1d2*(Ev -900d0)
      elseif (Ev.lt.1100) then
         f1 = read_histo(pi0_erec_am1000_1100,nrow,ncol,binw,icol,Erec)
         f2 = read_histo(pi0_erec_am1100_1100,nrow,ncol,binw,icol,Erec)
         f = f1 +(f2 -f1)/1d2*(Ev -1000d0)
      elseif (Ev.lt.1200) then
         f1 = read_histo(pi0_erec_am1100_1100,nrow,ncol,binw,icol,Erec)
         f2 = read_histo(pi0_erec_am1200_1100,nrow,ncol,binw,icol,Erec)
         f = f1 +(f2 -f1)/1d2*(Ev -1100d0)
      elseif (Ev.lt.1300) then
         f1 = read_histo(pi0_erec_am1200_1100,nrow,ncol,binw,icol,Erec)
         f2 = read_histo(pi0_erec_am1300_1100,nrow,ncol,binw,icol,Erec)
         f = f1 +(f2 -f1)/1d2*(Ev -1200d0)
      elseif (Ev.lt.1400) then
         f1 = read_histo(pi0_erec_am1300_1100,nrow,ncol,binw,icol,Erec)
         f2 = read_histo(pi0_erec_am1400_1100,nrow,ncol,binw,icol,Erec)
         f = f1 +(f2 -f1)/1d2*(Ev -1300d0)
      elseif (Ev.lt.1500) then
         f1 = read_histo(pi0_erec_am1400_1100,nrow,ncol,binw,icol,Erec)
         f2 = read_histo(pi0_erec_am1500_1100,nrow,ncol,binw,icol,Erec)
         f = f1 +(f2 -f1)/1d2*(Ev -1400d0)
      elseif (Ev.lt.1600) then
         f1 = read_histo(pi0_erec_am1500_1100,nrow,ncol,binw,icol,Erec)
         f2 = read_histo(pi0_erec_am1600_1100,nrow,ncol,binw,icol,Erec)
         f = f1 +(f2 -f1)/1d2*(Ev -1500d0)
      elseif (Ev.lt.1700) then
         f1 = read_histo(pi0_erec_am1600_1100,nrow,ncol,binw,icol,Erec)
         f2 = read_histo(pi0_erec_am1700_1100,nrow,ncol,binw,icol,Erec)
         f = f1 +(f2 -f1)/1d2*(Ev -1600d0)
      elseif (Ev.lt.1800) then
         f1 = read_histo(pi0_erec_am1700_1100,nrow,ncol,binw,icol,Erec)
         f2 = read_histo(pi0_erec_am1800_1100,nrow,ncol,binw,icol,Erec)
         f = f1 +(f2 -f1)/1d2*(Ev -1700d0)
      elseif (Ev.lt.1900) then
         f1 = read_histo(pi0_erec_am1800_1100,nrow,ncol,binw,icol,Erec)
         f2 = read_histo(pi0_erec_am1900_1100,nrow,ncol,binw,icol,Erec)
         f = f1 +(f2 -f1)/1d2*(Ev -1800d0)
      elseif (Ev.lt.2000) then
         f1 = read_histo(pi0_erec_am1900_1100,nrow,ncol,binw,icol,Erec)
         f2 = read_histo(pi0_erec_am2000_1100,nrow,ncol,binw,icol,Erec)
         f = f1 +(f2 -f1)/1d2*(Ev -1900d0)
      elseif (Ev.lt.2100) then
         f1 = read_histo(pi0_erec_am2000_1100,nrow,ncol,binw,icol,Erec)
         f2 = read_histo(pi0_erec_am2100_1100,nrow,ncol,binw,icol,Erec)
         f = f1 +(f2 -f1)/1d2*(Ev -2000d0)
      elseif (Ev.lt.2200) then
         f1 = read_histo(pi0_erec_am2100_1100,nrow,ncol,binw,icol,Erec)
         f2 = read_histo(pi0_erec_am2200_1100,nrow,ncol,binw,icol,Erec)
         f = f1 +(f2 -f1)/1d2*(Ev -2100d0)
      elseif (Ev.lt.2300) then
         f1 = read_histo(pi0_erec_am2200_1100,nrow,ncol,binw,icol,Erec)
         f2 = read_histo(pi0_erec_am2300_1100,nrow,ncol,binw,icol,Erec)
         f = f1 +(f2 -f1)/1d2*(Ev -2200d0)
      elseif (Ev.lt.2400) then
         f1 = read_histo(pi0_erec_am2300_1100,nrow,ncol,binw,icol,Erec)
         f2 = read_histo(pi0_erec_am2400_1100,nrow,ncol,binw,icol,Erec)
         f = f1 +(f2 -f1)/1d2*(Ev -2300d0)
      elseif (Ev.lt.2500) then
         f1 = read_histo(pi0_erec_am2400_1100,nrow,ncol,binw,icol,Erec)
         f2 = read_histo(pi0_erec_am2500_1100,nrow,ncol,binw,icol,Erec)
         f = f1 +(f2 -f1)/1d2*(Ev -2400d0)
      elseif (Ev.lt.2600) then
         f1 = read_histo(pi0_erec_am2500_1100,nrow,ncol,binw,icol,Erec)
         f2 = read_histo(pi0_erec_am2600_1100,nrow,ncol,binw,icol,Erec)
         f = f1 +(f2 -f1)/1d2*(Ev -2500d0)
      elseif (Ev.lt.2700) then
         f1 = read_histo(pi0_erec_am2600_1100,nrow,ncol,binw,icol,Erec)
         f2 = read_histo(pi0_erec_am2700_1100,nrow,ncol,binw,icol,Erec)
         f = f1 +(f2 -f1)/1d2*(Ev -2600d0)
      elseif (Ev.lt.2800) then
         f1 = read_histo(pi0_erec_am2700_1100,nrow,ncol,binw,icol,Erec)
         f2 = read_histo(pi0_erec_am2800_1100,nrow,ncol,binw,icol,Erec)
         f = f1 +(f2 -f1)/1d2*(Ev -2700d0)
      elseif (Ev.lt.2900) then
         f1 = read_histo(pi0_erec_am2800_1100,nrow,ncol,binw,icol,Erec)
         f2 = read_histo(pi0_erec_am2900_1100,nrow,ncol,binw,icol,Erec)
         f = f1 +(f2 -f1)/1d2*(Ev -2800d0)
      elseif (Ev.lt.3000) then
         f1 = read_histo(pi0_erec_am2900_1100,nrow,ncol,binw,icol,Erec)
         f2 = read_histo(pi0_erec_am3000_1100,nrow,ncol,binw,icol,Erec)
         f = f1 +(f2 -f1)/1d2*(Ev -2900d0)
      elseif (Ev.lt.3100) then
         f1 = read_histo(pi0_erec_am3000_1100,nrow,ncol,binw,icol,Erec)
         f2 = read_histo(pi0_erec_am3100_1100,nrow,ncol,binw,icol,Erec)
         f = f1 +(f2 -f1)/1d2*(Ev -3000d0)
      elseif (Ev.lt.3200) then
         f1 = read_histo(pi0_erec_am3100_1100,nrow,ncol,binw,icol,Erec)
         f2 = read_histo(pi0_erec_am3200_1100,nrow,ncol,binw,icol,Erec)
         f = f1 +(f2 -f1)/1d2*(Ev -3100d0)
      elseif (Ev.lt.3300) then
         f1 = read_histo(pi0_erec_am3200_1100,nrow,ncol,binw,icol,Erec)
         f2 = read_histo(pi0_erec_am3300_1100,nrow,ncol,binw,icol,Erec)
         f = f1 +(f2 -f1)/1d2*(Ev -3200d0)
      elseif (Ev.lt.3400) then
         f1 = read_histo(pi0_erec_am3300_1100,nrow,ncol,binw,icol,Erec)
         f2 = read_histo(pi0_erec_am3400_1100,nrow,ncol,binw,icol,Erec)
         f = f1 +(f2 -f1)/1d2*(Ev -3300d0)
      elseif (Ev.lt.3500) then
         f1 = read_histo(pi0_erec_am3400_1100,nrow,ncol,binw,icol,Erec)
         f2 = read_histo(pi0_erec_am3500_1100,nrow,ncol,binw,icol,Erec)
         f = f1 +(f2 -f1)/1d2*(Ev -3400d0)
      elseif (Ev.lt.3600) then
         f1 = read_histo(pi0_erec_am3500_1100,nrow,ncol,binw,icol,Erec)
         f2 = read_histo(pi0_erec_am3600_1100,nrow,ncol,binw,icol,Erec)
         f = f1 +(f2 -f1)/1d2*(Ev -3500d0)
      elseif (Ev.lt.3700) then
         f1 = read_histo(pi0_erec_am3600_1100,nrow,ncol,binw,icol,Erec)
         f2 = read_histo(pi0_erec_am3700_1100,nrow,ncol,binw,icol,Erec)
         f = f1 +(f2 -f1)/1d2*(Ev -3600d0)
      elseif (Ev.lt.3800) then
         f1 = read_histo(pi0_erec_am3700_1100,nrow,ncol,binw,icol,Erec)
         f2 = read_histo(pi0_erec_am3800_1100,nrow,ncol,binw,icol,Erec)
         f = f1 +(f2 -f1)/1d2*(Ev -3700d0)
      elseif (Ev.lt.3900) then
         f1 = read_histo(pi0_erec_am3800_1100,nrow,ncol,binw,icol,Erec)
         f2 = read_histo(pi0_erec_am3900_1100,nrow,ncol,binw,icol,Erec)
         f = f1 +(f2 -f1)/1d2*(Ev -3800d0)
      elseif (Ev.lt.4000) then
         f1 = read_histo(pi0_erec_am3900_1100,nrow,ncol,binw,icol,Erec)
         f2 = read_histo(pi0_erec_am4000_1100,nrow,ncol,binw,icol,Erec)
         f = f1 +(f2 -f1)/1d2*(Ev -3900d0)
      elseif (Ev.lt.4100) then
         f1 = read_histo(pi0_erec_am4000_1100,nrow,ncol,binw,icol,Erec)
         f2 = read_histo(pi0_erec_am4100_1100,nrow,ncol,binw,icol,Erec)
         f = f1 +(f2 -f1)/1d2*(Ev -4000d0)
      elseif (Ev.lt.4200) then
         f1 = read_histo(pi0_erec_am4100_1100,nrow,ncol,binw,icol,Erec)
         f2 = read_histo(pi0_erec_am4200_1100,nrow,ncol,binw,icol,Erec)
         f = f1 +(f2 -f1)/1d2*(Ev -4100d0)
      elseif (Ev.lt.4300) then
         f1 = read_histo(pi0_erec_am4200_1100,nrow,ncol,binw,icol,Erec)
         f2 = read_histo(pi0_erec_am4300_1100,nrow,ncol,binw,icol,Erec)
         f = f1 +(f2 -f1)/1d2*(Ev -4200d0)
      elseif (Ev.lt.4400) then
         f1 = read_histo(pi0_erec_am4300_1100,nrow,ncol,binw,icol,Erec)
         f2 = read_histo(pi0_erec_am4400_1100,nrow,ncol,binw,icol,Erec)
         f = f1 +(f2 -f1)/1d2*(Ev -4300d0)
      elseif (Ev.lt.4500) then
         f1 = read_histo(pi0_erec_am4400_1100,nrow,ncol,binw,icol,Erec)
         f2 = read_histo(pi0_erec_am4500_1100,nrow,ncol,binw,icol,Erec)
         f = f1 +(f2 -f1)/1d2*(Ev -4400d0)
      elseif (Ev.lt.4600) then
         f1 = read_histo(pi0_erec_am4500_1100,nrow,ncol,binw,icol,Erec)
         f2 = read_histo(pi0_erec_am4600_1100,nrow,ncol,binw,icol,Erec)
         f = f1 +(f2 -f1)/1d2*(Ev -4500d0)
      elseif (Ev.lt.4700) then
         f1 = read_histo(pi0_erec_am4600_1100,nrow,ncol,binw,icol,Erec)
         f2 = read_histo(pi0_erec_am4700_1100,nrow,ncol,binw,icol,Erec)
         f = f1 +(f2 -f1)/1d2*(Ev -4600d0)
      elseif (Ev.lt.4800) then
         f1 = read_histo(pi0_erec_am4700_1100,nrow,ncol,binw,icol,Erec)
         f2 = read_histo(pi0_erec_am4800_1100,nrow,ncol,binw,icol,Erec)
         f = f1 +(f2 -f1)/1d2*(Ev -4700d0)
      elseif (Ev.lt.4900) then
         f1 = read_histo(pi0_erec_am4800_1100,nrow,ncol,binw,icol,Erec)
         f2 = read_histo(pi0_erec_am4900_1100,nrow,ncol,binw,icol,Erec)
         f = f1 +(f2 -f1)/1d2*(Ev -4800d0)
      elseif (Ev.lt.5000) then
         f1 = read_histo(pi0_erec_am4900_1100,nrow,ncol,binw,icol,Erec)
         f2 = read_histo(pi0_erec_am5000_1100,nrow,ncol,binw,icol,Erec)
         f = f1 +(f2 -f1)/1d2*(Ev -4900d0)
      elseif (Ev.lt.5100) then
         f1 = read_histo(pi0_erec_am5000_1100,nrow,ncol,binw,icol,Erec)
         f2 = read_histo(pi0_erec_am5100_1100,nrow,ncol,binw,icol,Erec)
         f = f1 +(f2 -f1)/1d2*(Ev -5000d0)
      elseif (Ev.lt.5200) then
         f1 = read_histo(pi0_erec_am5100_1100,nrow,ncol,binw,icol,Erec)
         f2 = read_histo(pi0_erec_am5200_1100,nrow,ncol,binw,icol,Erec)
         f = f1 +(f2 -f1)/1d2*(Ev -5100d0)
      elseif (Ev.lt.5300) then
         f1 = read_histo(pi0_erec_am5200_1100,nrow,ncol,binw,icol,Erec)
         f2 = read_histo(pi0_erec_am5300_1100,nrow,ncol,binw,icol,Erec)
         f = f1 +(f2 -f1)/1d2*(Ev -5200d0)
      elseif (Ev.lt.5400) then
         f1 = read_histo(pi0_erec_am5300_1100,nrow,ncol,binw,icol,Erec)
         f2 = read_histo(pi0_erec_am5400_1100,nrow,ncol,binw,icol,Erec)
         f = f1 +(f2 -f1)/1d2*(Ev -5300d0)
      elseif (Ev.lt.5500) then
         f1 = read_histo(pi0_erec_am5400_1100,nrow,ncol,binw,icol,Erec)
         f2 = read_histo(pi0_erec_am5500_1100,nrow,ncol,binw,icol,Erec)
         f = f1 +(f2 -f1)/1d2*(Ev -5400d0)
      elseif (Ev.lt.5600) then
         f1 = read_histo(pi0_erec_am5500_1100,nrow,ncol,binw,icol,Erec)
         f2 = read_histo(pi0_erec_am5600_1100,nrow,ncol,binw,icol,Erec)
         f = f1 +(f2 -f1)/1d2*(Ev -5500d0)
      elseif (Ev.lt.5700) then
         f1 = read_histo(pi0_erec_am5600_1100,nrow,ncol,binw,icol,Erec)
         f2 = read_histo(pi0_erec_am5700_1100,nrow,ncol,binw,icol,Erec)
         f = f1 +(f2 -f1)/1d2*(Ev -5600d0)
      elseif (Ev.lt.5800) then
         f1 = read_histo(pi0_erec_am5700_1100,nrow,ncol,binw,icol,Erec)
         f2 = read_histo(pi0_erec_am5800_1100,nrow,ncol,binw,icol,Erec)
         f = f1 +(f2 -f1)/1d2*(Ev -5700d0)
      elseif (Ev.lt.5900) then
         f1 = read_histo(pi0_erec_am5800_1100,nrow,ncol,binw,icol,Erec)
         f2 = read_histo(pi0_erec_am5900_1100,nrow,ncol,binw,icol,Erec)
         f = f1 +(f2 -f1)/1d2*(Ev -5800d0)
      else
         f = 0d0
      endif
