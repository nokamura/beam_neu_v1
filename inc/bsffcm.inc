      integer ncall,itmx1,itmx2
      double precision acc1,acc2
      common /bsffcm/ acc1,acc2,ncall,itmx1,itmx2
      integer nkffbs,lkffbs
      parameter (nkffbs=2,lkffbs=20)
      double precision dvalbs(nkffbs)
      equivalence (dvalbs(1),acc1)
      character*80 ckffbs
      common /bsffc0/ ckffbs(nkffbs)
      character*4 ctagbs(nkffbs)
      common /bsffc1/ ctagbs
      save /bsffcm/,/bsffc0/,/bsffc1/
