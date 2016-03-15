      integer nkeys
      common /hmffdt/ nkeys
      integer nkymax
      parameter (nkymax = 200)
      character*20 ckey
      character*80 cval
      common /hmffd0/ ckey(nkymax),cval(nkymax)
      save /hmffdt/, /hmffd0/
