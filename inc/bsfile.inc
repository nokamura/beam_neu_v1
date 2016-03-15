      integer lubsf,luisrf,luhist,luevtf
      integer nospev,nevprt
      common /bsfile/
     &        lubsf,luisrf,luhist,luevtf,
     &        nospev,nevprt
      character*128 bsoutf,bsisrf,bshstf,spevtf
      common /bsfil0/ bsoutf,bsisrf,bshstf,spevtf
      save /bsfile/,/bsfil0/
