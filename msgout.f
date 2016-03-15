************************************************************************
*     ------------------------------------
      SUBROUTINE MSGOUT(CNAME,CTEXT,ICODE)
*     ------------------------------------
*
*      (Purpose)
*        The error message handling routine for VANAL.
*
*      (Input)
*        CNAME           : Character information of an error source.
*        CTEXT           : Detailed text.
*
*      (Slave Subroutines)
*         MSGOT1 MSGOT2
*
*      (Creation Date)
*        Revised from VERROR at Sep 2,1986 by T.Tsuboyama.
*
************************************************************************

      CHARACTER CNAME*6,CTEXT*(*)
      PARAMETER (NTABLE=20)
      PARAMETER (NPMAX =50)
      CHARACTER*8  NAMES(NTABLE)
      CHARACTER*60 TEXTS(NTABLE),CTEMP
      INTEGER COUNT(NTABLE)
      INTEGER IRUNF(NTABLE),IEVTF(NTABLE)
      INTEGER IRUNL(NTABLE),IEVTL(NTABLE)
      INTEGER IENTRY/0/
      LOGICAL FULL/.FALSE./
      INTEGER ICODE,ICODE1

      INCLUDE 'inc/jlched.h'

      INCLUDE 'inc/msglun.h'

      CTEMP=CTEXT
      DO 1 I=1,IENTRY
        IF ( NAMES(I).EQ.CNAME.AND.CTEMP.EQ.TEXTS(I) ) THEN
          COUNT(I)=COUNT(I)+1
          IRUNL(I)=NRUN
          IEVTL(I)=NEVT
          IF (COUNT(I).GT.NPMAX) RETURN
          IF (COUNT(I).EQ.NPMAX) WRITE(LUNMSG,*)
     >         'MSGOUT: The following message have been printed ',NPMAX,
     >         'Times. Further calls for this message will be ignored.'
          GOTO 2
        ENDIF
1     CONTINUE
      IF (.NOT.FULL) THEN
        IENTRY=IENTRY+1
        IF (IENTRY.GT.NTABLE) THEN
          FULL=.TRUE.
          WRITE(LUNMSG,*) 'MSGOUT:*** Error Table Full***'
          IOVRUN=NRUN
          IOVEVT=NEVT
         ELSE
          I=IENTRY
          NAMES(I)=CNAME
          TEXTS(I)=CTEMP
          IRUNF(I)=NRUN
          IRUNL(I)=NRUN
          IEVTF(I)=NEVT
          IEVTL(I)=NEVT
          COUNT(I)=1
        ENDIF
      ENDIF
2     CONTINUE
      WRITE(LUNMSG,'(A,A,I8,A,I8,A,I6,A,A)') ' ===>',CNAME,ICODE,
     >            ',Run#',NRUN, ', Event#',NEVT,'=>',CTEMP
      RETURN

      ENTRY MSGFIN
      WRITE(LUNMSG,'(/A)')' MSGFIN: Error summary for this run'
      IF (IENTRY.EQ.0) THEN
        WRITE(LUNMSG,'(A)') ' MSGFIN: No Errors'
       ELSE
        IF (IENTRY.GT.NTABLE) THEN
          WRITE(LUNMSG,'(A,I6,A,I6,A)')
     >   ' MSGFIN: Error table overflow at Run',IOVRUN,', Event',IOVEVT
          IENTRY=NTABLE
        ENDIF
* BUBBLE SORT
        DO 3 I=1,IENTRY
          DO 4 J=I+1,IENTRY
            IF (  LLT(NAMES(J),NAMES(I)) .OR.
     > (NAMES(J).EQ.NAMES(I).AND.LLT(TEXTS(J),TEXTS(I))) ) THEN
* SWAP ITEMS
              CALL MSGOT1(TEXTS(J),TEXTS(I))
              CALL MSGOT2(COUNT(J),COUNT(I))
              CALL MSGOT2(IRUNF(J),IRUNF(I))
              CALL MSGOT2(IRUNL(J),IRUNL(I))
              CALL MSGOT2(IEVTF(J),IEVTF(I))
              CALL MSGOT2(IEVTL(J),IEVTL(I))
              CALL MSGOT3(NAMES(J),NAMES(I))
            ENDIF
4         CONTINUE
3       CONTINUE
        WRITE(LUNMSG,'(A)')
     >'  Routine Count     Run Event     Run Event'//
C       123456781234561234567812345612345678123456
     >'                Detail Message'
        DO 6 I=1,IENTRY
          WRITE(LUNMSG,'(1X,A8,I6,2(I8,I6),4X,A)')
     >       NAMES(I),COUNT(I),
     >       IRUNF(I),IEVTF(I),IRUNL(I),IEVTL(I),TEXTS(I)
6       CONTINUE
        IENTRY=0
      ENDIF
      RETURN
      END
      SUBROUTINE MSGOT1(C1,C2)
* SWAP CHARACTER*80
      CHARACTER*60 C1,C2,CTEMP
      CTEMP=C1
      C1=C2
      C2=CTEMP
      RETURN
      END
      SUBROUTINE MSGOT2(I1,I2)
* SWAP INTEGER
      INTEGER I1,I2,ITEMP
      ITEMP=I1
      I1=I2
      I2=ITEMP
      RETURN
      END
      SUBROUTINE MSGOT3(C1,C2)
* SWAP CHARACTER*8
      CHARACTER*8 C1,C2,CTEMP
      CTEMP=C1
      C1=C2
      C2=CTEMP
      RETURN
      END

      BLOCK DATA
      INCLUDE 'inc/msglun.h'
      DATA LUNMSG/ 6/
      END
