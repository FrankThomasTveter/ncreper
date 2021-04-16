      SUBROUTINE CHINDEF(SATNAME,INSTRUMENT,NRCHAN,CHINDEX,IRC)
C     ***********************************************************************
C     WRITTEN/MODIFIED BY:                                                  *
C     --------------------------------------------------------------------- *
C     |    NAME      |   DATE   |                 REASON                  | *
C     --------------------------------------------------------------------- *
C     | F. TVETER    | 20/08/02 | NEW                                     | *
C     |              |          |                                         | *
C     --------------------------------------------------------------------- *
C     ***********************************************************************
C
C     SETS UP DEFAULT CHANNEL INDEXES FOR THE SATELLITES
C
      IMPLICIT NONE
      SAVE
C
      INTEGER  MAXCH
      PARAMETER (MAXCH = 35)
      CHARACTER*8 SATNAME
      CHARACTER*12 INSTRUMENT
      INTEGER  NRCHAN
      INTEGER  CHINDEX(MAXCH)
      INTEGER  IRC
C
      INTEGER  LENS, LENA,LENGTH
      EXTERNAL LENGTH
      CHARACTER*8 MYNAME
      DATA MYNAME /'CHINDEF'/
C
      INTEGER  NCHNOAA, NCHNOAAK, NCHNOAAHRS, INRCHAN,II
      PARAMETER (NCHNOAA=23, NCHNOAAK=15, NCHNOAAHRS=19)
      INTEGER  ICHNOAA(NCHNOAA),ICHNOAAK(NCHNOAAK),
     &     ICHNOAAKA(NCHNOAAK), ICHNOAAHRS(NCHNOAAHRS)
C
      DATA ICHNOAA/1,2,3,4,5,
     &     6, 7, 8, 9, 10,
     &     11,12,13,14,15,
     &     16,17,18,19,
     &     21,22,23,24/
      DATA ICHNOAAK/
     &     28,29,30,31,32,33,34,35,36,37,
     &     38,39,40,41,42/      ! amsu-a only
      DATA ICHNOAAHRS/1,2,3,4,5,
     &     6, 7, 8, 9, 10,
     &     11,12,13,14,15,
     &     16,17,18,19/         !hirs IR channels
C
      LENS=LENGTH(SATNAME,8,6)
      LENA=LENGTH(INSTRUMENT,12,6)
      IF (SATNAME(1:LENS).EQ.'noaa12'.AND.
     &     INSTRUMENT(1:LENA).EQ.'hirs+msu') THEN
         INRCHAN=NCHNOAA
         DO II=1,INRCHAN
            CHINDEX(II)=ICHNOAA(II)
         END DO
      ELSE IF (SATNAME(1:LENS).EQ.'noaa14'.AND.
     &     INSTRUMENT(1:LENA).EQ.'hirs+msu') THEN
         INRCHAN=NCHNOAA
         DO II=1,INRCHAN
            CHINDEX(II)=ICHNOAA(II)
         END DO
      ELSE IF (SATNAME(1:LENS).EQ.'noaa15'.AND.
     &     INSTRUMENT(1:LENA).EQ.'amsu-a') THEN
         INRCHAN=NCHNOAAK
         DO II=1,INRCHAN
            CHINDEX(II)=ICHNOAAK(II)
         END DO
      ELSE IF (SATNAME(1:LENS).EQ.'noaa16'.AND.
     &     INSTRUMENT(1:LENA).EQ.'amsu-a') THEN
         INRCHAN=NCHNOAAK
         DO II=1,INRCHAN
            CHINDEX(II)=ICHNOAAK(II)
         END DO
      ELSE IF (SATNAME(1:LENS).EQ.'noaa17'.AND.
     &     INSTRUMENT(1:LENA).EQ.'amsu-a') THEN
         INRCHAN=NCHNOAAK
         DO II=1,INRCHAN
            CHINDEX(II)=ICHNOAAK(II)
         END DO
      ELSE IF (SATNAME(1:LENS).EQ.'noaa18'.AND.
     &     INSTRUMENT(1:LENA).EQ.'amsu-a') THEN
         INRCHAN=NCHNOAAK
         DO II=1,INRCHAN
            CHINDEX(II)=ICHNOAAK(II)
         END DO
      ELSE IF (SATNAME(1:LENS).EQ.'metop-2'.AND.
     &     INSTRUMENT(1:LENA).EQ.'amsu-a') THEN
         INRCHAN=NCHNOAAK
         DO II=1,INRCHAN
            CHINDEX(II)=ICHNOAAK(II)
         END DO
      ELSE IF ((satname(1:lens) .eq. 'noaa14'          !hirs from noaa14-16
     &        .or. satname(1:lens) .eq. 'noaa15'
     &        .or. satname(1:lens) .eq. 'noaa16'
     &        .or. satname(1:lens) .eq. 'noaa17'
     &        .or. satname(1:lens) .eq. 'noaa18'
     &        ) .and.
     &        instrument(1:lena) .eq. 'hirs') THEN
         INRCHAN=NCHNOAAHRS
         DO II=1,INRCHAN
            CHINDEX(II)=ICHNOAAHRS(II)
         END DO
      ELSE
         WRITE(*,*) MYNAME,'Unknown satellite:',
     &        SATNAME(1:LENS)//
     &        ' ('//INSTRUMENT(1:LENA)//')'
         IRC=40
         RETURN
      END IF
C
      IF (NRCHAN.GT.0 .AND. NRCHAN.NE.INRCHAN) THEN
         WRITE(*,*) MYNAME,'DEFAULT CHANNEL '//
     &        'INDEX MISMATCH.',NRCHAN,INRCHAN
         IRC=41
         RETURN
      ELSE
         NRCHAN=INRCHAN
      END IF
C
      RETURN
      END
