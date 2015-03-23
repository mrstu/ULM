SUBROUTINE READ_SNOWBANDS()

  ! UW Land Surface Hydrology Group implementation of Noah model
  ! modified from NLDAS implementation
  ! author: Ted Bohn, tbohn@hydro.washington.edu

  ! READS VIC SNOWBANDS FILE (ascii)

  ! driverMod contains definitions for all global driver variables
  USE driverMod

  IMPLICIT NONE
      
  ! RCS Id string, for version control
  CHARACTER*60 RCSID
  DATA RCSID/"$Id: READ_SNOWBANDS.f90,v 1.1.1.1 2007/10/04 20:01:30 vicadmin Exp $"/

  ! Define local variables
  INTEGER :: I,J
  INTEGER :: tmp_cellid
  REAL    :: total_area, total_prec, avg_elev, band_elev_wgt,band_elev_sum,bandlapsed_band,bandlapsed_cell


  ! Get band parameters
  IF (nbands == 1) THEN

    DO I = 1, landlen
      band_area(I,1) = 1.0
      band_elev(I,1) = ELEV_2d(I)
      band_prec(I,1) = 1.0
    END DO

  ELSE

    OPEN(99,FILE=SNOWBAND_FILE,STATUS='OLD')
    DO I = 1, landlen
      READ (99,*)tmp_cellid,(band_area(I,J),J=1,nbands),(band_elev(I,J),J=1,nbands),(band_prec(I,J),J=1,nbands)
    END DO
    CLOSE(99)

  END IF


  ! Compute Tfactor, Pfactor
  DO I = 1, landlen

    ! Check that fractions add to 1
    total_area = 0
    total_prec = 0
    avg_elev = 0
    ! MRS20150310: start enhancement
    !   If avergage elevation from soil file and snowbands file don't match
    !   -> set to average elevation to weighted average of snowbands.
    !
    !   Assume snowbands is derived from higher resolution DEM
    !   Potential pitfall of covering up grid misalignment errors
    DO J = 1, nbands
      IF (band_area(I,J) > 0) THEN
        avg_elev = avg_elev + band_elev(I,J) * band_area(I,J)
      END IF
      total_area = total_area + band_area(I,J)
      total_prec = total_prec + band_prec(I,J)
    END DO
    !   If (total_area < 1) check not needed for avg_elev b/c it is normalized by total_area
    avg_elev = avg_elev / total_area
    ! arbitrary small difference threshold between cell/band sourced avg elevation, 5 meters.
    IF ( 0 == 1 ) THEN
!    IF (ABS(avg_elev - ELEV_2d(I)) > 0.1) THEN
        write(*,'("Error: cell elevation not equal to band averaged elevation")')
        write(*,'("Cell elevation, id:",f10.2,i2)')ELEV_2d(I),I
        write(*,'("Band average elevation:",f10.2)')avg_elev
        band_elev_sum = 0
        write(*,'("Band <2 bandnum>: <3 band_area>, <4 band_elev_wgt>, <5 band_elev>, <6-7 lapse band rel {bands,cell} average>, <8-9 height diff band-{band_avg, cell_avg}>)")')
        DO J = 1, nbands
          IF (band_area(I,J) > 0) THEN
            band_elev_wgt = band_area(I,J) / total_area
            band_elev_sum = band_elev_sum + band_elev(I,J) * band_elev_wgt
            bandlapsed_band = (avg_elev-band_elev(I,J)) / 1000. * LAPSE_RATE
            bandlapsed_cell = (ELEV_2d(I)-band_elev(I,J)) / 1000. * LAPSE_RATE
            write(*,'("Band",i2,": ",f5.2, f5.2, f10.1, f6.2, f6.2, f8.2, f8.2)')J,band_area(I,J),band_elev_wgt,band_elev(I,J),bandlapsed_band,bandlapsed_cell,(band_elev(I,J)-avg_elev),(band_elev(I,J)-ELEV_2d(I))
          END IF
        END DO
        ! write(*,'("Band ",i2,": ",f8.2,f8.2,f8.2,f8.2,f8.2)')J,band_elev(I,J),band_area(I,J),band_elev_wgt,band_elev_sum, band_elev_sum-ELEV_2d(I), band_elev_sum-ELEV_2d(I)
        ! write(*,'("Cell elevation (",f8.2,") not equal to weighted band average elevation (",f8.2,")")')ELEV_2d(I),avg_elev
        write(*,'("Now exiting with signal 12")')
        !STOP 12
    !    ELEV_2d(I)=avg_elev
    END IF

    ! MRS20150310: end enhancement

    DO J = 1, nbands
      IF (band_area(I,J) > 0) THEN
        Tfactor(I,J) = LAPSE_RATE * (ELEV_2d(I) - band_elev(I,J)) / 1000
      ELSE
        Tfactor(I,J) = 0
      END IF
      total_area = total_area + band_area(I,J)
      total_prec = total_prec + band_prec(I,J)
    END DO


    IF (total_area == 0) THEN
      write(*,*)'Error: sum of snowband area fractions is 0'
      stop
    END IF
    IF (total_prec == 0) THEN
      write(*,*)'Error: sum of snowband precip fractions is 0'
      stop
    END IF

    ! Adjust fractions to add to 1
    IF (total_area /= 1.0) THEN
      DO J = 1, nbands
        band_area(I,J) = band_area(I,J) / total_area
      END DO
    END IF
    IF (total_prec /= 1.0) THEN
      DO J = 1, nbands
        band_prec(I,J) = band_prec(I,J) / total_prec
      END DO
    END IF

    ! Scale prec fractions by area fractions
    DO J = 1, nbands
      IF (band_area(I,J) > 0) THEN
        Pfactor(I,J) = band_prec(I,J) / band_area(I,J)
      ELSE
        Pfactor(I,J) = 0
      END IF
!write(*,*)I,J,elev_2d(I),band_area(I,J),band_elev(I,J),band_prec(I,J),Tfactor(I,J),Pfactor(I,J)
    END DO

  END DO

END
