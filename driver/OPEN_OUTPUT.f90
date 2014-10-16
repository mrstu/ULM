SUBROUTINE OPEN_OUTPUT()

  ! UW Land Surface Hydrology Group implementation of Noah model
  ! modified from NLDAS implementation
  ! author: Ted Bohn, tbohn@hydro.washington.edu

  ! CREATES OUTPUT FILES

  ! Modifications:
  ! 2007-Nov-06 Added SnowSoilT (T12).						Ben Livneh
  ! 2008-May-05 Removed SnowSoilT (T12), for compatibility with NOAH 2.8.	Ben Livneh
  ! 2008-May-15 Prints error message and quits if can't open output file.	TJB
  ! 2008-Jul-24 Added SnowTProf and SliqFrac.					Ben Livneh
  ! 2011-May-2 Added the _FillValue attribute to compliment missing_value for nco tools B.L.


  ! driverMod contains definitions of all global driver variables
  USE driverMod

  IMPLICIT NONE

  ! RCS Id string, for version control
  CHARACTER*60 RCSID
  DATA RCSID/"$Id: OPEN_OUTPUT.f90,v 1.7 2008/07/30 22:33:03 vicadmin Exp $"/

  ! Define local variables
  INTEGER nvars,ngatts,natts,unlimited
  INTEGER xdimid,ydimid,landdimid,leveldimid,banddimid,tstepdimid
  INTEGER rowvarid,colvarid,latvarid,lonvarid,cellidvarid,landmaskvarid,varid
  CHARACTER*20 name
  INTEGER xtype, dimids(2),odim(2),outdim_2d(3),outdim_3d(4),outdim_3d_pe(4),ndims,ndims_2d,ndims_3d
  INTEGER I,J,K,L,NT
  INTEGER :: temp_int(landlen)
  REAL :: temp(landlen)
  INTEGER :: land_idx

  IF (COMP_OUTPUT) THEN
    ndims = 1
    ndims_2d = 2
    ndims_3d = 3
  ELSE
    ndims = 2
    ndims_2d = 3
    ndims_3d = 4
  END IF

  DO K = 1,7

!! ********* SW MONITOR HACK - only write the "wb" and "pe" files **********
!!    IF (K == 2 .or. K == 7) THEN
    IF (K == 2) THEN
!! ********* ULM HACK ****************************************************
!     IF (K == 2) THEN
!! *******************************************************
    
    status = NF_CREATE(OUTFILES(K),0,OUT_NCIDS(K))
    IF (status .ne. NF_NOERR) THEN
      WRITE(*,*)'ERROR: cannot open output file',OUTFILES(K)
      STOP
    END IF
    if (SUFFIX.EQ.".199505.nc")then
       write(*,*)'in open_output'
    endif
    status = NF_DEF_DIM(OUT_NCIDS(K),'tstep',NF_UNLIMITED,tstepdimid)
    IF (K == 2 .or. K == 4) THEN
      status = NF_DEF_DIM(OUT_NCIDS(K),'level',MAXNSOIL,leveldimid)
    ELSE IF (K == 7) THEN
      status = NF_DEF_DIM(OUT_NCIDS(K),'band',NBANDS,banddimid)
    END IF
    IF (COMP_OUTPUT) THEN
      status = NF_DEF_DIM(OUT_NCIDS(K),'land',landlen,landdimid)
    ELSE
      status = NF_DEF_DIM(OUT_NCIDS(K),'y',ylen,ydimid)
      status = NF_DEF_DIM(OUT_NCIDS(K),'x',xlen,xdimid)
    END IF

    IF (COMP_OUTPUT) THEN
      odim(1) = landdimid
      outdim_2d(1) = landdimid
      outdim_2d(2) = tstepdimid
      outdim_3d(1) = landdimid
      outdim_3d(2) = leveldimid
      outdim_3d(3) = tstepdimid
      outdim_3d_pe(1) = landdimid
      outdim_3d_pe(2) = banddimid
      outdim_3d_pe(3) = tstepdimid
    ELSE
      odim(1) = xdimid
      odim(2) = ydimid
      outdim_2d(1) = xdimid
      outdim_2d(2) = ydimid
      outdim_2d(3) = tstepdimid
      outdim_3d(1) = xdimid
      outdim_3d(2) = ydimid
      outdim_3d(3) = leveldimid
      outdim_3d(4) = tstepdimid
      outdim_3d_pe(1) = xdimid
      outdim_3d_pe(2) = ydimid
      outdim_3d_pe(3) = banddimid
      outdim_3d_pe(4) = tstepdimid
    END IF

    status = NF_PUT_ATT_TEXT(OUT_NCIDS(K),NF_GLOBAL,'Conventions',7,'GDT 1.3')
    status = NF_PUT_ATT_TEXT(OUT_NCIDS(K),NF_GLOBAL,'file_name',31,OUTFILES(K))

    status = NF_DEF_VAR(OUT_NCIDS(K),'row',NF_INT,ndims,odim,rowvarid)
    status = NF_PUT_ATT_TEXT(OUT_NCIDS(K),rowvarid,'long_name',8,'Grid Row')
    status = NF_PUT_ATT_INT(OUT_NCIDS(K),rowvarid,'missing_value', NF_INT,1,NODATA_INT)
    status = NF_PUT_ATT_INT(OUT_NCIDS(K),rowvarid,'_FillValue', NF_INT,1,NODATA_INT)

    status = NF_DEF_VAR(OUT_NCIDS(K),'col',NF_INT,ndims,odim,colvarid)
    status = NF_PUT_ATT_TEXT(OUT_NCIDS(K),colvarid,'long_name',11,'Grid Column')
    status = NF_PUT_ATT_INT(OUT_NCIDS(K),colvarid,'missing_value', NF_INT,1,NODATA_INT)
    status = NF_PUT_ATT_INT(OUT_NCIDS(K),colvarid,'_FillValue', NF_INT,1,NODATA_INT)

    status = NF_DEF_VAR(OUT_NCIDS(K),'nav_lat',NF_REAL,ndims,odim,latvarid)
    status = NF_PUT_ATT_TEXT(OUT_NCIDS(K),latvarid,'units',13,'degrees_north')
    status = NF_PUT_ATT_REAL(OUT_NCIDS(K),latvarid,'valid_min',NF_REAL,1,-90.)
    status = NF_PUT_ATT_REAL(OUT_NCIDS(K),latvarid,'valid_max',NF_REAL,1,90.)
    status = NF_PUT_ATT_TEXT(OUT_NCIDS(K),latvarid,'long_name',8,'Latitude')
    status = NF_PUT_ATT_REAL(OUT_NCIDS(K),latvarid,'missing_value', NF_REAL,1,NODATA)
    status = NF_PUT_ATT_REAL(OUT_NCIDS(K),latvarid,'_FillValue', NF_REAL,1,NODATA)

    status = NF_DEF_VAR(OUT_NCIDS(K),'nav_lon',NF_REAL,ndims,odim,lonvarid)
    status = NF_PUT_ATT_TEXT(OUT_NCIDS(K),lonvarid,'units',12,'degrees_east')
    status = NF_PUT_ATT_REAL(OUT_NCIDS(K),lonvarid,'valid_min',NF_REAL,1,-180.)
    status = NF_PUT_ATT_REAL(OUT_NCIDS(K),lonvarid,'valid_max',NF_REAL,1,180.)
    status = NF_PUT_ATT_TEXT(OUT_NCIDS(K),lonvarid,'long_name',9,'Longitude')
    status = NF_PUT_ATT_REAL(OUT_NCIDS(K),lonvarid,'missing_value', NF_REAL,1,NODATA)
    status = NF_PUT_ATT_REAL(OUT_NCIDS(K),lonvarid,'_FillValue', NF_REAL,1,NODATA)

    status = NF_DEF_VAR(OUT_NCIDS(K),'land',NF_INT,ndims,odim,landmaskvarid)
    status = NF_PUT_ATT_TEXT(OUT_NCIDS(K),landmaskvarid,'units',1,'-')
    status = NF_PUT_ATT_TEXT(OUT_NCIDS(K),landmaskvarid,'long_name',9,'Land Mask')
    status = NF_PUT_ATT_INT(OUT_NCIDS(K),landmaskvarid,'missing_value', NF_INT,1,NODATA_INT)
    status = NF_PUT_ATT_INT(OUT_NCIDS(K),landmaskvarid,'_FillValue', NF_INT,1,NODATA_INT)

    status = NF_DEF_VAR(OUT_NCIDS(K),'CellID',NF_INT,ndims,odim,cellidvarid)
    status = NF_PUT_ATT_TEXT(OUT_NCIDS(K),cellidvarid,'long_name',7,'Cell ID')
    status = NF_PUT_ATT_INT(OUT_NCIDS(K),cellidvarid,'missing_value', NF_INT,1,NODATA_INT)
    status = NF_PUT_ATT_INT(OUT_NCIDS(K),cellidvarid,'_FillValue', NF_INT,1,NODATA_INT)

    status = NF_DEF_VAR(OUT_NCIDS(K),'time',NF_REAL,1,tstepdimid,varid)
    status = NF_PUT_ATT_TEXT(OUT_NCIDS(K),varid,'units',7,'seconds')
    status = NF_PUT_ATT_TEXT(OUT_NCIDS(K),varid,'calendar',9,'gregorian')
    status = NF_PUT_ATT_TEXT(OUT_NCIDS(K),varid,'title',4,'Time')
    status = NF_PUT_ATT_TEXT(OUT_NCIDS(K),varid,'long_name',9,'Time axis')
    status = NF_PUT_ATT_TEXT(OUT_NCIDS(K),varid,'time_origin',19,MODEL_START_TIME)

    status = NF_DEF_VAR(OUT_NCIDS(K),'timestp',NF_INT,1,tstepdimid,varid)
    status = NF_PUT_ATT_TEXT(OUT_NCIDS(K),varid,'units',9,'timesteps')
    status = NF_PUT_ATT_TEXT(OUT_NCIDS(K),varid,'title',10,'Time steps')

    IF (K == 7) THEN
      status = NF_PUT_ATT_REAL(OUT_NCIDS(K),varid,'tstep_sec',NF_REAL,1,FORCING_DT_REAL)
    ELSE
      status = NF_PUT_ATT_REAL(OUT_NCIDS(K),varid,'tstep_sec',NF_REAL,1,OUTPUT_DT_REAL)
    END IF
    status = NF_PUT_ATT_TEXT(OUT_NCIDS(K),varid,'long_name',14,'Time step axis')
    status = NF_PUT_ATT_TEXT(OUT_NCIDS(K),varid,'time_origin',19,MODEL_START_TIME)

    ! Initialize var index L
    L = 1

    IF (K == 1) THEN
 
      ! Q.1) General energy balance components

      status = NF_DEF_VAR(OUT_NCIDS(K),'SWnet',NF_REAL,ndims_2d,outdim_2d,VARIDS(K,L))
      status = NF_PUT_ATT_TEXT(OUT_NCIDS(K),VARIDS(K,L),'units',5,'W/m^2')
      status = NF_PUT_ATT_REAL(OUT_NCIDS(K),VARIDS(K,L),'missing_value',NF_REAL,1,NODATA)
      status = NF_PUT_ATT_REAL(OUT_NCIDS(K),VARIDS(K,L),'_FillValue',NF_REAL,1,NODATA)

      L = L + 1
      status = NF_DEF_VAR(OUT_NCIDS(K),'LWnet',NF_REAL,ndims_2d,outdim_2d,VARIDS(K,L))
      status = NF_PUT_ATT_TEXT(OUT_NCIDS(K),VARIDS(K,L),'units',5,'W/m^2')
      status = NF_PUT_ATT_REAL(OUT_NCIDS(K),VARIDS(K,L),'missing_value',NF_REAL,1,NODATA)
      status = NF_PUT_ATT_REAL(OUT_NCIDS(K),VARIDS(K,L),'_FillValue',NF_REAL,1,NODATA)

      L = L + 1
      status = NF_DEF_VAR(OUT_NCIDS(K),'Qle',NF_REAL,ndims_2d,outdim_2d,VARIDS(K,L))
      status = NF_PUT_ATT_TEXT(OUT_NCIDS(K),VARIDS(K,L),'units',5,'W/m^2')
      status = NF_PUT_ATT_REAL(OUT_NCIDS(K),VARIDS(K,L),'missing_value',NF_REAL,1,NODATA)
      status = NF_PUT_ATT_REAL(OUT_NCIDS(K),VARIDS(K,L),'_FillValue',NF_REAL,1,NODATA)

      L = L + 1
      status = NF_DEF_VAR(OUT_NCIDS(K),'Qh',NF_REAL,ndims_2d,outdim_2d,VARIDS(K,L))
      status = NF_PUT_ATT_TEXT(OUT_NCIDS(K),VARIDS(K,L),'units',5,'W/m^2')
      status = NF_PUT_ATT_REAL(OUT_NCIDS(K),VARIDS(K,L),'missing_value',NF_REAL,1,NODATA)
      status = NF_PUT_ATT_REAL(OUT_NCIDS(K),VARIDS(K,L),'_FillValue',NF_REAL,1,NODATA)

      L = L + 1
      status = NF_DEF_VAR(OUT_NCIDS(K),'Qg',NF_REAL,ndims_2d,outdim_2d,VARIDS(K,L))
      status = NF_PUT_ATT_TEXT(OUT_NCIDS(K),VARIDS(K,L),'units',5,'W/m^2')
      status = NF_PUT_ATT_REAL(OUT_NCIDS(K),VARIDS(K,L),'missing_value',NF_REAL,1,NODATA)
      status = NF_PUT_ATT_REAL(OUT_NCIDS(K),VARIDS(K,L),'_FillValue',NF_REAL,1,NODATA)

      L = L + 1
      status = NF_DEF_VAR(OUT_NCIDS(K),'Qf',NF_REAL,ndims_2d,outdim_2d,VARIDS(K,L))
      status = NF_PUT_ATT_TEXT(OUT_NCIDS(K),VARIDS(K,L),'units',5,'W/m^2')
      status = NF_PUT_ATT_REAL(OUT_NCIDS(K),VARIDS(K,L),'missing_value',NF_REAL,1,NODATA)
      status = NF_PUT_ATT_REAL(OUT_NCIDS(K),VARIDS(K,L),'_FillValue',NF_REAL,1,NODATA)

      L = L + 1
      status = NF_DEF_VAR(OUT_NCIDS(K),'Qv',NF_REAL,ndims_2d,outdim_2d,VARIDS(K,L))
      status = NF_PUT_ATT_TEXT(OUT_NCIDS(K),VARIDS(K,L),'units',5,'W/m^2')
      status = NF_PUT_ATT_REAL(OUT_NCIDS(K),VARIDS(K,L),'missing_value',NF_REAL,1,NODATA)
      status = NF_PUT_ATT_REAL(OUT_NCIDS(K),VARIDS(K,L),'_FillValue',NF_REAL,1,NODATA)

      L = L + 1
      status = NF_DEF_VAR(OUT_NCIDS(K),'Qa',NF_REAL,ndims_2d,outdim_2d,VARIDS(K,L))
      status = NF_PUT_ATT_TEXT(OUT_NCIDS(K),VARIDS(K,L),'units',5,'W/m^2')
      status = NF_PUT_ATT_REAL(OUT_NCIDS(K),VARIDS(K,L),'missing_value', NF_REAL,1,NODATA)
      status = NF_PUT_ATT_REAL(OUT_NCIDS(K),VARIDS(K,L),'_FillValue', NF_REAL,1,NODATA)

      L = L + 1
      status = NF_DEF_VAR(OUT_NCIDS(K),'DelSurfHeat',NF_REAL,ndims_2d,outdim_2d,VARIDS(K,L))
      status = NF_PUT_ATT_TEXT(OUT_NCIDS(K),VARIDS(K,L),'units',5,'J/m^2')
      status = NF_PUT_ATT_REAL(OUT_NCIDS(K),VARIDS(K,L),'missing_value', NF_REAL,1,NODATA)
      status = NF_PUT_ATT_REAL(OUT_NCIDS(K),VARIDS(K,L),'_FillValue', NF_REAL,1,NODATA)

      L = L + 1
      status = NF_DEF_VAR(OUT_NCIDS(K),'DelColdCont',NF_REAL,ndims_2d,outdim_2d,VARIDS(K,L))
      status = NF_PUT_ATT_TEXT(OUT_NCIDS(K),VARIDS(K,L),'units',5,'J/m^2')
      status = NF_PUT_ATT_REAL(OUT_NCIDS(K),VARIDS(K,L),'missing_value', NF_REAL,1,NODATA)
      status = NF_PUT_ATT_REAL(OUT_NCIDS(K),VARIDS(K,L),'_FillValue', NF_REAL,1,NODATA)

    ELSE IF (K == 2) THEN

      ! Q.2) General water balance components

      status = NF_DEF_VAR(OUT_NCIDS(K),'Snowf',NF_REAL,ndims_2d,outdim_2d, VARIDS(K,L))
      status = NF_PUT_ATT_TEXT(OUT_NCIDS(K),VARIDS(K,L),'units',7,'kg/m^2s')
      status = NF_PUT_ATT_REAL(OUT_NCIDS(K),VARIDS(K,L),'missing_value', NF_REAL,1,NODATA)
      status = NF_PUT_ATT_REAL(OUT_NCIDS(K),VARIDS(K,L),'_FillValue', NF_REAL,1,NODATA)

      L = L + 1
      status = NF_DEF_VAR(OUT_NCIDS(K),'Rainf',NF_REAL,ndims_2d,outdim_2d, VARIDS(K,L))
      status = NF_PUT_ATT_TEXT(OUT_NCIDS(K),VARIDS(K,L),'units',7,'kg/m^2s')
      status = NF_PUT_ATT_REAL(OUT_NCIDS(K),VARIDS(K,L),'missing_value', NF_REAL,1,NODATA)
      status = NF_PUT_ATT_REAL(OUT_NCIDS(K),VARIDS(K,L),'_FillValue', NF_REAL,1,NODATA)

      L = L + 1
      status = NF_DEF_VAR(OUT_NCIDS(K),'Evap',NF_REAL,ndims_2d,outdim_2d,VARIDS(K,L))
      status = NF_PUT_ATT_TEXT(OUT_NCIDS(K),VARIDS(K,L),'units',7,'kg/m^2s')
      status = NF_PUT_ATT_REAL(OUT_NCIDS(K),VARIDS(K,L),'missing_value', NF_REAL,1,NODATA)
      status = NF_PUT_ATT_REAL(OUT_NCIDS(K),VARIDS(K,L),'_FillValue', NF_REAL,1,NODATA)

      L = L + 1
      status = NF_DEF_VAR(OUT_NCIDS(K),'Qs',NF_REAL,ndims_2d,outdim_2d,VARIDS(K,L))
      status = NF_PUT_ATT_TEXT(OUT_NCIDS(K),VARIDS(K,L),'units',7,'kg/m^2s')
      status = NF_PUT_ATT_REAL(OUT_NCIDS(K),VARIDS(K,L),'missing_value', NF_REAL,1,NODATA)
      status = NF_PUT_ATT_REAL(OUT_NCIDS(K),VARIDS(K,L),'_FillValue', NF_REAL,1,NODATA)

      L = L + 1
      status = NF_DEF_VAR(OUT_NCIDS(K),'Qsb',NF_REAL,ndims_2d,outdim_2d,VARIDS(K,L))
      status = NF_PUT_ATT_TEXT(OUT_NCIDS(K),VARIDS(K,L),'units',7,'kg/m^2s')
      status = NF_PUT_ATT_REAL(OUT_NCIDS(K),VARIDS(K,L),'missing_value', NF_REAL,1,NODATA)
      status = NF_PUT_ATT_REAL(OUT_NCIDS(K),VARIDS(K,L),'_FillValue', NF_REAL,1,NODATA)

      L = L + 1
      status = NF_DEF_VAR(OUT_NCIDS(K),'Qsm',NF_REAL,ndims_2d,outdim_2d,VARIDS(K,L))
      status = NF_PUT_ATT_TEXT(OUT_NCIDS(K),VARIDS(K,L),'units',7,'kg/m^2s')
      status = NF_PUT_ATT_REAL(OUT_NCIDS(K),VARIDS(K,L),'missing_value', NF_REAL,1,NODATA)
      status = NF_PUT_ATT_REAL(OUT_NCIDS(K),VARIDS(K,L),'_FillValue', NF_REAL,1,NODATA)

      L = L + 1
      status = NF_DEF_VAR(OUT_NCIDS(K),'Qfz',NF_REAL,ndims_2d,outdim_2d,VARIDS(K,L))
      status = NF_PUT_ATT_TEXT(OUT_NCIDS(K),VARIDS(K,L),'units',7,'kg/m^2s')
      status = NF_PUT_ATT_REAL(OUT_NCIDS(K),VARIDS(K,L),'missing_value', NF_REAL,1,NODATA)
      status = NF_PUT_ATT_REAL(OUT_NCIDS(K),VARIDS(K,L),'_FillValue', NF_REAL,1,NODATA)

      L = L + 1
      status = NF_DEF_VAR(OUT_NCIDS(K),'Qst',NF_REAL,ndims_2d,outdim_2d,VARIDS(K,L))
      status = NF_PUT_ATT_TEXT(OUT_NCIDS(K),VARIDS(K,L),'units',7,'kg/m^2s')
      status = NF_PUT_ATT_REAL(OUT_NCIDS(K),VARIDS(K,L),'missing_value', NF_REAL,1,NODATA)
      status = NF_PUT_ATT_REAL(OUT_NCIDS(K),VARIDS(K,L),'_FillValue', NF_REAL,1,NODATA)

      L = L + 1
      status = NF_DEF_VAR(OUT_NCIDS(K),'DelSoilMoist',NF_REAL,ndims_2d,outdim_2d, VARIDS(K,L))
      status = NF_PUT_ATT_TEXT(OUT_NCIDS(K),VARIDS(K,L),'units',6, 'kg/m^2')
      status = NF_PUT_ATT_REAL(OUT_NCIDS(K),VARIDS(K,L),'missing_value', NF_REAL,1,NODATA)
      status = NF_PUT_ATT_REAL(OUT_NCIDS(K),VARIDS(K,L),'_FillValue', NF_REAL,1,NODATA)

      L = L + 1
      status = NF_DEF_VAR(OUT_NCIDS(K),'DelSWE',NF_REAL,ndims_2d,outdim_2d,VARIDS(K,L))
      status = NF_PUT_ATT_TEXT(OUT_NCIDS(K),VARIDS(K,L),'units',6,'kg/m^2')
      status = NF_PUT_ATT_REAL(OUT_NCIDS(K),VARIDS(K,L),'missing_value', NF_REAL,1,NODATA)
      status = NF_PUT_ATT_REAL(OUT_NCIDS(K),VARIDS(K,L),'_FillValue', NF_REAL,1,NODATA)

      L = L + 1
      status = NF_DEF_VAR(OUT_NCIDS(K),'DelSurfStor',NF_REAL,ndims_2d,outdim_2d, VARIDS(K,L))
      status = NF_PUT_ATT_TEXT(OUT_NCIDS(K),VARIDS(K,L),'units',6, 'kg/m^2')
      status = NF_PUT_ATT_REAL(OUT_NCIDS(K),VARIDS(K,L),'missing_value', NF_REAL,1,NODATA)
      status = NF_PUT_ATT_REAL(OUT_NCIDS(K),VARIDS(K,L),'_FillValue', NF_REAL,1,NODATA)

      L = L + 1
      status = NF_DEF_VAR(OUT_NCIDS(K),'DelIntercept',NF_REAL,ndims_2d,outdim_2d, VARIDS(K,L))
      status = NF_PUT_ATT_TEXT(OUT_NCIDS(K),VARIDS(K,L),'units',6,'kg/m^2')
      status = NF_PUT_ATT_REAL(OUT_NCIDS(K),VARIDS(K,L),'missing_value', NF_REAL,1,NODATA)
      status = NF_PUT_ATT_REAL(OUT_NCIDS(K),VARIDS(K,L),'_FillValue', NF_REAL,1,NODATA)

! ************** SW MONITOR HACK - extra vars for "wb" file ***********************
      L = L + 1
      status = NF_DEF_VAR(OUT_NCIDS(K),'SWE',NF_REAL,ndims_2d,outdim_2d, VARIDS(K,L))
      status = NF_PUT_ATT_TEXT(OUT_NCIDS(K),VARIDS(K,L),'units',6,'kg/m^2')
      status = NF_PUT_ATT_REAL(OUT_NCIDS(K),VARIDS(K,L),'missing_value', NF_REAL,1,NODATA)
      status = NF_PUT_ATT_REAL(OUT_NCIDS(K),VARIDS(K,L),'_FillValue', NF_REAL,1,NODATA)

      L = L + 1
      status = NF_DEF_VAR(OUT_NCIDS(K),'SWEVeg',NF_REAL,ndims_2d,outdim_2d, VARIDS(K,L))
      status = NF_PUT_ATT_TEXT(OUT_NCIDS(K),VARIDS(K,L),'units',6,'kg/m^2')
      status = NF_PUT_ATT_REAL(OUT_NCIDS(K),VARIDS(K,L),'missing_value', NF_REAL,1,NODATA)
      status = NF_PUT_ATT_REAL(OUT_NCIDS(K),VARIDS(K,L),'_FillValue', NF_REAL,1,NODATA)

      L = L + 1
      status = NF_DEF_VAR(OUT_NCIDS(K),'SoilMoist',NF_REAL,ndims_3d,outdim_3d, VARIDS(K,L))
      status = NF_PUT_ATT_TEXT(OUT_NCIDS(K),VARIDS(K,L),'units',6,'kg/m^2')
      status = NF_PUT_ATT_REAL(OUT_NCIDS(K),VARIDS(K,L),'missing_value', NF_REAL,1,NODATA)
      status = NF_PUT_ATT_REAL(OUT_NCIDS(K),VARIDS(K,L),'_FillValue', NF_REAL,1,NODATA)

      L = L + 1
      status = NF_DEF_VAR(OUT_NCIDS(K),'SOILDEPTH',NF_REAL,ndims_3d,outdim_3d, VARIDS(K,L))
      status = NF_PUT_ATT_TEXT(OUT_NCIDS(K),VARIDS(K,L),'units',2,'mm')
      status = NF_PUT_ATT_REAL(OUT_NCIDS(K),VARIDS(K,L),'missing_value', NF_REAL,1,NODATA)
      status = NF_PUT_ATT_REAL(OUT_NCIDS(K),VARIDS(K,L),'_FillValue', NF_REAL,1,NODATA)

      L = L + 1
      status = NF_DEF_VAR(OUT_NCIDS(K),'SoilTemp',NF_REAL,ndims_3d,outdim_3d, VARIDS(K,L))
      status = NF_PUT_ATT_TEXT(OUT_NCIDS(K),VARIDS(K,L),'units',1,'K')
      status = NF_PUT_ATT_REAL(OUT_NCIDS(K),VARIDS(K,L),'missing_value', NF_REAL,1,NODATA)
      status = NF_PUT_ATT_REAL(OUT_NCIDS(K),VARIDS(K,L),'_FillValue', NF_REAL,1,NODATA)

      L = L + 1
      status = NF_DEF_VAR(OUT_NCIDS(K),'SMLiqFrac',NF_REAL,ndims_3d,outdim_3d, VARIDS(K,L))
      status = NF_PUT_ATT_TEXT(OUT_NCIDS(K),VARIDS(K,L),'units',1,'-')
      status = NF_PUT_ATT_REAL(OUT_NCIDS(K),VARIDS(K,L),'missing_value', NF_REAL,1,NODATA)
      status = NF_PUT_ATT_REAL(OUT_NCIDS(K),VARIDS(K,L),'_FillValue', NF_REAL,1,NODATA)

      L = L + 1
      status = NF_DEF_VAR(OUT_NCIDS(K),'SMFrozFrac',NF_REAL,ndims_3d,outdim_3d, VARIDS(K,L))
      status = NF_PUT_ATT_TEXT(OUT_NCIDS(K),VARIDS(K,L),'units',1,'-')
      status = NF_PUT_ATT_REAL(OUT_NCIDS(K),VARIDS(K,L),'missing_value', NF_REAL,1,NODATA)
      status = NF_PUT_ATT_REAL(OUT_NCIDS(K),VARIDS(K,L),'_FillValue', NF_REAL,1,NODATA)

      L = L + 1
      status = NF_DEF_VAR(OUT_NCIDS(K),'ECanop',NF_REAL,ndims_2d,outdim_2d, VARIDS(K,L))
      status = NF_PUT_ATT_TEXT(OUT_NCIDS(K),VARIDS(K,L),'units',7,'kg/m^2s')
      status = NF_PUT_ATT_REAL(OUT_NCIDS(K),VARIDS(K,L),'missing_value', NF_REAL,1,NODATA)
      status = NF_PUT_ATT_REAL(OUT_NCIDS(K),VARIDS(K,L),'_FillValue', NF_REAL,1,NODATA)

      L = L + 1
      status = NF_DEF_VAR(OUT_NCIDS(K),'TVeg',NF_REAL,ndims_2d,outdim_2d, VARIDS(K,L))
      status = NF_PUT_ATT_TEXT(OUT_NCIDS(K),VARIDS(K,L),'units',7,'kg/m^2s')
      status = NF_PUT_ATT_REAL(OUT_NCIDS(K),VARIDS(K,L),'missing_value', NF_REAL,1,NODATA)
      status = NF_PUT_ATT_REAL(OUT_NCIDS(K),VARIDS(K,L),'_FillValue', NF_REAL,1,NODATA)

      L = L + 1
      status = NF_DEF_VAR(OUT_NCIDS(K),'ESoil',NF_REAL,ndims_2d,outdim_2d, VARIDS(K,L))
      status = NF_PUT_ATT_TEXT(OUT_NCIDS(K),VARIDS(K,L),'units',7,'kg/m^2s')
      status = NF_PUT_ATT_REAL(OUT_NCIDS(K),VARIDS(K,L),'missing_value', NF_REAL,1,NODATA)
      status = NF_PUT_ATT_REAL(OUT_NCIDS(K),VARIDS(K,L),'_FillValue', NF_REAL,1,NODATA)

      L = L + 1
      status = NF_DEF_VAR(OUT_NCIDS(K),'SWnet',NF_REAL,ndims_2d,outdim_2d,VARIDS(K,L))
      status = NF_PUT_ATT_TEXT(OUT_NCIDS(K),VARIDS(K,L),'units',5,'W/m^2')
      status = NF_PUT_ATT_REAL(OUT_NCIDS(K),VARIDS(K,L),'missing_value',NF_REAL,1,NODATA)
      status = NF_PUT_ATT_REAL(OUT_NCIDS(K),VARIDS(K,L),'_FillValue',NF_REAL,1,NODATA)

      L = L + 1
      status = NF_DEF_VAR(OUT_NCIDS(K),'LWnet',NF_REAL,ndims_2d,outdim_2d,VARIDS(K,L))
      status = NF_PUT_ATT_TEXT(OUT_NCIDS(K),VARIDS(K,L),'units',5,'W/m^2')
      status = NF_PUT_ATT_REAL(OUT_NCIDS(K),VARIDS(K,L),'missing_value',NF_REAL,1,NODATA)
      status = NF_PUT_ATT_REAL(OUT_NCIDS(K),VARIDS(K,L),'_FillValue',NF_REAL,1,NODATA)

      L = L + 1
      status = NF_DEF_VAR(OUT_NCIDS(K),'Qle',NF_REAL,ndims_2d,outdim_2d,VARIDS(K,L))
      status = NF_PUT_ATT_TEXT(OUT_NCIDS(K),VARIDS(K,L),'units',5,'W/m^2')
      status = NF_PUT_ATT_REAL(OUT_NCIDS(K),VARIDS(K,L),'missing_value',NF_REAL,1,NODATA)
      status = NF_PUT_ATT_REAL(OUT_NCIDS(K),VARIDS(K,L),'_FillValue',NF_REAL,1,NODATA)

      L = L + 1
      status = NF_DEF_VAR(OUT_NCIDS(K),'Qh',NF_REAL,ndims_2d,outdim_2d,VARIDS(K,L))
      status = NF_PUT_ATT_TEXT(OUT_NCIDS(K),VARIDS(K,L),'units',5,'W/m^2')
      status = NF_PUT_ATT_REAL(OUT_NCIDS(K),VARIDS(K,L),'missing_value',NF_REAL,1,NODATA)
      status = NF_PUT_ATT_REAL(OUT_NCIDS(K),VARIDS(K,L),'_FillValue',NF_REAL,1,NODATA)

      L = L + 1
      status = NF_DEF_VAR(OUT_NCIDS(K),'Qg',NF_REAL,ndims_2d,outdim_2d,VARIDS(K,L))
      status = NF_PUT_ATT_TEXT(OUT_NCIDS(K),VARIDS(K,L),'units',5,'W/m^2')
      status = NF_PUT_ATT_REAL(OUT_NCIDS(K),VARIDS(K,L),'missing_value',NF_REAL,1,NODATA)
      status = NF_PUT_ATT_REAL(OUT_NCIDS(K),VARIDS(K,L),'_FillValue',NF_REAL,1,NODATA)

      L = L + 1
      status = NF_DEF_VAR(OUT_NCIDS(K),'PotEvap',NF_REAL,ndims_2d,outdim_2d, VARIDS(K,L))
      status = NF_PUT_ATT_TEXT(OUT_NCIDS(K),VARIDS(K,L),'units',7,'kg/m^2s')
      status = NF_PUT_ATT_REAL(OUT_NCIDS(K),VARIDS(K,L),'missing_value', NF_REAL,1,NODATA)
      status = NF_PUT_ATT_REAL(OUT_NCIDS(K),VARIDS(K,L),'_FillValue', NF_REAL,1,NODATA)

      L = L + 1
      status = NF_DEF_VAR(OUT_NCIDS(K),'AvgSurfT',NF_REAL,ndims_2d,outdim_2d, VARIDS(K,L))
      status = NF_PUT_ATT_TEXT(OUT_NCIDS(K),VARIDS(K,L),'units',1,'K')
      status = NF_PUT_ATT_REAL(OUT_NCIDS(K),VARIDS(K,L),'missing_value', NF_REAL,1,NODATA)
      status = NF_PUT_ATT_REAL(OUT_NCIDS(K),VARIDS(K,L),'_FillValue', NF_REAL,1,NODATA)
      
      L = L + 1
      status = NF_DEF_VAR(OUT_NCIDS(K),'SatSoil',NF_REAL,ndims_2d,outdim_2d, VARIDS(K,L))
      status = NF_PUT_ATT_TEXT(OUT_NCIDS(K),VARIDS(K,L),'units',1,'-')
      status = NF_PUT_ATT_REAL(OUT_NCIDS(K),VARIDS(K,L),'missing_value', NF_REAL,1,NODATA)
      status = NF_PUT_ATT_REAL(OUT_NCIDS(K),VARIDS(K,L),'_FillValue', NF_REAL,1,NODATA)

      L = L + 1
      status = NF_DEF_VAR(OUT_NCIDS(K),'WltSoil',NF_REAL,ndims_2d,outdim_2d, VARIDS(K,L))
      status = NF_PUT_ATT_TEXT(OUT_NCIDS(K),VARIDS(K,L),'units',1,'-')
      status = NF_PUT_ATT_REAL(OUT_NCIDS(K),VARIDS(K,L),'missing_value', NF_REAL,1,NODATA)
      status = NF_PUT_ATT_REAL(OUT_NCIDS(K),VARIDS(K,L),'_FillValue', NF_REAL,1,NODATA)

      L = L + 1
      status = NF_DEF_VAR(OUT_NCIDS(K),'RootMoist',NF_REAL,ndims_2d,outdim_2d, VARIDS(K,L))
      status = NF_PUT_ATT_TEXT(OUT_NCIDS(K),VARIDS(K,L),'units',6,'kg/m^2')
      status = NF_PUT_ATT_REAL(OUT_NCIDS(K),VARIDS(K,L),'missing_value', NF_REAL,1,NODATA)
      status = NF_PUT_ATT_REAL(OUT_NCIDS(K),VARIDS(K,L),'_FillValue', NF_REAL,1,NODATA)

      L = L + 1
      status = NF_DEF_VAR(OUT_NCIDS(K),'CanopInt',NF_REAL,ndims_2d,outdim_2d, VARIDS(K,L))
      status = NF_PUT_ATT_TEXT(OUT_NCIDS(K),VARIDS(K,L),'units',6,'kg/m^2')
      status = NF_PUT_ATT_REAL(OUT_NCIDS(K),VARIDS(K,L),'missing_value', NF_REAL,1,NODATA)
      status = NF_PUT_ATT_REAL(OUT_NCIDS(K),VARIDS(K,L),'_FillValue', NF_REAL,1,NODATA)

! ************** END SW MONITOR HACK - extra vars for "wb" file ***********************

    ELSE IF (K == 3) THEN

      ! Q.3) Surface state variables

      status = NF_DEF_VAR(OUT_NCIDS(K),'SnowT',NF_REAL,ndims_2d,outdim_2d, VARIDS(K,L))
      status = NF_PUT_ATT_TEXT(OUT_NCIDS(K),VARIDS(K,L),'units',1,'K')
      status = NF_PUT_ATT_REAL(OUT_NCIDS(K),VARIDS(K,L),'missing_value', NF_REAL,1,NODATA)
      status = NF_PUT_ATT_REAL(OUT_NCIDS(K),VARIDS(K,L),'_FillValue', NF_REAL,1,NODATA)

      L = L + 1
      status = NF_DEF_VAR(OUT_NCIDS(K),'SnowTProf',NF_REAL,ndims_2d,outdim_2d, VARIDS(K,L))
      status = NF_PUT_ATT_TEXT(OUT_NCIDS(K),VARIDS(K,L),'units',1,'K')
      status = NF_PUT_ATT_REAL(OUT_NCIDS(K),VARIDS(K,L),'missing_value', NF_REAL,1,NODATA)
      status = NF_PUT_ATT_REAL(OUT_NCIDS(K),VARIDS(K,L),'_FillValue', NF_REAL,1,NODATA)

      L = L + 1
      status = NF_DEF_VAR(OUT_NCIDS(K),'VegT',NF_REAL,ndims_2d,outdim_2d, VARIDS(K,L))
      status = NF_PUT_ATT_TEXT(OUT_NCIDS(K),VARIDS(K,L),'units',1,'K')
      status = NF_PUT_ATT_REAL(OUT_NCIDS(K),VARIDS(K,L),'missing_value', NF_REAL,1,NODATA)
      status = NF_PUT_ATT_REAL(OUT_NCIDS(K),VARIDS(K,L),'_FillValue', NF_REAL,1,NODATA)

      L = L + 1
      status = NF_DEF_VAR(OUT_NCIDS(K),'BaresoilT',NF_REAL,ndims_2d,outdim_2d, VARIDS(K,L))
      status = NF_PUT_ATT_TEXT(OUT_NCIDS(K),VARIDS(K,L),'units',1,'K')
      status = NF_PUT_ATT_REAL(OUT_NCIDS(K),VARIDS(K,L),'missing_value', NF_REAL,1,NODATA)
      status = NF_PUT_ATT_REAL(OUT_NCIDS(K),VARIDS(K,L),'_FillValue', NF_REAL,1,NODATA)

      L = L + 1
      status = NF_DEF_VAR(OUT_NCIDS(K),'AvgSurfT',NF_REAL,ndims_2d,outdim_2d, VARIDS(K,L))
      status = NF_PUT_ATT_TEXT(OUT_NCIDS(K),VARIDS(K,L),'units',1,'K')
      status = NF_PUT_ATT_REAL(OUT_NCIDS(K),VARIDS(K,L),'missing_value', NF_REAL,1,NODATA)
      status = NF_PUT_ATT_REAL(OUT_NCIDS(K),VARIDS(K,L),'_FillValue', NF_REAL,1,NODATA)

      L = L + 1
      status = NF_DEF_VAR(OUT_NCIDS(K),'RadT',NF_REAL,ndims_2d,outdim_2d, VARIDS(K,L))
      status = NF_PUT_ATT_TEXT(OUT_NCIDS(K),VARIDS(K,L),'units',1,'K')
      status = NF_PUT_ATT_REAL(OUT_NCIDS(K),VARIDS(K,L),'missing_value', NF_REAL,1,NODATA)
      status = NF_PUT_ATT_REAL(OUT_NCIDS(K),VARIDS(K,L),'_FillValue', NF_REAL,1,NODATA)

      L = L + 1
      status = NF_DEF_VAR(OUT_NCIDS(K),'Albedo',NF_REAL,ndims_2d,outdim_2d, VARIDS(K,L))
      status = NF_PUT_ATT_TEXT(OUT_NCIDS(K),VARIDS(K,L),'units',1,'-')
      status = NF_PUT_ATT_REAL(OUT_NCIDS(K),VARIDS(K,L),'missing_value', NF_REAL,1,NODATA)
      status = NF_PUT_ATT_REAL(OUT_NCIDS(K),VARIDS(K,L),'_FillValue', NF_REAL,1,NODATA)

      L = L + 1
      status = NF_DEF_VAR(OUT_NCIDS(K),'SWE',NF_REAL,ndims_2d,outdim_2d, VARIDS(K,L))
      status = NF_PUT_ATT_TEXT(OUT_NCIDS(K),VARIDS(K,L),'units',6,'kg/m^2')
      status = NF_PUT_ATT_REAL(OUT_NCIDS(K),VARIDS(K,L),'missing_value', NF_REAL,1,NODATA)
      status = NF_PUT_ATT_REAL(OUT_NCIDS(K),VARIDS(K,L),'_FillValue', NF_REAL,1,NODATA)

      L = L + 1
      status = NF_DEF_VAR(OUT_NCIDS(K),'SWEVeg',NF_REAL,ndims_2d,outdim_2d, VARIDS(K,L))
      status = NF_PUT_ATT_TEXT(OUT_NCIDS(K),VARIDS(K,L),'units',6,'kg/m^2')
      status = NF_PUT_ATT_REAL(OUT_NCIDS(K),VARIDS(K,L),'missing_value', NF_REAL,1,NODATA)
      status = NF_PUT_ATT_REAL(OUT_NCIDS(K),VARIDS(K,L),'_FillValue', NF_REAL,1,NODATA)

      L = L + 1
      status = NF_DEF_VAR(OUT_NCIDS(K),'SurfStor',NF_REAL,ndims_2d,outdim_2d, VARIDS(K,L))
      status = NF_PUT_ATT_TEXT(OUT_NCIDS(K),VARIDS(K,L),'units',6,'kg/m^2')
      status = NF_PUT_ATT_REAL(OUT_NCIDS(K),VARIDS(K,L),'missing_value', NF_REAL,1,NODATA)
      status = NF_PUT_ATT_REAL(OUT_NCIDS(K),VARIDS(K,L),'_FillValue', NF_REAL,1,NODATA)


    ELSE IF (K == 4) THEN

      ! Q.4) Subsurface state variables

      status = NF_DEF_VAR(OUT_NCIDS(K),'SoilMoist',NF_REAL,ndims_3d,outdim_3d, VARIDS(K,L))
      status = NF_PUT_ATT_TEXT(OUT_NCIDS(K),VARIDS(K,L),'units',6,'kg/m^2')
      status = NF_PUT_ATT_REAL(OUT_NCIDS(K),VARIDS(K,L),'missing_value', NF_REAL,1,NODATA)
      status = NF_PUT_ATT_REAL(OUT_NCIDS(K),VARIDS(K,L),'_FillValue', NF_REAL,1,NODATA)

      L = L + 1
      status = NF_DEF_VAR(OUT_NCIDS(K),'SoilTemp',NF_REAL,ndims_3d,outdim_3d, VARIDS(K,L))
      status = NF_PUT_ATT_TEXT(OUT_NCIDS(K),VARIDS(K,L),'units',6,'K')
      status = NF_PUT_ATT_REAL(OUT_NCIDS(K),VARIDS(K,L),'missing_value', NF_REAL,1,NODATA)
      status = NF_PUT_ATT_REAL(OUT_NCIDS(K),VARIDS(K,L),'_FillValue', NF_REAL,1,NODATA)

      L = L + 1
      status = NF_DEF_VAR(OUT_NCIDS(K),'SMLiqFrac',NF_REAL,ndims_3d,outdim_3d, VARIDS(K,L))
      status = NF_PUT_ATT_TEXT(OUT_NCIDS(K),VARIDS(K,L),'units',6,'-')
      status = NF_PUT_ATT_REAL(OUT_NCIDS(K),VARIDS(K,L),'missing_value', NF_REAL,1,NODATA)
      status = NF_PUT_ATT_REAL(OUT_NCIDS(K),VARIDS(K,L),'_FillValue', NF_REAL,1,NODATA)

      L = L + 1
      status = NF_DEF_VAR(OUT_NCIDS(K),'SMFrozFrac',NF_REAL,ndims_3d,outdim_3d, VARIDS(K,L))
      status = NF_PUT_ATT_TEXT(OUT_NCIDS(K),VARIDS(K,L),'units',6,'-')
      status = NF_PUT_ATT_REAL(OUT_NCIDS(K),VARIDS(K,L),'missing_value', NF_REAL,1,NODATA)
      status = NF_PUT_ATT_REAL(OUT_NCIDS(K),VARIDS(K,L),'_FillValue', NF_REAL,1,NODATA)

      L = L + 1
      status = NF_DEF_VAR(OUT_NCIDS(K),'SoilWet',NF_REAL,ndims_2d,outdim_2d, VARIDS(K,L))
      status = NF_PUT_ATT_TEXT(OUT_NCIDS(K),VARIDS(K,L),'units',1,'-')
      status = NF_PUT_ATT_REAL(OUT_NCIDS(K),VARIDS(K,L),'missing_value', NF_REAL,1,NODATA)
      status = NF_PUT_ATT_REAL(OUT_NCIDS(K),VARIDS(K,L),'_FillValue', NF_REAL,1,NODATA)

    ELSE IF (K == 5) THEN

      ! Q.5) Evaporation components

      status = NF_DEF_VAR(OUT_NCIDS(K),'PotEvap',NF_REAL,ndims_2d,outdim_2d, VARIDS(K,L))
      status = NF_PUT_ATT_TEXT(OUT_NCIDS(K),VARIDS(K,L),'units',7,'kg/m^2s')
      status = NF_PUT_ATT_REAL(OUT_NCIDS(K),VARIDS(K,L),'missing_value', NF_REAL,1,NODATA)
      status = NF_PUT_ATT_REAL(OUT_NCIDS(K),VARIDS(K,L),'_FillValue', NF_REAL,1,NODATA)

      L = L + 1
      status = NF_DEF_VAR(OUT_NCIDS(K),'ECanop',NF_REAL,ndims_2d,outdim_2d, VARIDS(K,L))
      status = NF_PUT_ATT_TEXT(OUT_NCIDS(K),VARIDS(K,L),'units',7,'kg/m^2s')
      status = NF_PUT_ATT_REAL(OUT_NCIDS(K),VARIDS(K,L),'missing_value', NF_REAL,1,NODATA)
      status = NF_PUT_ATT_REAL(OUT_NCIDS(K),VARIDS(K,L),'_FillValue', NF_REAL,1,NODATA)

      L = L + 1
      status = NF_DEF_VAR(OUT_NCIDS(K),'TVeg',NF_REAL,ndims_2d,outdim_2d, VARIDS(K,L))
      status = NF_PUT_ATT_TEXT(OUT_NCIDS(K),VARIDS(K,L),'units',7,'kg/m^2s')
      status = NF_PUT_ATT_REAL(OUT_NCIDS(K),VARIDS(K,L),'missing_value', NF_REAL,1,NODATA)
      status = NF_PUT_ATT_REAL(OUT_NCIDS(K),VARIDS(K,L),'_FillValue', NF_REAL,1,NODATA)

      L = L + 1
      status = NF_DEF_VAR(OUT_NCIDS(K),'ESoil',NF_REAL,ndims_2d,outdim_2d, VARIDS(K,L))
      status = NF_PUT_ATT_TEXT(OUT_NCIDS(K),VARIDS(K,L),'units',7,'kg/m^2s')
      status = NF_PUT_ATT_REAL(OUT_NCIDS(K),VARIDS(K,L),'missing_value', NF_REAL,1,NODATA)
      status = NF_PUT_ATT_REAL(OUT_NCIDS(K),VARIDS(K,L),'_FillValue', NF_REAL,1,NODATA)

      L = L + 1
      status = NF_DEF_VAR(OUT_NCIDS(K),'EWater',NF_REAL,ndims_2d,outdim_2d, VARIDS(K,L))
      status = NF_PUT_ATT_TEXT(OUT_NCIDS(K),VARIDS(K,L),'units',7,'kg/m^2s')
      status = NF_PUT_ATT_REAL(OUT_NCIDS(K),VARIDS(K,L),'missing_value', NF_REAL,1,NODATA)
      status = NF_PUT_ATT_REAL(OUT_NCIDS(K),VARIDS(K,L),'_FillValue', NF_REAL,1,NODATA)

      L = L + 1
      status = NF_DEF_VAR(OUT_NCIDS(K),'RootMoist',NF_REAL,ndims_2d,outdim_2d, VARIDS(K,L))
      status = NF_PUT_ATT_TEXT(OUT_NCIDS(K),VARIDS(K,L),'units',6,'kg/m^2')
      status = NF_PUT_ATT_REAL(OUT_NCIDS(K),VARIDS(K,L),'missing_value', NF_REAL,1,NODATA)
      status = NF_PUT_ATT_REAL(OUT_NCIDS(K),VARIDS(K,L),'_FillValue', NF_REAL,1,NODATA)

      L = L + 1
      status = NF_DEF_VAR(OUT_NCIDS(K),'CanopInt',NF_REAL,ndims_2d,outdim_2d, VARIDS(K,L))
      status = NF_PUT_ATT_TEXT(OUT_NCIDS(K),VARIDS(K,L),'units',6,'kg/m^2')
      status = NF_PUT_ATT_REAL(OUT_NCIDS(K),VARIDS(K,L),'missing_value', NF_REAL,1,NODATA)
      status = NF_PUT_ATT_REAL(OUT_NCIDS(K),VARIDS(K,L),'_FillValue', NF_REAL,1,NODATA)

      L = L + 1
      status = NF_DEF_VAR(OUT_NCIDS(K),'EvapSnow',NF_REAL,ndims_2d,outdim_2d, VARIDS(K,L))
      status = NF_PUT_ATT_TEXT(OUT_NCIDS(K),VARIDS(K,L),'units',7,'kg/m^2s')
      status = NF_PUT_ATT_REAL(OUT_NCIDS(K),VARIDS(K,L),'missing_value', NF_REAL,1,NODATA)
      status = NF_PUT_ATT_REAL(OUT_NCIDS(K),VARIDS(K,L),'_FillValue', NF_REAL,1,NODATA)

      L = L + 1
      status = NF_DEF_VAR(OUT_NCIDS(K),'SubSnow',NF_REAL,ndims_2d,outdim_2d, VARIDS(K,L))
      status = NF_PUT_ATT_TEXT(OUT_NCIDS(K),VARIDS(K,L),'units',7,'kg/m^2s')
      status = NF_PUT_ATT_REAL(OUT_NCIDS(K),VARIDS(K,L),'missing_value', NF_REAL,1,NODATA)
      status = NF_PUT_ATT_REAL(OUT_NCIDS(K),VARIDS(K,L),'_FillValue', NF_REAL,1,NODATA)

      L = L + 1
      status = NF_DEF_VAR(OUT_NCIDS(K),'SubSurf',NF_REAL,ndims_2d,outdim_2d, VARIDS(K,L))
      status = NF_PUT_ATT_TEXT(OUT_NCIDS(K),VARIDS(K,L),'units',7,'kg/m^2s')
      status = NF_PUT_ATT_REAL(OUT_NCIDS(K),VARIDS(K,L),'missing_value', NF_REAL,1,NODATA)
      status = NF_PUT_ATT_REAL(OUT_NCIDS(K),VARIDS(K,L),'_FillValue', NF_REAL,1,NODATA)

      L = L + 1
      status = NF_DEF_VAR(OUT_NCIDS(K),'ACond',NF_REAL,ndims_2d,outdim_2d, VARIDS(K,L))
      status = NF_PUT_ATT_TEXT(OUT_NCIDS(K),VARIDS(K,L),'units',3,'m/s')
      status = NF_PUT_ATT_REAL(OUT_NCIDS(K),VARIDS(K,L),'missing_value', NF_REAL,1,NODATA)
      status = NF_PUT_ATT_REAL(OUT_NCIDS(K),VARIDS(K,L),'_FillValue', NF_REAL,1,NODATA)

    ELSE IF (K == 6) THEN

      ! Q.7) Cold season processes

      status = NF_DEF_VAR(OUT_NCIDS(K),'SnowFrac',NF_REAL,ndims_2d,outdim_2d, VARIDS(K,L))
      status = NF_PUT_ATT_TEXT(OUT_NCIDS(K),VARIDS(K,L),'units',1,'-')
      status = NF_PUT_ATT_REAL(OUT_NCIDS(K),VARIDS(K,L),'missing_value', NF_REAL,1,NODATA)
      status = NF_PUT_ATT_REAL(OUT_NCIDS(K),VARIDS(K,L),'_FillValue', NF_REAL,1,NODATA)

      L = L + 1
      status = NF_DEF_VAR(OUT_NCIDS(K),'Fdepth',NF_REAL,ndims_2d,outdim_2d, VARIDS(K,L))
      status = NF_PUT_ATT_TEXT(OUT_NCIDS(K),VARIDS(K,L),'units',1,'m')
      status = NF_PUT_ATT_REAL(OUT_NCIDS(K),VARIDS(K,L),'missing_value', NF_REAL,1,NODATA)
      status = NF_PUT_ATT_REAL(OUT_NCIDS(K),VARIDS(K,L),'_FillValue', NF_REAL,1,NODATA)

      L = L + 1
      status = NF_DEF_VAR(OUT_NCIDS(K),'Tdepth',NF_REAL,ndims_2d,outdim_2d, VARIDS(K,L))
      status = NF_PUT_ATT_TEXT(OUT_NCIDS(K),VARIDS(K,L),'units',1,'m')
      status = NF_PUT_ATT_REAL(OUT_NCIDS(K),VARIDS(K,L),'missing_value', NF_REAL,1,NODATA)
      status = NF_PUT_ATT_REAL(OUT_NCIDS(K),VARIDS(K,L),'_FillValue', NF_REAL,1,NODATA)

      L = L + 1
      status = NF_DEF_VAR(OUT_NCIDS(K),'SAlbedo',NF_REAL,ndims_2d,outdim_2d, VARIDS(K,L))
      status = NF_PUT_ATT_TEXT(OUT_NCIDS(K),VARIDS(K,L),'units',1,'-')
      status = NF_PUT_ATT_REAL(OUT_NCIDS(K),VARIDS(K,L),'missing_value', NF_REAL,1,NODATA)
      status = NF_PUT_ATT_REAL(OUT_NCIDS(K),VARIDS(K,L),'_FillValue', NF_REAL,1,NODATA)

      L = L + 1
      status = NF_DEF_VAR(OUT_NCIDS(K),'SnowDepth',NF_REAL,ndims_2d,outdim_2d, VARIDS(K,L))
      status = NF_PUT_ATT_TEXT(OUT_NCIDS(K),VARIDS(K,L),'units',2,'mm')
      status = NF_PUT_ATT_REAL(OUT_NCIDS(K),VARIDS(K,L),'missing_value', NF_REAL,1,NODATA)
      status = NF_PUT_ATT_REAL(OUT_NCIDS(K),VARIDS(K,L),'_FillValue', NF_REAL,1,NODATA)

      L = L + 1
      status = NF_DEF_VAR(OUT_NCIDS(K),'SnowTProf',NF_REAL,ndims_2d,outdim_2d, VARIDS(K,L))
      status = NF_PUT_ATT_TEXT(OUT_NCIDS(K),VARIDS(K,L),'units',1,'K')
      status = NF_PUT_ATT_REAL(OUT_NCIDS(K),VARIDS(K,L),'missing_value', NF_REAL,1,NODATA)
      status = NF_PUT_ATT_REAL(OUT_NCIDS(K),VARIDS(K,L),'_FillValue', NF_REAL,1,NODATA)

      L = L + 1
      status = NF_DEF_VAR(OUT_NCIDS(K),'SliqFrac',NF_REAL,ndims_2d,outdim_2d, VARIDS(K,L))
      status = NF_PUT_ATT_TEXT(OUT_NCIDS(K),VARIDS(K,L),'units',1,'-')
      status = NF_PUT_ATT_REAL(OUT_NCIDS(K),VARIDS(K,L),'missing_value', NF_REAL,1,NODATA)
      status = NF_PUT_ATT_REAL(OUT_NCIDS(K),VARIDS(K,L),'_FillValue', NF_REAL,1,NODATA)

    ELSE IF (K == 7) THEN

      ! Potential Evap for use by SAC/SNOW17 model

      status = NF_DEF_VAR(OUT_NCIDS(K),'PotEvap',NF_REAL,ndims_3d,outdim_3d_pe, VARIDS(K,L))
      status = NF_PUT_ATT_TEXT(OUT_NCIDS(K),VARIDS(K,L),'units',7,'kg/m^2s')
      status = NF_PUT_ATT_REAL(OUT_NCIDS(K),VARIDS(K,L),'missing_value', NF_REAL,1,NODATA)
      status = NF_PUT_ATT_REAL(OUT_NCIDS(K),VARIDS(K,L),'_FillValue', NF_REAL,1,NODATA)

    END IF

    status = NF_ENDDEF(OUT_NCIDS(K))

    ! Write row, col, etc.

    IF (COMP_OUTPUT) THEN
      land_idx = 0
      DO I = 1,ylen
        DO J = 1,xlen
          IF (LANDMASK(J,I) == 1) THEN
            land_idx = land_idx + 1
            temp_int(land_idx) = ROW(J,I)
          END IF
        END DO
      END DO
      status = NF_PUT_VAR_INT(OUT_NCIDS(K),rowvarid,temp_int)
      land_idx = 0
      DO I = 1,ylen
        DO J = 1,xlen
          IF (LANDMASK(J,I) == 1) THEN
            land_idx = land_idx + 1
            temp_int(land_idx) = COL(J,I)
          END IF
        END DO
      END DO
      status = NF_PUT_VAR_INT(OUT_NCIDS(K),colvarid,temp_int)
      land_idx = 0
      DO I = 1,ylen
        DO J = 1,xlen
          IF (LANDMASK(J,I) == 1) THEN
            land_idx = land_idx + 1
            temp(land_idx) = LON(J,I)
          END IF
        END DO
      END DO
      status = NF_PUT_VAR_REAL(OUT_NCIDS(K),lonvarid,temp)
      land_idx = 0
      DO I = 1,ylen
        DO J = 1,xlen
          IF (LANDMASK(J,I) == 1) THEN
            land_idx = land_idx + 1
            temp(land_idx) = LAT(J,I)
          END IF
        END DO
      END DO
      status = NF_PUT_VAR_REAL(OUT_NCIDS(K),latvarid,temp)
      land_idx = 0
      DO I = 1,ylen
        DO J = 1,xlen
          IF (LANDMASK(J,I) == 1) THEN
            land_idx = land_idx + 1
            temp_int(land_idx) = LANDMASK(J,I)
          END IF
        END DO
      END DO
      status = NF_PUT_VAR_INT(OUT_NCIDS(K),landmaskvarid,temp_int)
      land_idx = 0
      DO I = 1,ylen
        DO J = 1,xlen
          IF (LANDMASK(J,I) == 1) THEN
            land_idx = land_idx + 1
            temp_int(land_idx) = CELLID(J,I)
          END IF
        END DO
      END DO
      status = NF_PUT_VAR_INT(OUT_NCIDS(K),cellidvarid,temp_int)
    ELSE
      status = NF_PUT_VAR_INT(OUT_NCIDS(K),rowvarid,ROW)
      status = NF_PUT_VAR_INT(OUT_NCIDS(K),colvarid,COL)
      status = NF_PUT_VAR_REAL(OUT_NCIDS(K),lonvarid,LON)
      status = NF_PUT_VAR_REAL(OUT_NCIDS(K),latvarid,LAT)
      status = NF_PUT_VAR_INT(OUT_NCIDS(K),landmaskvarid,LANDMASK)
      status = NF_PUT_VAR_INT(OUT_NCIDS(K),cellidvarid,CELLID)
    END IF

!! ********* END SW MONITOR HACK - only write the "wb" and "pe" files **********
!! ********* END ULM HACK *****************************************************  
 END IF
!! *******************************************************

END DO

END
