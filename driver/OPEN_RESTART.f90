SUBROUTINE OPEN_RESTART()

  ! UW Land Surface Hydrology Group implementation of Noah model
  ! modified from NLDAS implementation
  ! author: Ted Bohn, tbohn@hydro.washington.edu

  ! CREATES RESTART FILE

  ! Modifications:
  ! 2007-Nov-06 Added SnowSoilT (T12).						Ben Livneh
  ! 2007-Nov-12 Added LSTSNW.							TJB
  ! 2008-May-05 Removed T12, for compatibility with NOAH 2.8.			Ben Livneh
  ! 2008-May-15 Prints error message and quits if can't open output file.	TJB
  ! 2008-Jul-15 Added SnowSoilT (TPACK).					Ben Livneh
  ! 2008-Jul-24 Added PACH20.							Ben Livneh
  ! 2008-Oct-07 Included SAC parameters for use as unified model         Ben Livneh
  ! 2011-May-2 Added the _FillValue attribute to compliment missing_value for nco tools B.L.

  ! driverMod contains definitions of all global driver variables
  USE driverMod

  IMPLICIT NONE

  ! RCS Id string, for version control
  CHARACTER*60 RCSID
  DATA RCSID/"$Id: OPEN_RESTART.f90,v 1.10 2008/07/30 22:33:03 vicadmin Exp $"/

  ! Define local variables
  INTEGER ndims,nvars,ngatts,natts,unlimited
  INTEGER xdimid,ydimid,leveldimid,banddimid,sacdimid,frzdimid
  INTEGER rowvarid,colvarid,latvarid,lonvarid,cellidvarid,landmaskvarid,varid
  CHARACTER*20 name
  INTEGER xtype, dims2d(2),dims3d(3),dims4d(4),dims4dsac(4),dims4dfrz(4)
  INTEGER I,J,NT

  ! Create RESTART file
  ! Data will not be compressed
  write(*,*)'RESTARTFILE',RESTFILE

   status  = NF_CREATE(RESTFILE,0,RESTART_NCID)
write(*,*)'40,status,',status

  write(*,*)'OPEN RESTART',RESTFILE, status ,RESTART_NCID
write(*,*)'42,status,',status


  IF ( status  .ne. NF_NOERR) THEN
write(*,*)'44,status,',status
    WRITE(*,*)'ERROR: cannot open restart file',RESTFILE
    STOP
  END IF
  write(*,*)'maxnsoil restart',maxnsoil, status ,RESTART_NCID
write(*,*)'48,status,',status
   status  = NF_DEF_DIM(RESTART_NCID,'x',xlen,xdimid)
write(*,*)'49,status,',status
   status  = NF_DEF_DIM(RESTART_NCID,'y',ylen,ydimid)
write(*,*)'50,status,',status
   status  = NF_DEF_DIM(RESTART_NCID,'z',MAXNSOIL,leveldimid)
write(*,*)'51,status,',status
   status  = NF_DEF_DIM(RESTART_NCID,'band',NBANDS,banddimid)
write(*,*)'52,status,',status
   status  = NF_DEF_DIM(RESTART_NCID,'zone',6,sacdimid)
write(*,*)'53,status,',status
   status  = NF_DEF_DIM(RESTART_NCID,'frz',10,frzdimid)
!   status  = NF_DEF_DIM(RESTART_NCID,'--',10,frzdimid)
write(*,*)'54,status,',status
  dims2d(1) = xdimid
  dims2d(2) = ydimid
  dims3d(1) = xdimid
  dims3d(2) = ydimid
  dims3d(3) = banddimid
  dims4d(1) = xdimid
  dims4d(2) = ydimid
  dims4d(3) = banddimid
  dims4d(4) = leveldimid
  dims4dsac(1) = xdimid
  dims4dsac(2) = ydimid
  dims4dsac(3) = banddimid
  dims4dsac(4) = sacdimid
  dims4dfrz(1) = xdimid
  dims4dfrz(2) = ydimid
  dims4dfrz(3) = banddimid
  dims4dfrz(4) = frzdimid

   status  = NF_PUT_ATT_TEXT(RESTART_NCID,NF_GLOBAL,'Conventions',7,'GDT 1.3')
write(*,*)'73,status,',status
   status  = NF_PUT_ATT_TEXT(RESTART_NCID,NF_GLOBAL,'file_name',31,RESTART)
write(*,*)'74,status,',status

   status  = NF_DEF_VAR(RESTART_NCID,'col',NF_INT,2,dims2d,colvarid)
write(*,*)'76,status,',status
   status  = NF_PUT_ATT_TEXT(RESTART_NCID,colvarid,'units',1,'-')
write(*,*)'77,status,',status
   status  = NF_PUT_ATT_TEXT(RESTART_NCID,colvarid,'long_name',11,'Grid Column')
write(*,*)'78,status,',status
   status  = NF_PUT_ATT_INT(RESTART_NCID,colvarid,'missing_value',NF_INT,1, NODATA_INT)
write(*,*)'79,status,',status
   status  = NF_PUT_ATT_INT(RESTART_NCID,colvarid,'_FillValue',NF_INT,1, NODATA_INT)
write(*,*)'80,status,',status


   status  = NF_DEF_VAR(RESTART_NCID,'row',NF_INT,2,dims2d,rowvarid)
write(*,*)'83,status,',status
   status  = NF_PUT_ATT_TEXT(RESTART_NCID,rowvarid,'units',1,'-')
write(*,*)'84,status,',status
   status  = NF_PUT_ATT_TEXT(RESTART_NCID,rowvarid,'long_name',8,'Grid Row')
write(*,*)'85,status,',status
   status  = NF_PUT_ATT_INT(RESTART_NCID,rowvarid,'missing_value',NF_INT,1, NODATA_INT)
write(*,*)'86,status,',status
   status  = NF_PUT_ATT_INT(RESTART_NCID,rowvarid,'_FillValue',NF_INT,1, NODATA_INT)
write(*,*)'87,status,',status


   status  = NF_DEF_VAR(RESTART_NCID,'lon',NF_REAL,2,dims2d,lonvarid)
write(*,*)'90,status,',status
   status  = NF_PUT_ATT_TEXT(RESTART_NCID,lonvarid,'units',12,'Degrees East')
write(*,*)'91,status,',status
   status  = NF_PUT_ATT_REAL(RESTART_NCID,lonvarid,'valid_min',NF_REAL,1,-180.)
write(*,*)'92,status,',status
   status  = NF_PUT_ATT_REAL(RESTART_NCID,lonvarid,'valid_max',NF_REAL,1,180.)
write(*,*)'93,status,',status
   status  = NF_PUT_ATT_TEXT(RESTART_NCID,lonvarid,'long_name',9,'Longitude')
write(*,*)'94,status,',status
   status  = NF_PUT_ATT_REAL(RESTART_NCID,lonvarid,'missing_value',NF_REAL,1, NODATA)
write(*,*)'95,status,',status
   status  = NF_PUT_ATT_REAL(RESTART_NCID,lonvarid,'_FillValue',NF_REAL,1, NODATA)
write(*,*)'96,status,',status


   status  = NF_DEF_VAR(RESTART_NCID,'lat',NF_REAL,2,dims2d,latvarid)
write(*,*)'99,status,',status
   status  = NF_PUT_ATT_TEXT(RESTART_NCID,latvarid,'units',13,'Degrees North')
write(*,*)'100,status,',status
   status  = NF_PUT_ATT_REAL(RESTART_NCID,latvarid,'valid_min',NF_REAL,1,-90.)
write(*,*)'101,status,',status
   status  = NF_PUT_ATT_REAL(RESTART_NCID,latvarid,'valid_max',NF_REAL,1,90.)
write(*,*)'102,status,',status
   status  = NF_PUT_ATT_TEXT(RESTART_NCID,latvarid,'long_name',8,'Latitude')
write(*,*)'103,status,',status
   status  = NF_PUT_ATT_REAL(RESTART_NCID,latvarid,'missing_value',NF_REAL,1, NODATA)
write(*,*)'104,status,',status
   status  = NF_PUT_ATT_REAL(RESTART_NCID,latvarid,'_FillValue',NF_REAL,1, NODATA)
write(*,*)'105,status,',status


   status  = NF_DEF_VAR(RESTART_NCID,'land',NF_INT,2,dims2d,landmaskvarid)
write(*,*)'108,status,',status
   status  = NF_PUT_ATT_TEXT(RESTART_NCID,landmaskvarid,'units',1,'-')
write(*,*)'109,status,',status
   status  = NF_PUT_ATT_TEXT(RESTART_NCID,landmaskvarid,'long_name',8,'Landmask')
write(*,*)'110,status,',status
   status  = NF_PUT_ATT_INT(RESTART_NCID,landmaskvarid,'missing_value',NF_INT,1, NODATA_INT)
write(*,*)'111,status,',status
   status  = NF_PUT_ATT_INT(RESTART_NCID,landmaskvarid,'_FillValue',NF_INT,1, NODATA_INT)
write(*,*)'112,status,',status


   status  = NF_DEF_VAR(RESTART_NCID,'CellID',NF_INT,2,dims2d,cellidvarid)
write(*,*)'115,status,',status
   status  = NF_PUT_ATT_TEXT(RESTART_NCID,cellidvarid,'units',1,'-')
write(*,*)'116,status,',status
   status  = NF_PUT_ATT_TEXT(RESTART_NCID,cellidvarid,'long_name',7,'Cell ID')
write(*,*)'117,status,',status
   status  = NF_PUT_ATT_INT(RESTART_NCID,cellidvarid,'missing_value',NF_INT,1, NODATA_INT)
write(*,*)'118,status,',status
   status  = NF_PUT_ATT_INT(RESTART_NCID,cellidvarid,'_FillValue',NF_INT,1, NODATA_INT)
write(*,*)'119,status,',status


   status  = NF_DEF_VAR(RESTART_NCID,'SMC',NF_REAL,4,dims4d,varid)
write(*,*)'122,status,',status
   status  = NF_PUT_ATT_TEXT(RESTART_NCID,varid,'units',5,'m3/m3')
write(*,*)'123,status,',status
   status  = NF_PUT_ATT_TEXT(RESTART_NCID,varid,'axis',10,'z band y x')
write(*,*)'124,status,',status
   status  = NF_PUT_ATT_TEXT(RESTART_NCID,varid,'long_name',28,'Soil Layer Moisture Fraction')
write(*,*)'125,status,',status
   status  = NF_PUT_ATT_REAL(RESTART_NCID,varid,'missing_value',NF_REAL,1, NODATA)
write(*,*)'126,status,',status
   status  = NF_PUT_ATT_REAL(RESTART_NCID,varid,'_FillValue',NF_REAL,1, NODATA)
write(*,*)'127,status,',status
   status  = NF_PUT_ATT_TEXT(RESTART_NCID,varid,'associate',10,'z band y x')
write(*,*)'128,status,',status


   status  = NF_DEF_VAR(RESTART_NCID,'SH2O',NF_REAL,4,dims4d,varid)
write(*,*)'131,status,',status
   status  = NF_PUT_ATT_TEXT(RESTART_NCID,varid,'units',5,'m3/m3')
write(*,*)'132,status,',status
   status  = NF_PUT_ATT_TEXT(RESTART_NCID,varid,'axis',10,'z band y x')
write(*,*)'133,status,',status
   status  = NF_PUT_ATT_TEXT(RESTART_NCID,varid,'long_name',29,'Liquid Soil Moisture Fraction')
write(*,*)'134,status,',status
   status  = NF_PUT_ATT_REAL(RESTART_NCID,varid,'missing_value',NF_REAL,1, NODATA)
write(*,*)'135,status,',status
   status  = NF_PUT_ATT_REAL(RESTART_NCID,varid,'_FillValue',NF_REAL,1, NODATA)
write(*,*)'136,status,',status
   status  = NF_PUT_ATT_TEXT(RESTART_NCID,varid,'associate',10,'z band y x')
write(*,*)'137,status,',status


   status  = NF_DEF_VAR(RESTART_NCID,'STC',NF_REAL,4,dims4d,varid)
write(*,*)'140,status,',status
   status  = NF_PUT_ATT_TEXT(RESTART_NCID,varid,'units',1,'K')
write(*,*)'141,status,',status
   status  = NF_PUT_ATT_TEXT(RESTART_NCID,varid,'axis',10,'z band y x')
write(*,*)'142,status,',status
   status  = NF_PUT_ATT_TEXT(RESTART_NCID,varid,'long_name',22,'Soil Layer Temperature')
write(*,*)'143,status,',status
   status  = NF_PUT_ATT_REAL(RESTART_NCID,varid,'missing_value',NF_REAL,1, NODATA)
write(*,*)'144,status,',status
   status  = NF_PUT_ATT_REAL(RESTART_NCID,varid,'_FillValue',NF_REAL,1, NODATA)
write(*,*)'145,status,',status
   status  = NF_PUT_ATT_TEXT(RESTART_NCID,varid,'associate',10,'z band y x')
write(*,*)'146,status,',status


   status  = NF_DEF_VAR(RESTART_NCID,'SACST',NF_REAL,4,dims4dsac,varid)
write(*,*)'149,status,',status
   status  = NF_PUT_ATT_TEXT(RESTART_NCID,varid,'units',2,'mm')
write(*,*)'150,status,',status
   status  = NF_PUT_ATT_TEXT(RESTART_NCID,varid,'axis',13,'zone band y x')
write(*,*)'151,status,',status
   status  = NF_PUT_ATT_TEXT(RESTART_NCID,varid,'long_name',18,'SAC Moisture state')
write(*,*)'152,status,',status
   status  = NF_PUT_ATT_REAL(RESTART_NCID,varid,'missing_value',NF_REAL,1, NODATA)
write(*,*)'153,status,',status
   status  = NF_PUT_ATT_REAL(RESTART_NCID,varid,'_FillValue',NF_REAL,1, NODATA)
write(*,*)'154,status,',status
   status  = NF_PUT_ATT_TEXT(RESTART_NCID,varid,'associate',13,'zone band y x')
write(*,*)'155,status,',status


   status  = NF_DEF_VAR(RESTART_NCID,'FRZST',NF_REAL,4,dims4dfrz,varid)
write(*,*)'158,status,',status
   status  = NF_PUT_ATT_TEXT(RESTART_NCID,varid,'units',8,'K and mm')
write(*,*)'159,status,',status
   status  = NF_PUT_ATT_TEXT(RESTART_NCID,varid,'axis',11,'-- band y x')
write(*,*)'160,status,',status
   status  = NF_PUT_ATT_TEXT(RESTART_NCID,varid,'long_name',27,'Soil temp and frozen states')
write(*,*)'161,status,',status
   status  = NF_PUT_ATT_REAL(RESTART_NCID,varid,'missing_value',NF_REAL,1, NODATA)
write(*,*)'162,status,',status
   status  = NF_PUT_ATT_REAL(RESTART_NCID,varid,'_FillValue',NF_REAL,1, NODATA)
write(*,*)'163,status,',status
   status  = NF_PUT_ATT_TEXT(RESTART_NCID,varid,'associate',11,'-- band y x')
write(*,*)'164,status,',status


   status  = NF_DEF_VAR(RESTART_NCID,'T1',NF_REAL,3,dims3d,varid)
write(*,*)'167,status,',status
   status  = NF_PUT_ATT_TEXT(RESTART_NCID,varid,'units',1,'K')
write(*,*)'168,status,',status
   status  = NF_PUT_ATT_TEXT(RESTART_NCID,varid,'axis',8,'band y x')
write(*,*)'169,status,',status
   status  = NF_PUT_ATT_TEXT(RESTART_NCID,varid,'long_name',19,'Surface Temperature')
write(*,*)'170,status,',status
   status  = NF_PUT_ATT_REAL(RESTART_NCID,varid,'missing_value',NF_REAL,1, NODATA)
write(*,*)'171,status,',status
   status  = NF_PUT_ATT_REAL(RESTART_NCID,varid,'_FillValue',NF_REAL,1, NODATA)
write(*,*)'172,status,',status
   status  = NF_PUT_ATT_TEXT(RESTART_NCID,varid,'associate',8,'band y x')
write(*,*)'173,status,',status


   status  = NF_DEF_VAR(RESTART_NCID,'TPACK',NF_REAL,3,dims3d,varid)
write(*,*)'176,status,',status
   status  = NF_PUT_ATT_TEXT(RESTART_NCID,varid,'units',1,'K')
write(*,*)'177,status,',status
   status  = NF_PUT_ATT_TEXT(RESTART_NCID,varid,'axis',8,'band y x')
write(*,*)'178,status,',status
   status  = NF_PUT_ATT_TEXT(RESTART_NCID,varid,'long_name',20,'Snowpack Temperature')
write(*,*)'179,status,',status
   status  = NF_PUT_ATT_REAL(RESTART_NCID,varid,'missing_value',NF_REAL,1, NODATA)
write(*,*)'180,status,',status
   status  = NF_PUT_ATT_REAL(RESTART_NCID,varid,'_FillValue',NF_REAL,1, NODATA)
write(*,*)'181,status,',status
   status  = NF_PUT_ATT_TEXT(RESTART_NCID,varid,'associate',8,'band y x')
write(*,*)'182,status,',status


   status  = NF_DEF_VAR(RESTART_NCID,'PACH20',NF_REAL,3,dims3d,varid)
write(*,*)'185,status,',status
   status  = NF_PUT_ATT_TEXT(RESTART_NCID,varid,'units',2,'mm')
write(*,*)'186,status,',status
   status  = NF_PUT_ATT_TEXT(RESTART_NCID,varid,'axis',8,'band y x')
write(*,*)'187,status,',status
   status  = NF_PUT_ATT_TEXT(RESTART_NCID,varid,'long_name',21,'Snowpack Liquid Water')
write(*,*)'188,status,',status
   status  = NF_PUT_ATT_REAL(RESTART_NCID,varid,'missing_value',NF_REAL,1, NODATA)
write(*,*)'189,status,',status
   status  = NF_PUT_ATT_REAL(RESTART_NCID,varid,'_FillValue',NF_REAL,1, NODATA)
write(*,*)'190,status,',status
   status  = NF_PUT_ATT_TEXT(RESTART_NCID,varid,'associate',8,'band y x')
write(*,*)'191,status,',status


   status  = NF_DEF_VAR(RESTART_NCID,'CMC',NF_REAL,3,dims3d,varid)
write(*,*)'194,status,',status
   status  = NF_PUT_ATT_TEXT(RESTART_NCID,varid,'units',5,'kg/m2')
write(*,*)'195,status,',status
   status  = NF_PUT_ATT_TEXT(RESTART_NCID,varid,'axis',8,'band y x')
write(*,*)'196,status,',status
   status  = NF_PUT_ATT_TEXT(RESTART_NCID,varid,'long_name',19,'Canopy Interception')
write(*,*)'197,status,',status
   status  = NF_PUT_ATT_REAL(RESTART_NCID,varid,'missing_value',NF_REAL,1, NODATA)
write(*,*)'198,status,',status
   status  = NF_PUT_ATT_REAL(RESTART_NCID,varid,'_FillValue',NF_REAL,1, NODATA)
write(*,*)'199,status,',status
   status  = NF_PUT_ATT_TEXT(RESTART_NCID,varid,'associate',8,'band y x')
write(*,*)'200,status,',status


   status  = NF_DEF_VAR(RESTART_NCID,'UZTWC',NF_REAL,3,dims3d,varid)
write(*,*)'203,status,',status
   status  = NF_PUT_ATT_TEXT(RESTART_NCID,varid,'units',2,'mm')
write(*,*)'204,status,',status
   status  = NF_PUT_ATT_TEXT(RESTART_NCID,varid,'axis',5,'band y x')
write(*,*)'205,status,',status
   status  = NF_PUT_ATT_TEXT(RESTART_NCID,varid,'long_name',33,'Upper zone tension water contents')
write(*,*)'206,status,',status
   status  = NF_PUT_ATT_REAL(RESTART_NCID,varid,'missing_value',NF_REAL,1, NODATA)
write(*,*)'207,status,',status
   status  = NF_PUT_ATT_REAL(RESTART_NCID,varid,'_FillValue',NF_REAL,1, NODATA)
write(*,*)'208,status,',status
   status  = NF_PUT_ATT_TEXT(RESTART_NCID,varid,'associate',5,'band y x')
write(*,*)'209,status,',status


   status  = NF_DEF_VAR(RESTART_NCID,'UZFWC',NF_REAL,3,dims3d,varid)
write(*,*)'212,status,',status
   status  = NF_PUT_ATT_TEXT(RESTART_NCID,varid,'units',2,'mm')
write(*,*)'213,status,',status
   status  = NF_PUT_ATT_TEXT(RESTART_NCID,varid,'axis',5,'band y x')
write(*,*)'214,status,',status
   status  = NF_PUT_ATT_TEXT(RESTART_NCID,varid,'long_name',30,'Upper zone free water contents')
write(*,*)'215,status,',status
   status  = NF_PUT_ATT_REAL(RESTART_NCID,varid,'missing_value',NF_REAL,1, NODATA)
write(*,*)'216,status,',status
   status  = NF_PUT_ATT_REAL(RESTART_NCID,varid,'_FillValue',NF_REAL,1, NODATA)
write(*,*)'217,status,',status
   status  = NF_PUT_ATT_TEXT(RESTART_NCID,varid,'associate',5,'band y x')
write(*,*)'218,status,',status


   status  = NF_DEF_VAR(RESTART_NCID,'LZTWC',NF_REAL,3,dims3d,varid)
write(*,*)'221,status,',status
   status  = NF_PUT_ATT_TEXT(RESTART_NCID,varid,'units',2,'mm')
write(*,*)'222,status,',status
   status  = NF_PUT_ATT_TEXT(RESTART_NCID,varid,'axis',5,'band y x')
write(*,*)'223,status,',status
   status  = NF_PUT_ATT_TEXT(RESTART_NCID,varid,'long_name',33,'Lower zone tension water contents')
write(*,*)'224,status,',status
   status  = NF_PUT_ATT_REAL(RESTART_NCID,varid,'missing_value',NF_REAL,1, NODATA)
write(*,*)'225,status,',status
   status  = NF_PUT_ATT_REAL(RESTART_NCID,varid,'_FillValue',NF_REAL,1, NODATA)
write(*,*)'226,status,',status
   status  = NF_PUT_ATT_TEXT(RESTART_NCID,varid,'associate',5,'band y x')
write(*,*)'227,status,',status


   status  = NF_DEF_VAR(RESTART_NCID,'LZFSC',NF_REAL,3,dims3d,varid)
write(*,*)'230,status,',status
   status  = NF_PUT_ATT_TEXT(RESTART_NCID,varid,'units',2,'mm')
write(*,*)'231,status,',status
   status  = NF_PUT_ATT_TEXT(RESTART_NCID,varid,'axis',5,'band y x')
write(*,*)'232,status,',status
   status  = NF_PUT_ATT_TEXT(RESTART_NCID,varid,'long_name',37,'Lower zone free supplemental contents')
write(*,*)'233,status,',status
   status  = NF_PUT_ATT_REAL(RESTART_NCID,varid,'missing_value',NF_REAL,1, NODATA)
write(*,*)'234,status,',status
   status  = NF_PUT_ATT_REAL(RESTART_NCID,varid,'_FillValue',NF_REAL,1, NODATA)
write(*,*)'235,status,',status
   status  = NF_PUT_ATT_TEXT(RESTART_NCID,varid,'associate',5,'band y x')
write(*,*)'236,status,',status


   status  = NF_DEF_VAR(RESTART_NCID,'LZFPC',NF_REAL,3,dims3d,varid)
write(*,*)'239,status,',status
   status  = NF_PUT_ATT_TEXT(RESTART_NCID,varid,'units',2,'mm')
write(*,*)'240,status,',status
   status  = NF_PUT_ATT_TEXT(RESTART_NCID,varid,'axis',5,'band y x')
write(*,*)'241,status,',status
   status  = NF_PUT_ATT_TEXT(RESTART_NCID,varid,'long_name',32,'Lower zone free primary contents')
write(*,*)'242,status,',status
   status  = NF_PUT_ATT_REAL(RESTART_NCID,varid,'missing_value',NF_REAL,1, NODATA)
write(*,*)'243,status,',status
   status  = NF_PUT_ATT_REAL(RESTART_NCID,varid,'_FillValue',NF_REAL,1, NODATA)
write(*,*)'244,status,',status
   status  = NF_PUT_ATT_TEXT(RESTART_NCID,varid,'associate',5,'band y x')
write(*,*)'245,status,',status


   status  = NF_DEF_VAR(RESTART_NCID,'ADIMC',NF_REAL,3,dims3d,varid)
write(*,*)'248,status,',status
   status  = NF_PUT_ATT_TEXT(RESTART_NCID,varid,'units',2,'mm')
write(*,*)'249,status,',status
   status  = NF_PUT_ATT_TEXT(RESTART_NCID,varid,'axis',5,'band y x')
write(*,*)'250,status,',status
   status  = NF_PUT_ATT_TEXT(RESTART_NCID,varid,'long_name',40,'Tension water contents of the ADIMP area')
write(*,*)'251,status,',status
   status  = NF_PUT_ATT_REAL(RESTART_NCID,varid,'missing_value',NF_REAL,1, NODATA)
write(*,*)'252,status,',status
   status  = NF_PUT_ATT_REAL(RESTART_NCID,varid,'_FillValue',NF_REAL,1, NODATA)
write(*,*)'253,status,',status
   status  = NF_PUT_ATT_TEXT(RESTART_NCID,varid,'associate',5,'band y x')
write(*,*)'254,status,',status


   status  = NF_DEF_VAR(RESTART_NCID,'SNOWH',NF_REAL,3,dims3d,varid)
write(*,*)'257,status,',status
   status  = NF_PUT_ATT_TEXT(RESTART_NCID,varid,'units',1,'m')
write(*,*)'258,status,',status
   status  = NF_PUT_ATT_TEXT(RESTART_NCID,varid,'axis',8,'band y x')
write(*,*)'259,status,',status
   status  = NF_PUT_ATT_TEXT(RESTART_NCID,varid,'long_name',10,'Snow Depth')
write(*,*)'260,status,',status
   status  = NF_PUT_ATT_REAL(RESTART_NCID,varid,'missing_value',NF_REAL,1, NODATA)
write(*,*)'261,status,',status
   status  = NF_PUT_ATT_REAL(RESTART_NCID,varid,'_FillValue',NF_REAL,1, NODATA)
write(*,*)'262,status,',status
   status  = NF_PUT_ATT_TEXT(RESTART_NCID,varid,'associate',8,'band y x')
write(*,*)'263,status,',status


   status  = NF_DEF_VAR(RESTART_NCID,'SNEQV',NF_REAL,3,dims3d,varid)
write(*,*)'266,status,',status
   status  = NF_PUT_ATT_TEXT(RESTART_NCID,varid,'units',1,'m')
write(*,*)'267,status,',status
   status  = NF_PUT_ATT_TEXT(RESTART_NCID,varid,'axis',8,'band y x')
write(*,*)'268,status,',status
   status  = NF_PUT_ATT_TEXT(RESTART_NCID,varid,'long_name',19,'Snow Water Equivalent')
write(*,*)'269,status,',status
   status  = NF_PUT_ATT_REAL(RESTART_NCID,varid,'missing_value',NF_REAL,1, NODATA)
write(*,*)'270,status,',status
   status  = NF_PUT_ATT_REAL(RESTART_NCID,varid,'_FillValue',NF_REAL,1, NODATA)
write(*,*)'271,status,',status
   status  = NF_PUT_ATT_TEXT(RESTART_NCID,varid,'associate',8,'band y x')
write(*,*)'272,status,',status


   status  = NF_DEF_VAR(RESTART_NCID,'SNCOVR',NF_REAL,3,dims3d,varid)
write(*,*)'275,status,',status
   status  = NF_PUT_ATT_TEXT(RESTART_NCID,varid,'units',1,'-')
write(*,*)'276,status,',status
   status  = NF_PUT_ATT_TEXT(RESTART_NCID,varid,'axis',8,'band y x')
write(*,*)'277,status,',status
   status  = NF_PUT_ATT_TEXT(RESTART_NCID,varid,'long_name',17,'Snow Cover Extent')
write(*,*)'278,status,',status
   status  = NF_PUT_ATT_REAL(RESTART_NCID,varid,'missing_value',NF_REAL,1, NODATA)
write(*,*)'279,status,',status
   status  = NF_PUT_ATT_REAL(RESTART_NCID,varid,'_FillValue',NF_REAL,1, NODATA)
write(*,*)'280,status,',status
   status  = NF_PUT_ATT_TEXT(RESTART_NCID,varid,'associate',8,'band y x')
write(*,*)'281,status,',status


   status  = NF_DEF_VAR(RESTART_NCID,'LSTSNW',NF_INT,3,dims3d,varid)
write(*,*)'284,status,',status
   status  = NF_PUT_ATT_TEXT(RESTART_NCID,varid,'units',1,'-')
write(*,*)'285,status,',status
   status  = NF_PUT_ATT_TEXT(RESTART_NCID,varid,'axis',8,'band y x')
write(*,*)'286,status,',status
   status  = NF_PUT_ATT_TEXT(RESTART_NCID,varid,'long_name',17,'Last snow counter')
write(*,*)'287,status,',status
   status  = NF_PUT_ATT_INT(RESTART_NCID,varid,'missing_value',NF_INT,1, NODATA_INT)
write(*,*)'288,status,',status
   status  = NF_PUT_ATT_INT(RESTART_NCID,varid,'_FillValue',NF_INT,1, NODATA_INT)
write(*,*)'289,status,',status
   status  = NF_PUT_ATT_TEXT(RESTART_NCID,varid,'associate',8,'band y x')
write(*,*)'290,status,',status


   status  = NF_DEF_VAR(RESTART_NCID,'CH',NF_REAL,3,dims3d,varid)
write(*,*)'293,status,',status
   status  = NF_PUT_ATT_TEXT(RESTART_NCID,varid,'units',3,'m/s')
write(*,*)'294,status,',status
   status  = NF_PUT_ATT_TEXT(RESTART_NCID,varid,'axis',8,'band y x')
write(*,*)'295,status,',status
   status  = NF_PUT_ATT_TEXT(RESTART_NCID,varid,'long_name',19,'Heat Exchange Coeff')
write(*,*)'296,status,',status
   status  = NF_PUT_ATT_REAL(RESTART_NCID,varid,'missing_value',NF_REAL,1, NODATA)
write(*,*)'297,status,',status
   status  = NF_PUT_ATT_REAL(RESTART_NCID,varid,'_FillValue',NF_REAL,1, NODATA)
write(*,*)'298,status,',status
   status  = NF_PUT_ATT_TEXT(RESTART_NCID,varid,'associate',8,'band y x')
write(*,*)'299,status,',status


   status  = NF_DEF_VAR(RESTART_NCID,'CM',NF_REAL,3,dims3d,varid)
write(*,*)'302,status,',status
   status  = NF_PUT_ATT_TEXT(RESTART_NCID,varid,'units',3,'m/s')
write(*,*)'303,status,',status
   status  = NF_PUT_ATT_TEXT(RESTART_NCID,varid,'axis',8,'band y x')
write(*,*)'304,status,',status
   status  = NF_PUT_ATT_TEXT(RESTART_NCID,varid,'long_name',19,'Momentum Exchange Coeff')
write(*,*)'305,status,',status
   status  = NF_PUT_ATT_REAL(RESTART_NCID,varid,'missing_value',NF_REAL,1, NODATA)
write(*,*)'306,status,',status
   status  = NF_PUT_ATT_REAL(RESTART_NCID,varid,'_FillValue',NF_REAL,1, NODATA)
write(*,*)'307,status,',status
   status  = NF_PUT_ATT_TEXT(RESTART_NCID,varid,'associate',8,'band y x')
write(*,*)'308,status,',status


   status  = NF_ENDDEF(RESTART_NCID)
write(*,*)'311,status,',status
  write(*,*)'ENDDEF status',RESTFILE,status,RESTART_NCID

  ! Write row, col, etc.

   status  = NF_PUT_VAR_INT(RESTART_NCID,rowvarid,ROW)
write(*,*)'316,status,',status
   status  = NF_PUT_VAR_INT(RESTART_NCID,colvarid,COL)
write(*,*)'317,status,',status
   status  = NF_PUT_VAR_REAL(RESTART_NCID,lonvarid,LON)
write(*,*)'318,status,',status
   status  = NF_PUT_VAR_REAL(RESTART_NCID,latvarid,LAT)
write(*,*)'319,status,',status
   status  = NF_PUT_VAR_INT(RESTART_NCID,landmaskvarid,LANDMASK)
write(*,*)'320,status,',status
   status  = NF_PUT_VAR_INT(RESTART_NCID,cellidvarid,CELLID)
write(*,*)'321,status,',status
!hack
!   status  = NF_CLOSE(RESTART_NCID)
write(*,*)'323,status,',status
!  write(*,*)'CLOSE status',RESTFILE,status,RESTART_NCID
!  stop
END
