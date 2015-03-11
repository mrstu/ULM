cvk  1/2008 ------------------------------------------------------------------
cvk  1/2008 Added option to normalize soil moisture
cvk  1/2008 To select this option, change the statement in the subr. header to
cvk  1/2008         PARAMETER (NORMALIZED = 1)
cvk  1/2008 Default not normalized option is
cvk  1/2008         PARAMETER (NORMALIZED = 0)
cvk  1/2008 ------------------------------------------------------------------

C MEMBER FLAND2
C  (from old member FCFLAND1)
C
      SUBROUTINE FLAND2(PXV,EDMND,TA,DT,SACST,FRZST,SACPAR,
     +     FRZPAR,NSOIL,NUPL,NSAC,IVERS,SURF,GRND,TCI,TET,
     +     SMC,SH2O,SACST_PRV,id,prflag)

CBL Removed all the snow and frozen soil stuff from the former subroutine
CBL FLAND1, since this is done within the Noah portion of the code in ULM
CBL and hence this subroutine was renamed to avoid confusion FLAND2.
cvk  1/2008  Introduced option to generate normalized SM
cvk  1/2008  To select this option, change the next parameter from 0 to 1
cvk  1/2008
cvk Normalized soil moisture content
cc      PARAMETER (NORMALIZE = 1)
cvk Not normalized soil moisture content
cc      PARAMETER (NORMALIZE = 0)

c  DT is in days here
C  DTFRZ IN SEC., IDTFRZ IS # FRZ_STEPS
C.......................................
C     THIS SUBROUTINE EXECUTES THE 'SAC-SMA ' OPERATION FOR ONE TIME
C         PERIOD.
C.......................................
C     SUBROUTINE INITIALLY WRITTEN BY. . .
C            ERIC ANDERSON - HRL     APRIL 1979     VERSION 1
C.......................................

CVK  FROZEN GROUND CHANGES
CVK  UZTWC,UZFWC,LZTWC,LZFSC,LZFPC ARE TOTAL WATER STORAGES
CVK  UZTWH,UZFWH,LZTWH,LZFSH,LZFPH ARE UNFROZEN WATER STORAGES

      PARAMETER (T0 = 273.16)
      REAL SACPAR(*),FRZPAR(*),SACST(*),FRZST(*),SACST_PRV(*)
      real smc(*),sh2o(*)
c SACPAR() is array of original SAC parameters, and FRZPAR() is array
c of frozen ground parameters and calculated constants
c SACST() and FRZST() same for states
        
c  delited real FGCO(6),ZSOIL(6),TSOIL(8),FGPM(11)      
      REAL LZTWM,LZFSM,LZFPM,LZSK,LZPK,LZTWC,LZFSC,LZFPC
      REAL LZTWH,LZFSH,LZFPH
      INTEGER PRFLAG,delint

CVK----------------------------------------------------------------
CVK_02  NEW COMMON STATEMENT FOR DESIRED SOIL LAYERS
CVK     THIS VERSION HAS HARD CODED OUTPUT SOIL LAYERS
CVK     LATER ON IT SHOULD BE CHANGED TO MAKE THEM VARIABLE 
c      INTEGER NINT/5/,NINTW/5/
CBL      INTEGER NDSINT,NDINTW, NORMALIZE
CBL      REAL TSINT(*),SWINT(*),SWHINT(*)
ck      REAL DSINT(10)/0.075,0.15,0.35,0.75,1.5,0.0,0.0,0.,0.,0./
c      REAL DSINT(10)/0.10,0.40,0.6,0.75,1.5,0.0,0.0,0.,0.,0./      
ck      REAL DSINTW(10)/0.075,0.15,0.35,0.75,1.5,0.0,0.0,0.,0.,0./
CBL      REAL DSINT(*), DSINTW(*)
c      REAL DSINTW(10)/0.10,0.40,0.6,0.75,1.5,0.0,0.0,0.,0.,0./
CBL      REAL TSTMP(10),DSMOD(10),SWTMP(10),SWHTMP(10)
c      SAVE DSINT,DSINTW,NINT,NINTW
CVK----------------------------------------------------------------      

cc      DIMENSION EPDIST(24)
C     COMMON BLOCKS
cc      COMMON/FSMPM1/UZTWM,UZFWM,UZK,PCTIM,ADIMP,RIVA,ZPERC,REXP,LZTWM,
cc     1              LZFSM,LZFPM,LZSK,LZPK,PFREE,SIDE,SAVED,PAREA
CVK
CVK      COMMON/FPMFG1/FGPM(10)
c--      COMMON/FPMFG1/itta,FGPM(15),ivers,ifrze      

CVK_02  NEW COMMON BLOCK FOR INTERPOLATED SOIL TEMP AND SOIL PARAMETERS
c--      COMMON/TSLINT/TSINT(10),NINT,SWINT(10),SWHINT(10),NINTW
c--      COMMON/FRDSTFG/SMAX,PSISAT,BRT,SWLT,QUARTZ,STYPE,NUPL,NSAC,
c--     +               RTUZ,RTLZ,DZUP,DZLOW      
cc      COMMON/FRZCNST/ FRST_FACT,CKSOIL,ZBOT
c--      COMMON/FRZCNST/ FRST_FACT,ZBOT      
      
CVK
CVK  NEW FG VERSION PARAMETERS & SOIL LAYER DEFINITION:
CVK          FGPM(1) - SOIL TEXTURE CLASS
CVK          FGPM(2) - OPTIONAL, SOIL TEMPERATURE AT THE 3M DEPTH
CVK          FGPM(3) - OPTIONAL, POROSITY OF RESIDUE LAYER
CVK          PAR(18) [if no calb=FGPM(4)] - RUNOFF REDUCTION PARAMETER 1
CVK          PAR(19) [if no calb=FGPM(5)] - RUNOFF REDUCTION PARAMETER 2
CVK          PAR(20) [if no calb=FGPM(6)] - RUNOFF REDUCTION PARAMETER 3
CVK          FGPM(6) - RUNOFF REDUCTION PARAMETER 3 (FOR ERIC'S VERSION ONLY)
CVK          FGPM(7) - NUMBER OF SOIL LAYERS 
CVK          FGPM(8)-FGPM(15) - DEPTHS OF SOIL LAYERS (M), NEGATIVE.
CVK                             FIRST LAYER (RESIDUE) DEPTH=-0.03M

c--      COMMON/FSMCO1/UZTWC,UZFWC,LZTWC,LZFSC,LZFPC,ADIMC,FGCO(6),RSUM(7),
c--     1   PPE,PSC,PTA,PWE,PSH,TSOIL(8)     

c--      COMMON/FSUMS1/SROT,SIMPVT,SRODT,SROST,SINTFT,SGWFP,SGWFS,SRECHT,
c--     1              SETT,SE1,SE3,SE4,SE5

C
C    ================================= RCS keyword statements ==========
c--      CHARACTER*68     RCSKW1,RCSKW2
c--      DATA             RCSKW1,RCSKW2 /                                 '
c--     .$Source: /fs/hsmb5/hydro/CVS_root/gorms/sac/fland1.f,v $
c--     . $',                                                             '
c--     .$Id: fland1.f,v 2.5 2008/05/14 15:48:02 zcui Exp $
c--     . $' /
C    ===================================================================

CVK  ADDITION FOR FROZEN DEPTH ESTIMATION
c--      SAVE IIPREV,FRSTPREV,FRZPREV
c--      IF(FGCO(1) .EQ. 0.) THEN
c--       IIPREV=0
c--       FRSTPREV=0.
c--       FRZPREV=0.
c--      ENDIF
CVK--------------------------------------
       
C       write(*,*) ' FRZPAR:',(frzpar(ii),ii=1,6)
c define major parameters from the arra
      if (prflag == 1) then
         WRITE(*,*) '---------TOP OF FLAND2-----------'
         WRITE(*,*) 'FLAND2: SACST: ', (SACST(i), i=1, 6 )
         WRITE(*,*) 'FLAND2: FRZST: ', (FRZST(i), i= 6, 10 )
         WRITE(*,*) 'FLAND2: SACST_PRV: ', (SACST_PRV(i), i=1, 6 )
         WRITE(*,*) 'FLAND2: SACPAR: ', (SACPAR(i), i=1, 16)
         WRITE(*,*) 'FLAND2: FRZPAR: ', (FRZPAR(i), i=1, 13)
      endif

cbl2014 initialize some accounting variables in terms of what
cbl2014 is leaving and entering the storages      
      UZTWHminus=0.
      UZFWHminus=0.
      UZTWCminus=0.
      UZFWCminus=0.
      UZTWHplus=0.
      UZFWHplus=0.
      UZTWCplus=0.
      UZFWCplus=0.
cbl2014 Lower zone
      LZTWHminus=0.
      LZFPHminus=0.
      LZFSHminus=0.
      LZTWCminus=0.
      LZFPCminus=0.
      LZFSCminus=0.
      LZTWHplus=0.
      LZFPHplus=0.
      LZFSHplus=0.
      LZTWCplus=0.
      LZFPCplus=0.
      LZFSCplus=0.

      UZTWM=SACPAR(1)
      UZFWM=SACPAR(2)
      ADIMP=SACPAR(5)
      LZTWM=SACPAR(9)
      LZFSM=SACPAR(10)
      LZFPM=SACPAR(11)
      if(prflag==1) then
         write(*,*)'UZTWM=SACPAR(1)',UZTWM,SACPAR(1)
         write(*,*)'UZFWM=SACPAR(2)',UZFWM,SACPAR(2)
         write(*,*)'ADIMP=SACPAR(5)',ADIMP,SACPAR(5)
         write(*,*)'LZTWM=SACPAR(9)',LZTWM,SACPAR(9)
         write(*,*)'LZFSM=SACPAR(10)',LZFSM,SACPAR(10)
         write(*,*)'LZFPM=SACPAR(11)',LZFPM,SACPAR(11)
      endif
      PAREA=1.0-SACPAR(4)-ADIMP
      if(prflag==1) then
         write(*,*)'PAREA=1.0-SACPAR(4)-ADIMP'
         write(*,*)'PAREA SACPAR(4) ADIMP',PAREA,SACPAR(4),ADIMP
      endif

      IF(IVERS .NE. 0) CKSL=FRZPAR(4)
      if(prflag==1) then
         write(*,*)'IF(IVERS .NE. 0) CKSL=FRZPAR(4)'
         write(*,*)'IVERS CKSL FRZPAR(4)',IVERS,CKSL,FRZPAR(4)
      endif

c define states from the array
      UZTWC=SACST(1)
      UZFWC=SACST(2)
      LZTWC=SACST(3)
      LZFSC=SACST(4)
      LZFPC=SACST(5)
      ADIMC=SACST(6)

      if(prflag==1) then
         write(*,*)'UZTWC=SACST(1)',UZTWC,SACST(1)
         write(*,*)'UZFWC=SACST(2)',UZFWC,SACST(2)
          write(*,*)'LZTWC=SACST(3)',LZTWC,SACST(3)
         write(*,*)'LZFSC=SACST(4)',LZFSC,SACST(4)
         write(*,*)'LZFPC=SACST(5)',LZFPC,SACST(5)
         write(*,*)'ADIMC=SACST(6)',ADIMC,SACST(6)
      endif

      if(prflag==1) then
         write(*,*)'UZTW',UZTWH,UZTWC,UZTWM
         write(*,*)'UZFW',UZFWH,UZFWC,UZFWM
         write(*,*)'LZTW',LZTWH,LZTWC,LZTWM
         write(*,*)'LZFP',LZFPH,LZFPC,LZFPM
         write(*,*)'LZFS',LZFSH,LZFSC,LZFSM
      endif
      if(prflag==1) then
cbl      write(*,*) 'pars - ',UZTWM,UZFWM,SACPAR(3),SACPAR(4),SACPAR(5),
cbl     &        SACPAR(6),SACPAR(7),SACPAR(8),LZTWM,LZFSM,LZFPM,
cbl     &        SACPAR(12),SACPAR(13),SACPAR(14),SACPAR(15),SACPAR(16)
cbl      write(*,*) 'start sac1 - states ', UZTWC,UZFWC,LZTWC,LZFSC,LZFPC,
cbl     &           ADIMC
cbl      write(*,*) '           - runoff ', ROIMP,SDRO,SSUR,SIF,BFS,BFP
cbl      write(*,*) '           - ET ', E1,E2,E3,E4,E5,TET
      endif
c      WRITE(*,*) 'FLAND2: UZFWC= ', UZFWC
      IF(IVERS .EQ. 0) THEN
CVK  OLD FROZEN GROUND VERSION: KEEP UNFROZEN WATER = TOTAL
C--       RUZICE=UZK
C--       RLZICE=LZSK
C--       RUZPERC=1.0
       UZTWH=UZTWC
       UZFWH=UZFWC
       LZTWH=LZTWC
       LZFSH=LZFSC
       LZFPH=LZFPC
      if(prflag==1) then
         write(*,*)'IF(IVERS .EQ. 0) THEN'
         write(*,*)'###UZTWH=UZTWC',UZTWH,UZTWC
         write(*,*)'###UZFWH=UZFWC',UZFWH,UZFWC
         write(*,*)'LZTWH=LZTWC',LZTWH,LZTWC
         write(*,*)'LZFSH=LZFSC',LZFSH,LZFSC
         write(*,*)'LZFPH=LZFPC',LZFPH,LZFPC
      endif
      ELSE        
CVK  NEW FROZEN GROUND VERSION: USE ESTIMATED UNFROZEN WATER
CVK  REDEFINE UNFROZEN WATER VARIABLES IF THEY ARE NEGATIVE
C       
       UZTWH=FRZST(6)
       IF(UZTWH .LT. 0.) UZTWH=0.
       UZFWH=FRZST(7)
       IF(UZFWH .LT. 0.) UZFWH=0.
       LZTWH=FRZST(8)
       IF(LZTWH .LT. 0.) LZTWH=0.
       LZFSH=FRZST(9)
       IF(LZFSH .LT. 0.) LZFSH=0.
       LZFPH=FRZST(10)
       IF(LZFPH .LT. 0.) LZFPH=0.
      if(prflag==1) then
         write(*,*)'ELSE'
         write(*,*)'###UZTWH=FRZST(6)',UZTWH,FRZST(6)
         write(*,*)'###UZFWH=FRZST(7)',UZFWH,FRZST(7)
         write(*,*)'LZTWH=FRZST(8)',LZTWH,FRZST(8)
         write(*,*)'LZFSH=FRZST(9)',LZFSH,FRZST(9)
         write(*,*)'LZFPH=FRZST(10)',LZFPH,FRZST(10)
      endif
CVK  RUZICE & RLZICE ARE REDUCTION OF FREE WATER MOVEMENT 
CVK  BASED ON KULIK'S THEORY: Kfrz = Kunfrz/(1+FGPM(4)*ICE)**2 
       ZUP=FRZPAR(9+NUPL)
      if (prflag==1) then
         write(*,*)'ZUP=FRZPAR(9+NUPL)'
         delint=(9+NUPL)
         write(*,*)'ZUP FRZPAR(9+NUPL)',ZUP,FRZPAR(delint)
      endif
CBL    ZLW=FRZPAR(9+NSAC) 
CBL    Changed NSAC to NSOIL since these are not equal here
       ZLW=FRZPAR(9+NSOIL)
      if (prflag==1) then
         write(*,*)'ZLW=FRZPAR(9+NSOIL)'
         delint=(9+NSOIL)
         write(*,*)'ZLW FRZPAR(9+NSOIL)',ZLW,FRZPAR(delint)
      endif
CBL       RUZICE=0.001*FRZPAR(6)*(UZTWC-UZTWH+UZFWC-UZFWH)/
CBL     +        (FRZPAR(10)-ZUP)
CBL  Changed FRZPAR(10) to 0 since surface layer not ignored here
cbl2014 change back to the above?
       RUZICE=0.001*FRZPAR(6)*(UZTWC-UZTWH+UZFWC-UZFWH)/
     +        (0.0-ZUP)
      if (prflag==1) then
         write(*,*)'RUZICE=0.001*FRZPAR(6)*(UZTWC-UZTWH+UZFWC-UZFWH)/'
         write(*,*)'+        (0.0-ZUP)'
         write(*,*)'RUZICE FRZPAR(6) UZTWC',RUZICE,FRZPAR(6),UZTWC
         write(*,*)'UZTWH UZFWC UZFWH ZUP',UZTWH,UZFWC,UZFWH,ZUP
      endif
       RLZICE=0.001*FRZPAR(7)*(LZTWC-LZTWH+LZFSC-LZFSH)/
     +        (ZUP-ZLW)
      if (prflag==1) then
         write(*,*)'RLZICE=0.001*FRZPAR(7)*(LZTWC-LZTWH+LZFSC-LZFSH)/'
         write(*,*)'+        (ZUP-ZLW)'
         write(*,*)'RLZICE FRZPAR(7) LZTWC',RLZICE,FRZPAR(7),LZTWC
         write(*,*)'LZTWH LZFSC LZFSH ZUP',LZTWH,LZFSC,LZFSH,ZUP
      endif
       RUZPERC=1.0
       if (prflag==1)write(*,*)'RUZPERC=1.0',RUZPERC
       IF(RUZICE .EQ. 0.) THEN
        RUZICE = SACPAR(3)
        if (prflag==1) then
           write(*,*)'IF(RUZICE .EQ. 0.) THEN'
           write(*,*)'RUZICE = SACPAR(3)',RUZICE,SACPAR(3)
        endif
       ELSE 
        RUZPERC=1.0/((1.+CKSL*RUZICE)**2)
        if (prflag==1) then
           write(*,*)'ELSE'
           delme=(1.+CKSL*RUZICE)
           write(*,*)'RUZPERC product',RUZPERC,delme
           write(*,*)'CKSL RUZICE',CKSL,RUZICE
        endif
        RUZICE=1.-EXP(LOG(1.-SACPAR(3))*RUZPERC)
        if (prflag==1) then
           write(*,*)'RUZICE=1.-EXP(LOG(1.-SACPAR(3))*RUZPERC)'
           delme=1.-EXP(LOG(1.-SACPAR(3)))
           write(*,*)'RUZICE explog ',RUZICE,delme
           write(*,*)'SACPAR(3) RUZPERC',SACPAR(3),RUZPERC
        endif
       ENDIF
       IF(RLZICE .EQ. 0.) THEN
        RLZICE = SACPAR(12)
        if (prflag==1) then
           write(*,*)'IF(RLZICE .EQ. 0.) THEN'
           write(*,*)'RLZICE = SACPAR(12)',RLZICE,SACPAR(12)
        endif
       ELSE  
        RLZICE=1.0/((1.+CKSL*RLZICE)**2)
        if (prflag==1) then
           write(*,*)'ELSE'
           write(*,*)'RLZICE=1.0/((1.+CKSL*RLZICE)**2)'
           write(*,*)'RLZICE CKSL RLZICE',RLZICE,CKSL,RLZICE
        endif
        RLZICE=1.-EXP(LOG(1.-SACPAR(12))*RLZICE)
        if (prflag==1) then
           write(*,*)'RLZICE=1.-EXP(LOG(1.-SACPAR(12))*RLZICE)'
           delme=1.-EXP(LOG(1.-SACPAR(12)))
           write(*,*)'RLZICE explog',RLZICE,delme
           write(*,*)'SACPAR(12) RLZICE',SACPAR(12),RLZICE
        endif
       ENDIF 
      ENDIF
c      if(uztwc .ne. uztwh) WRITE(*,*) 'ST1=',uztwc,uztwh
c      if(uzfwc .ne. uzfwh) WRITE(*,*) 'ST2=',uzfwc,uzfwh
c      if(lztwc .ne. lztwh) WRITE(*,*) 'ST3=',lztwc,lztwh
c      if(lzfsc .ne. lzfsh) WRITE(*,*) 'ST4=',lzfsc,lzfsh
c      if(lzfpc .ne. lzfph) WRITE(*,*) 'ST5=',lzfpc,lzfph

C.......................................
C     COMPUTE EVAPOTRANSPIRATION LOSS FOR THE TIME INTERVAL.
C        EDMND IS THE ET-DEMAND FOR THE TIME INTERVAL
cc      EDMND=EP*EPDIST(KINT)
cVK ADJUST EDMND FOR EFFECT OF SNOW & FOREST COVER.
c from EX1 OFS subroutine
CBL    Removed the adjustment of PET for snow and forest cover
CBL    since these are handled elsewhere in the model
CBL      EDMND=(1.-(1.0-SACPAR(17))*AESC)*EDMND
C
C     COMPUTE ET FROM UPPER ZONE.
CVK      E1=EDMND*(UZTWC/UZTWM)
CVK  ONLY UNFROZEN WATER CAN BE EVAPORATED
      E1=EDMND*(UZTWH/UZTWM)
      if (prflag==1) then
         write(*,*)'E1=EDMND*(UZTWH/UZTWM)'
         write(*,*)'E1 EDMND',E1,EDMND
         write(*,*)'UZTWH UZTWM',UZTWH,UZTWM
      endif
      RED=EDMND-E1
      if (prflag==1)then
         write(*,*)'RED=EDMND-E1'
         write(*,*)'RED EDMND E1',RED,EDMND,E1
      endif
C     RED IS RESIDUAL EVAP DEMAND
CVK      UZTWC=UZTWC-E1
      UZTWH=UZTWH-E1
      if (prflag==1)then
         write(*,*)'###UZTWH=UZTWH-E1'
         write(*,*)'###UZTWH E1',UZTWH,E1
         UZTWHminus=UZTWHminus-E1
      endif
      E2=0.0
      if (prflag==1)write(*,*)'E2=0.0',E2
CV.K      IF(UZTWC.GE.0.) THEN
      IF(UZTWH.GE.0.) THEN
CV.K    SUBTRACT ET FROM TOTAL WATER STORAGE
       UZTWC=UZTWC-E1
       if (prflag==1)then
          write(*,*)'IF(UZTWH.GE.0.) THEN'
          write(*,*)'###UZTWC=UZTWC-E1'
          write(*,*)'###UZTWC E1',UZTWC,E1
          UZTWCminus=UZTWCminus-E1
          write(*,*)'goto 220'
       endif
       GO TO 220
      ENDIF 

C     E1 CAN NOT EXCEED UZTWC
CV.K      E1=E1+UZTWC
CV.K      UZTWC=0.0
      E1=E1+UZTWH
      if (prflag==1)then
         write(*,*)'E1=E1+UZTWH'
         write(*,*)'E1 UZTWH',E1,UZTWH
      endif
      UZTWH=0.0
      if (prflag==1)write(*,*)'###UZTWH=0.0',UZTWH
CV.K   REDUCE TOTAL TENSION WATER BY ACTUAL E1
      UZTWC=UZTWC-E1
      if (prflag==1)then
         write(*,*)'###UZTWC=UZTWC-E1'
         write(*,*)'###UZTWC E1',UZTWC,E1
         UZTWCminus=UZTWCminus-E1
      endif
      IF(UZTWC .LT. 0.0) UZTWC=0.0
      if (prflag==1)then
         write(*,*)'###IF(UZTWC .LT. 0.0) UZTWC=0.0'
         write(*,*)'###UZTWC',UZTWC
      endif
      RED=EDMND-E1
      if (prflag==1)then
         write(*,*)'RED=EDMND-E1'
         write(*,*)'RED EDMND E1',RED,EDMND,E1
      endif
CV.K      IF(UZFWC.GE.RED) GO TO 221
      if (prflag==1)then
         write(*,*)'IF(UZFWH.GE.RED) GO TO 221'
         write(*,*)'UZFWH RED',UZFWH,RED
      endif
      IF(UZFWH.GE.RED) GO TO 221

C     E2 IS EVAP FROM UZFWC.
CV.K      E2=UZFWC
CV.K      UZFWC=0.0
      E2=UZFWH
      if (prflag==1)then
         write(*,*)'E2=UZFWH'
         write(*,*)'E2 UZFWH',E2,UZFWH
      endif
      UZFWH=0.0
      if (prflag==1)write(*,*)'###UZFWH=0.0',UZFWH
CV.K   REDUCE TOTAL FREE WATER BY ACTUAL E2
      UZFWC=UZFWC-E2
      if (prflag==1)then
         write(*,*)'UZFWC=UZFWC-E2'
         write(*,*)'###UZFWC E2',UZFWC,E2
         UZFWCminus=UZFWCminus-E1
      endif
      IF(UZFWC .LT. 0.0) UZFWC=0.0     
      if (prflag==1)then
         write(*,*)'###IF(UZFWC .LT. 0.0) UZFWC=0.0',UZFWC
      endif     
      RED=RED-E2
      if (prflag==1)write(*,*)'RED=RED-E2',RED,E2
      if (prflag==1)write(*,*)'GO TO 225'
      GO TO 225
  221 E2=RED
      if (prflag==1)write(*,*)'221 E2=RED',E2,RED
CVK   SUBTRACT E2 FROM TOTAL & UNFROZEN FREE WATER STORAGES
      UZFWC=UZFWC-E2
      if (prflag==1)write(*,*)'###UZFWC=UZFWC-E2',UZFWC,E2
      UZFWCminus=UZFWCminus-E2
      UZFWH=UZFWH-E2
      if (prflag==1)write(*,*)'###UZFWH=UZFWH-E2',UZFWH,E2
      UZFWHminus=UZFWHminus-E2
      RED=0.0
      if (prflag==1)write(*,*)'RED=0.0',RED
      if (prflag==1)then
         write(*,*)'220 IF((UZTWC/UZTWM).GE.(UZFWC/UZFWM)) GO TO 225'
         delme=UZTWC/UZTWM
         delme2=UZFWC/UZFWM
         write(*,*)'(UZTWC/UZTWM) (UZFWC/UZFWM',delme,delme2
      endif     
  220 IF((UZTWC/UZTWM).GE.(UZFWC/UZFWM)) GO TO 225
C     UPPER ZONE FREE WATER RATIO EXCEEDS UPPER ZONE
C     TENSION WATER RATIO, THUS TRANSFER FREE WATER TO TENSION
      UZRAT=(UZTWC+UZFWC)/(UZTWM+UZFWM)
      if (prflag==1)then
         write(*,*)'UZRAT=(UZTWC+UZFWC)/(UZTWM+UZFWM)'
         delme=UZTWC+UZFWC
         delme2=UZTWM+UZFWM
         write(*,*)'UZRAT (UZTWC+UZFWC) (UZTWM+UZFWM)'
         write(*,*)UZRAT,delme,delme2
      endif     
CV.K  ACCOUNT FOR RATIO OF UNFROZEN WATER ONLY
CV.K  AND ADJUST FOUR SOIL STATES 
CV.K      UZTWC=UZTWM*UZRAT
CV.K      UZFWC=UZFWM*UZRAT
      DUZTWC=UZTWM*UZRAT-UZTWC
      if (prflag==1)then
         write(*,*)'DUZTWC=UZTWM*UZRAT-UZTWC'
         write(*,*)'DUZTWC UZTWM UZRAT UZTWC'
         write(*,*)DUZTWC,UZTWM,UZRAT,UZTWC
      endif     
      IF(DUZTWC .GT. UZFWH) DUZTWC=UZFWH 
      if (prflag==1)then
         write(*,*)'IF(DUZTWC .GT. UZFWH) DUZTWC=UZFWH '
         write(*,*)'DUZTWC UZFWH',DUZTWC,UZFWH
         write(*,*)DUZTWC,UZTWM,UZRAT,UZTWC
      endif     
CV.K  TRANSFERED WATER CAN NOT EXCEED UNFROZEN FREE WATER
      UZTWC=UZTWC+DUZTWC
      if (prflag==1)write(*,*)'###UZTWC=UZTWC+DUZTWC',UZTWC,DUZTWC
      UZTWCplus=UZTWCplus+DUZTWC
      UZTWH=UZTWH+DUZTWC
      if (prflag==1)write(*,*)'###UZTWH=UZTWH+DUZTWC',UZTWH,DUZTWC
      UZTWHplus=UZTWHplus+DUZTWC
      UZFWC=UZFWC-DUZTWC
      if (prflag==1)write(*,*)'###UZFWC=UZFWC-DUZTWC',UZFWC,DUZTWC
      UZFWCminus=UZFWCminus-DUZTWC
      UZFWH=UZFWH-DUZTWC
      if (prflag==1)write(*,*)'###UZFWH=UZFWH-DUZTWC',UZFWH,DUZTWC
      UZFWHminus=UZFWHminus-DUZTWC
CV.K  CHECK UNFROZEN WATER STORAGES TOO
      if (prflag==1)then
         write(*,*)'###225 IF (UZTWC.LT.0.00001) THEN'
         write(*,*)'###UZTWC',UZTWC
      endif     
  225 IF (UZTWC.LT.0.00001) THEN
       UZTWC=0.0
       if (prflag==1)write(*,*)'###UZTWC=0.0',UZTWC
       UZTWH=0.0
       if (prflag==1)write(*,*)'###UZTWH=0.0',UZTWH
      ENDIF 
      if (prflag==1)write(*,*)'IF (UZFWC.LT.0.00001) THEN',UZFWC
      IF (UZFWC.LT.0.00001) THEN
         UZFWC=0.0
         if (prflag==1)write(*,*)'###UZFWC=0.0',UZFWC
         UZFWH=0.0
         if (prflag==1)write(*,*)'###UZFWH=0.0',UZFWH
      ENDIF 
C
C     COMPUTE ET FROM THE LOWER ZONE.
C     COMPUTE ET FROM LZTWC (E3)
CV.K      E3=RED*(LZTWC/(UZTWM+LZTWM))
CV.K      LZTWC=LZTWC-E3
CV.K      IF(LZTWC.GE.0.0) THEN
CV.K  ONLY UNFROZEN WATER CAN BE EVAPORATED
cbl do not allow soil moisture to be removed from lower zone
cbl via 'direct soil' evaporation. Rather, this 'deep' soil
cbl moisture will be removed via root water uptake from
cbl transpiration.
cbl      E3=RED*(LZTWH/(UZTWM+LZTWM))
      E3=0.0
      if (prflag==1)write(*,*)'E3=0.0',E3
      LZTWH=LZTWH-E3
      if(prflag==1)write(*,*)'LZTWH=LZTWH-E3',LZTWH,E3
      if(prflag==1)write(*,*)'IF(LZTWH.GE.0.0) THEN',LZTWH
      IF(LZTWH.GE.0.0) THEN
         LZTWC=LZTWC-E3
         if (prflag==1)write(*,*)'LZTWC=LZTWC-E3',LZTWC,E3
         if (prflag==1)write(*,*)'GO TO 226'
         GO TO 226
      ENDIF
       
C     E3 CAN NOT EXCEED LZTWC
CV.K      E3=E3+LZTWC
CV.K      LZTWC=0.0
      E3=E3+LZTWH
      if (prflag==1)then
         write(*,*)'E3=E3+LZTWH'
         write(*,*)'E3 LZTWH',E3,LZTWH
      endif     
      LZTWH=0.0
      if (prflag==1)write(*,*)'LZTWH=0.0',LZTWH
CV.K   REDUCE TOTAL TENSION WATER BY E3
       LZTWC=LZTWC-E3
      if (prflag==1)then
         write(*,*)'LZTWC=LZTWC-E3'
         write(*,*)'LZTWC E3',LZTWC,E3
      endif     
  226 RATLZT=LZTWC/LZTWM
      if (prflag==1)then
         write(*,*)'226 RATLZT=LZTWC/LZTWM'
         write(*,*)'RATLZT LZTWC LZTWM',RATLZT,LZTWC,LZTWM
      endif     
      RATLZ=(LZTWC+LZFPC+LZFSC-SACPAR(16))/(LZTWM+LZFPM+LZFSM
     +       -SACPAR(16))
      if (prflag==1)then
         write(*,*)'RATLZ=(LZTWC+LZFPC+LZFSC-SACPAR(16))/'
         write(*,*)'(LZTWM+LZFPM+LZFSM-SACPAR(16))'
         delme=(LZTWC+LZFPC+LZFSC-SACPAR(16))
         delme2=(LZTWM+LZFPM+LZFSM-SACPAR(16))
         write(*,*)'RATLZ numerator denominator',RATLZ,delme,delme2
         write(*,*)'LZTWC LZFPC LZFSC',LZTWC,LZFPC,LZFSC
         write(*,*)'SACPAR(16)',SACPAR(16)
         write(*,*)'LZTWM LZFPM LZFSM',LZTWM,LZFPM,LZFSM
         write(*,*)'SACPAR(16)',SACPAR(16)
      endif     
      if (prflag==1)then
         write(*,*)'IF(RATLZT.GE.RATLZ) GO TO 230'
         write(*,*)'RATLZT RATLZ',RATLZT,RATLZ
      endif     
      IF(RATLZT.GE.RATLZ) GO TO 230
C     RESUPPLY LOWER ZONE TENSION WATER FROM LOWER
C     ZONE FREE WATER IF MORE WATER AVAILABLE THERE.
      DEL=(RATLZ-RATLZT)*LZTWM
      if (prflag==1)then
         write(*,*)'DEL=(RATLZ-RATLZT)*LZTWM'
         write(*,*)'DEL RATLZ RATLZT LZTWM',DEL,RATLZ,RATLZT,LZTWM
      endif     
CV.K  ONLY UNFROZEN WATER CAN BE TRANSFERED
c       if(lzfsc .ne. lzfsh) write(*,*) 'BST4=',lzfsc,lzfsh
      SFH=LZFSH+LZFPH
      if (prflag==1)then
         write(*,*)'SFH=LZFSH+LZFPH'
         write(*,*)'SFH LZFSH LZFPH',SFH,LZFSH,LZFPH
      endif     
      IF(DEL .GT. SFH) DEL=SFH
      if (prflag==1)then
         write(*,*)'IF(DEL .GT. SFH) DEL=SFH'
         write(*,*)'DEL SFH',DEL,SFH
      endif     
      LZFSH=LZFSH-DEL
      if (prflag==1)then
         write(*,*)'LZFSH=LZFSH-DEL'
         write(*,*)'LZFSH LZFSH DEL',LZFSH,LZFSH,DEL
      endif     
      IF(LZFSH .GE. 0.0) THEN
      if (prflag==1)then
         write(*,*)'IF(LZFSH .GE. 0.0) THEN'
         write(*,*)'LZFSH',LZFSH
      endif     
C     TRANSFER FROM LZFSC TO LZTWC.      
       LZFSC=LZFSC-DEL
       if (prflag==1)then
          write(*,*)'LZFSC=LZFSC-DEL'
          write(*,*)'LZFSC LZFSC DEL',LZFSC,DEL
       endif     
c         if(lzfsc .lt. lzfsh) then
c          write(*,*) ' lzfsc1: ',lzfsc,lzfsh,del
c          stop
c         endif 
      ELSE
C     IF TRANSFER EXCEEDS LZFSC THEN REMAINDER COMES FROM LZFPC
       LZFPC=LZFPC+LZFSH
       if (prflag==1)then
          write(*,*)'LZFPC=LZFPC+LZFSH'
          write(*,*)'LZFPC LZFPC LZFSH',LZFPC,LZFPC,LZFSH
       endif     
       LZFPH=LZFPH+LZFSH
       if (prflag==1)then
          write(*,*)'LZFPH=LZFPH+LZFSH'
          write(*,*)'LZFPH LZFPH LZFSH',LZFPH,LZFSH
       endif     
       xx=LZFSH+DEL
       if (prflag==1)then
          write(*,*)'xx=LZFSH+DEL'
          write(*,*)'xx LZFSH DEL',xx,LZFSH,DEL
       endif     
       LZFSC=LZFSC-xx
       if (prflag==1)then
          write(*,*)'LZFSC=LZFSC-xx'
          write(*,*)'LZFSC LZFSC xx',LZFSC,LZFSC,xx
       endif     
c         if(lzfsc .lt. lzfsh) then
c          write(*,*) ' lzfsc2: ',lzfsc,lzfsh,del,xx
c          stop
c         endif 
       LZFSH=0.0
      if (prflag==1)write(*,*)'LZFSH=0.0',LZFSH
      ENDIF
      LZTWC=LZTWC+DEL
       if (prflag==1)then
          write(*,*)'LZTWC=LZTWC+DEL'
          write(*,*)'LZTWC LZTWC DEL',LZTWC,LZTWC,DEL
       endif     
      LZTWH=LZTWH+DEL
       if (prflag==1)then
          write(*,*)'LZTWH=LZTWH+DEL'
          write(*,*)'LZTWH LZTWH DEL',LZTWH,LZTWH,DEL
       endif     

CV.K      LZTWC=LZTWC+DEL
CV.K      LZFSC=LZFSC-DEL
CV.K      IF(LZFSC.GE.0.0) GO TO 230
CV.K      LZFPC=LZFPC+LZFSC
CV.K      LZFSC=0.0

CV.K  CHECK UNFROZEN WATER STORAGE
 230   IF (LZTWC.LT.0.00001) THEN
          LZTWC=0.0
          if (prflag==1)write(*,*)'230 IF (LZTWC.LT.0.00001) THEN'
             if (prflag==1)write(*,*)'LZTWC=0.0',LZTWC
             LZTWH=0.0 
             if (prflag==1)write(*,*)'LZTWH=0.0 ',LZTWH
          ENDIF 
C
C     COMPUTE ET FROM ADIMP AREA.-E5
      E5=E1+(RED+E2)*((ADIMC-E1-UZTWC)/(UZTWM+LZTWM))
       if (prflag==1)then
          write(*,*)'E5=E1+(RED+E2)*((ADIMC-E1-UZTWC)/(UZTWM+LZTWM))'
          write(*,*)'E5 E1 (RED+E2)',E5,E1,(RED+E2)
          delme=(ADIMC-E1-UZTWC)
          delme2=(UZTWM+LZTWM)
          delme3=((ADIMC-E1-UZTWC)/(UZTWM+LZTWM))
          write(*,*)'(ADIMC-E1-UZTWC) (UZTWM+LZTWM)',delme,delme2
          write(*,*)'((ADIMC-E1-UZTWC)/(UZTWM+LZTWM))',delme3
       endif     
C      ADJUST ADIMC,ADDITIONAL IMPERVIOUS AREA STORAGE, FOR EVAPORATION.
      ADIMC=ADIMC-E5
       if (prflag==1)then
          write(*,*)'ADIMC=ADIMC-E5'
          write(*,*)'ADIMC ADIMC E5',ADIMC,ADIMC,E5
       endif     
       if (prflag==1)write(*,*)'IF(ADIMC.GE.0.0) GO TO 231',ADIMC
      IF(ADIMC.GE.0.0) GO TO 231
C     E5 CAN NOT EXCEED ADIMC.
      E5=E5+ADIMC
       if (prflag==1)then
          write(*,*)'E5=E5+ADIMC'
          write(*,*)'E5 ADIMC',E5,ADIMC
       endif     
      ADIMC=0.0
      if (prflag==1)write(*,*)'ADIMC 0.0',ADIMC
  231 E5=E5*ADIMP
       if (prflag==1)then
          write(*,*)'231 E5=E5*ADIMP'
          write(*,*)'231 E5 ADIMP',E5,ADIMP
       endif     
C     E5 IS ET FROM THE AREA ADIMP.
C.......................................
C     COMPUTE PERCOLATION AND RUNOFF AMOUNTS.
      TWX=PXV+UZTWC-UZTWM  
       if (prflag==1)then
          write(*,*)'TWX=PXV+UZTWC-UZTWM'
          write(*,*)'TWX PXV UZTWC UZTWM',TWX,PXV,UZTWC,UZTWM
       endif          
C     TWX IS THE TIME INTERVAL AVAILABLE MOISTURE IN EXCESS
C     OF UZTW REQUIREMENTS.
       if (prflag==1)write(*,*)'IF(TWX.GE.0.0) GO TO 232',TWX
      IF(TWX.GE.0.0) GO TO 232
C     ALL MOISTURE HELD IN UZTW--NO EXCESS.
      UZTWC=UZTWC+PXV
       if (prflag==1)then
          write(*,*)'###UZTWC=UZTWC+PXV'
          write(*,*)'###UZTWC UZTWC PXV',UZTWC,UZTWC,PXV
          UZTWCplus=UZTWCplus+PXV
       endif          
CV.K  ADJUST UNFROZEN TENSION WATER
      UZTWH=UZTWH+PXV      
      if (prflag==1)then
         write(*,*)'###UZTWH=UZTWH+PXV'
         write(*,*)'###UZTWH PXV',UZTWH,PXV
          UZTWHplus=UZTWHplus+PXV
      endif          
      TWX=0.0
      if (prflag==1)write(*,*)'TWX=0.0',TWX
      if (prflag==1)write(*,*)'GO TO 233'
      GO TO 233
C      MOISTURE AVAILABLE IN EXCESS OF UZTWC STORAGE.
CV.K  232 UZTWC=UZTWM
  232 UZTWH=UZTWH+(UZTWM-UZTWC)
      if (prflag==1)then
         write(*,*)'###232 UZTWH=UZTWH+(UZTWM-UZTWC)'
         write(*,*)'###UZTWH UZTWH UZTWM UZTWC',UZTWH,UZTWH,UZTWM,UZTWC
         UZTWHplus=UZTWHplus+(UZTWM-UZTWC)
      endif          
      UZTWC=UZTWM
      if (prflag==1)write(*,*)'###UZTWC=UZTWM',UZTWC
  233 ADIMC=ADIMC+PXV-TWX
      if (prflag==1)then
         write(*,*)'233 ADIMC=ADIMC+PXV-TWX'
         write(*,*)'ADIMC ADIMC PXV TWX',ADIMC,ADIMC,PXV,TWX
      endif          
C
C     COMPUTE IMPERVIOUS AREA RUNOFF.
      ROIMP=PXV*SACPAR(4)
      if (prflag==1)then
         write(*,*)'ROIMP=PXV*SACPAR(4)'
         write(*,*)'ROIMP PXV SACPAR(4)',ROIMP,PXV,SACPAR(4)
      endif          
C      ROIMP IS RUNOFF FROM THE MINIMUM IMPERVIOUS AREA.
      SIMPVT=SIMPVT+ROIMP
      if (prflag==1)then
         write(*,*)'SIMPVT=SIMPVT+ROIMP'
         write(*,*)'SIMPVT SIMPVT ROIMP',SIMPVT,SIMPVT,ROIMP
      endif          
C
C     INITIALIZE TIME INTERVAL SUMS.
      SBF=0.0
      if (prflag==1)write(*,*)'SBF=0.0',SBF
      SSUR=0.0
      if (prflag==1)write(*,*)'SSUR=0.0',SSUR
      SIF=0.0
      if (prflag==1)write(*,*)'SIF=0.0',SIF
      SPERC=0.0
      if (prflag==1)write(*,*)'SPERC=0.0',SPREC
      SDRO=0.0
      if (prflag==1)write(*,*)'SDRO=0.0',SDRO
      SPBF=0.0
      if (prflag==1)write(*,*)'SPBF=0.0',SPBF
C
C     DETERMINE COMPUTATIONAL TIME INCREMENTS FOR THE BASIC TIME
C     INTERVAL
CV.K      NINC=1.0+0.2*(UZFWC+TWX)
CV.K  PERCOLATE UNFROZEN WATER ONLY
      NINC=1.0+0.2*(UZFWH+TWX)
      if (prflag==1)then
         write(*,*)'NINC=1.0+'
         delme=(0.2*(UZFWH+TWX))
         write(*,*)'NINC 1.0 0.2*(UZFWH+TWX)',NINC,1.0,delme
      endif          
C     NINC=NUMBER OF TIME INCREMENTS THAT THE TIME INTERVAL
C     IS DIVIDED INTO FOR FURTHER
C     SOIL-MOISTURE ACCOUNTING.  NO ONE INCREMENT
C     WILL EXCEED 5.0 MILLIMETERS OF UZFWC+PAV
      DINC=(1.0/NINC)*DT
      if (prflag==1)then
         write(*,*)'DINC=(1.0/NINC)*DT'
         delme=(1.0/NINC)
         write(*,*)'DINC (1.0/NINC) DT',DINC,delme,DT
      endif          
C     DINC=LENGTH OF EACH INCREMENT IN DAYS.
      PINC=TWX/NINC
      if (prflag==1)then
         write(*,*)'PINC=TWX/NINC'
         write(*,*)'PINC TWX NINC',PINC,TWX,NINC
      endif          
C     PINC=AMOUNT OF AVAILABLE MOISTURE FOR EACH INCREMENT.
C      COMPUTE FREE WATER DEPLETION FRACTIONS FOR
C     THE TIME INCREMENT BEING USED-BASIC DEPLETIONS
C      ARE FOR ONE DAY
CVK INTRODUCED REDUCTION (RUZICE & RLZICE) DUE FROZEN GROUND
CVK HOWEVER, PRIMARY RUNOFF IS UNCHANGED
CVK      DUZ=1.0-((1.0-UZK)**DINC)
CVK      DLZS=1.0-((1.0-LZSK)**DINC)
CVK  Linear transformation for frozen ground
cc      DUZ=1.0-((1.0-UZK*RUZICE)**DINC)
cc      DLZS=1.0-((1.0-LZSK*RLZICE)**DINC)
CVK  Non-linear (correct) transformation for frozen ground
      IF(IVERS .EQ. 0) THEN
       DUZ =1.0-((1.0-SACPAR(3))**DINC)
       DLZS=1.0-((1.0-SACPAR(12))**DINC)
      if (prflag==1)then
         write(*,*)'IF(IVERS .EQ. 0) THEN'
         write(*,*)'DUZ =1.0-((1.0-SACPAR(3))**DINC)'
         delme=(1.0-SACPAR(3))
         write(*,*)'DUZ ((1.0-SACPAR(3)) DINC',DUZ,delme,DINC
         write(*,*)'DLZS=1.0-((1.0-SACPAR(12))**DINC)'
         delme=(1.0-SACPAR(12))
         write(*,*)'DLZS (1.0-SACPAR(12)) DINC',DLZS,delme,DINC
      endif          
      ELSE        
       DUZ=1.0-((1.0-RUZICE)**DINC)
       DLZS=1.0-((1.0-RLZICE)**DINC)
       if (prflag==1)then
          write(*,*)'ELSE'
          write(*,*)' DUZ=1.0-((1.0-RUZICE)**DINC)'
          delme=(1.0-RUZICE)
          write(*,*)'DUZ ((1.0-RUZICE) RUZICE DINC',DUZ,delme,DINC,RUZICE
          write(*,*)'DLZS=1.0-((1.0-RLZICE)**DINC)'
          delme=(1.0-RLZICE)
          write(*,*)'DLZS (1.0-RLZICE) RLZICE DINC',DLZS,delme,DINC,RLZICE
       endif
      ENDIF 
      DLZP=1.0-((1.0-SACPAR(13))**DINC)
      if (prflag==1)then
         write(*,*)'DLZP=1.0-((1.0-SACPAR(13))**DINC)'
         delme=(1.0-SACPAR(13))
         write(*,*)'DLZP ((1.0-SACPAR(13)) SACPAR(13) DINC',DLZP,delme,
     +        SACPAR(13),DINC
      endif          
c      write(*,*)'dlzp lzpk dinc',dlzp,sacpar(13),dinc
C
C CVK  ADJUSTMENT TO DEPLETIONS DUE TO FROZEN WATER
         
C.......................................
C     START INCREMENTAL DO LOOP FOR THE TIME INTERVAL.
      if (prflag==1)write(*,*)'DO 240 I=1,NINC',NINC
      DO 240 I=1,NINC
      ADSUR=0.0
      if (prflag==1)write(*,*)'ADSUR=0.0',ADSUR
C     COMPUTE DIRECT RUNOFF (FROM ADIMP AREA).
      RATIO=(ADIMC-UZTWC)/LZTWM
      if (prflag==1)then
         write(*,*)'RATIO=(ADIMC-UZTWC)/LZTWM'
         write(*,*)'RATIO ADIMC UZTWC LZTWM',RATIO,ADIMC,UZTWC,LZTWM
      endif          
      IF (RATIO.LT.0.0) RATIO=0.0
      if (prflag==1)write(*,*)'IF (RATIO.LT.0.0) RATIO=0.0',RATIO
      ADDRO=PINC*(RATIO**2)
      if (prflag==1)then
         write(*,*)'ADDRO=PINC*(RATIO**2)'
         write(*,*)'ADDRO PINC RATIO',ADDRO,PINC,RATIO
      endif          
C     ADDRO IS THE AMOUNT OF DIRECT RUNOFF FROM THE AREA ADIMP.
C
C     COMPUTE BASEFLOW AND KEEP TRACK OF TIME INTERVAL SUM.
CV.K      BF=LZFPC*DLZP
CV.K      LZFPC=LZFPC-BF
CV.K      IF (LZFPC.GT.0.0001) GO TO 234
CV.K      BF=BF+LZFPC
CV.K      LZFPC=0.0
CV.K  BASEFLOW FROM UNFROZEN WATER ONLY   
      BF=LZFPH*DLZP
      if (prflag==1)then
         write(*,*)'BF=LZFPH*DLZP'
         write(*,*)'BF LZFPH DLZP',BF,LZFPH,DLZP
      endif          
      LZFPH=LZFPH-BF
      if (prflag==1)then
         write(*,*)'LZFPH=LZFPH-BF'
         write(*,*)'LZFPH LZFPH BF',LZFPH,LZFPH,BF
      endif          
      IF (LZFPH.GT.0.0001) THEN
       LZFPC=LZFPC-BF
      if (prflag==1)then
         write(*,*)'IF (LZFPH.GT.0.0001) THEN'
         write(*,*)'LZFPC=LZFPC-BF'
         write(*,*)'LZFPC LZFPC BF',LZFPC,LZFPC,BF
         write(*,*)'GO TO 234'
      endif          
       GO TO 234
      ENDIF
      BF=BF+LZFPH
      if (prflag==1)write(*,*)'BF=BF+LZFPH',BF,LZFPH
      LZFPH=0.0
      if (prflag==1)write(*,*)'LZFPH=0.0',LZFPH
      LZFPC=LZFPC-BF
      if (prflag==1)write(*,*)'LZFPC=LZFPC-BF',LZFPC,BF
      IF(LZFPC .LE. 0.0001) LZFPC=0.0
      if (prflag==1)write(*,*)'IF(LZFPC .LE. 0.0001) LZFPC=0.0',LZFPC
CV.K-------------------------------------
C      
  234 SBF=SBF+BF
      if (prflag==1)write(*,*)'234 SBF=SBF+BF',SBF,BF
      SPBF=SPBF+BF
      if (prflag==1)write(*,*)'SPBF=SPBF+BF',SPBF,BF
CV.K  SUPPLAMENTAL FLOW FROM UNFROZEN WATER ONLY (NOTE, DLZS
CV.K  NOTE, DLZS IS REDUCED DUE FROZEN GROUND
CV.K      BF=LZFSC*DLZS
CV.K      LZFSC=LZFSC-BF
CV.K      IF(LZFSC.GT.0.0001) GO TO 235
CV.K      BF=BF+LZFSC
CV.K      LZFSC=0.0
      BF=LZFSH*DLZS
      if (prflag==1)then
         write(*,*)'BF=LZFSH*DLZS'
         write(*,*)'BF LZFSH DLZS',BF,LZFSH,DLZS
      endif          
      LZFSH=LZFSH-BF
      if (prflag==1)write(*,*)'LZFSH=LZFSH-BF',LZFSH,BFBF
      IF(LZFSH.GT.0.0001) THEN
cc?      IF(LZFSH.GT.0.0) THEN      
       LZFSC=LZFSC-BF
       if (prflag==1)then
          write(*,*)'IF(LZFSH.GT.0.0001) THEN'
          write(*,*)'LZFSC=LZFSC-BF'
          write(*,*)'LZFSC LZFSC BF',LZFSC,LZFSC,BF
          write(*,*)'GO TO 235'
       endif          
c         if(abs(lzfsc-lzfsh) .gt. 0.000001) then
c         if(abs(lzfsc-lzfsh) .gt. 0.000001) then
c          write(*,*) ' lzfsc3: ',lzfsc,lzfsh,bf
c         endif 
       GO TO 235
      ENDIF
      BF=BF+LZFSH
      if (prflag==1)write(*,*)'BF=BF+LZFSH',BF,LZFSH
      LZFSH=0.0
      if (prflag==1)write(*,*)'LZFSH=0.0',LZFSH
      LZFSC=LZFSC-BF
      if (prflag==1)write(*,*)'LZFSC=LZFSC-BF',LZFSC,BF
      IF(LZFSC .LE. 0.0001) LZFSC=0.0   
      if (prflag==1)write(*,*)'IF(LZFSC .LE. 0.0001) LZFSC=0.0',LZFSC
CV.K--------------------------------------------
C       
  235 SBF=SBF+BF
      if (prflag==1)write(*,*)'235 SBF=SBF+BF',SBF,BF
C
C      COMPUTE PERCOLATION-IF NO WATER AVAILABLE THEN SKIP
ccvk      IF((PINC+UZFWC).GT.0.01) GO TO 251
      xx1=PINC+UZFWH
      if (prflag==1)then
         write(*,*)'xx1=PINC+UZFWH'
         write(*,*)'xx1 PINC UZFWH',xx1,PINC,UZFWH
      endif          
      if (prflag==1)write(*,*)'IF(xx1.GT.0.01) GO TO 251',xx1
      IF(xx1.GT.0.01) GO TO 251
      UZFWC=UZFWC+PINC
      if (prflag==1)write(*,*)'###UZFWC=UZFWC+PINC',UZFWC,PINC
CV.K  ADD TO UNFROZEN WATER ALSO
      UZFWH=UZFWH+PINC
      if (prflag==1)write(*,*)'###UZFWH=UZFWH+PINC',UZFWH,PINC
      if (prflag==1)write(*,*)'GO TO 249'
      GO TO 249
  251 PERCM=LZFPM*DLZP+LZFSM*DLZS
      if (prflag==1)then
         write(*,*)'251 PERCM=LZFPM*DLZP+LZFSM*DLZS'
         write(*,*)'PERCM LZFPM DLZP LZFSM DLZS',PERCM,LZFPM,DLZP
         write(*,*)'LZFSM DLZS',LZFSM,DLZS
      endif          
c      write(*,*)'percm lzfpm',percm,lzfpm
c      write(*,*)'dlzp lzfsm dlzs',dlzp,lzfsm,dlzs
CVK      PERC=PERCM*(UZFWC/UZFWM)
CV.K  USE ONLY UNFROZEN WATER RATIOS 
ccvk  new change: PERCOLATION REDUCED BY RUZPERC 
CC       PERC=PERCM*(UZFWH/UZFWM)*RUZICE
      PERC=PERCM*(UZFWH/UZFWM)
      if (prflag==1)then
         write(*,*)'PERC=PERCM*(UZFWH/UZFWM)'
         delme=(UZFWH/UZFWM)
         write(*,*)'PERC PERCM (UZFWH/UZFWM)',PERC,PERCM,delme
      endif          
c      write(*,*)'perc uzfwc uzfwm',perc,uzfwc,uzfwm
      IF(IVERS .NE. 0) PERC=PERC*RUZPERC
      if (prflag==1)then
         write(*,*)'IF(IVERS .NE. 0) PERC=PERC*RUZPERC'
         write(*,*)'PERC RUZPERC',PERC,RUZPERC
      endif          
C--      PERC=PERCM*(UZFWH/UZFWM)*RUZPERC
c      write(*,*)'perc ruzperc',perc,ruzperc
CV.K      DEFR=1.0-((LZTWC+LZFPC+LZFSC)/(LZTWM+LZFPM+LZFSM))
cvk 6/22/00      DEFR=1.0-((LZTWH+LZFPH+LZFSH)/(LZTWM+LZFPM+LZFSM))
cvk  better to keep original definition of DEFR using total water
      DEFR=1.0-((LZTWC+LZFPC+LZFSC)/(LZTWM+LZFPM+LZFSM))
      if (prflag==1)then
         write(*,*)'DEFR=1.0-((LZTWC+LZFPC+LZFSC)/(LZTWM+LZFPM+LZFSM))'
         delme=(LZTWC+LZFPC+LZFSC)
         delme2=(LZTWM+LZFPM+LZFSM)
         delme3=(LZTWC+LZFPC+LZFSC)/(LZTWM+LZFPM+LZFSM)
         write(*,*)'DEFR LZTWC LZFPC LZFSC',DEFR,LZTWC,LZFPC,LZFSC
         write(*,*)'LZTWM LZFPM LZFSM',LZTWM,LZFPM,LZFSM
         write(*,*)'ratio numerator denominator',delme3,delme,delme2
      endif          
c      write(*,*)'defr',defr      
C     DEFR IS THE LOWER ZONE MOISTURE DEFICIENCY RATIO
c--      FR=1.0
C     FR IS THE CHANGE IN PERCOLATION WITHDRAWAL DUE TO FROZEN GROUND.
c--      FI=1.0
C     FI IS THE CHANGE IN INTERFLOW WITHDRAWAL DUE TO FROZEN GROUND.
c--      IF (IFRZE.EQ.0) GO TO 239
c--       UZDEFR=1.0-((UZTWC+UZFWC)/(UZTWM+UZFWM))
CVK
CVK     CALL FGFR1(DEFR,FR,UZDEFR,FI)
CVK      IF( IVERS .EQ. 1) THEN
CVK  IF IVERS=1, OLD VERSION; IF IVERS=2, NEW VERS. FROST INDEX,
CVK  BUT OLD VERS. OF PERCOLAT. AND INTERFLOW REDUCTION
c--      IF( IVERS .LE. 2) CALL FGFR1(DEFR,FR,UZDEFR,FI)
      
c--      IF(IVERS .EQ. 3 .AND. FGPM(5) .GT. 0.) THEN
CVK  OPTIONAL VERSION TO ACCOUNT FOR ADDITIONAL IMPERVIOUS
CVK  AREAS EFFECTS DUE FROZEN GROUND
c--       FR=1-SURFRZ1(FGCO(1),FGPM(6),FGPM(5))
c--       FI=FR
c--      ENDIF 
      
c--  239 PERC=PERC*(1.0+ZPERC*(DEFR**REXP))*FR
  239 PERC=PERC*(1.0+SACPAR(7)*(DEFR**SACPAR(8)))
      if (prflag==1)then
         write(*,*)'DEFR=1.0-((LZTWC+LZFPC+LZFSC)/(LZTWM+LZFPM+LZFSM))'
         delme=(LZTWC+LZFPC+LZFSC)
         delme2=(LZTWM+LZFPM+LZFSM)
         delme3=(LZTWC+LZFPC+LZFSC)/(LZTWM+LZFPM+LZFSM)
         write(*,*)'DEFR LZTWC LZFPC LZFSC',DEFR,LZTWC,LZFPC,LZFSC
         write(*,*)'LZTWM LZFPM LZFSM',LZTWM,LZFPM,LZFSM
         write(*,*)'ratio numerator denominator',delme3,delme,delme2
      endif          
C     NOTE...PERCOLATION OCCURS FROM UZFWC BEFORE PAV IS ADDED.
CV.K      IF(PERC.LT.UZFWC) GO TO 241
      if (prflag==1)write(*,*)'IF(PERC.LT.UZFWH) GO TO 241',PERC,UZFWH
      IF(PERC.LT.UZFWH) GO TO 241
C      PERCOLATION RATE EXCEEDS UZFWH.
CV.K      PERC=UZFWC
      PERC=UZFWH
      if (prflag==1)write(*,*)'PERC=UZFWH',PERC,UZFWH
C     PERCOLATION RATE IS LESS THAN UZFWH.
  241 UZFWC=UZFWC-PERC
      if (prflag==1)write(*,*)'###241 UZFWC=UZFWC-PERC',UZFWC,PERC
CV.K  ADJUST UNFROZEN STORAGE ALSO  
      UZFWH=UZFWH-PERC    
      if (prflag==1)write(*,*)'###UZFWH=UZFWH-PERC',UZFWH,PERC

C     CHECK TO SEE IF PERCOLATION EXCEEDS LOWER ZONE DEFICIENCY.
      CHECK=LZTWC+LZFPC+LZFSC+PERC-LZTWM-LZFPM-LZFSM
      if (prflag==1)then
         write(*,*)'CHECK=LZTWC+LZFPC+LZFSC+PERC-LZTWM-LZFPM-LZFSM'
         write(*,*)'CHECK LZTWC LZFPC LZFSC',CHECK,LZTWC,LZFPC,LZFSC
         write(*,*)'PERC LZTWM LZFPM LZFSM',PERC,LZTWM,LZFPM,LZFSM
      endif          
      if (prflag==1)write(*,*)'IF(CHECK.LE.0.0) GO TO 242',CHECK
      IF(CHECK.LE.0.0) GO TO 242
      PERC=PERC-CHECK
      if (prflag==1)write(*,*)'PERC=PERC-CHECK',PERC,CHECK
      UZFWC=UZFWC+CHECK
      if (prflag==1)write(*,*)'**UZFWC=UZFWC+CHECK',UZFWC,CHECK
CV.K  ADJUST UNFROZEN STARAGE ALSO
      UZFWH=UZFWH+CHECK 
      if (prflag==1)write(*,*)'###UZFWH=UZFWH+CHECK',UZFWH,CHECK
  242 SPERC=SPERC+PERC
      if (prflag==1)write(*,*)'242 SPERC=SPERC+PERC',SPERC,PERC
C     SPERC IS THE TIME INTERVAL SUMMATION OF PERC
C
C     COMPUTE INTERFLOW AND KEEP TRACK OF TIME INTERVAL SUM.
C     NOTE...PINC HAS NOT YET BEEN ADDED
CV.K      DEL=UZFWC*DUZ*FI
CVK  INTERFLOW ALSO REDUCED DUE FROFEN GROUND (DUZ REDUCED BY RUZICE)
CVK  ADDITIONAL REDUCTION DUE IMPERVIOUS FROZEN AREAS (FI) IS OPTIONAL
CVK  IN THE NEW VERSION. BASIC OPTION IS FI=1
c--      DEL=UZFWH*DUZ*FI
      DEL=UZFWH*DUZ
      if (prflag==1)then
         write(*,*)'DEL=UZFWH*DUZ'
         write(*,*)'DEL UZFWH DUZ',DEL,UZFWH,DUZ
      endif          
      SIF=SIF+DEL
      if (prflag==1)write(*,*)'SIF=SIF+DEL',SIF,DEL
      UZFWC=UZFWC-DEL
      if (prflag==1)write(*,*)'###UZFWC=UZFWC-DEL',UZFWC,UZFWC,DEL
CV.K  ADJUST UNFROZEN STORAGE ALSO
      UZFWH=UZFWH-DEL      
      if (prflag==1)write(*,*)'###UZFWH=UZFWH-DEL',UZFWH,UZFWH,DEL
C     DISTRIBE PERCOLATED WATER INTO THE LOWER ZONES
C     TENSION WATER MUST BE FILLED FIRST EXCEPT FOR THE PFREE AREA.
C     PERCT IS PERCOLATION TO TENSION WATER AND PERCF IS PERCOLATION
C         GOING TO FREE WATER.
      PERCT=PERC*(1.0-SACPAR(14))
      if (prflag==1)then
         write(*,*)'PERCT=PERC*(1.0-SACPAR(14))'
         write(*,*)'PERCT PERC SACPAR(14)',PERCT,PERC,SACPAR(14)
      endif          
      xx1=PERCT+LZTWC
      if (prflag==1)then
         write(*,*)'xx1=PERCT+LZTWC'
         write(*,*)'xx1 PERCT LZTWC',xx1,PERCT,LZTWC
      endif          
      if (prflag==1)write(*,*)'IF (xx1.GT.LZTWM) GO TO 243',xx1,LZTWM
      IF (xx1.GT.LZTWM) GO TO 243
      LZTWC=LZTWC+PERCT
      if (prflag==1)write(*,*)'LZTWC=LZTWC+PERCT',LZTWC,PERCT
CV.K  ADJUST UNFROZEN STORAGE ALSO
      LZTWH=LZTWH+PERCT      
      if (prflag==1)write(*,*)'LZTWH=LZTWH+PERCT',LZTWH,PERCT
      PERCF=0.0
      if (prflag==1)write(*,*)'PERCF=0.0',PERCF
      if (prflag==1)write(*,*)'GO TO 244'
      GO TO 244
  243 PERCF=PERCT+LZTWC-LZTWM
      if (prflag==1)then
         write(*,*)'243 PERCF=PERCT+LZTWC-LZTWM'
         write(*,*)'PERCF PERCT LZTWC LZTWM',PERCF,PERCT,LZTWC,LZTWM
      endif          
CV.K  CHANGE UNFROZEN WATER STORAGE
      LZTWH=LZTWH+LZTWM-LZTWC  
      if (prflag==1)then
         write(*,*)'LZTWH=LZTWH+LZTWM-LZTWC'
         write(*,*)'LZTWH LZTWH LZTWM LZTWC',LZTWH,LZTWH,LZTWM,LZTWC
      endif          
      LZTWC=LZTWM
      if (prflag==1)write(*,*)'LZTWC=LZTWM',LZTWC,LZTWM
C
C      DISTRIBUTE PERCOLATION IN EXCESS OF TENSION
C      REQUIREMENTS AMONG THE FREE WATER STORAGES.
  244 PERCF=PERCF+PERC*SACPAR(14)
      if (prflag==1)then
         write(*,*)'244 PERCF=PERCF+PERC*SACPAR(14)'
         delme=PERC*SACPAR(14)
         write(*,*)'PERCF PERC SACPAR(14)',PERCF,PERC,SACPAR(14)
         write(*,*)'PERC*SACPAR(14)',delme
      endif          
      if (prflag==1)write(*,*)'IF(PERCF.EQ.0.0) GO TO 245',PERCF
      IF(PERCF.EQ.0.0) GO TO 245
      HPL=LZFPM/(LZFPM+LZFSM)
      if (prflag==1)then
         write(*,*)'HPL=LZFPM/(LZFPM+LZFSM)'
         delme=(LZFPM+LZFSM)
         write(*,*)'HPL LZFPM LZFSM',HPL,LZFPM,LZFPM,LZFSM
         write(*,*)'((LZFPM+LZFSM)',delme
      endif          
C     HPL IS THE RELATIVE SIZE OF THE PRIMARY STORAGE
C     AS COMPARED WITH TOTAL LOWER ZONE FREE WATER STORAGE.

c VK changed to account for ZERO MAX storage
      if(LZFPM .ne. 0.) then
       RATLP=LZFPC/LZFPM
      if (prflag==1)then
         write(*,*)'if(LZFPM .ne. 0.) then'
         write(*,*)'RATLP=LZFPC/LZFPM'
         write(*,*)'RATLP LZFPC LZFPM',RATLP,LZFPC,LZFPM
      endif          
      else
       RATLP = 1.
       if (prflag==1)write(*,*)'else RATLP=1.',RATLP
      endif
      if(LZFSM .ne. 0.) then
       RATLS=LZFSC/LZFSM
      if (prflag==1)then
         write(*,*)'if(LZFSM .ne. 0.) then'
         write(*,*)'RATLS=LZFSC/LZFSM'
         write(*,*)'RATLS LZFSC LZFSM',RATLS,LZFSC,LZFSM
      endif          
      else
       RATLS = 1.
       if (prflag==1)write(*,*)'else RATLS=1.',RATLS
      endif
        
C     RATLP AND RATLS ARE CONTENT TO CAPACITY RATIOS, OR
C     IN OTHER WORDS, THE RELATIVE FULLNESS OF EACH STORAGE
      FRACP=(HPL*2.0*(1.0-RATLP))/((1.0-RATLP)+(1.0-RATLS))
      if (prflag==1)then
         write(*,*)'FRACP=(HPL*2.0*(1.0-RATLP))/'
         write(*,*)'((1.0-RATLP)+(1.0-RATLS))'
         write(*,*)'FRACP HPL RATLP',FRACP,HPL,RATLP
         write(*,*)'RATLP RATLS',RATLP,RATLS

      endif          
C     FRACP IS THE FRACTION GOING TO PRIMARY.
      IF (FRACP.GT.1.0) FRACP=1.0
       if (prflag==1)write(*,*)'IF (FRACP.GT.1.0) FRACP=1.0',FRACP
      PERCP=PERCF*FRACP
       if (prflag==1)write(*,*)'PERCP=PERCF*FRACP',PERCP,FRACP
      PERCS=PERCF-PERCP
       if (prflag==1)write(*,*)'PERCS=PERCF-PERCP',PERCS,PERCP
C     PERCP AND PERCS ARE THE AMOUNT OF THE EXCESS
C     PERCOLATION GOING TO PRIMARY AND SUPPLEMENTAL
C      STORGES,RESPECTIVELY.
      LZFSC=LZFSC+PERCS
       if (prflag==1)write(*,*)'LZFSC LZFSC PERCS',LZFSC,PERCS
CV.K      IF(LZFSC.LE.LZFSM) GO TO 246
      IF(LZFSC.LE.LZFSM) THEN
       LZFSH=LZFSH+PERCS
      if (prflag==1)then
         write(*,*)'LZFSH=LZFSH+PERCS',LZFSH,PERCS
         write(*,*)'GO TO 246'
      endif          
       GO TO 246
      ENDIF
       
      PERCS=PERCS-LZFSC+LZFSM
      if (prflag==1)then
         write(*,*)'PERCS=PERCS-LZFSC+LZFSM'
         write(*,*)'PERCS PERCS LZFSC LZFSM',PERCS,PERCS,LZFSC,LZFSM
      endif          
CV.K  ADJUST UNFROZEN STORAGE ALSO
      LZFSH=LZFSH+PERCS
      if (prflag==1)write(*,*)'LZFSH=LZFSH+PERCS',LZFSH,PERCS
      LZFSC=LZFSM
       if (prflag==1)write(*,*)'LZFSC=LZFSM',LZFSC,LZFSM
  246 LZFPC=LZFPC+(PERCF-PERCS)
      if (prflag==1)then
         write(*,*)'246 LZFPC=LZFPC+(PERCF-PERCS)'
         write(*,*)'LZFPC LZFPC PERCF PERCS',LZFPC,LZFPC,PERCF,PERCS
      endif          
C     CHECK TO MAKE SURE LZFPC DOES NOT EXCEED LZFPM.
CV.K      IF (LZFPC.LE.LZFPM) GO TO 245
      IF (LZFPC.LE.LZFPM) THEN
       LZFPH=LZFPH+(PERCF-PERCS)
      if (prflag==1)then
         write(*,*)'IF (LZFPC.LE.LZFPM) THEN'
         write(*,*)'LZFPH=LZFPH+(PERCF-PERCS)'
         write(*,*)'LZFPH PERCF PERCS',LZFPH,PERCF,PERCS
         write(*,*)'GO TO 245'
      endif          
       GO TO 245
      ENDIF
       
      EXCESS=LZFPC-LZFPM
      if (prflag==1)then
         write(*,*)'EXCESS=LZFPC-LZFPM'
         write(*,*)'EXCESS LZFPC LZFPM',EXCESS,LZFPC,LZFPM
      endif          
      LZTWC=LZTWC+EXCESS
      if (prflag==1)write(*,*)'LZTWC=LZTWC+EXCESS',LZTWC,EXCESS
CV.K  ADJUST UNFROZEN STORAGES ALSO
      LZTWH=LZTWH+EXCESS
      if (prflag==1)write(*,*)'LZTWH=LZTWH+EXCESS',LZTWH,EXCESS
      LZFPH=LZFPH+(PERCF-PERCS)-EXCESS
      if (prflag==1)then
         write(*,*)'LZFPH=LZFPH+(PERCF-PERCS)-EXCESS'
         write(*,*)'LZFPH PERCF PERCS EXCESS',LZFPH,PERCF,PERCS,EXCESS
      endif          
      LZFPC=LZFPM
      if (prflag==1)write(*,*)'LZFPC LZFPM',LZFPC,LZFPM
C
C     DISTRIBUTE PINC BETWEEN UZFWC AND SURFACE RUNOFF.
      if (prflag==1)write(*,*)'245 IF(PINC.EQ.0.0) GO TO 249',PINC
  245 IF(PINC.EQ.0.0) GO TO 249
C     CHECK IF PINC EXCEEDS UZFWM
      xx1=PINC+UZFWC
      if (prflag==1)then
         write(*,*)'xx1=PINC+UZFWC'
         write(*,*)'xx1 PINC UZFWC',xx1,PINC,UZFWC
      endif          
      if (prflag==1)write(*,*)'IF(xx1.GT.UZFWM) GO TO 248',xx1,UZFWM
      IF(xx1.GT.UZFWM) GO TO 248
C     NO SURFACE RUNOFF
      UZFWC=UZFWC+PINC
      if (prflag==1)write(*,*)'###UZFWC=UZFWC+PINC',UZFWC,PINC
CV.K  ADJUST UNFROZEN STORAGE ALSO
      UZFWH=UZFWH+PINC
      if (prflag==1)write(*,*)'###UZFWH=UZFWH+PINC',UZFWH,PINC
      if (prflag==1)write(*,*)'GO TO 249'
      GO TO 249
C
C     COMPUTE SURFACE RUNOFF (SUR) AND KEEP TRACK OF TIME INTERVAL SUM.
  248 SUR=PINC+UZFWC-UZFWM
      if (prflag==1)then
         write(*,*)'248 SUR=PINC+UZFWC-UZFWM'
         write(*,*)'SUR PINC UZFWC UZFWM',SUR,PINC,UZFWC,UZFWM
      endif          
      UZFWC=UZFWM
      if (prflag==1)write(*,*)'###UZFWC=UZFWM',UZFWC,UZFWM
CV.K  ADJUST UNFROZEN STORAGE ALSO
      UZFWH=UZFWH+PINC-SUR
      if (prflag==1)then
         write(*,*)'###UZFWH=UZFWH+PINC-SUR'
         write(*,*)'###UZFWH PINC SUR',UZFWH,PINC,SUR
      endif          
      SSUR=SSUR+SUR*PAREA
      if (prflag==1)then
         write(*,*)'SSUR=SSUR+SUR*PAREA'
         write(*,*)'SSUR SUR PAREA',SSUR,SUR,PAREA
      endif          
      ADSUR=SUR*(1.0-ADDRO/PINC)
      if (prflag==1)then
         write(*,*)'ADSUR=SUR*(1.0-ADDRO/PINC)'
         write(*,*)'ADSUR SUR ADDRO PINC)',ADSUR,SUR,ADDRO,PINC
      endif          
C     ADSUR IS THE AMOUNT OF SURFACE RUNOFF WHICH COMES
C     FROM THAT PORTION OF ADIMP WHICH IS NOT
C     CURRENTLY GENERATING DIRECT RUNOFF.  ADDRO/PINC
C     IS THE FRACTION OF ADIMP CURRENTLY GENERATING
C     DIRECT RUNOFF.
      SSUR=SSUR+ADSUR*ADIMP
      if (prflag==1)then
         write(*,*)'SSUR=SSUR+ADSUR*ADIMP'
         write(*,*)'SSUR ADSUR ADIMP',SSUR,ADSUR,ADIMP
      endif          
C
C     ADIMP AREA WATER BALANCE -- SDRO IS THE 6 HR SUM OF
C          DIRECT RUNOFF.
  249 ADIMC=ADIMC+PINC-ADDRO-ADSUR  
      if (prflag==1)then
         write(*,*)'249 ADIMC=ADIMC+PINC-ADDRO-ADSUR  '
         write(*,*)'ADIMC ADIMC PINC ADDRO ADSUR',ADIMC,PINC,ADDRO,ADSUR
      endif          
      xx1=UZTWM+LZTWM
      if (prflag==1)then
         write(*,*)'xx1=UZTWM+LZTWM'
         write(*,*)'xx1 UZTWM LZTWM',xx1,UZTWM,LZTWM
      endif          
      if (prflag==1)write(*,*)'IF (ADIMC.LE.xx1) GO TO 247',ADIMC,xx1
      IF (ADIMC.LE.xx1) GO TO 247
      ADDRO=ADDRO+ADIMC-xx1
      if (prflag==1)then
         write(*,*)'ADDRO=ADDRO+ADIMC-xx1'
         write(*,*)'ADDRO ADIMC xx1',ADDRO,ADIMC,xx1
      endif          
      ADIMC=xx1
      if (prflag==1)write(*,*)'ADIMC=xx1',ADIMC,xx1
  247 SDRO=SDRO+ADDRO*ADIMP
      if (prflag==1)then
         write(*,*)'247 SDRO=SDRO+ADDRO*ADIMP'
         write(*,*)'SDRO ADDRO ADIMP',SDRO,ADDRO,ADIMP
      endif          
      IF (ADIMC.LT.0.00001) ADIMC=0.0
      if (prflag==1)write(*,*)'IF (ADIMC.LT.0.00001) ADIMC=0.0',ADIMC
  240 CONTINUE
      if (prflag==1)write(*,*)'240 CONTINUE'

C.......................................
C     END OF INCREMENTAL DO LOOP.
C.......................................

C     COMPUTE SUMS AND ADJUST RUNOFF AMOUNTS BY THE AREA OVER
C     WHICH THEY ARE GENERATED.
      EUSED=E1+E2+E3
      if (prflag==1)then
         write(*,*)'EUSED=E1+E2+E3'
         write(*,*)'EUSED E1 E2 E3',EUSED,E1,E2,E3
      endif          
C     EUSED IS THE ET FROM PAREA WHICH IS 1.0-ADIMP-PCTIM
      SIF=SIF*PAREA
      if (prflag==1)write(*,*)'SIF=SIF*PAREA',SIF,PAREA
C
C     SEPARATE CHANNEL COMPONENT OF BASEFLOW
C     FROM THE NON-CHANNEL COMPONENT
      TBF=SBF*PAREA
      if (prflag==1)write(*,*)'TBF=SBF*PAREA',TBF,SBF,PAREA
C     TBF IS TOTAL BASEFLOW
      BFCC=TBF*(1.0/(1.0+SACPAR(15)))
      if (prflag==1)then
         write(*,*)'BFCC=TBF*(1.0/(1.0+SACPAR(15)))'
         delme=(1.0/(1.0+SACPAR(15)))
         write(*,*)'BFCC ratio TBF SACPAR(15)',BFCC,delme,TBF,SACPAR(15)
      endif          
C     BFCC IS BASEFLOW, CHANNEL COMPONENT
      BFP=SPBF*PAREA/(1.0+SACPAR(15))
      if (prflag==1)then
         write(*,*)'BFP=SPBF*PAREA/(1.0+SACPAR(15))'
         delme=SPBF*PAREA
         write(*,*)'BFP product SPBF',BFP,delme,SPBF
         write(*,*)'PAREA SACPAR(15)',PAREA,SACPAR(15)
      endif          
      BFS=BFCC-BFP
      if (prflag==1)write(*,*)'BFS=BFCC-BFP',BFS,BFCC,BFP
      IF(BFS.LT.0.0)BFS=0.0
      if (prflag==1)write(*,*)'IF(BFS.LT.0.0)BFS=0.0',BFS
      BFNCC=TBF-BFCC
      if (prflag==1)write(*,*)'BFNCC=Baseflow, non-channel component'
      if (prflag==1)write(*,*)'BFNCC=TBF-BFCC',BFNCC,TBF,BFCC
C     BFNCC IS BASEFLOW,NON-CHANNEL COMPONENT
C
C     ADD TO MONTHLY SUMS.
c--      SINTFT=SINTFT+SIF
c--      SGWFP=SGWFP+BFP
c--      SGWFS=SGWFS+BFS
c--      SRECHT=SRECHT+BFNCC
c--      SROST=SROST+SSUR
c--      SRODT=SRODT+SDRO
C
C     COMPUTE TOTAL CHANNEL INFLOW FOR THE TIME INTERVAL.
      TCI=ROIMP+SDRO+SSUR+SIF+BFCC
      if (prflag==1)then
         write(*,*)'TCI=TOTAL CHANNEL INFLOW'
         write(*,*)'TCI=ROIMP+SDRO+SSUR+SIF+BFCC'
         delme=SPBF*PAREA
         write(*,*)'TCI ROIMP SDRO',TCI,ROIMP,SDRO
         write(*,*)'SSUR SIF BFCC',SSUR,SIF,BFCC
      endif          
        GRND = SIF + BFCC   ! interflow is part of ground flow
        if (prflag==1)write(*,*)'GRND = SIF + BFCC'
        if (prflag==1)write(*,*)'GRND SIF BFCC',GRND,SIF,BFCC
CC	GRND = BFCC         ! interflow is part of surface flow
	SURF = TCI - GRND
        if (prflag==1)write(*,*)'SURF = TCI - GRND'
        if (prflag==1)write(*,*)'SURF TCI GRND',SURF,TCI,GRND
C
C     COMPUTE E4-ET FROM RIPARIAN VEGETATION.
	E4=(EDMND-EUSED)*SACPAR(6)
      if (prflag==1)then
         write(*,*)'E4=(EDMND-EUSED)*SACPAR(6)'
         write(*,*)'E4 EDMND EUSED SACPAR(6)',E4,EDMND,EUSED,SACPAR(6)
      endif          
C
C     SUBTRACT E4 FROM CHANNEL INFLOW
	TCI=TCI-E4
        if (prflag==1)write(*,*)'TCI=TCI-E4',TCI,E4
        if (prflag==1)write(*,*)'IF(TCI.GE.0.0) GO TO 250',TCI
	IF(TCI.GE.0.0) GO TO 250
	E4=E4+TCI
        if (prflag==1)write(*,*)'E4=E4+TCI',E4,TCI
	TCI=0.0
        if (prflag==1)write(*,*)'TCI=0.0',TCI
cc  250 SROT=SROT+TCI
250	CONTINUE
        if (prflag==1)write(*,*)'250	CONTINUE'
	GRND = GRND - E4
        if (prflag==1)write(*,*)'GRND=GRND-E4',GRND,E4
	IF (GRND .LT. 0.) THEN
	   SURF = SURF + GRND
	   GRND = 0.
           if (prflag==1)then
              write(*,*)'IF (GRND .LT. 0.) THEN'
              write(*,*)'SURF = SURF + GRND',SURF,GRND
              write(*,*)'GRND = 0.',GRND
           endif          
	 IF (SURF .LT. 0.) SURF = 0.
        if (prflag==1)write(*,*)'IF (SURF .LT. 0.) SURF = 0.',SURF
	END IF
C
C     COMPUTE TOTAL EVAPOTRANSPIRATION-TET
      EUSED=EUSED*PAREA
      if (prflag==1)write(*,*)'EUSED=EUSED*PAREA',EUSED,PAREA
      TET=EUSED+E5+E4
      if (prflag==1)then
         write(*,*)'TET=EUSED+E5+E4'
         write(*,*)'TET EUSED E5 E4',TET,EUSED,E5,E4
      endif          
c--      SETT=SETT+TET
c--      SE1=SE1+E1*PAREA
c--      SE3=SE3+E3*PAREA
c--      SE4=SE4+E4
c--      SE5=SE5+E5
C     CHECK THAT ADIMC.GE.UZTWC
      IF (ADIMC.LT.UZTWC) ADIMC=UZTWC
      if(prflag==1)write(*,*)'IF(ADIMC.LT.UZTWC)ADIMC=UZTWC',ADIMC,UZTWC
C
c  Return back SAC states
      SACST(1)=UZTWC
      SACST(2)=UZFWC      
      SACST(3)=LZTWC
      SACST(4)=LZFSC
      SACST(5)=LZFPC
      SACST(6)=ADIMC

c new change: check negative states
      do i=1,6
       if(sacst(i) .lt. -1.0) then
        write(*,*) ' SAC state#',i,'<-1.',sacst(i)
        stop
       endif
       if(sacst(i) .lt. 0.0) sacst(i)=0.0
      enddo
      if(uztwh .lt. 0.0) uztwh=0.0
      if(uzfwh .lt. 0.0) uzfwh=0.0
      if(lztwh .lt. 0.0) lztwh=0.0
      if(lzfsh .lt. 0.0) lzfsh=0.0
      if(lzfph .lt. 0.0) lzfph=0.0
      if(sacst(1) .lt. uztwh) uztwh=sacst(1)
      if(sacst(2) .lt. uzfwh) uzfwh=sacst(2)
      if(sacst(3) .lt. lztwh) lztwh=sacst(3)
      if(sacst(4) .lt. lzfsh) lzfsh=sacst(4)
      if(sacst(5) .lt. lzfph) lzfph=sacst(5)
c new change        
       
c      WRITE(*,*) 'FLAND2: end, UZFWC = ', UZFWC
CVK  NEW VERSION OF FROST INDEX  ------------------------------
       IF (IVERS .NE. 0) THEN
        IF(FRZST(6) .LT. 0.) THEN
         FRZST(6)=FRZST(6)+UZTWH
        ELSE
         FRZST(6)=UZTWH
        ENDIF
        IF(FRZST(7) .LT. 0.) THEN
         FRZST(7)=FRZST(7)+UZFWH
        ELSE
         FRZST(7)=UZFWH
        ENDIF
        IF(FRZST(8) .LT. 0.) THEN
         FRZST(8)=FRZST(8)+LZTWH
        ELSE
         FRZST(8)=LZTWH
        ENDIF
        IF(FRZST(9) .LT. 0.) THEN
         FRZST(9)=FRZST(9)+LZFSH
        ELSE
         FRZST(9)=LZFSH
        ENDIF
        IF(FRZST(10) .LT. 0.) THEN
         FRZST(10)=FRZST(10)+LZFPH
        ELSE
         FRZST(10)=LZFPH
        ENDIF
CBL Need to look at the above "IF" statement and may wish to include
CBL above code, since frozen ground is being considered, just not here
      END IF

CBL2014 Check water balance
      if(prflag==1) then
cbl total water
         sacupper0=( SACST_PRV(1) + SACST_PRV(2) )
         sacupper=( SACST(1) + SACST(2) )
         write(*,*)'Cont:sacupper0 sacupper',sacupper0,sacupper
         saclower0=( SACST_PRV(3) + SACST_PRV(4) + SACST_PRV(5) )
         saclower=( SACST(3) + SACST(4) + SACST(5) )
         write(*,*)'Cont:saclower0 saclower',saclower0,saclower
         sacupperdif=sacupper-sacupper0
         write(*,*)'sacupper-sacupper0',sacupperdif
         saclowerdif=saclower-saclower0
         write(*,*)'saclower-saclower0',saclowerdif
         adimcdif=( SACST(6) - SACST_PRV(6) )
         write(*,*)'ADIMCDIF',adimcdif,SACST(6),SACST_PRV(6)
         do i=1,6
         delme=SACST(i)-SACST_PRV(i)
         write(*,*)i,'SACST(i)-SACST_PRV(i)',delme,SACST(i),SACST_PRV(i)
         enddo
         sactot=sacupper+saclower
         sactot0=sacupper0+saclower0
         sactotdif=sactot-sactot0
         write(*,*)'sactotdif sactot sactot0',sactotdif,sactot,sactot0
cbl liquid water
cbl         sacupper0=( FRZST_PRV(6) + FRZST_PRV(7) )
cbl         sacupper=( FRZST(6) + FRZST(7) )
cbl         write(*,*)'Cont:sacupper0 sacupper',sacupper0,sacupper
cbl         saclower0=( FRZST_PRV(8) + FRZST_PRV(9) + FRZST_PRV(10) )
cbl         saclower=( FRZST(8) + FRZST(9) + FRZST(10) )
cbl         write(*,*)'Cont:saclower0 saclower',saclower0,saclower
cbl         frztot=sacupper+saclower
cbl         frztot0=sacupper0+saclower0
cbl         frztotdif=frztot-frztot0
cbl         write(*,*)'frztotdif frztot frztot0',frztotdif,frztot,frztot0
cbl DS should be negative, since all other fluxes are 'leaving'
cbl the soil column         
         DS = sactotdif
         BAL = PXV-TET-SURF-GRND-DS
         write(*,*)'WITHIN FLAND wb_error BAL',BAL
         write(*,*)'PXV-TET-SURF-GRND-DS'
         write(*,"(5(xf13.8))")PXV,TET,SURF,GRND,DS

         write(*,*)'LOWER ZONE WATER BALANCE'
         delinput=PXV-sacupperdif-TET-SURF
         write(*,*)'delinput PXV sacupperdif TET SURF'
         write(*,*)delinput,PXV,sacupperdif,TET,SURF
         output=GRND
         delstorage=saclowerdif
         bal=delinput-output-delstorage
         write(*,*)'balance',bal
         write(*,*)'delinput output ds',delinput,output,delstorage

      endif
      if(prflag==1) then
         write(*,*)'UZTW',UZTWH,UZTWC,UZTWM
         write(*,*)'UZFW',UZFWH,UZFWC,UZFWM
         write(*,*)'LZTW',LZTWH,LZTWC,LZTWM
         write(*,*)'LZFP',LZFPH,LZFPC,LZFPM
         write(*,*)'LZFS',LZFSH,LZFSC,LZFSM
      endif

CBL        CALL FROST2_1(PXV,TA,WE,AESC,SH,FRZPAR,SACPAR,FRZST,SACST,
CBL     +        SACST_PRV,SMC,SH2O,DTFRZ,IDTFRZ,NSOIL,NUPL,NSAC,IVERS,
CBL     +        FRZDUP,FRZDBT,FROST, SMAX)
CBL
CBLCVK_02  NEW OPTION TO INTERPOLATE MODEL SOIL LAYER TEMP. INTO DESIRED LAYERS
CBL      NMOD=NSOIL+1
CBL      DSMOD(1)=0.
CBL      DSMOD(1)=-0.5*FRZPAR(10)
CBL      TSTMP(1)=FRZCBLST(1)
CBL      SWTMP(1CBLCBL)=SMC(2)
CBL      SWHTMP(1)=SH2O(2)
CBL      DSMOD(NMOD)=-FRZPAR(CBL5)
CBL      TSTMP(NMOD)=FRZPAR(2)-T0
CBL      SWTMP(NMOD)=SMAX
CBL      SWHTMP(NMOD)=SMAX
CBL      do i=2,nmod-1
CBL       DSMOD(I)=-0.5*(FRZPAR(I+8)+FRZPAR(I+9))
CBL       TSTMP(I)=FRZST(I)
CBL       SWTMP(I)=SMC(I)
CBL       SWHTMP(I)=SH2O(I)
CBL      ENDDO 
CBLcc-      do ii = 1, nmod
CBLcc-        WRITE(*,*) SWTMP(ii), SWHTMP(ii), DSMOD(ii)
CBLcc-      ENDDO
CBL
CBL      CALL SOIL_INT1(TSTMP,NMOD,DSMOD,DSINT,NDSINT,TSINT)
CBL      CALL SOIL_INT1(SWTMP,NMOD,DSMOD,DSINTW,NDINTW,SWINT)
CBL      CALL SOIL_INT1(SWHTMP,NMOD,DSMOD,DSINTW,NDINTW,SWHINT)
CBL
CBLcvk  1/2008 Option to generate normalized soil moisture content (SR)
CBL      if(NORMALIZE .eq. 1) then
CBL       DO I=1,NDINTW
CBL        SWINT(I) =  (SWINT(I) -  FRZPAR(9))/(SMAX-FRZPAR(9))
CBL        SWHINT(I) = (SWHINT(I) - FRZPAR(9))/(SMAX-FRZPAR(9CBL))
CBL       ENDDO
CBL     endif 
CBLcvk  1/2008  end soil moisture normalization
CBL       
CBLC--      DO I=1,NINTW
CBLC--       IF(I .EQ. 1) THEN
CBLC--        SWINT(I)=SWINT(I)*DSINTW(I)*1000.
CBLC--        SWHINT(I)=SWHINT(I)*DSINTW(I)*1000
CBLC--       ELSE	
CBLC--        SWINT(I)=SWINT(I)*(DSINTW(I)-DSINTW(I-1))*1000.
CBLC--        SWHINT(I)=SWHINT(I)*(DSINTW(I)-DSINTW(I-1))*1000
CBLC--       ENDIF
CBLC--      ENDDO 	
CBL
CBLc        WRITE (*,905) (sacst(ii),ii=1,6),
CBLc     1  SPERC,ROIMP,SDRO,SSUR,SIF,BFS,BFP,TCI,
CBLc     2  EDMND,TET,PXV,(frzst(II),II=6,10),WE,SH,AESC,TA,frost,
CBLc     3  (frzst(II),II=1,nsoil),frzdup,frzdbt             
CBLc        write(*,977) smax,(swint(ii),swhint(ii),tstmp(ii),ii=1,nintw)
CBL       ELSE
CBLcc        WRITE (*,977) (sacst(ii),ii=1,6),
CBLcc     1  SPERC,ROIMP,SDRO,SSUR,SIF,BFS,BFP,TCI,
CBLcc     2  EDMND,TET,PXV,UZTWH,UZFWH,LZTWH,LZFSH,LZFPH,WE,SH,AESC,TA
CBL      ENDIF
CBL  905 FORMAT (1H ,14x,F7.2,F7.3,F7.2,F7.3,2F7.2,7F7.3,2F8.3,F7.3,
CBL     +        F9.4,5f7.2,2f6.1,F6.3,F9.4,f7.2,8f7.2,20F7.1)
CBL  977 FORMAT (1H ,14x,f7.2,6(3F7.2))
CBLc             ,F7.3,F7.2,F7.3,2F7.2,7F7.3,2F8.3,F7.3,
CBLc     +        F9.4,5F7.2,2f6.1,F6.3,F9.4,8f7.2,20F7.1)
CBL
C.......................................

      RETURN
      END
