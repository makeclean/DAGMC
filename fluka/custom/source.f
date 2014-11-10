*$ CREATE SOURCE.FOR
*COPY SOURCE
*
*=== source ===========================================================*
*
      SUBROUTINE SOURCE ( NOMORE )

      INCLUDE '(DBLPRC)'
      INCLUDE '(DIMPAR)'
      INCLUDE '(IOUNIT)'
*
*----------------------------------------------------------------------*
*                                                                      *
*     Copyright (C) 1990-2010      by    Alfredo Ferrari & Paola Sala  *
*     All Rights Reserved.                                             *
*                                                                      *
*                                                                      *
*     New source for FLUKA9x-FLUKA20xy:                                *
*                                                                      *
*     Created on 07 January 1990   by    Alfredo Ferrari & Paola Sala  *
*                                                   Infn - Milan       *
*                                                                      *
*     Last change on  17-Oct-10    by    Alfredo Ferrari               *
*                                                                      *
*  This is just an example of a possible user written source routine.  *
*  note that the beam card still has some meaning - in the scoring the *
*  maximum momentum used in deciding the binning is taken from the     *
*  beam momentum.  Other beam card parameters are obsolete.            *
*                                                                      *
*       Output variables:                                              *
*                                                                      *
*              Nomore = if > 0 the run will be terminated              *
*                                                                      *
*----------------------------------------------------------------------*
*
      INCLUDE '(BEAMCM)'
      INCLUDE '(FHEAVY)'
      INCLUDE '(FLKSTK)'
      INCLUDE '(IOIOCM)'
      INCLUDE '(LTCLCM)'
      INCLUDE '(PAPROP)'
      INCLUDE '(SOURCM)'
      INCLUDE '(SUMCOU)'
*
      LOGICAL LFIRST
*
      SAVE LFIRST
      DATA LFIRST / .TRUE. /
      DOUBLE PRECISION RANDS(6)
      
* MIMICS THE MCNP IDUM/RDUM
      DOUBLE PRECISION RDUM(50)
      INTEGER IDUM(50)
      INTEGER BIN_IDX
      DOUBLE PRECISION BIN_WIDTH
      DOUBLE PRECISION ERG,XXX,YYY
      DOUBLE PRECISION R,Z

      SAVE RDUM
      SAVE IDUM
      SAVE BIN_WIDTH

*======================================================================*
*                                                                      *
*                 BASIC VERSION                                        *
*                                                                      *
*======================================================================*
      NOMORE = 0
*  +-------------------------------------------------------------------*
*  |  First call initializations:
      IF ( LFIRST ) THEN
*  |  *** The following 3 cards are mandatory ***
         TKESUM = ZERZER
         LFIRST = .FALSE.
         LUSSRC = .TRUE.
*     |  *** User initialization ***
*         CALL SETUP()

         RDUM(1)=1.D20
         RDUM(2)=3.E19
         RDUM(3)=1.0D20
         RDUM(4)=6.09D0
         RDUM(5)=0.1D0
         RDUM(6)=45.9D0
         RDUM(7)=0.7D0
         RDUM(8)=1.0E0
         RDUM(9)=8.06D0
         RDUM(10)=0.94D0
         RDUM(11)=1.682D0
         RDUM(12)=2.85D0
         RDUM(13)=0.53D0
         RDUM(14)=0.35D0
         
         IDUM(1)=0
         IDUM(2)=40

         CALL SETUP_PARAMETRIC(RDUM(1),RDUM(2),RDUM(3),RDUM(4), 
     &                         RDUM(5),RDUM(6),RDUM(7),RDUM(8),                                
     &                         RDUM(9),RDUM(10),RDUM(11),RDUM(12),                                
     &                         RDUM(13),RDUM(14),' ',IDUM(1),IDUM(2))                                
*        calculate the width of the radial bin
         BIN_WIDTH = RDUM(10)/REAL(IDUM(2))

         
      END IF
*  |
*  +-------------------------------------------------------------------*
*  Push one source particle to the stack. Note that you could as well
*  push many but this way we reserve a maximum amount of space in the
*  stack for the secondaries to be generated
*  Npflka is the stack counter: of course any time source is called it
*  must be =0
      NPFLKA = NPFLKA + 1

*  Particle type (1=proton.....). Ijbeam is the type set by the BEAM
*  card
*  +-------------------------------------------------------------------*
*  |  (Radioactive) isotope:
      IF ( IJBEAM .EQ. -2 .AND. LRDBEA ) THEN
         IARES  = IPROA
         IZRES  = IPROZ
         IISRES = IPROM
         CALL STISBM ( IARES, IZRES, IISRES )
         IJHION = IPROZ  * 1000 + IPROA
         IJHION = IJHION * 100 + KXHEAV
         IONID  = IJHION
         CALL DCDION ( IONID )
         CALL SETION ( IONID )
*  |
*  +-------------------------------------------------------------------*
*  |  Heavy ion:
      ELSE IF ( IJBEAM .EQ. -2 ) THEN
         IJHION = IPROZ  * 1000 + IPROA
         IJHION = IJHION * 100 + KXHEAV
         IONID  = IJHION
         CALL DCDION ( IONID )
         CALL SETION ( IONID )
         ILOFLK (NPFLKA) = IJHION
*  |  Flag this is prompt radiation
         LRADDC (NPFLKA) = .FALSE.
*  |  Group number for "low" energy neutrons, set to 0 anyway
         IGROUP (NPFLKA) = 0
*  |
*  +-------------------------------------------------------------------*
*  |  Normal hadron:
      ELSE
         IONID = IJBEAM
         ILOFLK (NPFLKA) = IJBEAM
*  |  Flag this is prompt radiation
         LRADDC (NPFLKA) = .FALSE.
*  |  Group number for "low" energy neutrons, set to 0 anyway
         IGROUP (NPFLKA) = 0
      END IF
*  |
*  +-------------------------------------------------------------------*
*  From this point .....
*  Particle generation (1 for primaries)
      LOFLK  (NPFLKA) = 1
*  User dependent flag:
      LOUSE  (NPFLKA) = 0
*  No channeling:
      LCHFLK (NPFLKA) = .FALSE.
      DCHFLK (NPFLKA) = ZERZER
*  User dependent spare variables:
      DO 100 ISPR = 1, MKBMX1
         SPAREK (ISPR,NPFLKA) = ZERZER
 100  CONTINUE
*  User dependent spare flags:
      DO 200 ISPR = 1, MKBMX2
         ISPARK (ISPR,NPFLKA) = 0
 200  CONTINUE

* set the rands for this particle
      RANDS(1)=FLRNDM(XDUMMY)
      RANDS(2)=FLRNDM(XDUMMY)
      RANDS(3)=FLRNDM(XDUMMY)
      RANDS(4)=FLRNDM(XDUMMY)
      RANDS(5)=FLRNDM(XDUMMY)
      RANDS(6)=FLRNDM(XDUMMY)
* particle direction cosines
      CALL RACO(TXFLK(NPFLKA),TYFLK(NPFLKA),TZFLK(NPFLKA))
      WRITE(*,*) TXFLK(NPFLKA),TYFLK(NPFLKA),TZFLK(NPFLKA)
*  sample
*      CALL SAMPLE(RANDS,XFLK(NPFLKA),YFLK(NPFLKA),ZFLK(NPFLKA),
*     &            TKEFLK(NPFLKA),WTFLK(NPFLKA))
      CALL SAMPLE_SOURCE_PARAMETRIC(IDUM(2),BIN_WIDTH,RAD_SMP,
     &                             RAD_IDX,RANDS(1),RANDS(2))
*     Convert radial bin to RZ coords      
      CALL CONVERT_RAD_TO_RZ(RDUM(11),RDUM(10),RDUM(12),RDUM(13)
     &                      ,RDUM(14),RAD_SMP,RANDS(3),R,Z)
*     Convert rz to xyz
      CALL CONVERT_RZ_TO_XYZ(R,RANDS(4),XXX,YYY,RDUM(15),
     &                       RDUM(16),360.0D0)

      XFLK(NPFLKA) = XXX*100.0
      YFLK(NPFLKA) = YYY*100.0
      ZFLK(NPFLKA) = Z*100.0

*     Sample Energy
      CALL SAMPLE_ENERGY_PARAMETRIC(RAD_IDX,RANDS(5),RANDS(6),ERG)

*     Convert Kinteic energy to GeV
      TKEFLK(NPFLKA) = ERG/1000.0 

      WRITE(*,*) XFLK(NPFLKA),YFLK(NPFLKA),ZFLK(NPFLKA),
     &            TKEFLK(NPFLKA),WTFLK(NPFLKA)

*  Wt is the weight of the particle
      WEIPRI = WEIPRI + WTFLK (NPFLKA)

*  Save the track number of the stack particle:
      ISPARK (MKBMX2,NPFLKA) = NPFLKA
      NPARMA = NPARMA + 1
      NUMPAR (NPFLKA) = NPARMA
      NEVENT (NPFLKA) = 0
      DFNEAR (NPFLKA) = +ZERZER
*  ... to this point: don't change anything
*  Particle age (s)
      AGESTK (NPFLKA) = +ZERZER
      AKNSHR (NPFLKA) = -TWOTWO
*  Kinetic energy of the particle (GeV)
      TKEFLK (NPFLKA) = TKEFLK(NPFLKA)/1000.
*  Particle momentum
      PMOFLK (NPFLKA) = SQRT ( TKEFLK (NPFLKA) * ( TKEFLK (NPFLKA)
     &                       + TWOTWO * AM (IONID) ) )

*     TZFLK  (NPFLKA) = SQRT ( ONEONE - TXFLK (NPFLKA)**2
*    &                       - TYFLK (NPFLKA)**2 )
*  Polarization cosines:
      TXPOL  (NPFLKA) = -TWOTWO
      TYPOL  (NPFLKA) = +ZERZER
      TZPOL  (NPFLKA) = +ZERZER

*  Calculate the total kinetic energy of the primaries: don't change
      IF ( ILOFLK (NPFLKA) .EQ. -2 .OR. ILOFLK (NPFLKA) .GT. 100000 )
     &   THEN
         TKESUM = TKESUM + TKEFLK (NPFLKA) * WTFLK (NPFLKA)
      ELSE IF ( ILOFLK (NPFLKA) .NE. 0 ) THEN
         TKESUM = TKESUM + ( TKEFLK (NPFLKA) + AMDISC (ILOFLK(NPFLKA)) )
     &          * WTFLK (NPFLKA)
      ELSE
         TKESUM = TKESUM + TKEFLK (NPFLKA) * WTFLK (NPFLKA)
      END IF
      RADDLY (NPFLKA) = ZERZER
*  Here we ask for the region number of the hitting point.
*     NREG (NPFLKA) = ...
*  The following line makes the starting region search much more
*  robust if particles are starting very close to a boundary:
      CALL GEOCRS ( TXFLK (NPFLKA), TYFLK (NPFLKA), TZFLK (NPFLKA) )
      CALL GEOREG ( XFLK  (NPFLKA), YFLK  (NPFLKA), ZFLK  (NPFLKA),
     &              NRGFLK(NPFLKA), IDISC )
*  Do not change these cards:
      CALL GEOHSM ( NHSPNT (NPFLKA), 1, -11, MLATTC )
      NLATTC (NPFLKA) = MLATTC
      CMPATH (NPFLKA) = ZERZER
      CALL SOEVSV
      RETURN
*=== End of subroutine Source =========================================*
      END

