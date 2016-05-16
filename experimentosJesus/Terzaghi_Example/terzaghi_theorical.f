!     ************************************************************
!     *                                                          *
!     *                 T  E  R  Z  A  G  H  I                   *
!     *                                                          *
!     *    PROGRAM TO COMPUTE ANALYTICAL SOLUTION OF             * 
!     *    TERZAGHI PROBLEM WITH BOUNDARY AND INTIAL CONDITIONS  *
!     *    GIVEN IN INPUT'S FILE FROM 'TERZAGHI_BEN' PATH        *
!     *                                                          *
!     *----------------------------------------------------------*
!     *                                                          *
!     *               JESUS ALEXEI   - UFF/LNCC                  *
!     *                                                          *
!     *----------------------------------------------------------*
!     *                                                          *
!     ************************************************************
!
      IMPLICIT NONE
!
      INTEGER I,J,K,L
      INTEGER IPRSR, IOUT
      real(8) :: zero, PT5, one, two, three, four, pi
      real(8) :: EE, KAPPA, GAMMAW, CV, HH, SERIE
      REAL(8) :: XSTAR, T, tempo, cof, xm, xm2, p0
!
      REAL(8), DIMENSION(50) :: Z
!
      CHARACTER*30 NIPRSR, NIOUT 
      CHARACTER*3 ASTEP
!
      IPRSR = 20
!
      Z(1) = 0.1D0
!
      DO 5 I=2,50
         Z(I) = Z(I-1)+0.2D0
 5    CONTINUE
!
      P0 = 1.0D+1
!
      PI    = 4.0D0*DATAN(1.0D0)
!
      ZERO  = 0.0D0
      PT5   = 0.5D0
      ONE   = 1.0D0
      TWO   = 2.0D0
      THREE = 3.0D0
      FOUR  = 4.0D0
!
!.. SOLID PARAMETERS
!.... YOUNG MODULUS: EE=1.0 [kN/m2]
!  
      EE=1.0E+2
!
!.... POISSON RATIO: XNU=0.0
!
!      XNU=0.0D0
!
!.... BULK MODULUS: KBULK=EE/[3(1-2*XNU)] [kN/m2]
!
!      BULK=EE/(THREE*(ONE-TWO*XNU))
!
!.... SHEAR MODULUS: GSHEAR=EE/[2(1+XNU)] [kN/m2]
!
!      GSHEAR=EE/(TWO*(ONE+XNU))
!
!.... CONSTRAINED MODULUS: D= K+4G/3 [kN/m2]
!
!      D=BULK+FOUR*GSHEAR/THREE
!
!.... COMPRESSIBILITY COEFICIENT: MV=1/(K+4G/3) [m2/kN]
!
!        XMV=ONE/D
!
!.. POROUS PARAMETERS
!.... HYDRAULIC CONDUCTIVITY: KAPPA=1.0 [m/s] {pg.48Verruijt}
!
      KAPPA = 1.0D0
!
!.... FLUID PARAMETERS
!..... VOLUMETRIC WEIGHT OF WATER GAMMAW=10.0D0  [kN/m3]
!
      GAMMAW=1.0D0
!
!..... COEFICIENTE DE CONSOLIDACCAO (VERRUIJT) CV [m2/s]
!
      CV=KAPPA*EE/GAMMAW
!
!.... TAMANHO DA COLUNA
!
      HH=10.0D0
!
!.... STEP TIME
!
      DO 50 I=1,10
         TEMPO = DFLOAT(I)/1000  
         WRITE(*,3001) I,TEMPO
         T = CV*TEMPO/(HH*HH)
         WRITE(ASTEP,'(I3.3)') I
!                123456789+123456789+    123    4567
         NIPRSR='analytical/prsr_step'//ASTEP//'.dat'
         OPEN(UNIT=IPRSR, FILE=NIPRSR) 
!
         DO 30 J=1,50
            XSTAR = Z(J)/HH
            SERIE = 0.0D0
            DO 20 L=1,500
               XM  = PI*(TWO*DFLOAT(L)-ONE)*PT5
               XM2 = XM**2
               COF= (-ONE)**(L-1)
               SERIE=SERIE+(TWO*COF*DCOS(XM*XSTAR)*DEXP(-XM2*T))/XM
 20       CONTINUE
!
          WRITE(IPRSR,3001) J, Z(J), P0*SERIE
!
 30      CONTINUE
!
         CLOSE(IPRSR)
!
 50   CONTINUE
!
 3001 FORMAT(I5,2X,40(1PE15.8,2X)) 
 4000 FORMAT(2X,40(1PE15.8,2X)) 
 4001 FORMAT(I5,2X,I2,2X,40(1PE15.8,2X)) 
!
! ----------------------------------------------------------------- 
!
      END
