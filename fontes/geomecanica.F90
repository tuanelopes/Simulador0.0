!
!     *******************************************************************
!     *                                                                 *
!     *            * * *    G E O M E C H A N I C      * * *            *
!     *                                                                 *
!     *                                                                 *
!     *           A  FORTRAN 77  PROGRAM FOR SIMULATION                 *
!     *                                                                 *
!     *       GEOMECHANICAL BEHAVIOUR WITH HETEROGENEOUS YOUNG          *
!     *                                                                 *
!     *       MODULUS OF ELASTIC SATURATED POROUS MEDIA                 *
!     *                                                                 *
!     *       FRAMEWORK: FINITE ELEMENT AND B--BAR METHODS              *
!     *                                                                 *
!     *       BASED ON LINEAR CODE OF MARCIO MURAD - LNCC               *
!     *                                                                 *
!     *    --------------------------------------------------------     *
!     *                                                                 *
!     *       NON-LINEAR AND B-BAR FEM IMPLEMENTATION                   * 
!     *       FORMULATION JESUS ALEXEI - UFF-LNCC                       *
!     *                                                                 *
!     *                      MODIFY ON 30/11/2015                       *
!     *                                                                 *
!**** new module ********************************************************
!
      module mGeomecanica
!
      use mGlobaisEscalares, only: novaMalha
!
! ESCALARES
      REAL(8) :: KGYNG, RHOYNG ! YOUNG's GEOMETRIC MEAN, STRENGHT
      INTEGER :: NDAUX, NDT
      REAL*8  :: TMANDEL, GEOTIME
      INTEGER :: NROWB, NINTD, NESD, NED, NED2, NSTR
      INTEGER :: NEESQ2, NEESQ, NEE, NEE2, IOPT, IBBAR, NLOOPS
!
      INTEGER, ALLOCATABLE :: IDDIS(:,:)
      REAL(8), ALLOCATABLE :: GEOPRSR(:)
      REAL(8), ALLOCATABLE :: SIGMAT(:),SIGMA0(:),DIVU(:),DIVU0(:)
      REAL(8), ALLOCATABLE :: DIS(:,:), DIS0(:,:), STRSS(:,:,:), VDP(:,:)
! CREEP
      REAL*8, ALLOCATABLE, DIMENSION(:,:)   :: DTRL, AVCREP, AVSTRS
      REAL*8, ALLOCATABLE, DIMENSION(:,:)   :: STRSS0
      REAL*8, ALLOCATABLE, DIMENSION(:,:,:) :: HMTTG, ECREEP, AUXM
      INTEGER :: NROWB2, NWTNITER
      REAL*8  :: RESMAX, RESIDUAL
!
      real*8,  allocatable :: fDesloc(:,:)
!
      CHARACTER(len=128) :: YNG_IN
      real(8) :: kg,kgphi    ! media geometrica
      real(8) :: rho,rhophi  ! coeficiente da variancia (strenght)
!
      contains
!
!**** new *********************************************************************
!
      subroutine montarSistEqAlgGeo(TASK,satElem)
!
      use mAlgMatricial,     only: idDesloc, ALHSD, BRHSD, CLHSD, NEQD
      use mAlgMatricial,     only: LMD, IDIAGD, NALHSD, LOAD, FTOD, FACTOR
      use mMalha,            only: conecNodaisElem, listaDosElemsPorNo
      use mMalha,            only: x, numnp, nen, nsd, XC, numelReserv
      use mGlobaisEscalares, only: NDOFD, NLVECTD, TypeProcess
      use mGlobaisEscalares, only: optSolver
      use mHidroDinamicaRT,  only: pressaoElem
!
      implicit none
!
      integer :: i
      character(len=16) :: TASK
      REAL(8), DIMENSION(numelReserv) :: satElem
! 
!.... ALLOCATE MEMORY FOR GLOBAL EQUATION SYSTEM 
! 
      print*, " *** TASK FOR MATRIX CONFIGURATION = ", TASK
!
      IF (TASK.EQ.'right_hand_reset') THEN
!
         IF(.not.allocated(BRHSD)) allocate(BRHSD(NEQD))
!
!.... CLEAR RIGHT HAND SIDE FOR THE ELASTIC  PROBLEM 
! 
         BRHSD = 0.0D0
!
!.... ACCOUNT THE NODAL FORCES IN THE R.H.S. 
!
         IF (NLVECTD.GT.0) &
            CALL LOAD(idDesloc,fDesloc,brhsd,NDOFD,NUMNP,NLVECTD)
!
         IF (NLVECTD.GT.0) &
            CALL FTOD(idDesloc,DIS,fDesloc,NDOFD,NUMNP,NLVECTD)
!
         RETURN
!
      ENDIF
!
      IF (TASK.EQ.'bbarmatrix_elast') THEN
!
!.... MOUNT ONLY LEFT HAND SIDE AT ELEMENT LEVEL 
! 
         CALL BBARMTRX_ELAST(x, conecNodaisElem, alhsd, brhsd, idiagD, lmD)
!
         RETURN
!
      ENDIF
!
      IF (TASK.EQ.'bbarmatrix_creep') THEN
!
         if(.not.allocated(ALHSD)) allocate(ALHSD(NALHSD))
         if(.not.allocated(BRHSD)) allocate(BRHSD(NEQD))
!
!.... CLEAR LEFT AND RIGHT HAND SIDE FOR THE ELASTIC  PROBLEM 
! 
         ALHSD = 0.0D00
         BRHSD = 0.0D00
!
!.... ACCOUNT THE NODAL FORCES IN THE R.H.S. 
!
         IF (NLVECTD.GT.0) &
            CALL LOAD(idDesloc,fDesloc,brhsd,NDOFD,NUMNP,NLVECTD)
!
         IF (NLVECTD.GT.0) &
            CALL FTOD(idDesloc,DIS,fDesloc,NDOFD,NUMNP,NLVECTD)
!
!.... MOUNT LEFT AND RIGHT SIDE AT ELEMENT LEVEL 
!
         CALL BBARMTRX_CREEP(X,pressaoElem,satElem,conecNodaisElem, &
              alhsd,brhsd,idiagD,lmD)
!
         RETURN
!
      ENDIF
!
      IF (TASK.EQ.'right_hand_elast') THEN
!
!.... MOUNT VECTOR SOURCE OF RIGHT HAND SIDE IN ELASTIC PROBLEM 
!
         CALL VECTOR_SOURC(satElem, pressaoElem, GEOPRSR, brhsd, lmD)
!
         RETURN
!
      ENDIF
!
      END SUBROUTINE
!
!**** new ********************************************************************
!
      SUBROUTINE BBARMTRX_CREEP(x,p,satElem,conecNodaisElem,alhsd,brhsd,idiagD,lmD)
!
!.... PROGRAM TO CALCULATE STIFNESS MATRIX AND FORCE ARRAY FOR THE 
!        STOKE'S DISPLACEMENT  ELEMENT AND 
!        ASSEMBLE INTO THE GLOBAL LEFT-HAND-SIDE MATRIX 
!        AND RIGHT-HAND SIDE VECTOR 
! 
      use mGlobaisEscalares, only: ndofD,ndofP,nrowsh,nnp,S3DIM,IBBAR,optSolver
      use mGlobaisArranjos,  only: grav, c
      use mPropGeoFisica,    only: YOUNG, GRAINBLK, RHODVECT, PORELAGR
      use mPropGeoFisica,    only: RHOW, RHOO, BULKWATER, BULKOIL
      use mPropGeoFisica,    only: GEOFORM, GEOINDIC, BULK
      use mPropGeoFisica,    only: POISVECT,YUNGVECT,RHODVECT
      use mPropGeoFisica,    only: PORE
!
      use mAlgMatricial,     only: rowdot, coldot
      use mAlgMatricial,     only: neqD, nalhsD
      use mSolverGaussSkyline, only: addrhs, addlhs
!
      use mFuncoesDeForma,   only: shgq, shlq
      use mMalha,            only: local, numnp, multab
      use mMalha,            only: nsd, numel, numelReserv, nen
      use mSolverPardiso, only: addlhsCRS, ApGeo, AiGeo
      use mSolverHypre
!
      use mGlobaisEscalares, only: PRESSAOREF, COTAREF
!
      implicit none
!
      real*8,  intent(in)    :: x(nsd,numnp), p(1,numelReserv)
      real*8,  intent(in)    :: satelem(numelReserv)
      integer, intent(in)    :: conecNodaisElem(nen,numel)
      real(8), intent(inout) :: alhsd(nalhsD), brhsd(neqD)
      integer, intent(in)    :: idiagD(*)
      integer, intent(in)    :: lmD(ned2,nen,numel)
!
      real(8) :: xl(nesd,nen), disl(ned2,nen)
      real(8) :: ELRESFD(nee2), ELEFFMD(nee2,nee2)
!
!.... REMOVE ABOVE CARD FOR SINGLE-PRECISION OPERATION  
!  
      LOGICAL DIAG,QUAD,ZERODL,LSYM
!
!.... LOCAL VECTORS AND MATRIZES
!
      REAL(8), DIMENSION(NROWB,NESD)    :: BBARJ, BBARI, QIXIBBAR
      REAL(8), DIMENSION(NROWB,NROWB)   :: QIXI
      REAL(8), DIMENSION(NROWB)         :: TENSAO, PRESSURE
!
      REAL(8) :: SHGD(NROWSH,NEN,NINTD), SHLD(NROWSH,NEN,NINTD)
      REAL(8) :: SHGBR(NROWSH,NEN)
      REAL(8) :: DETD(NINTD), R(NINTD), WD(NINTD)
!
      REAL(8) :: C1, BULKROCK, BIOTCOEF, POROSITY, ONEMPORE
      REAL(8) :: RHOMAT, DENSOLID, DENFLUID, SATURAT, ONEMSAT
      REAL(8) :: BULKGRAIN, BULKFLUD, POISSON, RHOWLIN, RHOOLIN
!
      INTEGER :: NEL, I, J, L, K
!
      integer :: tid, omp_get_thread_num,omp_get_num_threads, numThreads
      integer :: numPrimeiroElemento, numUltimoElemento, inicioSol, fimSol
!
!.... GENERATION OF LOCAL SHAPE FUNCTIONS AND WEIGHT VALUES 
! 
      stop "BBARMTRX_CREEP" 
      CALL SHLQ(SHLD,WD,NINTD,NEN)
!
!      CONSISTENT MATRIX
! 
      DIAG       = .FALSE. 
      TENSAO     = 0.0D0
      tid        = 1
      numThreads = 1

! !$OMP PARALLEL FIRSTPRIVATE(tid) &
! !$OMP PRIVATE (numPrimeiroElemento, numUltimoElemento, inicioSol,fimSol) &
! !$OMP PRIVATE (I,J,L,K,NEL,BBARI,BBARJ,QIXIBBAR,QIXI,PRESSURE,XL,DISL,TENSAO,C1,IOPT,QUAD,LSYM) &
! !$OMP PRIVATE (R,SHGD,SHGBR,DETD) &
! !$OMP REDUCTION(+:ELEFFMD,ELRESFD) !,BRHSD,ALHSD)

! 
! #ifdef withOMP
!        tid=tid+omp_get_thread_num()
!        numThreads=omp_get_num_threads()
! #endif
!
      numPrimeiroElemento = 1
      numUltimoElemento   = numel
      call dividirTrabalho(numPrimeiroElemento, numUltimoElemento, numThreads, tid-1, inicioSol, fimSol)
!
! !$OMP PARALLEL DO ORDERED SCHEDULE(DYNAMIC)
!
      DO 500 NEL=inicioSol,fimSol
!      DO 500 NEL=1,numel
!
         BBARI    = 0.0D0
         BBARJ    = 0.0D0
         QIXIBBAR = 0.0D0
         QIXI     = 0.0D0
         R        = 0.0D0
!
!.... SETUP DENSITY ELEMENT
!
         POROSITY = GEOINDIC('POROSTY',GEOFORM(NEL))
         DENSOLID = GEOINDIC('ROKDENS',GEOFORM(NEL))
         DENFLUID = GEOINDIC('FLUDENS',GEOFORM(NEL))
         BULKFLUD = GEOINDIC('BLKFLUD',GEOFORM(NEL))
         DENFLUID = DENFLUID*(1.0D0+BULKFLUD*(GEOPRSR(NEL)-PRESSAOREF))
!
         IF (GEOFORM(NEL).EQ.'RESERVATORIO') THEN
            RHOWLIN  = RHOW*(1.0D0+(P(1,NEL)-PRESSAOREF)/BULKWATER)
            RHOOLIN  = RHOO*(1.0D0+(P(1,NEL)-PRESSAOREF)/BULKOIL)
            SATURAT  = satElem(NEL)
            ONEMSAT  = 1.0D0-satElem(NEL)
            BULKFLUD = SATURAT/BULKWATER + ONEMSAT/BULKOIL
            DENFLUID = SATURAT*RHOWLIN   + ONEMSAT*RHOOLIN
            DENFLUID = DENFLUID*(1.0D0+BULKFLUD*(P(1,NEL)-PRESSAOREF))
            POROSITY = PORE(NEL)
         ENDIF
!
         RHOMAT    = POROSITY*DENFLUID + (1.0D0-POROSITY)*DENSOLID
!
!         write(2025,*) 'AQUI ====>', nel, RHOMAT, DENFLUID, densolid ! rhomat, poisson, bulkgrain
!
!.... SETUP BIOT COEFICIENT ALSO CALLED BIOT'S ALPHA
!
         POISSON   = GEOINDIC('POISSON',GEOFORM(NEL))
         BULKGRAIN = GEOINDIC('BLKGRIN',GEOFORM(NEL))
         BULKROCK  = BULK(YOUNG(NEL), POISSON, S3DIM) 
         BIOTCOEF  = 1.0D0 - BULKROCK/BULKGRAIN
!
!.... SETUP ELEMENT PRESSURE VECTOR
!
         PRESSURE(1) = BIOTCOEF*GEOPRSR(NEL)
         PRESSURE(2) = PRESSURE(1)
         PRESSURE(3) = 0.0D0
         PRESSURE(4) = PRESSURE(1)
!
!LINE TO TEST EQUIVALENT GRAVITY:                  PRESSURE    = 0.0D0
!LINE TO TEST EQUIVALENT GRAVITY:         PRESSURE(2) = GEOPRSR(NEL)
!
         IF (GEOFORM(NEL).EQ.'RESERVATORIO') THEN
            PRESSURE(1) = BIOTCOEF*P(1,NEL)
            PRESSURE(2) = PRESSURE(1)
            PRESSURE(4) = PRESSURE(1)
!
!LINE TO TEST EQUIVALENT GRAVITY:                 PRESSURE    = 0.0D0
!LINE TO TEST EQUIVALENT GRAVITY:        PRESSURE(2) = P(1,NEL)/BIOTCOEF
!
         ENDIF
!
!.... CLEAR STIFFNESS MATRIX AND FORCE ARRAY 
! 
      ELEFFMD=0.0D0 
      ELRESFD=0.0D0 
! 
!.... LOCALIZE COORDINATES AND DIRICHLET B.C. 
! 
      CALL LOCAL(conecNodaisElem(1,NEL),X,XL,NEN,NSD,NESD) 
      CALL LOCAL(conecNodaisElem(1,NEL),DIS,DISL,NEN,NDOFD,NED2) 
!
      QUAD = .TRUE. 
      IF (conecNodaisElem(3,NEL).EQ.conecNodaisElem(4,NEL)) QUAD = .FALSE.  
      CALL SHGQ(XL,DETD,SHLD,SHGD,NINTD,NEL,QUAD,NEN) 
!
!.... SETUP FOR AXISYMMETRIC OPTION
!
      IF (IOPT.EQ.2) THEN
         DO 100 L=1,NINTD
            R(L)    = ROWDOT(SHGD(NROWSH,1,L),XL,NROWSH,NESD,NEN)
            DETD(L) = DETD(L)*R(L) 
 100     CONTINUE
      ENDIF
!
!.... MOUNT STIFFNESS MATRIX 
!
!.... CALCULATE MEAN VALUES OF SHAPE FUNCTION GLOBAL DERIVATIVES
!.... FOR MEAN-DILATATIONAL B-BAR FORMULATION
!
      CALL XMEANSH(SHGBR,WD,DETD,R,SHGD,NEN,NINTD,IOPT,NESD,NROWSH)
!
!.... LOOP OVER INTEGRATIONN POINTS
!
      DO 400 L=1,NINTD 
!
!... ... SETUP TANGENT MATRIX QIXI: ORDER 4X4 FOR MULTIPLICATION
!
         CALL TANG2QX(HMTTG(NEL,L,1:16),QIXI)
!
!...... SETUP ELEMENT DETERMINANT AND GAUSS WEIGHTS
!
         C1=DETD(L)*WD(L)
!
!.... ...UPLOAD B-BAR MATRIX AT NODE J
!
         DO 200 J=1,NEN
            CALL SETBB(BBARJ,SHGD(1:NROWSH,J,L),  &
     &           SHGBR(1:NROWSH,J),R(L),NROWSH,NROWB,IOPT,IBBAR)
!
!.... .... MULTIPLY QIXI*BBARJ ===>QIXIBBAR
!
            QIXIBBAR=matmul(QIXI,BBARJ)
!             CALL MULTAB(QIXI,BBARJ,QIXIBBAR,4,4,4,4,4,2,1)
!
!.... .... UPLOAD B-BAR MATRIX AT NODE I
!
            DO 200 I=1,NEN
               CALL SETBB(BBARI,SHGD(1:NROWSH,I,L),  &
     &              SHGBR(1:NROWSH,I),R(L),NROWSH,NROWB,IOPT,IBBAR)
!
!.... .......  MOUNT ELEMNT STIFFNESS NODAL MATRIX:
!.... .......  K^E_IJ= MULTIPLY BBAR^T_I*(QIXI*BBAR_J)
!
               ELEFFMD(NED2*I-1,NED2*J-1)= ELEFFMD(NED2*I-1,NED2*J-1)  &
     &             +COLDOT(BBARI(1:4,1),QIXIBBAR(1:4,1),4)*C1
!
               ELEFFMD(NED2*I-1,NED2*J)= ELEFFMD(NED2*I-1,NED2*J)      &
     &             +COLDOT(BBARI(1:4,1),QIXIBBAR(1:4,2),4)*C1
!
               ELEFFMD(NED2*I,NED2*J-1)= ELEFFMD(NED2*I,NED2*J-1)      &
     &             +COLDOT(BBARI(1:4,2),QIXIBBAR(1:4,1),4)*C1
!
               ELEFFMD(NED2*I,NED2*J)= ELEFFMD(NED2*I,NED2*J)          &
     &             +COLDOT(BBARI(1:4,2),QIXIBBAR(1:4,2),4)*C1
 200        CONTINUE
      BBARI=0.0D0
!
      DO 350 I=1,NEN
         CALL SETBB(BBARI,SHGD(1:NROWSH,I,L), &
     &              SHGBR(1:NROWSH,I),R(L),NROWSH,NROWB,IOPT,IBBAR)
         DO 340 K=1,NROWB
            TENSAO(K) = STRSS(NEL,L,K)-PRESSURE(K)
 340     CONTINUE
         ELRESFD(NED2*I-1)= ELRESFD(NED2*I-1)  &
     &                      - COLDOT(BBARI(1:4,1),TENSAO,4)*C1 
!
         ELRESFD(NED2*I)  = ELRESFD(NED2*I)    &
     &                      - COLDOT(BBARI(1:4,2),TENSAO,4)*C1 &
     &                      + RHOMAT*GRAV(2)*SHGD(3,I,L)*C1
 350  CONTINUE
!
 400  CONTINUE
!  
!.... COMPUTATION OF DIRICHLET B.C. CONTRIBUTION 
! 
      CALL ZTEST(DISL,NEE2,ZERODL) 
!
      IF(.NOT.ZERODL)CALL KDBCGEO(ELEFFMD,ELRESFD,DISL,NEE2,LMD(1,1,NEL),NEL) 
!
!.... ASSEMBLE ELEMENT STIFFNESS MATRIX AND FORCE ARRAY INTO GLOBAL 
!....     LEFT-HAND-SIDE MATRIX AND RIGHT-HAND SIDE VECTOR 
! 
      LSYM=.TRUE.

! $OMP ORDERED

     if (optSolver=='skyline')   then
        CALL ADDLHS(ALHSD,ELEFFMD,LMD(1,1,NEL),IDIAGD,NEE2,DIAG,LSYM) 
     endif     

     if (optSolver=='pardiso')   then
           call addlhsCRS(ALHSD,ELEFFMD,LMD(1,1,nel),ApGeo,AiGeo,NEE2) 
     endif
     
     if (optSolver=='hypre')   then 
           call addnslHYPRE(A_HYPRE_G, ELEFFMD, LMD(1,1,nel), NEE2, LSYM)
     endif
!
     CALL ADDRHS(BRHSD,ELRESFD,LMD(1,1,NEL),NEE2) 

! $OMP END ORDERED
       
 500   CONTINUE 
! OMP END DO


      if (optSolver=='hypre')   then     ! creep
       do i = 1, neqD     ! creep
           rows_G(i) = i-1      ! creep
       end do     ! creep
      call adicionarValoresVetor_HYPRE(b_HYPRE_G, 1, neqD, rows_G, BRHSD)

      endif  ! creep


! $OMP END PARALLEL

!
      RETURN
!
 2000 FORMAT(5(1PE15.8,2X))
 2222  FORMAT('ELEMENTO (NEL)=',I5,2X,' GAUSS POINT (L)=',I2/5X   &
     &' C11, C21, C31, C41 =',4(1PE9.2,2X)/5X,                    &
     &' C12, C22, C32, C42 =',4(1PE9.2,2X)/5X,                    &
     &' C13, C23, C33, C43 =',4(1PE9.2,2X)/5X,                    &
     &' C14, C24, C34, C44 =',4(1PE9.2,2X)//)
!
      END SUBROUTINE
!
!**** new *******************************************************************
!
      subroutine BBARMTRX_ELAST(x, conecNodaisElem,alhsd, brhsd, idiagD, lmD)
!
!.... PROGRAM TO CALCULATE STIFNESS MATRIX FOR THE ELASTIC DISPLACEMENT    
!        ELEMENT AND ASSEMBLE INTO THE GLOBAL LEFT-HAND-SIDE MATRIX 
!
      use mGlobaisEscalares, only: ndofD, nrowsh, IBBAR, optSolver
      use mGlobaisArranjos,  only: grav, c
      use mPropGeoFisica,    only: YOUNG, GEOFORM, GEOINDIC
      use mPropGeoFisica,    only: POISVECT, YUNGVECT, RHODVECT
      use mAlgMatricial,     only: rowdot, coldot
      use mAlgMatricial,     only: neqD, nalhsD
      use mSolverGaussSkyline, only: addrhs, addlhs
      use mFuncoesDeForma,   only: shgq, shlq
      use mMalha,            only: local, numnp, multab
      use mMalha,            only: nsd, numel, nen, numelReserv
      use mSolverPardiso, only: addlhsCRS, ApGeo, AiGeo
      use mSolverHypre
!
      implicit none
!
      real*8,  intent(in)    :: x(nsd,numnp)
      integer, intent(in)    :: conecNodaisElem(nen,numel)
      real(8), intent(inout) :: alhsd(nalhsD), brhsd(neqD)
      integer, intent(in)    :: idiagD(*)
      integer, intent(in)    :: lmD(ned2,nen,numel)
!
      real(8) :: xl(nesd,nen), disl(ned2,nen)
      real(8) :: ELRESFD(nee2), ELEFFMD(nee2,nee2)
!
!.... REMOVE ABOVE CARD FOR SINGLE-PRECISION OPERATION  
!  
      LOGICAL DIAG,QUAD,ZERODL,LSYM
!
!.... LOCAL VECTORS AND MATRIZES
!
      REAL(8) :: BBARJ(NROWB,NESD), BBARI(NROWB,NESD)
      REAL(8) :: QIXI(NROWB,NROWB), QIXIBBAR(NROWB,NESD)
      INTEGER :: NEL, I, J, L, II, JJ
      REAL(8) :: DETD(NINTD), R(NINTD), WD(NINTD)
      REAL(8) :: SHGD(NROWSH,NEN,NINTD), SHLD(NROWSH,NEN,NINTD)
      REAL(8) :: SHGBR(NROWSH,NEN)
      REAL(8) :: CBBAR(NROWB, NROWB)
      REAL(8) :: C1, POISSON
!
      character(len=11) :: nomeA
      character(len=11) :: nomeB
      integer :: ierr
      nomeA = "Alhs00.out."
      nomeB = "brhs00.out."

!
!.... GENERATION OF LOCAL SHAPE FUNCTIONS AND WEIGHT VALUES 
! 
      CALL SHLQ(SHLD,WD,NINTD,NEN)
!
!.... CONSISTENT MATRIX
! 
      DIAG = .FALSE.
!
      BBARI = 0.0D0
      BBARJ = 0.0D0
      R     = 0.0D0
!
      DO 500 NEL=1,NUMEL
! 
!.... CLEAR STIFFNESS MATRIX AND FORCE ARRAY 
! 
      ELEFFMD=0.D0
      ELRESFD=0.D0
! 
!.... LOCALIZE COORDINATES AND DIRICHLET B.C. 
! 
      CALL LOCAL(conecNodaisElem(1,NEL),X,XL,NEN,NSD,NESD)
      CALL LOCAL(conecNodaisElem(1,NEL),DIS,DISL,NEN,NDOFD,NED2)
      QUAD = .TRUE.
      IF (conecNodaisElem(3,NEL).EQ.conecNodaisElem(4,NEL)) QUAD = .FALSE.
! 
      CALL SHGQ(XL,DETD,SHLD,SHGD,NINTD,NEL,QUAD,NEN)
!
!.... SETUP FOR AXISYMMETRIC OPTION
!
      IF (IOPT.EQ.2) THEN
         DO 100 L=1,NINTD
            R(L)    = ROWDOT(SHGD(NROWSH,1,L),XL,NROWSH,NESD,NEN)
            DETD(L) = DETD(L)*R(L)
 100    CONTINUE
      ENDIF
!
      POISSON = GEOINDIC('POISSON',GEOFORM(NEL))
!
!.... SETUP STOCHASTIC ELASTICITY TENSOR FOR BBAR METHOD 
! 
      CALL SETUPC(CBBAR,YOUNG(NEL),POISSON,NROWB,IOPT)
!
!.... FORM STIFFNESS MATRIX 
!
!.... .. CALCULATE MEAN VALUES OF SHAPE FUNCTION GLOBAL DERIVATIVES
!.... .. FOR MEAN-DILATATIONAL B-BAR FORMULATION
!
         CALL XMEANSH(SHGBR,WD,DETD,R,SHGD,NEN,NINTD,IOPT,NESD,NROWSH)
!
!.... .. LOOP OVER INTEGRATIONN POINTS
!
         DO 400  L=1,NINTD
!
!.... ...SETUP ELEMENT DETERMINANT AND GAUSS WEIGHTS
!
          C1=DETD(L)*WD(L)
!
!.... ... UPLOAD B-BAR MATRIX AT NODE J
!
          DO 200 J=1,NEN
!
            CALL SETBB(BBARJ,SHGD(1:NROWSH,J,L), &
     &            SHGBR(1:NROWSH,J),R(L),NROWSH,NROWB,IOPT,IBBAR)
!
!.... .... MULTIPLY CBBAR*BBARJ ===>QIXIBBAR
!
            CALL MULTAB(CBBAR,BBARJ,QIXIBBAR,4,4,4,4,4,2,1)
!
!.... .... UPLOAD B-BAR MATRIX AT NODE I
!
            DO 200 I=1,NEN
!
            CALL SETBB(BBARI,SHGD(1:NROWSH,I,L), &
     &            SHGBR(1:NROWSH,I),R(L),NROWSH,NROWB,IOPT,IBBAR)
!
!.... .... MOUNT ELEMNT STIFFNESS NODAL MATRIX:
!.... ...        K^E_IJ= MULTIPLY BBAR^T_I*(QIXI*BBAR_J)
!
             ELEFFMD(NED2*I-1,NED2*J-1)= ELEFFMD(NED2*I-1,NED2*J-1) &
     &             +COLDOT(BBARI(1:4,1),QIXIBBAR(1:4,1),4)*C1
!
             ELEFFMD(NED2*I-1,NED2*J)= ELEFFMD(NED2*I-1,NED2*J) &
     &             +COLDOT(BBARI(1:4,1),QIXIBBAR(1:4,2),4)*C1
!
             ELEFFMD(NED2*I,NED2*J-1)= ELEFFMD(NED2*I,NED2*J-1) &
     &             +COLDOT(BBARI(1:4,2),QIXIBBAR(1:4,1),4)*C1
!
             ELEFFMD(NED2*I,NED2*J)= ELEFFMD(NED2*I,NED2*J)  &
     &             +COLDOT(BBARI(1:4,2),QIXIBBAR(1:4,2),4)*C1
!
 200      CONTINUE
 400  CONTINUE
!  
!      COMPUTATION OF DIRICHLET B.C. CONTRIBUTION 
! 
      CALL ZTEST(DISL,NEE2,ZERODL)
!
      IF(.NOT.ZERODL)CALL KDBCGEO(ELEFFMD,ELRESFD,DISL,NEE2,LMD(1,1,NEL),NEL)
!
!.... ASSEMBLE ELEMENT STIFFNESS MATRIX AND FORCE ARRAY INTO GLOBAL 
!        LEFT-HAND-SIDE MATRIX AND RIGHT-HAND SIDE VECTOR 
! 
      LSYM=.TRUE.
!
      if (optSolver=='skyline')   then
         CALL ADDLHS(ALHSD,ELEFFMD,LMD(1,1,NEL),IDIAGD,NEE2,DIAG,LSYM) 
      endif

      if (optSolver=='pardiso')   then      
          call addlhsCRS(ALHSD,ELEFFMD,LMD(1,1,nel),ApGeo,AiGeo,NEE2) 
      endif

     if (optSolver=='hypre')   then ! elast
           !write(*,*) "  call addnslHYPRE(A_HYPRE_G, ELEFFMD" ! elast
           !write(*,*) nel, "  lmD =" , LMD(:,:,nel) 
           call addnslHYPRE(A_HYPRE_G, ELEFFMD, LMD(1,1,nel), NEE2, LSYM) ! elast
     endif ! elast
!
      CALL ADDRHS(BRHSD,ELRESFD,LMD(1,1,NEL),NEE2)
!
        DO 450 II=1,NEE2
        DO 450 JJ=1,NEE2
           AUXM(NEL,II,JJ)=ELEFFMD(II,JJ)
 450   CONTINUE
!
 500   CONTINUE

       if (optSolver=='hypre')   then   ! elast
          do i = 1, neqD   ! elast
             rows_G(i) = i-1    ! elast
          end do   ! elast
          write(*,*) " brhsD =", brhsd(1:5)
          write(*,*) " brhsD =", brhsd(neqd-4:neqd)
       
          call adicionarValoresVetor_HYPRE(b_HYPRE_G, 1, neqD, rows_G, BRHSD)   ! elast
                                                                   ! elast
        ! call fecharMatriz_HYPRE            (A_HYPRE_G, parcsr_A_G)   ! elast
          call fecharVetor_HYPRE             (b_HYPRE_G, par_b_G   )   ! elast
          call fecharVetor_HYPRE             (u_HYPRE_G, par_u_G   )   ! elast

       end if   ! elast
!
      RETURN
!
 2000 FORMAT(5(1PE15.8,2X))
!
      END SUBROUTINE
!
!*** NEW *** FOR STOCHASTIC YOUNG MODULUS ******************************* 
!
      SUBROUTINE SETUPC(CBBAR,YOUNG,POISSON,NROWB,IOPT)
! 
!     PROGRAM TO SETUP ELASTICITY TENSOR 
! 
!      IMPLICIT REAL*8 (A-H,O-Z) 
      IMPLICIT NONE
! 
      REAL*8  :: YOUNG,POISSON
      INTEGER :: NROWB,IOPT
      REAL(8) :: CBBAR(NROWB, NROWB)
!
      REAL(8) :: AMU2,ALAM

!     REMOVE ABOVE CARD FOR SINGLE PRECISION OPERATION 
!
!..... SET MATERIAL PARAMETERS FOR OUT-OF-PLANE COMPONENTS
!
!      C(5,M) == YOUNG   e   C(6,M)==POISSON
!          AMU2 = YOUNG/(ONE+POISSON)
!          ALAM = AMU2*POISSON/(ONE-TWO*POISSON)
!
          CBBAR = 0.0
          AMU2 = YOUNG/(1.0D0+POISSON)
          ALAM = AMU2*POISSON/(1.0D0-2.0D0*POISSON)
!
!..... COLUMN MATRIX D3 
!
          CBBAR(1,4) = ALAM
          CBBAR(2,4) = ALAM
          CBBAR(3,4) = 0.0D0
          CBBAR(4,4) = ALAM + AMU2
!
!..... TRANSPOSE OF D3
!
          CBBAR(4,1) = CBBAR(1,4)
          CBBAR(4,2) = CBBAR(2,4)
          CBBAR(4,3) = CBBAR(3,4)
!
!..... SET MATERIAL PARAMETERS FOR IN-PLANE COMPONENTS
!
!..... .. SET STRESS PLANE CONDITION IOPT.EQ.0
!
          IF (IOPT.EQ.0) ALAM = ALAM*AMU2/(ALAM + AMU2)
!
!..... .. MATRIX D_33: UP TRIANGULAR
!
          CBBAR(1,1) = ALAM + AMU2
          CBBAR(1,2) = ALAM
          CBBAR(2,2) = CBBAR(1,1)
          CBBAR(1,3) = 0.0D0
          CBBAR(2,3) = 0.0D0
          CBBAR(3,3) = 0.5D0*AMU2
!
!..... .. MATRIX D_33: LOW TRIANGULAR
! 
          CBBAR(2,1) = CBBAR(1,2)
          CBBAR(3,1) = CBBAR(1,3)
          CBBAR(3,2) = CBBAR(2,3)
!
 100   CONTINUE
!
       RETURN
!
 2000 FORMAT(5(1PE15.8,2X))
!
       END SUBROUTINE
!
!**** FROM  HUGHES: THE FINITE ELEMENT METHOD ** PAG 760 *************** 
!
      SUBROUTINE XMEANSH(SHGBAR,W,DET,R,SHG,NEN,NINT,IOPT,NESD,NROWSH)
!
!.... PROGRAM TO CALCULATE MEAN VALUES OF SHAPE FUNCTION
!        GLOBAL DERIVATIVES FOR B-BAR METHOD
!
!        NOTE: IF IOPT.EQ.2, DET(L)=DET(L)*R(L) UPON ENTRY
!
      use mAlgMatricial, only: coldot
!
      IMPLICIT NONE
!
      INTEGER :: NEN, NROWSH, NINT, IOPT, NESD
      REAL(8) :: SHGBAR(3,NEN),W(*),DET(*),R(*),SHG(NROWSH,NEN,*)
!
      REAL(8) :: VOLINV, TEMP1, TEMP2
      INTEGER :: L, I, J
!
      CALL CLEAR(SHGBAR,3*NEN)
!
      VOLINV = 1.0D0/COLDOT(W,DET,NINT)
!
      DO 300 L=1,NINT
         TEMP1 = W(L)*DET(L)*VOLINV
         IF (IOPT.EQ.2) TEMP2 = TEMP1/R(L)
!
         DO 200 J=1,NEN
            DO 100 I=1,NESD
               SHGBAR(I,J) = SHGBAR(I,J) + TEMP1*SHG(I,J,L)
 100        CONTINUE
 200     CONTINUE
         IF (IOPT.EQ.2) SHGBAR(3,J)=SHGBAR(3,J)+TEMP2*SHG(3,J,L)
 300  CONTINUE
!
      RETURN
!
      END SUBROUTINE
!
!**** FROM  HUGHES: THE FINITE ELEMENT METHOD ** PAG 780 **************** 
!
      SUBROUTINE SETBB(BBAR,SHG,SHGBAR,R,NROWSH,NROWB,IOPT,IBBAR) 
!
!..... PROGRAM TO SET UP THE STRAIN-DISPLACEMENT MATRIX "B" FOR
!         TWO-DIMENSIONAL CONTINUUM ELEMENTS
!
!         IBBAR = 0,  STANDARD B-MATRIX
!
!         IBBAR = 1,  MEAN-DILATATIONAL B-MATRIX
!
      IMPLICIT NONE
!
!.... DEACTIVATE ABOVE CARD(S) FOR SINGLE-PRECISION OPERATION
!
      integer :: NROWSH,NROWB,IOPT,IBBAR
      real*8 ::  SHG(NROWSH), SHGBAR(NROWSH), BBAR(NROWB,2), R
!
      real*8 :: temp1, temp2, temp3, constb
!
       BBAR(1,1) = SHG(1)
       BBAR(2,1) = 0.0D0
       BBAR(3,1) = SHG(2)
       BBAR(4,1) = 0.0D0
!
       BBAR(1,2) = 0.0D0 
       BBAR(2,2) = SHG(2)
       BBAR(3,2) = SHG(1)
       BBAR(4,2) = 0.0D0
!
!.... AXISYMMETRIC CASE
!
      IF (IOPT.EQ.2) THEN 
        BBAR(4,1) = SHG(3)/R
        BBAR(4,2) = 0.0D0
      ENDIF
!
      IF (IBBAR.EQ.0) RETURN
!
      CONSTB = 1.0D0/3.0D0
!
!.... ADD CONTRIBUTIONS TO FORM B-BAR
!
!.... DEFINE VALUES FOR ROW=4 CASES: PLAIN STRAIN OR AXISYMMETRIC
!
      IF (IOPT.EQ.2) THEN
         TEMP3 = CONSTB*(SHGBAR(3)-SHG(3)/R)
         BBAR(1,1) = BBAR(1,1)+TEMP3
         BBAR(2,1) = BBAR(2,1)+TEMP3
         BBAR(4,1) = BBAR(4,1)+TEMP3
       ELSE
         BBAR(4,1) = 0.0D0
         BBAR(4,1) = 0.0D0
      ENDIF
!
      TEMP1 = CONSTB*(SHGBAR(1)-SHG(1))
      TEMP2 = CONSTB*(SHGBAR(2)-SHG(2))
!
      BBAR(1,1) = BBAR(1,1) + TEMP1
      BBAR(2,1) = BBAR(2,1) + TEMP1
      BBAR(4,1) = BBAR(4,1) + TEMP1
!
      BBAR(1,2) = BBAR(1,2) + TEMP2
      BBAR(2,2) = BBAR(2,2) + TEMP2
      BBAR(4,2) = BBAR(4,2) + TEMP2
!
      RETURN
!
      END SUBROUTINE
!
!**** NEW **************************************************************** 
!
      SUBROUTINE KDBCGEO(ELEFFM,ELRESF,DL,NEE,LM,NEL)
! 
!.... PROGRAM TO ADJUST LOAD VECTOR FOR PRESCRIBED DISPLACEMENT 
!     BOUNDARY CONDITION 
! 
      USE mGlobaisEscalares
!
      IMPLICIT NONE
!
      INTEGER :: NEE, NEL
      REAL*8  :: ELEFFM(NEE,*),ELRESF(*),DL(*)
      INTEGER :: LM(*)
!
      INTEGER :: I,J,L
      REAL(8) :: VAL
! 
!    THIS VERSION OF KDBC IS ONLY VALID FOR THERMOELASTIC CONS. 
! 
      DO 200 J=1,NEE
         L   = LM(J) 
         VAL = DL(J)     
         IF(L.GT.0)      GO TO 200
         IF(VAL.EQ.0.0D0) GO TO 200 
         DO 100 I=1,NEE
            ELRESF(I)=ELRESF(I)-ELEFFM(I,J)*VAL
100      CONTINUE
! 
200   CONTINUE
! 
      RETURN
!
      END SUBROUTINE
!
!**** NEW **** FOR INCOMPRESSIBILITY *************************************** 
!
      SUBROUTINE POS4ITER(x, conecNodaisElem, STRESS, DIVU)
!
!.... PROGRAM TO COMPUTE ELEMENT STRESS AND VOLUMETRIC DEFORMATION 
!....            FROM GEOMECHANIC ELASTIC MODEL
!
      use mMalha,            only: nsd, numnp, numel, numelReserv
      use mMalha,            only: XC, nen, LOCAL
      use mGlobaisEscalares, only: NNP, NDOFD, nrowsh
      use mAlgMatricial,     only: rowdot, coldot
      use mFuncoesDeForma,   only: shgq, shlq
      use mPropGeoFisica,    only: YOUNG, GEOFORM, GEOINDIC 
!
      IMPLICIT NONE
!
      REAL*8,  intent(in)              :: X(NSD,NUMNP)
      INTEGER, intent(in)              :: conecNodaisElem(NEN,NUMEL)
      REAL(8), DIMENSION(NROWB,NUMEL)  :: STRESS
      REAL(8), DIMENSION(NUMEL)        :: DIVU
!
      LOGICAL QUAD
!
      INTEGER :: I,J,K,L,NEL,NODE
!
      REAL(8), DIMENSION(NESD,NEN)          :: XL
      REAL(8), DIMENSION(NED2,NEN)          :: DISL
      REAL(8), DIMENSION(NINTD)             :: WD, DETD, R
      REAL(8), DIMENSION(NROWSH,NEN,NINTD)  :: SHLD,SHGD
      REAL(8), DIMENSION(NROWSH,NEN)        :: SHGBR
      REAL(8), DIMENSION(NROWB,NROWB)       :: CBBAR
!
!.... LOCAL VECTORS AND MATRIZES
!
      REAL(8), DIMENSION(NROWB,NESD)    :: BBARJ
      REAL(8), DIMENSION(NROWB)         :: STRAIN
      REAL(8), DIMENSION(4)             :: UNITVEC
!
      REAL(8) :: POISSON, AREA, C1
!
      UNITVEC(1)= 1.0D0
      UNITVEC(2)= 1.0D0 
      UNITVEC(3)= 0.0D0 
      UNITVEC(4)= 1.0D0
!
!.... GENERATION OF LOCAL SHAPE FUNCTIONS AND WEIGHT VALUES 
! 
      CALL SHLQ(SHLD,WD,NINTD,NEN)
!
      DO 500 NEL=1,NUMEL
!
!.... ..LOCALIZE GEOMECHANICAL REGION TO SETUP ELASTIC PARAMETERS
!
         POISSON = GEOINDIC('POISSON',GEOFORM(NEL))
!
!.... ..SETUP STOCHASTIC ELASTICITY TENSOR FOR BBAR METHOD 
! 
         CALL SETUPC(CBBAR,YOUNG(NEL),POISSON,NROWB,IOPT)
!
!...  ..LOCALIZE COORDINATES AND DIRICHLET B.C.
!
         CALL LOCAL(conecNodaisElem(1,NEL),X,XL,NEN,NSD,NESD)
         CALL LOCAL(conecNodaisElem(1,NEL),DIS,DISL,NEN,NDOFD,NED2)
!
         QUAD = .TRUE.
         IF (conecNodaisElem(3,NEL).EQ.conecNodaisElem(4,NEL)) QUAD = .FALSE.
!
         CALL SHGQ(XL,DETD,SHLD,SHGD,NINTD,NEL,QUAD,NEN)
!
         DIVU(NEL) = 0.0D0
!
!...  ..CLEAR INITIAL STRAIN
!
         STRAIN = 0.0D0 
!
!.... ..DEFINE ELEMENT AREA
!
         AREA = 0.0D0
!
!.... ..SETUP FOR AXISYMMETRIC OPTION
!
         IF (IOPT.EQ.2) THEN
            DO 150 L=1,NINTD
               R(L)    = ROWDOT(SHGD(NROWSH,1,L),XL,NROWSH,NESD,NEN)
               DETD(L) = DETD(L)*R(L) 
 150        CONTINUE
         ENDIF
!
!.... ..CALCULATE MEAN VALUES OF SHAPE FUNCTION GLOBAL DERIVATIVES
!.... ..FOR MEAN-DILATATIONAL B-BAR FORMULATION
!
         CALL XMEANSH(SHGBR,WD,DETD,R,SHGD,NEN,NINTD,IOPT,NESD,NROWSH)
!
!.... ..LOOP OVER INTEGRATIONN POINTS
!
         DO 300 L=1,NINTD
            C1 = WD(L)*DETD(L)
            AREA = AREA + C1
!
            DO 200 J=1,NEN  
!.... ..UPLOAD B-BAR MATRIX AT NODE J
               CALL SETBB(BBARJ,SHGD(1:NROWSH,J,L),SHGBR(1:NROWSH,J), &
     &                 R(L),NROWSH,NROWB,IOPT,IBBAR)
! 
!.... ..COMPUTE STRAINS WITHIN INTRINSIC B-BAR FORMULATION 
!
               DO 200 K=1,NROWB
                  STRAIN(K)=STRAIN(K)+COLDOT(BBARJ(K,1:2),DISL(1:2,J),2)*C1
 200        CONTINUE
 300     CONTINUE
!
!.... ..COMPUTE MEAN DEFORMATION OVER ELEMENT 
!
         DO 350 K=1,NROWB 
            STRAIN(K)=STRAIN(K)/AREA
 350     CONTINUE
!..
         DIVU(NEL) = COLDOT(UNITVEC,STRAIN,NROWB)
!
!.... ..COMPUTE MEAN VOLUMETRIC STRESS 
!
         DO 420 K=1,NROWB
            STRESS(K,NEL)=COLDOT(CBBAR(K,1:4),STRAIN(1:4),4)+STRSS0(K,NEL)
 420     CONTINUE
!
 500  CONTINUE
!
      RETURN 
!
 4000 FORMAT(2X,40(1PE15.8,2X)) 
 4500 FORMAT(I8,X,40(1PE15.8,2X))
 4600 FORMAT(A12,X,40(1PE15.8,2X))
!
      END SUBROUTINE
!
!**** NEW **** FOR INCOMPRESSIBILITY *************************************** 
!
      SUBROUTINE POS4CREEP(x, conecNodaisElem)
!
!.... PROGRAM TO UPDATE STRESS FOR NON-LINEAR CREEP MODEL 
! 
      use mMalha,            only: nsd, numnp, numel, nen, LOCAL
      use mMalha,            only: numelReserv, multab
      use mGlobaisEscalares, only: ndofp, ndofD, nrowsh, nnp, NCREEP
      use mAlgMatricial,     only: rowdot, coldot
      use mFuncoesDeForma,   only: shgq, shlq
      use mPropGeoFisica,    only: YOUNG, PORE, PORE0
      use mPropGeoFisica,    only: GEOFORM, GEOINDIC, FNCMECLAW
!
      IMPLICIT NONE
!
      REAL*8,  intent(in)    :: x(nsd,numnp)
      INTEGER, intent(in)    :: conecNodaisElem(nen,numel)
      CHARACTER(5)           :: MECLAW
!
      LOGICAL QUAD
!
      INTEGER :: J,K,L,NEL      
!
!.... INPUT/OUTPUT VECTORS AND MATRIZES OF SUBROUTINE
!
      REAL(8), DIMENSION(NESD,NEN)      :: XL
      REAL(8), DIMENSION(NED2,NEN)      :: DLTRL
!
!.... LOCAL VECTORS AND MATRIZES
!
      REAL(8), DIMENSION(NROWB,NESD)    :: BBARJ
      REAL(8), DIMENSION(NROWB,NROWB)   :: QIXI, CBBAR
      REAL(8), DIMENSION(NROWB)         :: DEVSTRS,TENSAO
      REAL(8), DIMENSION(4)             :: STRAIN
!
      REAL(8) :: SHLD(NROWSH,NEN,NINTD), SHGD(NROWSH,NEN,NINTD)
      REAL(8) :: SHGBR(NROWSH,NEN)
      REAL(8) :: DETD(NINTD), R(NINTD), WD(NINTD)

      REAL(8) :: POISSON, xpoisson
      REAL(8) :: GTIMES2, GTIMES3, BULKROCK, TRCSTRS
      REAL(8) :: GAMMA,FX,DERFX,DELTAG, QTRIAL, ROOT3D2 
!
      DATA ROOT3D2/1.224744871391589D0/
!
!.... GENERATION OF LOCAL SHAPE FUNCTIONS AND WEIGHT VALUES 
! 
      CALL SHLQ(SHLD,WD,NINTD,NEN)
!
      TENSAO = 0.0D0
!
      DO 500 NEL=1,NUMEL
!
         MECLAW  = FNCMECLAW(GEOFORM(NEL))
!
         POISSON = GEOINDIC('POISSON',GEOFORM(NEL))
         GTIMES2 = YOUNG(NEL)/(1.0D0+POISSON)
         GTIMES3 = 1.5D0*GTIMES2
!
         BULKROCK = YOUNG(NEL)/(3.0D0*(1.0D0-2.0D0*POISSON))    
!
!.... SETUP STOCHASTIC ELASTICITY TENSOR FOR BBAR METHOD 
! 
         CALL SETUPC(CBBAR,YOUNG(NEL),POISSON,NROWB,IOPT)
!
!..... ..LOCALIZE COORDINATES AND DIRICHLET B.C.
!
         CALL LOCAL(conecNodaisElem(1,NEL),X,XL,NEN,NSD,NESD)
         CALL LOCAL(conecNodaisElem(1,NEL),DTRL,DLTRL,NEN,NDOFD,NED2)
!
         QUAD = .TRUE.
         IF (conecNodaisElem(3,NEL).EQ.conecNodaisElem(4,NEL)) QUAD = .FALSE.
!
         CALL SHGQ(XL,DETD,SHLD,SHGD,NINTD,NEL,QUAD,NEN)
!..
!..... ..SETUP FOR AXISYMMETRIC OPTION
!..
         IF (IOPT.EQ.2) THEN
            DO 100 L=1,NINTD
               R(L)    = ROWDOT(SHGD(NROWSH,1,L),XL,NROWSH,NESD,NEN)
               DETD(L) = DETD(L)*R(L) 
 100        CONTINUE
         ENDIF
!
!.... .. CALCULATE MEAN VALUES OF SHAPE FUNCTION GLOBAL DERIVATIVES
!....     FOR MEAN-DILATATIONAL B-BAR FORMULATION
!
         CALL XMEANSH(SHGBR,WD,DETD,R,SHGD,NEN,NINTD,IOPT,NESD,NROWSH)
!..
!.... ...LOOP OVER INTEGRATION POINTS
!..
         DO 400 L=1,NINTD
!..
!..... ..CLEAR INITIAL STRAIN
!..
            STRAIN=0.0D0
!
            DO 200 J=1,NEN  
!..
!.... ..... UPLOAD B-BAR MATRIX AT NODE J
!..
               CALL SETBB(BBARJ,SHGD(1:NROWSH,J,L),SHGBR(1:NROWSH,J), &
     &                    R(L),NROWSH,NROWB,IOPT,IBBAR)
!.. 
!.... ..... COMPUTE STRAINS WITHIN INTRINSIC B-BAR FORMULATION 
!..
               DO 200 K=1,NROWB
                  STRAIN(K)=STRAIN(K)+ &
     &                    COLDOT(BBARJ(K,1:2),DLTRL(1:2,J),2)
 200        CONTINUE
!..
!.... ..... COMPUTE STRESS 
!..
!            CALL MULTAB(CBBAR,STRAIN,TENSAO,4,4,4,4,4,1,1)
            TENSAO=matmul(CBBAR,STRAIN)
!..
!.... ..... COMPUTE DEVIATOR STRESS TENSOR
!..
            CALL DEVTENSOR(DEVSTRS,TRCSTRS,TENSAO,NROWB)
!
!.... ..... COMPUTE EFFECTIVE TRIAL STRESS
!
            QTRIAL=ROOT3D2*DSQRT(DEVNORM2(DEVSTRS,NROWB))
!..
!... SOLVE NON LINEAR EQUATION FOR GAMMA FRAMEWORK: NEWTON METHOD
!..
            GAMMA=0.0D0
!
!            IF (MOD(NNP,NCREEP).EQ.0) THEN
               IF (MECLAW.EQ.'CREEP') THEN 
 220              CONTINUE
                  FX =  FNOTLIN(QTRIAL,GTIMES3,1,GAMMA)
                  DERFX = FNOTLIN(QTRIAL,GTIMES3,2,GAMMA) 
                  DELTAG = -FX/DERFX
                  GAMMA = GAMMA+DELTAG
                  IF (DABS(DELTAG).GT.1.0D-6) GOTO 220
               ENDIF
!            ENDIF
!..
!.... ...  UPDATE DEVIATOR STRESS
!..
            DO 230 K=1,NROWB
               DEVSTRS(K)=(1.0D0-GAMMA*GTIMES3/QTRIAL)*DEVSTRS(K)
 230        CONTINUE
!
            STRSS(NEL,L,1) = DEVSTRS(1)+TRCSTRS
            STRSS(NEL,L,2) = DEVSTRS(2)+TRCSTRS
            STRSS(NEL,L,3) = DEVSTRS(3)
            STRSS(NEL,L,4) = DEVSTRS(4)+TRCSTRS
!..
!.... .... UPDATE LOCAL TANGENT MATRIX
!..
            CALL SETTGMX(QIXI,DEVSTRS,GTIMES2,GTIMES3,BULKROCK, &
     &                   QTRIAL,GAMMA,MECLAW)
!..
!.... .... TRANSFER 4X4-ORDER MATRIX TO GLOBAL TANGENT ARRAY
!..
            CALL QX2TANG(QIXI,HMTTG(NEL,L,1:16))
!
 400     CONTINUE
!
 500  CONTINUE
!
      RETURN 
!
2222  FORMAT('ELEMENTO (NEL)=',I5,2X,' GAUSS POINT (L)=',I2/5X  &
     &' C11, C21, C31, C41 =',4(1PE9.2,2X)/5X,                  &
     &' C12, C22, C32, C42 =',4(1PE9.2,2X)/5X,                  &
     &' C13, C23, C33, C43 =',4(1PE9.2,2X)/5X,                  &
     &' C14, C24, C34, C44 =',4(1PE9.2,2X)//)
 4000 FORMAT(2X,40(1PE15.8,2X)) 
 5000 FORMAT(I4,2X,I1,2X,40(1PE15.8,2X)) 
!
      END SUBROUTINE
!
!**** NEW **** FOR STOCHASTIC AND NON-LINEAR FORMULATION ***************** 
!
      SUBROUTINE GEOSETUP(YOUNG,GEOFORM, &
     &                     NUMEL,NROWB,NINTD,NROWB2,IOPT)
!
      use mPropGeoFisica,    only: GEOINDIC
!
      IMPLICIT NONE
!
!.... PROGRAM TO SETUP INITIAL INELASTIC TANGENT MATRIX 
! 
      INTEGER :: NUMEL,NROWB,NINTD,NROWB2,IOPT
!
!.... NEXT LINES: GLOBAL VECTORS AND MATRIZES
!
      CHARACTER*12, DIMENSION(NUMEL)         :: GEOFORM
      REAL(8),      DIMENSION(NUMEL)         :: YOUNG
      REAL(8) :: CBBAR(NROWB, NROWB)
! 
      REAL(8) :: POISSON
      integer :: NEL, L
!
      DO 500 NEL=1,NUMEL 
         POISSON = GEOINDIC('POISSON',GEOFORM(NEL))
!
!.... SETUP STOCHASTIC ELASTICITY TENSOR FOR BBAR METHOD 
!
         CALL SETUPC(CBBAR,YOUNG(NEL),POISSON,NROWB,IOPT)
!
!.... SETUP INITIAL TANGENT MATRIX AT ELEMENT GAUSS POINT 
!
         DO 100 L=1,NINTD
            CALL QX2TANG(CBBAR,HMTTG(NEL,L,1:16))
 100     CONTINUE
 500  CONTINUE
!
      RETURN
!
 4000 FORMAT(2X,40(1PE15.8,2X)) 
!
      END SUBROUTINE
!
!**** NEW ****************************************************************
!
      FUNCTION POWERLAW(X,NDERHO)
!
      USE mPropGeoFisica, only: DTCREEP, CREEPZERO, POWERN, SIGMAREF
!
!.....PROGRAM TO COMPUTE POWER LAW LIKE FUCNTION 
!
      IMPLICIT NONE
!
      REAL*8 :: X
      INTEGER NDERHO
!
      REAL*8 :: POWERLAW, XIMAGE
!
      XIMAGE=(DTCREEP*CREEPZERO)*(X/SIGMAREF)**POWERN
!
      GOTO(100,200) NDERHO
!
!.... FUNCTION ONLY
!
 100  CONTINUE
        POWERLAW = XIMAGE
      RETURN
!..
!.... FIRST DERIVATIVE OF POWER LAW 
!..
 200  CONTINUE
        POWERLAW = POWERN*XIMAGE/X
      RETURN
!
      END FUNCTION
!
!**** NEW ****************************************************************
!
      FUNCTION FNOTLIN(QTRIAL,THREEG,NDERHO,XINPUT)
!
!
!.....PROGRAM TO COMPUTE NON LINEAR FUNCTION FROM POWER LAW
!
      IMPLICIT NONE
!
      INTEGER :: NDERHO
      REAL(8) :: QTRIAL, THREEG, X, XINPUT, FNOTLIN
!
      X = QTRIAL-THREEG*XINPUT
!
      GOTO(100,200) NDERHO
!
!... FUNCTION
!
 100  CONTINUE
        FNOTLIN = XINPUT-POWERLAW(X,NDERHO)
      RETURN
!
!... FIRST DERIVATIVE FUNCTION
!
 200  CONTINUE
        FNOTLIN = 1.0D0+POWERLAW(X,NDERHO)*THREEG
      RETURN
!
      END FUNCTION
!
!**** NEW ****************************************************************
!
      FUNCTION DEVNORM2(X,NROWB)
!
!.....PROGRAM TO COMPUTE SQUARE NORM OF DEVIATOR TENSOR (2D MODEL)
!
      IMPLICIT NONE
!
      INTEGER :: NROWB
      REAL*8  :: X(NROWB)
!
      REAL*8  :: DEVNORM2
!
      DEVNORM2 = X(1)*X(1)+X(2)*X(2)+X(4)*X(4)+2.0D0*X(3)*X(3)
!
      RETURN
!
      END FUNCTION
!      
!**** NEW ****************************************************************
!
      SUBROUTINE COMPTRACE(P,AVSTRS,SIGMAT,NROWB,NUMEL,numelReserv)
!
!.... PROGRAM TO COMPUTE TRACE OF TOTAL STRESS
!.... TOTAL STRESS = SOLID_STRSS_TRACE+FLUID_PRESSURE
!
      use mGlobaisEscalares, only: S3DIM
      use mPropGeoFisica,    only: YUNGVECT, POISVECT, GRAINBLK
      use mPropGeoFisica,    only: YOUNG, GEOFORM, BULK
!
      IMPLICIT NONE
!
      INTEGER :: NEL,NUMEL,NROWB,numelReserv
!
      REAL(8), DIMENSION(numelReserv) :: P, SIGMAT
      REAL(8), DIMENSION(NROWB,NUMEL) :: AVSTRS
!      
      REAL(8) :: BULKROCK, BIOTCOEF
!
      DO 500 NEL=1,numelReserv
         BULKROCK = BULK(YOUNG(NEL),POISVECT(1),S3DIM)
         BIOTCOEF = 1.0D0 - BULKROCK/GRAINBLK(1)
!
         IF (GEOFORM(NEL).EQ.'RESERVATORIO') THEN
         SIGMAT(NEL)=AVSTRS(1,NEL)+AVSTRS(2,NEL)+AVSTRS(4,NEL) &
    &                  -S3DIM*BIOTCOEF*P(NEL)
         ENDIF
500   CONTINUE
!
      RETURN
!
      END SUBROUTINE
!
!**** NEW ***** FOR CREEP MODELING *************************************
!
      SUBROUTINE DEVTENSOR(DEVIATOR,TRACED3,TENSORIN,NROWB)
!..
!.... PROGRAM TO COMPUTE DEVIATOR FOR PLANE DEFORMATIONS STATE  
!..
      IMPLICIT NONE
!
      INTEGER :: NROWB
      REAL*8 :: DEVIATOR(NROWB), TENSORIN(NROWB), TRACED3
!
      TRACED3 = (TENSORIN(1)+TENSORIN(2)+TENSORIN(4))/3.0D0
!
      DEVIATOR(1)= TENSORIN(1)-TRACED3
      DEVIATOR(2)= TENSORIN(2)-TRACED3
      DEVIATOR(3)= TENSORIN(3)
      DEVIATOR(4)= TENSORIN(4)-TRACED3
!
      RETURN
!
      END SUBROUTINE
!
!*** NEW *** FOR STOCHASTIC YOUNG MODULUS ******************************* 
!
      SUBROUTINE SETTGMX(QIXI,DEVSTRS,GTIMES2,GTIMES3,BULKROCK,& 
     &                   QTRIAL,GAMMA, MECLAW)
! 
!..... PROGRAM TO SETUP TANGENT MATRIX 
!
      USE mGlobaisEscalares, only: nnp, NCREEP
! 
      IMPLICIT NONE
! 
      REAL*8  :: GTIMES2, GTIMES3, BULKROCK, QTRIAL, GAMMA
      CHARACTER(5) :: MECLAW
!
      REAL*8 :: UNITVECT(4),UNITTENS(4,4),DEVPROJ(4,4),QIXI(4,4)
      REAL*8 :: DEVSTRS(4), VECTNORM, X
      REAL*8 :: BLOCOA, BLOCO2G, BLOCO6G2, BLOCOG2
      INTEGER :: I, J
!
      UNITVECT(1)=1.0D0
      UNITVECT(2)=1.0D0
      UNITVECT(3)=0.0D0
      UNITVECT(4)=1.0D0
!
      UNITTENS(1,1)=1.0D0
      UNITTENS(1,2)=0.0D0
      UNITTENS(1,3)=0.0D0
      UNITTENS(1,4)=0.0D0
      UNITTENS(2,1)=0.0D0
      UNITTENS(2,2)=1.0D0
      UNITTENS(2,3)=0.0D0
      UNITTENS(2,4)=0.0D0
      UNITTENS(3,1)=0.0D0
      UNITTENS(3,2)=0.0D0
      UNITTENS(3,3)=0.5D0
      UNITTENS(3,4)=0.0D0 
      UNITTENS(4,1)=0.0D0 
      UNITTENS(4,2)=0.0D0 
      UNITTENS(4,3)=0.0D0 
      UNITTENS(4,4)=1.0D0 
!..
!... MOUNT DEVIATOR PROJETOR
!..
      DO 20 I=1,4
         DO 10 J=1,4
            DEVPROJ(I,J)=UNITTENS(I,J)-UNITVECT(I)*UNITVECT(J)/3.0D0
 10      CONTINUE
 20   CONTINUE
!
!      WRITE(102,2222) 1,1, ((DEVPROJ(I,J),I=1,4),J=1,4)
!..
!.... COMPUTE DEVIATOR NORM
!
       VECTNORM = DEVNORM2(DEVSTRS,4)
!..
!.... COMPUTE FACTORS THAT MULTIPLY FOURTH-ORDER MATRICES
!..
       X = QTRIAL-GTIMES3*GAMMA
       BLOCOA   = POWERLAW(X,2)/(1.0D0+GTIMES3*POWERLAW(X,2))
!
       IF (MECLAW.EQ.'CREEP') THEN 
            BLOCO2G  = GTIMES2*(1.0D0-GAMMA*GTIMES3/QTRIAL)
            BLOCO6G2 = GTIMES2*GTIMES3*(GAMMA/QTRIAL-BLOCOA)/VECTNORM
          ELSE
            BLOCO6G2 = 0.0D0
            BLOCO2G  = GTIMES2
       ENDIF
!
       DO 400 I=1,4
        DO 300 J=1,4
           QIXI(I,J)=  BLOCO2G*DEVPROJ(I,J) &
     &               + BULKROCK*UNITVECT(I)*UNITVECT(J) &
     &               + BLOCO6G2*DEVSTRS(I)*DEVSTRS(J)
 300    CONTINUE
 400   CONTINUE
!
       RETURN
!
 2000 FORMAT(5(1PE15.8,2X))
!
       END SUBROUTINE
!
!*** NEW *** FOR STOCHASTIC YOUNG MODULUS ******************************* 
!
      SUBROUTINE QX2TANG(QMATR4X4,TANGENT)
! 
!     PROGRAM TO TRANSFER FROM QIXI 4X4 MATRIX TO TANGENT ARRAY
! 
      IMPLICIT NONE
! 
!     REMOVE ABOVE CARD FOR SINGLE PRECISION OPERATION 
! 
      REAL*8 :: TANGENT(16),QMATR4X4(4,4)
!
            TANGENT(1)  = QMATR4X4(1,1)
            TANGENT(2)  = QMATR4X4(1,2)
            TANGENT(3)  = QMATR4X4(1,3)
            TANGENT(4)  = QMATR4X4(1,4)
            TANGENT(5)  = QMATR4X4(2,1)
            TANGENT(6)  = QMATR4X4(2,2)
            TANGENT(7)  = QMATR4X4(2,3)
            TANGENT(8)  = QMATR4X4(2,4)
            TANGENT(9)  = QMATR4X4(3,1)
            TANGENT(10) = QMATR4X4(3,2)
            TANGENT(11) = QMATR4X4(3,3)
            TANGENT(12) = QMATR4X4(3,4)
            TANGENT(13) = QMATR4X4(4,1)
            TANGENT(14) = QMATR4X4(4,2)
            TANGENT(15) = QMATR4X4(4,3)
            TANGENT(16) = QMATR4X4(4,4)
!
       RETURN
!
 2000 FORMAT(5(1PE15.8,2X))
!
      END SUBROUTINE
!
!**** NEW *** FOR PLASTICITY ********************************************* 
!
      SUBROUTINE UPDATEINCR(IDDIS,BRHSD,NDOFD,NUMNP)
! 
!.... UPDATE TRIAL DISPLACEMENT ARRAY 
! 
      IMPLICIT NONE
! 
!.... REMOVE ABOVE CARD FOR SINGLE-PRECISION OPERATION 
! 
      INTEGER :: NDOFD,NUMNP
      INTEGER :: IDDIS(NDOFD,*)
      REAL*8  :: BRHSD(*)
!
      INTEGER :: I, J, K
! 
!      write(32,*) 'ENTRANDO update iter =',nwtniter
!      call xxwrite(dtrl,ndof2,numnp)
!      call xxwrite(brhsd,1,10000)

      DO 200 I=1,NDOFD
! 
        DO 100 J=1,NUMNP
           K = IDDIS(I,J)
           IF (K.GT.0) then
           DTRL(I,J) = DTRL(I,J) + BRHSD(K)
           ENDIF
  100   CONTINUE
! 
  200 CONTINUE
!
!      write(32,*) 'SAINDO update iter =',nwtniter
!      call xxwrite(dtrl,ndof2,numnp)
!      call xxwrite(brhsd,1,10000)
      RETURN
!
      END SUBROUTINE
!
!*** NEW *** FOR STOCHASTIC YOUNG MODULUS ******************************* 
!
      SUBROUTINE TANG2QX(TANGENT,QMATR4X4)
! 
!     PROGRAM TO TRANSFER FROM TANGENT ARRAY TO QIXI 4X4 MATRIX 
!
      IMPLICIT NONE
! 
!     REMOVE ABOVE CARD FOR SINGLE PRECISION OPERATION 
! 
      real*8 :: TANGENT(16),QMATR4X4(4,4)
!
         QMATR4X4(1,1) = TANGENT(1)
         QMATR4X4(1,2) = TANGENT(2)
         QMATR4X4(1,3) = TANGENT(3)
         QMATR4X4(1,4) = TANGENT(4)
         QMATR4X4(2,1) = TANGENT(5)
         QMATR4X4(2,2) = TANGENT(6)
         QMATR4X4(2,3) = TANGENT(7)
         QMATR4X4(2,4) = TANGENT(8)
         QMATR4X4(3,1) = TANGENT(9)
         QMATR4X4(3,2) = TANGENT(10)
         QMATR4X4(3,3) = TANGENT(11)
         QMATR4X4(3,4) = TANGENT(12)
         QMATR4X4(4,1) = TANGENT(13)
         QMATR4X4(4,2) = TANGENT(14)
         QMATR4X4(4,3) = TANGENT(15)
         QMATR4X4(4,4) = TANGENT(16)
!
       RETURN
!
 2000 FORMAT(5(1PE15.8,2X))
!
      END SUBROUTINE
!
!**** NEW **** FOR INCOMPRESSIBILITY *************************************** 
!
      SUBROUTINE PRINT_DX(DIS ,PORE  ,P     ,S     ,MASCN ,STRSS0, &
     &                    VC  ,AVSTRS,NDOF2 ,NUMEL ,NROWB ,NUMNP)
!
      use mMalha,            only: XC
      USE mGlobaisEscalares, only: NUMDX, NNP, PATHDX
      use mLeituraEscrita,   only: PRINT_DXINFO, IFEDX
      use mPropGeoFisica,    only: GEOFORM, GEOINDIC
!
!.... PROGRAM TO PRINT DATA FROM GEOMECHANIC MODEL
! 
      IMPLICIT NONE
!  
!.... REMOVE ABOVE CARD FOR SINGLE-PRECISION OPERATION  
!
      CHARACTER*30 NIDISP,NIPRSR,NIPORE,NICREP,NISATR,NIVELT
      CHARACTER*30 NISIGX,NISIGY,NISGTA,NISIGZ
      CHARACTER*30 NIS2S1,NIMASC
!
      CHARACTER*3 ASTEP
!
!.... INPUT/OUTPUT VECTORS AND MATRIZES OF SUBROUTINE
!
      INTEGER :: I,J,K,NEL,NDOF2,NUMEL,NROWB,NUMNP
      INTEGER :: IDISP, IPRSR, IPORE, ICREP, ISATR, IVELT
      INTEGER :: IMASC, ISIGX, ISIGY, ISGTA, ISIGZ, IS2S1
!
      REAL(8), DIMENSION(NDOF2,NUMNP) :: DIS
      REAL(8), DIMENSION(NDOF2,*)     :: VC
      REAL(8), DIMENSION(NUMEL)       :: PORE,P,S,MASCN
      REAL(8), DIMENSION(NROWB,NUMEL) :: AVSTRS 
      REAL(8), DIMENSION(NROWB,NUMEL) :: STRSS0 
!
      REAL(8), DIMENSION(NROWB)       :: TENSAO, DEVSTRS
!
      REAL(8) :: ROOT3D2,TRCSTRS,QTRIAL,STRIAL,SPORE,SBLNC
      REAL(8) :: PHIDRO, ROOT, S1, S2, S3, TAU, SIGMAN
!
      DATA ROOT3D2/1.224744871391589D0/
!
!... OUTPUT DATA FILES 
!
      IF (NUMDX.EQ.0) RETURN
!
      IDISP = 631
      IPRSR = 632
      IPORE = 633
      ICREP = 634
      ISATR = 635
      IVELT = 636
      IMASC = 637
      ISIGX = 638
      ISIGY = 639
      ISGTA = 640
      ISIGZ = 641
      IS2S1 = 642
!
!.... SETUP FILES COUNTER
!
!      WRITE(ASTEP,'(I3.3)') IDINT(TPRT_PHI/DTPRT_PHI)
      WRITE(ASTEP,'(I3.3)') NNP/NUMDX
!
!.... OUT-PUT FILES FOR OPEN-DX VISUALIZATION
!.... DATA VALUES AT NODAL FEM POINTS
!
      NIDISP = PATHDX//'/disp'//ASTEP//'.stoc'
      NIPRSR = PATHDX//'/prsr'//ASTEP//'.stoc'
      NIPORE = PATHDX//'/pore'//ASTEP//'.stoc'
      NICREP = PATHDX//'/crep'//ASTEP//'.stoc'
      NISATR = PATHDX//'/satr'//ASTEP//'.stoc'
      NIVELT = PATHDX//'/velt'//ASTEP//'.stoc'
      NIMASC = PATHDX//'/masc'//ASTEP//'.stoc'
      NISIGX = PATHDX//'/sigx'//ASTEP//'.stoc'
      NISIGY = PATHDX//'/sigy'//ASTEP//'.stoc'
      NISGTA = PATHDX//'/sigt'//ASTEP//'.stoc'
      NISIGZ = PATHDX//'/sigz'//ASTEP//'.stoc'
      NIS2S1 = PATHDX//'/s2s1'//ASTEP//'.stoc'
!
!.....OPEN NODAL DATA FILES
!
      OPEN(UNIT=IDISP, FILE= NIDISP)
      OPEN(UNIT=IPRSR, FILE= NIPRSR)
      OPEN(UNIT=IPORE, FILE= NIPORE)
      OPEN(UNIT=ICREP, FILE= NICREP)
      OPEN(UNIT=ISATR, FILE= NISATR)
      OPEN(UNIT=IVELT, FILE= NIVELT)
      OPEN(UNIT=IMASC, FILE= NIMASC)
      OPEN(UNIT=ISIGX, FILE= NISIGX)
      OPEN(UNIT=ISIGY, FILE= NISIGY)
      OPEN(UNIT=ISGTA, FILE= NISGTA)
      OPEN(UNIT=ISIGZ, FILE= NISIGZ)
      OPEN(UNIT=IS2S1, FILE= NIS2S1)
!
!.... PRINT DISPLACEMENTS
!
       DO 30 I=1,NUMNP
          WRITE(IDISP,4000) (DIS(J,I),J=1,NDOF2) 
 30    CONTINUE
!
      DO 500 NEL=1,NUMEL 
!
!.... PRINT PRESSURE
!
      IF (GEOFORM(NEL).EQ.'RESERVATORIO') SPORE = P(NEL)
!      IF (GEOFORM(NEL).NE.'RESERVATORIO') SPORE = -10.0D0
      IF (GEOFORM(NEL).NE.'RESERVATORIO') SPORE = -GEOPRSR(NEL)
!      WRITE(IPRSR,4000) XC(1,NEL),XC(2,NEL),SPORE
      WRITE(IPRSR,4000) SPORE
!
!.... PRINT GEOMECHANIC POROSITY: "PORE" 
! 
      IF (GEOFORM(NEL).EQ.'RESERVATORIO') SPORE = PORE(NEL)
      IF (GEOFORM(NEL).NE.'RESERVATORIO') SPORE = -10.0D0
      WRITE(IPORE,4000) SPORE
!
!.... PRINT STRESS CREEP VALUE:
! 
      DO 100 K=1,NROWB
!         TENSAO(KK) = STRSS(NEL,1,KK)
         TENSAO(K) = AVSTRS(K,NEL)  ! +STRSS0(K,NEL)
 100  CONTINUE
!
!      IF (GEOFORM(NEL).EQ.'RESERVATORIO') THEN
!            QTRIAL = -10.0D0
!         ELSE
!.... ... .. COMPUTE DEVIATOR STRESS TENSOR
      CALL DEVTENSOR(DEVSTRS,TRCSTRS,TENSAO,NROWB)
!.... ... .. COMPUTE EFFECTIVE TRIAL STRESS
      QTRIAL=ROOT3D2*DSQRT(DEVNORM2(DEVSTRS,NROWB))
!      ENDIF
!
      WRITE(ICREP,4000) QTRIAL
!
!.... PRINT SATURATION
!
      IF (GEOFORM(NEL).EQ.'RESERVATORIO') STRIAL =  S(NEL)
      IF (GEOFORM(NEL).NE.'RESERVATORIO') STRIAL = -10.0D0
      WRITE(ISATR,3000) STRIAL
!
!.... PRINT MASS CONTENT 
!
      IF (GEOFORM(NEL).EQ.'RESERVATORIO') STRIAL = MASCN(NEL)
      IF (GEOFORM(NEL).NE.'RESERVATORIO') STRIAL = -10.0D0
      WRITE(IMASC,4000) STRIAL
!
!.... PRINT TOTAL VELOCITY 
!
      IF (GEOFORM(NEL).EQ.'RESERVATORIO') THEN 
             STRIAL = VC(1,NEL)
             SPORE  = VC(2,NEL)
         ELSE
             STRIAL = 0.0D0
             SPORE  = 0.0D0
      ENDIF
!       WRITE(IVELT,4000) (VC(J,NEL),J=1,NDOF2)
      WRITE(IVELT,4000)  STRIAL, SPORE
!
!.... PRINT STRESS_X COMPONENT
! 
      WRITE(ISIGX,4000) TENSAO(1)
!
!.... PRINT STRESS_Y COMPONENT
!
      WRITE(ISIGY,4000) TENSAO(2)
!
!.... PRINT STRESS_XY COMPONENT 
!
      WRITE(ISGTA,4000) TENSAO(3)
!      WRITE(ISGTA,4000) (TENSAO(1)+tensao(2))*GEOINDIC('POISSON',GEOFORM(NEL)), tensao(4)
!.... COMPUTE PRINCIPAL STRESS
!
      PHIDRO = 0.5D0*(TENSAO(1)+TENSAO(2))
      ROOT=DSQRT((0.5D0*(TENSAO(1)-TENSAO(2)))**2+TENSAO(3)**2)
      S1=PHIDRO+ROOT
      S2=PHIDRO-ROOT
!
!      WRITE(ISGTA,4000) S1
!
!      S3=GEOINDIC('POISSON',GEOFORM(NEL))*(s1+s2)
!
!      TAU = (S1-S2)**2+(S2-S3)**2+(S3-S1)**2
!      TAU = DSQRT(TAU)/3.0D0
!      TAU = ROOT
!      SIGMAN = DABS(S1+S2+S3)/3.0D0
!      SIGMAN = DABS(PHIDRO)

!      WRITE(ISGTA,4000) s3- tensao(4)
!.... PRINT STRESS_Z COMPONENT 
!
      WRITE(ISIGZ,4000) TENSAO(4)
!
!.... DILATANCY EXPERIMENT
!
!      IF (GEOFORM(NEL).EQ.'RESERVATORIO') THEN 
!              STRIAL = 0.0D0
!           ELSE
!              STRIAL = TAU-0.8996*SIGMAN+0.01697*SIGMAN**2
!      ENDIF
!
!      WRITE(ISIGZ,4000) STRIAL
!
!.... PRINT YOUNG MODULUS 
!
!      IF (NNP.EQ.0) S1=1.0D0

      WRITE(IS2S1,4000) s2  !/s1 ! tensao(2)-spore
!
 500  CONTINUE
! 
!       IF (NNP.EQ.NVEL) THEN
!       PI=4.0D0*DATAN(1.0D0)
!       DO 600 NEL=51,NUMEL,100
!          radio=DSQRT(XC(1,NEL)**2+XC(2,NEL)**2)
!          angle=DATAN(XC(2,NEL)/XC(1,NEL))
! 
!          sgxxcos2 = (dcos(angle))**2*STRSS(NEL,1,1)
!          sgyysin2 = (dsin(angle))**2*STRSS(NEL,1,2)
!          sgxysin2 = (dsin(2.0d0*angle))*STRSS(NEL,1,3)
!          sgxxcos2 = (dcos(angle))**2*AVSTRS(1,NEL)
!          sgyysin2 = (dsin(angle))**2*AVSTRS(2,NEL)
!          sgxysin2 = (dsin(2.0d0*angle))*AVSTRS(3,NEL)
! 
!          sigmar = sgxxcos2+sgyysin2+sgxysin2 
!
!         WRITE(101,4000) RADIO,ANGLE*180.0D0/PI,SIGMAR
! 600  CONTINUE
!      ENDIF
      CLOSE(IDISP)
      CLOSE(IPRSR)
      CLOSE(IPORE)
      CLOSE(ICREP)
      CLOSE(ISATR)
      CLOSE(IMASC)
      CLOSE(ISIGX)
      CLOSE(ISIGY)
      CLOSE(ISIGZ)
      CLOSE(IS2S1)
      CLOSE(IVELT)
!
!.... PRINT INFORMATION ON NODAL DX FILE "nodestoc.dx"
!
      CALL PRINT_DXINFO('WRITE_FEDX_DATA',IFEDX,NUMNP,NUMEL)
!
      RETURN 
!
 2000 FORMAT(I5,2X,40(1PE15.8,2X)) 
 3000 FORMAT(2X,5(F25.15,2x)) 
 4000 FORMAT(2X,40(1PE15.8,2X)) 
!
      END SUBROUTINE
!
!**** NEW **** FORCES APPLIED ON IRREGULAR MESH ******************************
!
      SUBROUTINE InSeaLoad(FDIS,NDOF2,NUMNP,NLVECT,XTERLOAD)
!
!**** PROGRAM TO MULTIPLY VERTICAL NODAL FORCES ON TOP POST-SALT BY SEALOAD 
!
      IMPLICIT NONE
!
      INTEGER  :: NODE, NDOF2, NUMNP, NLVECT
      REAL(8), DIMENSION(NDOF2,NUMNP,NLVECT) :: FDIS
      REAL(8) :: XTERLOAD
!
      DO 120 NODE=1,NUMNP
!         IF (FDIS(1,NODE,1).NE.0.0D0) FDIS(1,NODE,1) = XTERLOAD*FDIS(1,NODE,1)
         IF (FDIS(2,NODE,1).NE.0.0D0) FDIS(2,NODE,1) = XTERLOAD*FDIS(2,NODE,1)
!         IF (FDIS(3,NODE,1).NE.0.0D0) FDIS(3,NODE,1) = XTERLOAD*FDIS(3,NODE,1)
 120  CONTINUE
!
      RETURN
!
      END SUBROUTINE  
!      
!**** NEW ****************************************************************
!
      SUBROUTINE PRSRINIT(GEOPRSR,P,XC,NUMEL,numelReserv)
!
!.... PROGRAM TO SETUP INITIAL HIDROSTATIC PRESSURE 
! 
      use mGlobaisEscalares, only: PATHDX, TypeProcess, SOLIDONLY
      use mGlobaisArranjos,  only: GRAV
      use mMalha,            only: NEN, NSD, LEFTLINE, RGHTLINE 
      use mPropGeoFisica,    only: GEOFORM, FNCMECLAW, RHOW, RHOO
      use mPropGeoFisica,    only: PWELL, nelXReserv, XTERLOAD
!
      IMPLICIT NONE
!
      INTEGER ::  NEL, NUMEL, NUMELRESERV
!
      REAL(8), DIMENSION(NUMEL)         :: GEOPRSR
      REAL(8), DIMENSION(1,numelReserv) :: P
      REAL*8,  DIMENSION(NSD,NUMEL)     :: XC
      REAL*8  :: THREEA, BSKEMPT, UNDRAINU, PLOAD
      
!
!
!test lines for translation use ymark(1) plus ymark(4) for plane
!      nel = numnp-16*76+1 
!      DO 145 I=nel,nel+75
!         xref = DOMESURF(X(1,I),YMARK(1),IDome) ! YMARK(1)  ! 
!         WRITE(2025,1000) X(1,I),X(2,I),XREF !+YMARK(1)
! 145  CONTINUE
!
!      NISTSV = PATHDX//'/stress0.dat'
!      OPEN(UNIT=ISTSV, FILE= NISTSV)
!

      IF (TRIM(TypeProcess).EQ.'TERZAGHI') THEN 
         DO 100 NEL=1,numel
            GEOPRSR(NEL) = 0.0D0 
            IF (GEOFORM(NEL).EQ.'RESERVATORIO') P(1,NEL) = XTERLOAD
100      CONTINUE
         RETURN
      ENDIF
!
      IF (TRIM(TypeProcess).EQ.'MANDEL') THEN 
         THREEA   = 3.0D0*(RGHTLINE-LEFTLINE)
         BSKEMPT  = 1.0D0
         UNDRAINU = 0.5D0
         PLOAD    = XTERLOAD*BSKEMPT*(1.0D0+UNDRAINU)/THREEA
         DO 110 NEL=1,NUMEL
            GEOPRSR(NEL) = 0.0D0 
            IF (GEOFORM(NEL).EQ.'RESERVATORIO') P(1,NEL) = PLOAD 
110      CONTINUE
         RETURN
      ENDIF
!
      IF (SOLIDONLY) THEN
         DO 120 NEL=1,NUMEL
            GEOPRSR(NEL) = 0.0D0
            IF (GEOFORM(NEL).EQ.'RESERVATORIO') P(1,NEL) = GEOPRSR(NEL)
120      CONTINUE       
         RETURN
      ENDIF
!
      PLOAD = 0.0D0
      IF (GRAV(2).EQ.0.0D0) PLOAD = XTERLOAD
      DO 200 NEL=1,NUMEL
         GEOPRSR(NEL) = SUMPRSR(XC(1,NEL),XC(2,NEL),GRAV(2),'LINEAR',NEL)
         GEOPRSR(NEL) = GEOPRSR(NEL) + PLOAD
!         GEOPRSR(NEL) = SUMPRSR(XC(1,NEL),XC(2,NEL),GRAV(2),'INCOMP',NEL)
!         GEOPRSR(NEL) = SUMPRSR(XC(1,NEL),XC(2,NEL),GRAV(2),'EXACTL',NEL)
         IF (GEOFORM(NEL).EQ.'RESERVATORIO') P(1,NEL) = GEOPRSR(NEL)

200   CONTINUE
!
!... projection of reservoir pressure to boundary condition for test 
! 
      DO 300 NEL=nelXReserv,numelReserv,nelXReserv
         PWELL(NEL/nelXReserv) = P(1,NEL) 
300   CONTINUE
!
!      CLOSE(ISTSV)
!
      RETURN
 1000 FORMAT(2X,40(1PE15.8,2X)) 
 4500 FORMAT(I8,X,40(1PE15.8,2X))
!
      END SUBROUTINE
!
!**** NEW **** FOR INITIAL STRESS ***************************************** 
!
      FUNCTION SUMPRSR(XC, YC, G, TASK,ELEMNT)
!
      use mMalha,         only: LEFTLINE, RGHTLINE
      use mPropGeoFisica, only: MEANDENS, MEANBULK
!
!.... PROGRAM TO COMPUTE HIDROSTATIC FLUID PRESSURE LOAD 
!....    TASK = 'INCOMP'  ---> INCOMPRESSIBLE FLUID MODEL
!....    TASK = 'EXACTL'  ---> COMPRESSIBLE FLUID MODEL
!....    TASK = 'LINEAR'  ---> LINEAR APPROXIMATION MODEL 
!....                          FOR COMPRESSIBLE FLUID
!
      IMPLICIT NONE
!
      INTEGER      :: ELEMNT
      REAL(8)      :: XC, YC, SUMPRSR, DENSRSRV, BULKRSRV, G
      CHARACTER(6) :: TASK
!
      IF (XC.LT.LEFTLINE) THEN 
         DENSRSRV = MEANDENS(5)
         BULKRSRV = MEANBULK(5)
      ENDIF
!
      IF ((XC.GE.LEFTLINE).AND.(XC.LE.RGHTLINE)) THEN
         DENSRSRV = MEANDENS(1)
         BULKRSRV = MEANBULK(1) 
      ENDIF
!
      IF (XC.GT.RGHTLINE) THEN 
         DENSRSRV = MEANDENS(3)
         BULKRSRV = MEANBULK(3)
      ENDIF
!
      SUMPRSR = DABS(ADDPRSR(XC,YC,DENSRSRV,BULKRSRV,G,TASK,ELEMNT))
!
      RETURN
!
      END FUNCTION
!
!**** NEW **** FOR INITIAL STRESS ***************************************** 
!
      FUNCTION ADDPRSR(XC,YC,DENSRSRV,BULKRSRV,G,TASK,ELEMNT)
!
      use mGlobaisEscalares, only: SALTCREEP, S3DIM
      use mMalha,          only: POSTLINE, SALTLINE, RTOPLINE
      use mMalha,          only: RBTTLINE, RIFTLINE, IDOME
      use mPropGeoFisica,  only: RHOW, BULKWATER, SEADEPTH, GRAINBLK
      use mPropGeoFisica,  only: DOMESURF, PORELAGR, YOUNG, POISVECT
      use mPropGeoFisica,  only: MEANDENS, MEANBULK, RHODVECT, BULK
      use mPropGeoFisica,  only: MEANRHOW, MEANBLKW
!
      use mGlobaisEscalares, only: PRESSAOREF, COTAREF
!
!.... FUNCTION TO ADD DENSITIES ON VERTICAL DIRECTION  
!
      IMPLICIT NONE
!
      INTEGER :: ELEMNT
      REAL(8) :: XC, YC, ADDPRSR, DENSRSRV, BULKRSRV, WEIGHT
      REAL(8) :: SOMAPARC, SOMA, YMIN, YMAX, G, RHOEQ, BULKEQ
      REAL(8) :: BULKROCK, BIOTCOEF, POROSITY, ONEMPORE
      REAL(8), DIMENSION(4) :: YMARK
      CHARACTER(6) :: TASK
!
      SOMA     = 0.0D0
      SOMAPARC = 0.0D0
      RHOEQ    = MEANRHOW(6)
      BULKEQ   = 1.0D0/MEANBLKW(6)
      SOMAPARC = DABS(HIDROSTAT(SEADEPTH,COTAREF,RHOEQ,BULKEQ,G,TASK))
      SOMAPARC = SOMAPARC + PRESSAOREF
      WEIGHT   = 0.0D0
!
!.... SET VERTICAL COORDINATES LIMITS FOR POST-SALT REGION 
!
      YMIN     = POSTLINE
      YMAX     = DOMESURF(XC,SALTLINE,IDome) ! YMARK(1)  ! 
      RHOEQ    = MEANDENS(6)
      BULKEQ   = MEANBULK(6)
!
!.... NEXT LINES: LOAD INTO POS-SALT REGION
!
      IF ((YC.LT.YMIN).AND.(YC.GT.YMAX)) THEN 
         SOMA    = DABS(HIDROSTAT(DABS(YC),DABS(YMIN),RHOEQ,BULKEQ, &
     &                       G,TASK))
         ADDPRSR = SOMA + SOMAPARC
         RETURN
      ENDIF
!
!.... NEXT LINE: COMPUTE PARTIAL LOAD FOR POS_SALT REGION
!
      SOMAPARC = DABS(HIDROSTAT(DABS(YMAX),DABS(YMIN),RHOEQ,BULKEQ, & 
     &                     G,TASK))+SOMAPARC
!
      IF ((PORELAGR(4).EQ.0.0D0).OR.(SALTCREEP)) THEN 
         RHOEQ   = RHODVECT(6)
         BULKEQ  = 1.0D0
         POROSITY= PORELAGR(6)
         ONEMPORE= 1.0D0-POROSITY
         WEIGHT  = POROSITY*SOMAPARC + &
     &             ONEMPORE*DABS(HIDROSTAT(DABS(YMAX),DABS(YMIN),RHOEQ, &
     &             BULKEQ, G,'INCOMP')) + WEIGHT 
      ENDIF
!
!.... RE-SET VERTICAL COORDINATES LIMITS FOR SALT REGION 
!
      YMIN   = YMAX
      YMAX   = RTOPLINE
      RHOEQ  = MEANDENS(4) 
      BULKEQ = MEANBULK(4)
!
!.... NEXT LINES: LOAD INTO SALT REGION
!
      IF ((YC.LT.YMIN).AND.(YC.GT.YMAX)) THEN 
         SOMA    = DABS(HIDROSTAT(DABS(YC),DABS(YMIN),RHOEQ,BULKEQ, &
     &                       G,TASK))
         ADDPRSR = SOMA + SOMAPARC
         IF ((PORELAGR(4).EQ.0.0D0).OR.(SALTCREEP)) ADDPRSR = 0.0D0 
         RETURN
      ENDIF

!
!.... NEXT LINE: COMPUTE PARTIAL LOAD FOR SALT+POS_SALT REGION
!
      SOMAPARC = DABS(HIDROSTAT(DABS(YMAX),DABS(YMIN),RHOEQ,BULKEQ, & 
     &                     G,TASK)) + SOMAPARC
!
      IF ((PORELAGR(4).EQ.0.0D0).OR.(SALTCREEP)) THEN 
         RHOEQ   = RHODVECT(4)
         BULKEQ  = 1.0D0
         POROSITY= PORELAGR(4)
         ONEMPORE= 1.0D0-POROSITY
         WEIGHT  = POROSITY*SOMAPARC + &
     &             ONEMPORE*DABS(HIDROSTAT(DABS(YMAX),DABS(YMIN),RHOEQ, &
     &             BULKEQ, G,'INCOMP')) + WEIGHT 
!
!.... SETUP BIOT COEFICIENT ALSO CALLED BIOT'S ALPHA
!
         BULKROCK  = BULK(YOUNG(ELEMNT), POISVECT(1), S3DIM) 
         BIOTCOEF  = 1.0D0 - BULKROCK/GRAINBLK(1)
!
        SOMAPARC = BIOTCOEF*WEIGHT
! 
      ENDIF
!
!.... RE-SET VERTICAL COORDINATES LIMITS FOR RESERVOIR REGION 
!
      YMIN   = YMAX
      YMAX   = RBTTLINE
      POROSITY = 1.0D0 
      IF ((PORELAGR(4).EQ.0.0D0).OR.(SALTCREEP)) POROSITY = PORELAGR(1) 
      ONEMPORE =  1.0D0-POROSITY
      RHOEQ  = POROSITY*DENSRSRV + ONEMPORE*RHODVECT(1)
      BULKEQ = POROSITY*BULKRSRV + ONEMPORE/GRAINBLK(1)
!
!.... NEXT LINES: LOAD INTO RESERVOIR REGION
!
      IF ((YC.LT.YMIN).AND.(YC.GT.YMAX)) THEN
         SOMA    = DABS(HIDROSTAT(DABS(YC),DABS(YMIN),RHOEQ,BULKEQ, & 
     &             G, TASK))
         ADDPRSR = SOMA + SOMAPARC
         RETURN
      ENDIF
!
!.... NEXT LINE: COMPUTE PARTIAL LOAD FOR: 
!                RESERVOIR+SALT+POS_SALT REGION
!
      SOMAPARC = DABS(HIDROSTAT(DABS(YMAX),DABS(YMIN),RHOEQ,BULKEQ, & 
     &                     G,TASK)) + SOMAPARC
!
!.... RE-SET VERTICAL COORDINATES LIMITS FOR RIFT REGION 
!
      YMIN   = YMAX
      YMAX   = RIFTLINE
      RHOEQ  = MEANDENS(2) 
      BULKEQ = MEANBULK(2)
!
!.... NEXT LINES: LOAD INTO RIFT REGION
!
      ADDPRSR = SOMAPARC
!
      IF ((YC.LT.YMIN).AND.(YC.GT.YMAX)) THEN 
          SOMA    = DABS(HIDROSTAT(DABS(YC),DABS(YMIN),RHOEQ,BULKEQ, &
     &                        G,TASK))
          ADDPRSR = SOMA + SOMAPARC
          RETURN
      ENDIF
!
      RETURN
!
      END FUNCTION
!
!**** NEW **** FOR INITIAL STRESS ***************************************** 
!
      FUNCTION HIDROSTAT(Z1, Z2, RHOEQ, BULKEQ, G, TASK)
!
!.... PROGRAM TO COMPUTE HIDROSTATIC PRESSURE FIELD 
! 
      IMPLICIT NONE
!
      REAL(8) :: Z1, Z2, HIDROSTAT, RHOEQ, BULKEQ, G
      CHARACTER(6) TASK
!                  123456
      IF (TASK.EQ.'INCOMP') THEN
          HIDROSTAT = PRSRINCOMP(Z1,Z2,RHOEQ,G)
          RETURN
      ENDIF
!
      IF (TASK.EQ.'EXACTL') THEN
          HIDROSTAT = PRSREXACT(Z1,Z2,RHOEQ,BULKEQ,G)
          RETURN
      ENDIF
!
      IF (TASK.EQ.'LINEAR') THEN
          HIDROSTAT = PRSRLINEAR(Z1,Z2,RHOEQ,BULKEQ,G)
          RETURN
      ENDIF
!
      RETURN
!
      END FUNCTION
!
!**** NEW **** FOR INITIAL STRESS ***************************************** 
!
      FUNCTION PRSRINCOMP(Z1, Z2, RHOEQ, G)
!
!.... PROGRAM TO COMPUTE PRESSURE OF INCOMPRESSIBLE FLUID
! 
      IMPLICIT NONE
!
      REAL(8) :: Z1, Z2, PRSRINCOMP, RHOEQ, G
!
      PRSRINCOMP = RHOEQ*G*(Z1-Z2)
!
      RETURN
!
      END FUNCTION
!
!**** NEW **** FOR INITIAL STRESS ***************************************** 
!
      FUNCTION PRSREXACT(Z1, Z2, RHOEQ, BULKEQ, G)
!
!.... PROGRAM TO COMPUTE PRESSURE OF COMPRESSIBLE FLUID: EXACT FORMULATION
! 
      IMPLICIT NONE
!
      REAL(8) :: Z1, Z2, PRSREXACT, RHOEQ, BULKEQ, G
!
      PRSREXACT = PRSRLOG(Z1,RHOEQ,BULKEQ,G)-PRSRLOG(Z2,RHOEQ,BULKEQ,G)
!
      RETURN
!
      END FUNCTION
!
!
!**** NEW **** FOR INITIAL STRESS ***************************************** 
!
      FUNCTION PRSRLOG(Z, RHOEQ, BULKEQ, G)
!
!.... PROGRAM TO COMPUTE PRESSURE OF COMPRESSIBLE FLUID: EXACT FORMULATION
! 
      IMPLICIT NONE
!
      REAL(8) :: Z, PRSRLOG, RHOEQ, BULKEQ, G
!
      PRSRLOG = -(DLOG(1.0D0-RHOEQ*G*BULKEQ*Z))/BULKEQ
!
      RETURN
!
      END FUNCTION
!
!**** NEW **** FOR INITIAL STRESS ***************************************** 
!
      FUNCTION PRSRLINEAR(Z1, Z2, RHOEQ, BULKEQ, G)
!
!.... PROGRAM TO COMPUTE PRESSURE OF COMPRESSIBLE FLUID: LINEAR APROXIMATION
! 
      IMPLICIT NONE
!
      REAL(8) :: Z1, Z2, PRSRLINEAR, RHOEQ, BULKEQ, G
!
      PRSRLINEAR = PRSREXP(Z1,RHOEQ,BULKEQ,G)-PRSREXP(Z2,RHOEQ,BULKEQ,G)
!
      RETURN
!
      END FUNCTION
!
!**** NEW **** FOR INITIAL STRESS ***************************************** 
!
      FUNCTION PRSREXP(Z, RHOEQ, BULKEQ, G)
!
!.... PROGRAM TO COMPUTE PRESSURE OF COMPRESSIBLE FLUID: LINEAR APROXIMATION
! 
      IMPLICIT NONE
!
      REAL(8) :: Z, PRSREXP, RHOEQ, BULKEQ, G
!
      PRSREXP = (DEXP((RHOEQ*BULKEQ)*(G*Z))-1.0D0)/BULKEQ
!
      RETURN
!
      END FUNCTION
!
!**** NEW ************************************************************** 
!
      SUBROUTINE VECTOR_SOURC(satElem, p, GEOPRSR, brhsd, lmD)
!
      use mGlobaisEscalares, only: ndofD, nrowsh, S3DIM, optSolver
      use mGlobaisArranjos,  only: grav
      use mAlgMatricial,     only: rowdot, coldot,  neqD
      use mSolverGaussSkyline, only: addrhs
      use mFuncoesDeForma,   only: shgq,shlq
      use mMalha,            only: local, multab
      use mMalha,            only: X, conecNodaisElem,xc
      use mMalha,            only: numnp, numel, nen, nsd, numelReserv
      use mPropGeoFisica,    only: YOUNG, GRAINBLK, RHODVECT, PORELAGR
      use mPropGeoFisica,    only: RHOW, RHOO, BULKWATER, BULKOIL
      use mPropGeoFisica,    only: GEOFORM, GEOINDIC, BULK
      use mPropGeoFisica,    only: PORE
!
      use mGlobaisEscalares, only: PRESSAOREF, COTAREF
      use mSolverHypre
!
!.... PROGRAM TO COMPUTE PLASTIC DEFORMATION OVER THE MACRO DOMAIN 
! 
      IMPLICIT NONE
!  
!.... REMOVE ABOVE CARD FOR SINGLE-PRECISION OPERATION  
!
      LOGICAL DIAG,QUAD,ZERODL
!
!.... INPUT/OUTPUT VECTORS AND MATRIZES OF SUBROUTINE
!
      REAL*8,  intent(in)    :: P(1,numelReserv)
      REAL*8,  intent(in)    :: satElem(numelReserv)
      REAL*8,  intent(in)    :: GEOPRSR(NUMEL)
      REAL(8), intent(inout) :: BRHSD(neqD)
      INTEGER, intent(in)    :: LMD(NED2,NEN,NUMEL)
!
      real(8) :: xl(nesd,nen), disl(ned2,nen)
      REAL*8 :: ELEFFMD(NEE2,NEE2),ELRESFD(NEE2)
!
!.... LOCAL VECTORS AND MATRIZES
!
      REAL(8) :: BBARI(NROWB,NESD), PRESSURE(NROWB)
      REAL(8) :: SHGD(NROWSH,NEN,NINTD), SHLD(NROWSH,NEN,NINTD)
      REAL(8) :: SHGBR(NROWSH,NEN)
      REAL(8) :: DETD(NINTD), R(NINTD), WD(NINTD)
!
      REAL(8) :: C1, BULKROCK, BIOTCOEF, POROSITY, ONEMPORE
      REAL(8) :: RHOMAT, DENSOLID, DENFLUID, SATURAT, ONEMSAT
      REAL(8) :: BULKGRAIN, BULKFLUD, POISSON, RHOWLIN, RHOOLIN
      INTEGER :: NEL, I, L, II, JJ
!
!.... GENERATION OF LOCAL SHAPE FUNCTIONS AND WEIGHT VALUES 
! 
      CALL SHLQ(SHLD,WD,NINTD,NEN)
!
!... CONSISTENT MATRIX
!
      DIAG = .FALSE.
!
      DO 500 NEL=1,NUMEL
!
!.... SETUP DENSITY ELEMENT
!
         POROSITY = GEOINDIC('POROSTY',GEOFORM(NEL))
         DENSOLID = GEOINDIC('ROKDENS',GEOFORM(NEL))
         DENFLUID = GEOINDIC('FLUDENS',GEOFORM(NEL))
         BULKFLUD = GEOINDIC('BLKFLUD',GEOFORM(NEL))
         DENFLUID = DENFLUID*(1.0D0+BULKFLUD*(GEOPRSR(NEL)-PRESSAOREF))
!
         IF (GEOFORM(NEL).EQ.'RESERVATORIO') THEN
            RHOWLIN  = RHOW*(1.0D0+(P(1,NEL)-PRESSAOREF)/BULKWATER)
            RHOOLIN  = RHOO*(1.0D0+(P(1,NEL)-PRESSAOREF)/BULKOIL)
            SATURAT  = satElem(NEL)
            ONEMSAT  = 1.0D0-satElem(NEL)
            BULKFLUD = SATURAT/BULKWATER + ONEMSAT/BULKOIL
            DENFLUID = SATURAT*RHOWLIN   + ONEMSAT*RHOOLIN
            DENFLUID = DENFLUID*(1.0D0+BULKFLUD*(P(1,NEL)-PRESSAOREF))
            POROSITY = PORE(NEL)
         ENDIF
!
         RHOMAT = POROSITY*DENFLUID + (1.0D0-POROSITY)*DENSOLID
!
!.... SETUP BIOT COEFICIENT ALSO CALLED BIOT'S ALPHA
!
         POISSON   = GEOINDIC('POISSON',GEOFORM(NEL))
         BULKGRAIN = GEOINDIC('BLKGRIN',GEOFORM(NEL))
         BULKROCK  = BULK(YOUNG(NEL), POISSON, S3DIM) 
         BIOTCOEF  = 1.0D0 - BULKROCK/BULKGRAIN
!
!.... SETUP ELEMENT PRESSURE VECTOR
!
         PRESSURE(1) = BIOTCOEF*GEOPRSR(NEL)
         PRESSURE(2) = PRESSURE(1)
         PRESSURE(3) = 0.0D0
         PRESSURE(4) = PRESSURE(1)
!
!LINE TO TEST EQUIVALENT GRAVITY:                  PRESSURE    = 0.0D0
!LINE TO TEST EQUIVALENT GRAVITY:         PRESSURE(2) = GEOPRSR(NEL)
!
         IF (GEOFORM(NEL).EQ.'RESERVATORIO') THEN
            PRESSURE(1) = BIOTCOEF*P(1,NEL)
            PRESSURE(2) = PRESSURE(1)
            PRESSURE(4) = PRESSURE(1)
!
!LINE TO TEST EQUIVALENT GRAVITY:                 PRESSURE    = 0.0D0
!LINE TO TEST EQUIVALENT GRAVITY:        PRESSURE(2) = P(1,NEL)/BIOTCOEF
!
         ENDIF
!
!         pressure(2) = (xc(2,numel)-xc(2,nel))*grav(2)*rhomat
!         pressure(2) = pressure(2)+0.5d0*grav(2)*rhomat
!         pressure(1) = pressure(2)
!         pressure(3) = 0.0d0
!         pressure(4) = pressure(2)
!         pressure = 0.0d0
!
!         write(2025,4500) nel,xc(2,nel),pressure(2) 
!RHOMAT, grav(2) !DENFLUID, densolid ! rhomat, poisson, bulkgrain
!
!... CLEAR STIFFNESS MATRIX AND FORCE ARRAY
!
         ELEFFMD = 0.0D0
         ELRESFD = 0.0D0
!
!... LOCALIZE COORDINATES AND DIRICHLET B.C.
!
         CALL LOCAL(conecNodaisElem(1,NEL),X,XL,NEN,NSD,NESD)
         CALL LOCAL(conecNodaisElem(1,NEL),DIS,DISL,NEN,NDOFD,NED2)
!
         QUAD = .TRUE.
         IF (conecNodaisElem(3,NEL).EQ.conecNodaisElem(4,NEL)) QUAD = .FALSE.
!
         CALL SHGQ(XL,DETD,SHLD,SHGD,NINTD,NEL,QUAD,NEN)
!
!.... SETUP FOR AXISYMMETRIC OPTION
!
         IF (IOPT.EQ.2) THEN
            DO 100 L=1,NINTD
               R(L)   = ROWDOT(SHGD(NROWSH,1,L),XL,NROWSH,NESD,NEN)
               DETD(L) = DETD(L)*R(L)
 100        CONTINUE
         ENDIF
!
!.... FORM STIFFNESS MATRIX 
!
!.... .. CALCULATE MEAN VALUES OF SHAPE FUNCTION GLOBAL DERIVATIVES
!.... .. FOR MEAN-DILATATIONAL B-BAR FORMULATION
!
         CALL XMEANSH(SHGBR,WD,DETD,R,SHGD,NEN,NINTD,IOPT,NESD,NROWSH)
!
!.... .. LOOP OVER INTEGRATIONN POINTS
!
         DO 400 L=1,NINTD
!
            C1 = DETD(L)*WD(L)
!
!**** **** MOUNT FORCE VECTOR ****************** 
!
            CALL CLEAR(BBARI,NROWB*NESD)
!
            DO 300 I=1,NEN
!
            CALL SETBB(BBARI,SHGD(1:NROWSH,I,L),&
     &              SHGBR(1:NROWSH,I),R(L),NROWSH,NROWB,IOPT,IBBAR)
!
                 ELRESFD(NED2*I-1)= ELRESFD(NED2*I-1) &
     &                      + COLDOT(BBARI(1:4,1),PRESSURE,4)*C1
!
                 ELRESFD(NED2*I)  = ELRESFD(NED2*I) &
     &                      + COLDOT(BBARI(1:4,2),PRESSURE,4)*C1 &
     &                      + RHOMAT*GRAV(2)*SHGD(3,I,L)*C1
!
  300       CONTINUE
  400    CONTINUE
!
         DO 450 II=1,NEE2
            DO 450 JJ=1,NEE2
               ELEFFMD(II,JJ) = AUXM(NEL,II,JJ)
 450     CONTINUE
!
!...     COMPUTATION OF DIRICHLET B.C. CONTRIBUTION
!   
         CALL ZTEST(DISL,NEE2,ZERODL)
!
         IF (.NOT.ZERODL) &
     &      CALL KDBCGEO(ELEFFMD,ELRESFD,DISL,NEE2,LMD(1,1,NEL),NEL)
!
!.... ASSEMBLE ELEMENT FORCE ARRAY INTO GLOBAL RIGHT-HAND SIDE VECTOR
!
         CALL ADDRHS(BRHSD,ELRESFD,LMD(1,1,NEL),NEE2)
!
 500  CONTINUE
!
!      do 660 i=1,neqD
!          write(3030,*) i, brhsd(i)
!660   continue

      if (optSolver=='hypre')   then
         call adicionarValoresVetor_HYPRE(b_HYPRE_G, 1, neqD, rows_G, BRHSD)
   !      call fecharMatriz_HYPRE            (A_HYPRE_G, parcsr_A_G)  ! creep
         call fecharVetor_HYPRE             (b_HYPRE_G, par_b_G   )  ! creep
   !      call fecharVetor_HYPRE             (u_HYPRE_G, par_u_G   )  ! creep
      endif

      RETURN
 4500 FORMAT(I7,X,40(1PE15.8,2X))
!
      END SUBROUTINE
!
!**** NEW *************************************************************** 
!
      SUBROUTINE POS4_MASSCNT(DIVU,DIVU0, P, P0, &
     &         PORE,PORE0,YOUNG,MASCN,MASCN0,PHIEULER,NUMEL,NUMELRESERV)
!
!.... PROGRAM TO COMPUTE POROSITY FROM GEOMECHANIC MODEL
! 
      use mLeituraEscrita,   only: CODERROR
      use mGlobaisEscalares, only: S3DIM
      use mPropGeoFisica,    only: POISVECT, GRAINBLK, BULKWATER
      use mPropGeoFisica,    only: GEOFORM, BULK
!
      IMPLICIT NONE
!
      LOGICAL QUAD
!
      INTEGER :: I,J,K,L,NEL, NUMEL, numelReserv
!
!.... INPUT/OUTPUT VECTORS AND MATRIZES OF SUBROUTINE
!
      REAL*8,  INTENT(in)             :: P(1,numelReserv), P0(numelReserv)
      REAL(8), DIMENSION(numelReserv) :: PORE, PORE0, MASCN, MASCN0, PHIEULER
      REAL(8), DIMENSION(NUMEL)       :: DIVU, DIVU0, YOUNG
!
!.... ..LOCAL VECTORS AND MATRIZES
!
      REAL(8) :: POISSON, DIFFDIVU, DIFFPRES
!
      REAL(8) :: BULKROCK, BIOTCOEF
      REAL(8) :: BULKFLUID, DEFNM1, DEFMM1, JACOBIAN
!
      DO 500 NEL=1,NUMELRESERV
!
         IF (GEOFORM(NEL).EQ.'RESERVATORIO') THEN
!
           BULKROCK  = BULK(YOUNG(NEL),POISVECT(1),S3DIM)
           BIOTCOEF  = 1.0D0 - BULKROCK/GRAINBLK(1)   !also called ALPHA
!
!... FIRST DEFINE 1/N from Coussy 2.Edt   EQ. 4.35 
!
           DEFNM1    = (BIOTCOEF-PORE0(NEL))/GRAINBLK(1)
!
!... SECOND DEFINE 1/M FROM Coussy 2.Edt  EQ. 4.61
!
           DEFMM1 = DEFNM1 + PORE0(NEL)/BULKWATER 
!
!... DEFINITION OF JACOBIAN OF INFINITESSIMAL TRANSFORMATION COUSSY EQ 1.27 
!
           JACOBIAN  = 1.0D0 + DIVU(NEL)
!
           DIFFDIVU  = DIVU(NEL)-DIVU0(NEL)
           DIFFPRES  = P(1,NEL)-P0(NEL)
!
! ..COMPUTE LAGRANGIAN POROSITY
!           PORE(NEL) =  1.0D0-(1.0D0-PORE0(NEL))*DEXP(-diffdivu) 
!
!new line: linearized See: Coussy 2.Edt. Eqs. 4.19
           PORE(NEL) = PORE0(NEL)+BIOTCOEF*DIFFDIVU+DEFNM1*DIFFPRES

!
!...COMPUTE EULERAIN POROSITY
           PHIEULER(NEL) = PORE(NEL)/JACOBIAN
!
!...COMPUTE MASS CONTENT: linearized form Ref: Coussy 2.Edt. Eqs. 4.62
!
           MASCN(NEL)= MASCN0(NEL)+(BIOTCOEF*DIFFDIVU+DEFMM1*DIFFPRES)

           IF ((MASCN(NEL).LT.0.0D0).OR.(PORE(NEL).LT.0.0D0)) THEN
               write(2020,*) 'nel=',nel
               write(2020,*) 'alfa=',BIOTCOEF
               write(2020,*) 'divu=',divu(nel)
               write(2020,*) 'divu0=',divu0(nel)
               write(2020,*) 'biotmod=',defmm1
               write(2020,*) 'diffpress=',diffpres
               write(2020,*) '1/m(p-p0)=',defmm1*DIFFPRES
               write(2020,*) 'masscont=',mascn(nel), mascn0(nel) 
!                               123456789+12345678
               CALL CODERROR(4,'ALSO SEE fort.2020')
           ENDIF
!
         ENDIF
!
 500  CONTINUE
!
      RETURN 
!
 4000 FORMAT(2X,40(1PE15.8,2X)) 
!
      END SUBROUTINE
!
!**** NEW **** FOR VISCOELASTICITY *************************************** 
!
      SUBROUTINE POS4STRS(x, conecNodaisElem, STRESS, DIVU)
! 
      use mLeituraEscrita,   only: CODERROR
      use mMalha,            only: nen, LOCAL
      use mMalha,            only: nsd, numnp, numel, numelReserv
      use mGlobaisEscalares, only: ndofp, ndofD, nrowsh
      use mAlgMatricial,     only: rowdot, coldot
      use mFuncoesDeForma,   only: shgq, shlq
      use mPropGeoFisica,    only: YOUNG, PORE, PORE0, GRAINBLK
      use mPropGeoFisica,    only: GEOFORM, GEOINDIC
!
!.... PROGRAM TO COMPUTE SOLID ELEMENT: STRESS AND VOLUMETRIC DEFORMATION 
!....  ......  .....    COMPUTED FROM CREEP MODEL
!
      IMPLICIT NONE
!
      INTEGER, intent(in)             :: conecNodaisElem(NEN,NUMEL)
      REAL(8), intent(in)             :: X(NSD,NUMNP)
      REAL(8), DIMENSION(NROWB,NUMEL) :: STRESS
      REAL(8), DIMENSION(NUMEL)       :: DIVU
!
      LOGICAL QUAD
      REAL(8), DIMENSION(NESD,NEN)         :: XL
      REAL(8), DIMENSION(NED2,NEN)         :: DISL
      REAL(8), DIMENSION(NINTD)            :: WD, DETD, R
      REAL(8), DIMENSION(NROWSH,NEN,NINTD) :: SHLD,SHGD
      REAL(8), DIMENSION(NROWSH,NEN)       :: SHGBR
      REAL(8), DIMENSION(NROWB,NROWB)      :: CBBAR
!
!.... LOCAL VECTORS AND MATRIZES
!
      REAL(8), DIMENSION(NROWB,NESD)  :: BBARJ
      REAL(8), DIMENSION(NROWB)       :: STRAIN, DEVSTRS, TENSAO
      REAL(8), DIMENSION(4)           :: UNITVEC
      REAL(8), DIMENSION(NUMEL)       :: DIVULOC
!
      REAL*8  :: POISSON, AREA, C1, ROOT3D2, TRCSTRS, QVM 
      INTEGER :: NEL, L, J, K
!
      UNITVEC(1) = 1.0D0
      UNITVEC(2) = 1.0D0 
      UNITVEC(3) = 0.0D0 
      UNITVEC(4) = 1.0D0  
      ROOT3D2    = 1.224744871391589D0
!
!.... GENERATION OF LOCAL SHAPE FUNCTIONS AND WEIGHT VALUES 
! 
      CALL SHLQ(SHLD,WD,NINTD,NEN)
!
      DO 500 NEL=1,NUMEL
!
         POISSON = GEOINDIC('POISSON',GEOFORM(NEL))
!
!.... ..SETUP STOCHASTIC ELASTICITY TENSOR FOR BBAR METHOD 
! 
         CALL SETUPC(CBBAR, YOUNG(NEL),POISSON,NROWB,IOPT)
!
!...  ..LOCALIZE COORDINATES AND DIRICHLET B.C.
!
         CALL LOCAL(conecNodaisElem(1,NEL),X,XL,NEN,NSD,NESD)
         CALL LOCAL(conecNodaisElem(1,NEL),DIS,DISL,NEN,NDOFD,NED2)
!
         QUAD = .TRUE.
         IF (conecNodaisElem(3,NEL).EQ.conecNodaisElem(4,NEL)) QUAD = .FALSE.
!
         CALL SHGQ(XL,DETD,SHLD,SHGD,NINTD,NEL,QUAD,NEN)
!
         DIVULOC(NEL) = 0.0D0
!..
!.... ...CLEAR INITIAL STRAIN
!..
         CALL CLEAR(STRAIN,NROWB) 
         CALL CLEAR(TENSAO,NROWB) 
!..
!.... ...DEFINE ELEMENT AREA
!..
         AREA = 0.0D0
!..
!.... ...SETUP FOR AXISYMMETRIC OPTION
!..
         IF (IOPT.EQ.2) THEN
           DO 10 L=1,NINTD
             R(L)    = ROWDOT(SHGD(NROWSH,1,L),XL,NROWSH,NESD,NEN)
             DETD(L) = DETD(L)*R(L) 
 10        CONTINUE
         ENDIF
!
!.... ..CALCULATE MEAN VALUES OF SHAPE FUNCTION GLOBAL DERIVATIVES
!.... ..FOR MEAN-DILATATIONAL B-BAR FORMULATION
!
        CALL XMEANSH(SHGBR,WD,DETD,R,SHGD,NEN,NINTD,IOPT,NESD,NROWSH)
!
!.... ..LOOP OVER INTEGRATIONN POINTS
!
        DO 300 L=1,NINTD
!
        C1=WD(L)*DETD(L)
        AREA = AREA + C1
!
        DO 200 J=1,NEN  
!
!.... ..UPLOAD B-BAR MATRIX AT NODE J
!
        CALL SETBB(BBARJ,SHGD(1:NROWSH,J,L),SHGBR(1:NROWSH,J), &
     &             R(L),NROWSH,NROWB,IOPT,IBBAR)
! 
!.... ..COMPUTE STRAINS WITHIN INTRINSIC B-BAR FORMULATION 
!
        DO 200 K=1,NROWB
           STRAIN(K)=STRAIN(K)+                                &
     &                COLDOT(BBARJ(K,1:2),DISL(1:2,J),2)*C1
 200    CONTINUE
!..
!.... ..COMPUTE STRESS WITHIN INTRINSIC B-BAR FORMULATION 
!..
!.... ..TO COMPUTE MEAN VOLUMETRIC CREEP 
!.... ....FIRST: COMPUTE DEVIATOR STRESS TENSOR
!..
            CALL DEVTENSOR(DEVSTRS,TRCSTRS,STRSS(NEL,L,1:4),NROWB)
!
!.... ....SECOND: COMPUTE VON MISES EFFECTIVE STRESS
!
            QVM=ROOT3D2*DSQRT(DEVNORM2(DEVSTRS,NROWB))
!
        DO 220 K=1,NROWB
           ECREEP(NEL,L,K)=ECREEP(NEL,L,K)+                    &
     &                     ROOT3D2*POWERLAW(QVM,1)*DEVSTRS(K)/QVM
 220    CONTINUE
!
        DO 250 K=1,NROWB
           TENSAO(K) = TENSAO(K)+C1*STRSS(NEL,L,K)
 250    CONTINUE
!
 300    CONTINUE
!
        DO 350 K=1,NROWB       
           STRAIN(K)=STRAIN(K)/AREA
           STRESS(K,NEL)=TENSAO(K)/AREA! - STRSS0(K,NEL)
 350    CONTINUE
!
!.... ..COMPUTE MEAN DEFORMATION AND MEAN STRESS OVER ELEMENT 
!
        DIVULOC(NEL) = COLDOT(UNITVEC,STRAIN,NROWB)
!
 500  CONTINUE
!..
      DO 700 NEL=1,NUMEL
         DIVU(NEL) = DIVULOC(NEL)
700   CONTINUE
!
      RETURN 
!
 4000 FORMAT(2X,40(1PE15.8,2X)) 
!
      END SUBROUTINE
!
!**** NEW **** FOR INCOMPRESSIBILITY *************************************** 
!
      subroutine POS4TRNS2(x, conecNodaisElem, p, p0)
! 
      use mLeituraEscrita,   only: CODERROR
      use mMalha,            only: nen, LOCAL
      use mMalha,            only: nsd, numnp, numel, numelReserv
      use mGlobaisEscalares, only: ndofp, ndofD, nrowsh, S3DIM
      use mAlgMatricial,     only: rowdot, coldot
      use mFuncoesDeForma,   only: shgq, shlq
      use mPropGeoFisica,    only: YOUNG, PORE, PORE0, PHIEULER
      use mPropGeoFisica,    only: GRAINBLK, BULKWATER, MASCN, MASCN0
      use mPropGeoFisica,    only: POISVECT, GEOFORM, GEOINDIC, BULK
!
!.... PROGRAM TO COMPUTE POROSITY FROM GEOMECHANIC MODEL
!
      IMPLICIT NONE
!
      real*8,  intent(in)    :: x(nsd,numnp), p(1,numelReserv), p0(numelReserv)
      integer, intent(in)    :: conecNodaisElem(nen,numel)
!  
!.... REMOVE ABOVE CARD FOR SINGLE-PRECISION OPERATION  
!
      LOGICAL QUAD
      real(8) :: xl(nesd,nen), disl(ned2,nen), AREA, xdivu, c1
      REAL*8  :: POISSON
      integer :: nel, l, j, k
!
!.... ..LOCAL VECTORS AND MATRIZES
!
      REAL(8) :: DIVULOC(NUMEL)
      REAL(8) :: UNITVEC(NROWB),BBARJ(NROWB,NESD)
      REAL(8) :: CBBAR(NROWB, NROWB)
      REAL(8) :: SHLD(NROWSH,NEN,NINTD), SHGD(NROWSH,NEN,NINTD)
      REAL(8) :: SHGBR(NROWSH,NEN)
      REAL(8) :: DETD(NINTD), R(NINTD),WD(NINTD)
!
      REAL*8, DIMENSION(NROWB) :: TENSAO
!
      REAL*8 :: STRAIN(NROWB),DEVSTRS(NROWB)
      REAL*8  :: ROOT3D2, TRCSTRS, QVM, JACOBIAN
      REAL(8) :: BIOTMOD, BULKROCK, BIOTCOEF
      REAL(8) :: DEFNM1, DEFMM1, DIFFDIVU, DIFFPRES

      UNITVEC(1)= 1.0D0
      UNITVEC(2)= 1.0D0 
      UNITVEC(3)= 0.0D0 
      UNITVEC(4)= 1.0D0  
      ROOT3D2=1.224744871391589D0
!
!.... GENERATION OF LOCAL SHAPE FUNCTIONS AND WEIGHT VALUES 
! 
      CALL SHLQ(SHLD,WD,NINTD,NEN)
!
      DO 500 NEL=1,NUMEL
!
         POISSON = GEOINDIC('POISSON',GEOFORM(NEL))
!
!.... ..SETUP STOCHASTIC ELASTICITY TENSOR FOR BBAR METHOD 
! 
        CALL SETUPC(CBBAR, YOUNG(NEL),POISSON,NROWB,IOPT)
!
!...  ..LOCALIZE COORDINATES AND DIRICHLET B.C.
!
        CALL LOCAL(conecNodaisElem(1,NEL),X,XL,NEN,NSD,NESD)
        CALL LOCAL(conecNodaisElem(1,NEL),DIS,DISL,NEN,NDOFD,NED2)
!
        QUAD = .TRUE.
        IF (conecNodaisElem(3,NEL).EQ.conecNodaisElem(4,NEL)) QUAD = .FALSE.
!
        CALL SHGQ(XL,DETD,SHLD,SHGD,NINTD,NEL,QUAD,NEN)
!
        DIVU(NEL) = 0.0D0
!..
!.... ..CLEAR INITIAL STRAIN
!..
        CALL CLEAR(STRAIN,NROWB) 
        CALL CLEAR(TENSAO,NROWB) 
!..
!.... ..DEFINE ELEMENT AREA
!..
        AREA = 0.0D0
!..
!.... ..SETUP FOR AXISYMMETRIC OPTION
!..
        IF (IOPT.EQ.2) THEN
           DO 10 L=1,NINTD
              R(L)    = ROWDOT(SHGD(NROWSH,1,L),XL,NROWSH,NESD,NEN)
              DETD(L) = DETD(L)*R(L) 
 10        CONTINUE
        ENDIF
!
!.... ..CALCULATE MEAN VALUES OF SHAPE FUNCTION GLOBAL DERIVATIVES
!.... ..FOR MEAN-DILATATIONAL B-BAR FORMULATION
!
        CALL XMEANSH(SHGBR,WD,DETD,R,SHGD,NEN,NINTD,IOPT,NESD,NROWSH)
!
!.... ..LOOP OVER INTEGRATIONN POINTS
!
        DO 300 L=1,NINTD
!
           C1 = WD(L)*DETD(L)
           AREA = AREA + C1
           DO 200 J=1,NEN  
!
!.... ..UPLOAD B-BAR MATRIX AT NODE J
!
              CALL SETBB(BBARJ,SHGD(1:NROWSH,J,L),SHGBR(1:NROWSH,J), &
     &                R(L),NROWSH,NROWB,IOPT,IBBAR)
! 
!.... ..COMPUTE STRAINS WITHIN INTRINSIC B-BAR FORMULATION 
!
              DO 200 K=1,NROWB
              STRAIN(K)=STRAIN(K)+                                &
     &                  COLDOT(BBARJ(K,1:2),DISL(1:2,J),2)*C1
 200       CONTINUE
!..
!.... ..COMPUTE STRESS WITHIN INTRINSIC B-BAR FORMULATION 
!..
!.... ..TO COMPUTE MEAN VOLUMETRIC CREEP 
!.... ....FIRST: COMPUTE DEVIATOR STRESS TENSOR
!..
           CALL DEVTENSOR(DEVSTRS,TRCSTRS,STRSS(NEL,L,1:4),NROWB)
!
!.... ....SECOND: COMPUTE VON MISES EFFECTIVE STRESS
!
           QVM=ROOT3D2*DSQRT(DEVNORM2(DEVSTRS,NROWB))
!
           DO 220 K=1,NROWB
              ECREEP(NEL,L,K)=ECREEP(NEL,L,K)+                    &
     &                     ROOT3D2*POWERLAW(QVM,1)*DEVSTRS(K)/QVM
 220       CONTINUE
!
           DO 250 K=1,NROWB
              TENSAO(K) = TENSAO(K)+C1*STRSS(NEL,L,K)
 250       CONTINUE
!
 300    CONTINUE
!
        DO 350 K=1,NROWB       
           STRAIN(K)=STRAIN(K)/AREA
           AVSTRS(K,NEL)=TENSAO(K)/AREA
 350    CONTINUE
!
!.... ..COMPUTE MEAN DEFORMATION AND MEAN STRESS OVER ELEMENT 
        DIVU(NEL) = COLDOT(UNITVEC,STRAIN,NROWB)

!        write(888,*) divuloc(nel)
 500  CONTINUE
!..
!      DO 700 NEL=1,NUMEL
!         DIVU(NEL)=DIVULOC(NEL)
!700   CONTINUE

      DO 800 NEL=1,NUMELRESERV
!.... ..COMPUTE MEAN VOLUMETRIC DEFORMATION FOR RESERVOIR
!..
         BULKROCK  = BULK(YOUNG(NEL),POISVECT(1),S3DIM)
         BIOTCOEF  = 1.0D0 - BULKROCK/GRAINBLK(1)   !also called ALPHA
!
!... FIRST DEFINE 1/N from Coussy 2.Edt   EQ. 4.35 
!
         DEFNM1    = (BIOTCOEF-PORE0(NEL))/GRAINBLK(1)
!
!... SECOND DEFINE 1/M FROM Coussy 2.Edt  EQ. 4.61
!
         DEFMM1 = DEFNM1 + PORE0(NEL)/BULKWATER 
!
!... DEFINITION OF JACOBIAN OF INFINITESSIMAL TRANSFORMATION COUSSY EQ 1.27 
!
         JACOBIAN  = 1.0D0 + DIVU(NEL)
         DIFFDIVU  = DIVU(NEL)-DIVU0(NEL)
         DIFFPRES  = P(1,NEL)-P0(NEL)
!
! ..COMPUTE LAGRANGIAN POROSITY
!           PORE(NEL) =  1.0D0-(1.0D0-PORE0(NEL))*DEXP(-diffdivu) 
!
!new line: linearized See: Coussy 2.Edt. Eqs. 4.19
         PORE(NEL) = PORE0(NEL)+BIOTCOEF*DIFFDIVU+DEFNM1*DIFFPRES
!
!...COMPUTE EULERAIN POROSITY
         PHIEULER(NEL) = PORE(NEL)/JACOBIAN
!
!...COMPUTE MASS CONTENT: linearized form Ref: Coussy 2.Edt. Eqs. 4.62
!
         MASCN(NEL)= MASCN0(NEL)+(BIOTCOEF*DIFFDIVU+DEFMM1*DIFFPRES)

           IF ((MASCN(NEL).LT.0.0D0).OR.(PORE(NEL).LT.0.0D0)) THEN
               write(2020,*) 'nel=',nel
               write(2020,*) 'alfa=',BIOTCOEF
               write(2020,*) 'divu=',divu(nel)
               write(2020,*) 'divu0=',divu0(nel)
               write(2020,*) 'biotmod=',defmm1
               write(2020,*) 'diffpress=',diffpres
               write(2020,*) '1/m(p-p0)=',defmm1*DIFFPRES
               write(2020,*) 'masscont=',mascn(nel), mascn0(nel) 
!                               123456789+12345678
               CALL CODERROR(4,'ALSO SEE fort.2020')
           ENDIF
!..
800   CONTINUE
!
!      DO 550 NEL=1,NELX1*NELY1
!         DO 550 JJ=1,NROWB
!            STRSS(JJ,NEL)=0.0D0
! 550  CONTINUE
!
      RETURN 
!
 4000 FORMAT(2X,40(1PE15.8,2X)) 
!
      END SUBROUTINE
!
!**** NEW ************************************************************** 
!
      SUBROUTINE FTODMANDL(F,NDOF,NUMNP,NLVECT,XTERLOAD,TEMPO)
! 
!.... PROGRAM TO SETUP DISPLACMENTS FOR MANDEL PROBLEM 
!
      IMPLICIT NONE
!
!.... REMOVE ABOCE CARD FOR SINGLE-PRECISION OPERATION 
! 
      INTEGER :: I,NDOF, NUMNP, NLVECT
      REAL(8) :: XTERLOAD,TEMPO
      REAL(8), DIMENSION(NDOF,NUMNP,NLVECT)  :: F
!
      DO 200 I=NUMNP-100,NUMNP
         F(2,I,1) = SOLTEORM2(XTERLOAD,TEMPO) 
 200  CONTINUE
!
      RETURN 
!
      END SUBROUTINE
!
!**** NEW ***** FOR MANDEL PROBLEM **********************************
!
      FUNCTION SOLTEORM2(XTERLOAD,TEMPO)
!
      use mMalha,         only: RTOPLINE, RBTTLINE
      use mMalha,         only: LEFTLINE, RGHTLINE
      use mPropGeoFisica, only: POISVECT, YUNGVECT, PERMINICIAL
!
!.....PROGRAM TO COMPUTE EFFECTIVE VON MISES ENERGY 
!
      IMPLICIT NONE
!
      INTEGER :: I, IINPUT
      REAL(8) :: PI, A, B, EE, XNU, GSHEAR, UXNU, XKAPPA, XTERLOAD
      REAL(8) :: GAMMAW, C, SKEM, FORCE, FSA, TWOFSA, SERIE
      REAL(8) :: COEFUY, CONSTUY, COEFPRES, COEFSIGM, TEMPO
      REAL(8) :: BETA, TSTAR, SINCOS, BMSINCOS, UY, SOLTEORM2
      CHARACTER*30 NIINPUT
!
!.... NUMBER FOR INPUT FILE
!
      IINPUT  = 741
!
      NIINPUT = 'beta4mandel.dat'
!
      OPEN(UNIT=IINPUT, FILE=NIINPUT, STATUS='OLD')
!
      PI=4.0D0*DATAN(1.0D0)
!
!.... GEOMETRIC DATA: 
!
!....  HORIZONTAL LENGTH -> A=100 M
!....  VERTICAL LENGTH   -> B= 10 M 
!
      A = RGHTLINE - LEFTLINE
!
      B = RTOPLINE - RBTTLINE
!
!.... SOLID PARAMETERS
!
!.... YOUNG MODULUS: EE=1.0E8 [Pa]
!  
      EE  = YUNGVECT(1)
!
!..... DRAINED POISSON RATIO:
!
      XNU = POISVECT(1)
!
!.... SHEAR MODULUS: GSHEAR=EE/[2(1+XNU)] [Pa]
!
      GSHEAR=EE/(2.0D0*(1.0D0+XNU))
!
!.... UNDRAINED POISSON RATIO:
!
      UXNU = 0.5D0
!
!.... FLUID PARAMETERS
!
!.... PERMEABILITY XKAPPA: 
!       1.0 [mD {miliDarcy}]=1.0D-12 [m2 {square meters}]
!      XKAPPA = 1.0D-10
!
      XKAPPA = PERMINICIAL
!
!.... VOLUMETRIC WEIGHT OF WATER GAMMAW=1.0D0  [kN/m3]
!
      GAMMAW = 1.0D0
!
!.... COEFICIENTE DE CONSOLIDACCAO (VERRUIJT) C [m2/s]
!
!       CV=KAPPA*EE/GAMMAW
!
!.... C GENERALIZED CONSOLIDATION CONSTANT 
!
      C = 2.0D0*XKAPPA*GSHEAR*(1.0D0-XNU)
!
      C = C/(GAMMAW*(1.0D0-2.0D0*XNU))
!
!.... SKEMPTON
!
      SKEM = 1.0D0
!
!.... FORCE 
!
      FORCE = XTERLOAD
!
!.... CONSTANT OF SOLUTIONS
!
       FSA = FORCE/A
!
       TWOFSA = 2.0D0*FSA
!
       COEFUY = FSA*(1.0D0-UXNU)/GSHEAR
!
       CONSTUY = -0.5D0*FSA *(1.0D0-XNU)/GSHEAR
!
       COEFPRES = TWOFSA*SKEM*(1.0D0+UXNU)/(3.0D0)
!
       COEFSIGM = TWOFSA*(UXNU-XNU)/(1.0D0-XNU) 
!
       SERIE    = 0.0D0
!
       DO 50 I=1,20
!
         READ(IINPUT,4000) BETA
!
         TSTAR = -BETA*BETA*C*TEMPO/(A*A)
!
         SINCOS = DSIN(BETA)*DCOS(BETA)
!
         BMSINCOS = BETA-SINCOS
!
         SERIE = SERIE + SINCOS*DEXP(TSTAR)/BMSINCOS
!
 50   CONTINUE
!
      REWIND(IINPUT)
!
      UY = (CONSTUY + COEFUY*SERIE)*B
!
 100  CONTINUE
!
      CLOSE(IINPUT)
!
      SOLTEORM2 = UY
!
      RETURN  
!
 4000 FORMAT(2X,40(1PE15.8,2X)) 
!
      END FUNCTION
!
!***** %%% ***** %%% ***** %%% ***** %%% ***** %%% ***** %%% ***** 


     END MODULE MGEOMECANICA
!
!*********************** ************** **************************
