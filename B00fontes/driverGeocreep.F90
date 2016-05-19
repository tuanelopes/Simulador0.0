!  
!         programa de elementos finitos em fortran 90 
!         baseado em: The Finite Element Method, Hughes, T. J. R., (2003)
!
!         Eduardo Garcia e Tuane Lopes
!         bidu@lncc.br, tuane@lncc.br
!
!         LNCC/MCT
!         Petropolis, 07.2013
! 
!     ************************************************************
!     *                                                          *
!     *                                                          *
!     *         A LINEAR STATIC FINITE ELEMENT PROGRAM FOR       *
!     *                                                          *
!     *                 GALERKIN METHOD                          *
!     *                                                          *
!     *                                                          *
!     ************************************************************
!
  program reservoirSimulator
!
      use mGlobaisEscalares, only: exec, SALTCREEP
      use mLeituraEscrita,   only: fecharArquivosBase, abrirArquivosInformacoesMalha
      use mLeituraEscritaSimHidroGeoMec,   only: fecharArquivosSimHidroGeoMec
      use mGeomecanica, only: SOLIDONLY
      use mSolverHypre,      only: inicializarMPI, finalizarMPI
      use mSolverHypre,      only: myid, num_procs, mpi_comm
!
      implicit none
!
      CHARACTER*128    :: FLAG 
      real*8 :: t1, t2, t3, t4
      
      call timing(t3)
      
!-----------------------------------------------------------------------
#ifdef withHYPRE
      print*, "inicializando MPI"
      call inicializarMPI                 (myid, num_procs, mpi_comm)
#endif
!
!.... initialization phase
!
      call abrirArquivosInformacoesMalha ()
      print*, ""
      print*, "PREPROCESSAMENTO" 
      call timing(t1)
      call preprocessador_DS()
      call timing(t2)
#ifdef mostrarTempos
      write(*,*) "Tempo de preprocessamento=", t2-t1 
#endif
!
!.... solution phase
!
      IF (EXEC==1) THEN
         print*, ""
         print*, "Iniciando o PROCESSAMENTO..."
         CALL STRESS_INIT()
      
         IF (.NOT.SOLIDONLY) THEN 
            IF (SALTCREEP) THEN 
                 CALL processamento_creep()
              ELSE                 
                 CALL processamento_elast()
            ENDIF
         ENDIF 
      ENDIF
!
      call fecharArquivosBase()
      call fecharArquivosSimHidroGeoMec()
!
#ifdef withHYPRE
      call finalizarMPI()
#endif
      call timing(t4)
      write(*,*) "TEMPO DE PAREDE TOTAL DE MEDIDO=", t4 - t3, " segundos " 
!
end program reservoirSimulator
!
!**** NEW ********************************************************************
!
      subroutine preprocessador_DS()
!
      use mGlobaisArranjos,  only: title, availableSolvers
      use mGlobaisEscalares
      use mGeomecanica, only: ndofD, nlvectD , CYLINDER, SOLIDONLY
      use mHidrodinamicaRT, only: ndofV, nlvectV, nlvectP, ndofP 
!
      use mGeomecanica,      only: idDesloc, idiagD, neqD, nalhsD
      use mHidroDinamicaRT,  only: neqV, nalhsV, idiagV, idVeloc,simetriaVel
!
      use mMalha,            only: nen, nsd, numel
      use mMalha,            only: numLadosElem, numLadosReserv
      use mMalha,            only: numnp, x, numelReserv, numnpReserv
      use mMalha,            only: conecLadaisElem,conecNodaisElem
      use mMalha,            only: listaDosElemsPorFace
      use mMalha,            only: IrregMesh, Dirichlet, Neumann
      use mMalha,            only: IDome, I4SeaLoad
      use mMalha,            only: FNCPROCESS, FNCINIT, FNCSOLID
!
      use mLeituraEscrita,   only: iin, iecho, icoords, echo, iflag_tipoPrint, LEITuraGERAcaoCOORdenadas_DS
      use mLeituraEscritaSimHidroGeoMec,   only: leituraGeoformations_DS, ifdata
      use mLeituraEscritaSimHidroGeoMec,   only: inittime,lerDataIn_DS,lerNumericParam_DS
      use mLeituraEscritaSimHidroGeoMec,   only: lerCreepParam_DS, lerGeoMechParam_DS
      use mLeituraEscritaSimHidroGeoMec,   only: SETUPDX, abrirArquivosResultados, lerSimulatorParam_DS, lerRandfilesIn_DS
      use mLeituraEscritaSimHidroGeoMec,   only: readSetupPhaseDS
!
      use mPropGeoFisica,    only: lerPropriedadesFisicas
      use mPropGeoFisica,    only: nelx, nely, nelz
      use mPropGeoFisica,    only: nr, nrand, hx, hy, hz
      use mPropGeoFisica,    only: nelxReserv, nelyReserv, nelzReserv
      use mPropGeoFisica,    only: XTERLOAD, RHOW, RHOO, GEOMECLAW 
!
      use mHidrodinamicaRT,  only: fVeloc, pressaoElem, pressaoElemAnt,leituraCoordenadasPoco
      use mHidroDinamicaRT,  only: NCONDP,PCONDP, lerParametrosHidrodinamica_DS, optSolverV
      use mTransporte,       only: satElem, satinit
      use mGeomecanica,      only: fDesloc, InSeaLoad, optSolverD, simetriaGeo
      use mInputReader,      only: readInputFileDS, LEIturaGERacaoCOORdenadasDS
      use mInputReader,      only: leituraCodigosCondContornoDS,leituraValoresCondContornoDS

      use mMCMC,             only: INITMCMC
!
      implicit none
!
      real*8            :: t1, t2
      character(len=21) :: label
      integer           :: i, neq
      
      character(len=50) :: keyword_name
!     
      tempoTotalVelocidade  = 0.d00
      tempoTotalPressao     = 0.d00
      tempoTotalTransporte  = 0.d00
      tempoTotalGeomecanica = 0.d00
      tempoMontagemGeo=0.d00; tempoSolverGeo=0.d00
      tempoMontagemVel=0.d00; tempoSolverVel=0.d00
!
!.... Tipo Malha
!
      novaMalha = .false.
!
!.... input phase
!
      INITS3    = .TRUE.
      CYLINDER  = .FALSE.
      SALTCREEP = .FALSE.
      SOLIDONLY = .FALSE.
      availableSolvers=.false.
      
      optSolverV='skyline'
      optSolverD='skyline'
      
      open(unit=ifdata, file= 'inputDS.dat'  )
      call readInputFileDS(ifdata)
      call readSetupPhaseDS(nlvectV, nlvectP, nlvectD, optSolverV, optSolverD)
    
      availableSolvers(1)=.true. !skyline
#ifdef withPardiso
      availableSolvers(2)=.true. !pardiso
#endif
#ifdef withHYPRE
      availableSolvers(3)=.true. !hypre
#endif

     if(optSolverV=='pardiso'.or.optSolverD=='pardiso')then
        if(availableSolvers(2).eqv..false.) then
           print*, "O Solver escolhido (Pardiso) não está disponível"
           print*, "Solvers disponiveis:"
           print*, "                      SKYLINE"
           if(availableSolvers(3).eqv..true.) print*, "                      HYPRE"
           stop
        endif
     endif
     if(optSolverV=='hypre'.or.optSolverD=='hypre')then
        if(availableSolvers(3).eqv..false.) then
           print*, "O Solver escolhido (HYPRE) não está disponível"
           print*, "Solvers disponiveis:"
           print*, "                      SKYLINE"
           if(availableSolvers(2).eqv..true.) print*, "                      PARDISO"
           stop
        endif
     endif

      write(*,9001)  optSolverV, optSolverD

      if(optSolverV=="pardiso") simetriaVel=.true.
      if(optSolverD=="pardiso") simetriaGeo=.true.
      if(optSolverV=="hypre") simetriaVel=.false.
      if(optSolverD=="hypre") simetriaGeo=.false.
     
!
!.... STRESS DIMENSION PHISICS: S3DIM
!
      S3DIM = FNCPROCESS(TypeProcess)
!
      INITS3 = FNCINIT(TypeProcess)
!
      SOLIDONLY = FNCSOLID(TypeProcess)
!
      IF (SOLIDONLY) INITS3 = .TRUE.
!
!.... Logical novaMalha
      IF (IrregMesh.eq.1) novaMalha=.true.
!
      write(iecho,1000) title 
      write(iecho,3000) exec, iprtin, nsd
      write(iecho,4000) numnp, numLadosReserv, ndofP, ndofV, ndofD, &
                        nlvectP, nlvectV, nlvectD
      WRITE(IECHO,5000) nelx,nely,nelz,numnp, & 
                        nelxReserv,nelyReserv,nelzReserv,numnpReserv
      WRITE(IECHO,6000) IrregMesh, Dirichlet, Neumann, I4SeaLoad
!
      numLadosReserv = (nelxReserv+1)*(nelyReserv+1)*2
      numLadosReserv = numLadosReserv-((nelxReserv+1)+(nelyReserv+1)) 
      print*, "numLadosReserv=", numLadosReserv
!
!.... inicializa os parametros da simulacao
!
      call lerDataIn_DS
      call lerGeoMechParam_DS
      call lerNumericParam_DS
      call lerCreepParam_DS
      call lerSimulatorParam_DS
      call lerParametrosHidrodinamica_DS
      call abrirArquivosResultados
      call lerRandfilesIn_DS
      if(iflag_tipoPrint.eq.3)  CALL SETUPDX(NUMDX,  NITGEO, NITHIDRO, SALTCREEP)
!
      IF ((SOLIDONLY).AND.(SALTCREEP)) THEN 
          CYLINDER = .TRUE.
          GEOMECLAW(1) = 'CREEP'
      ENDIF
!
!.... initialization phase
!
      ndofP = 1
      ndofV = 1
!
      numLadosElem = 2*NSD
      nen          = 2**NSD
      ndofD        = NSD
!
      nr    = 1
      nrand = 1
!
!.... input coordinate data well
!      
      call leituraCoordenadasPoco(NCONDP,PCONDP)
      
!  
      call alocarMemoria()
!
!.... input coordinate data
!
      call LEITuraGERAcaoCOORdenadas_DS(x,nsd,numnp, iin, icoords, iprtin)
      
      call leituraGeoformations_DS(x, nsd, numnp)
!
!.... input boundary condition data and establish equation numbers
!
      keyword_name = "codigos_cond_contorno_veloc"
      call leituraCodigosCondContornoDS(keyword_name, idVeloc,ndofV,numLadosReserv,neqV,iecho,iprtin)    
!
      allocate(idiagV(neqV));  idiagV=0
!
!.... INPUT GEOMECHANIC DIRICHLET CONDITION DATA AND ESTABLISH EQUATION NUMBERS
!
      keyword_name = "codigos_cond_contorno_desloc"
      call leituraCodigosCondContornoDS(keyword_name,idDesloc,ndofD,numnp,neqD,iecho,iprtin)  

      allocate(idiagD(neqD));  idiagD=0
!
      print*, "ndofV=", ndofV, "neqV=", neqV
      print*, "ndofD=", ndofD, "neqD=", neqD
!
!.... input nodal force and prescribed kinematic boundary-value data
!
      keyword_name = "valores_cond_contorno_veloc"
      if (nlvectV.gt.0) call leituraValoresCondContornoDS(keyword_name,fVeloc,ndofV,numLadosReserv,1,nlvectV,iprtin)
!
!.... INPUT GEOMECHANIC DIRICHLET CONDITION DATA AND ESTABLISH EQUATION NUMBERS
!
      keyword_name = "valores_cond_contorno_desloc"
      if(geomech==1.and.nlvectD.gt.0) then
         call leituraValoresCondContornoDS(keyword_name, fDesloc, ndofD, numnp, 1_4, nlvectD, iprtin)
      endif     
!
!.... NEXT MULTIPLY SEALOAD ON Y DIRECTION 2-D MODEL NEUMANN CONDITIONS
!
      IF (I4SeaLoad.EQ.1) CALL InSeaLoad(FDESLOC,NDOFD,NUMNP, & 
     &                         NLVECTD,XTERLOAD)
!
!.... input element data
!
      call TOPologiaMALhaSistEQUAcoesDS(NALHSV, NEQV, NALHSD, NEQD) 
!
!.... inicializa os tempos de impressao
!
      call inittime
!
!.... estabelece a condicao inicial para a saturacao
! 
      call satinit(nsd,numelReserv,satElem)
! 
!.... inicializa as variaveis para o MCMC
!
      CALL INITMCMC(NCONDP,PCONDP)
!
!
 1000 format(20a4)
 3000 format(//&
     ' e x e c u t i o n   c o n t r o l   i n f o r m a t i o n '//5x,&
     ' execution code  . . . . . . . . . . . . . . (exec ) = ',i10//5x,&
     '    eq. 0, data check                                  ',   /5x,&
     '    eq. 1, execution                                   ',  //5x,&
     ' input data print code . . . . . . . . . . .(iprtin) = ',i10//5x,&
     '    eq. 0, print nodal and element input data          ',   /5x,&
     '    eq. 1, do not print nodal and element input data   ',   /5x, &
     ' number of space dimensions  . . . . . . . .(nsd   ) = ',i10)
 4000 format(5x,&
     ' number of nodal points  . . . . . . . . .  (numnp ) = ',i10//5x,&
     ' number of Lados         . . . . . . . . .(numLados) = ',i10//5x,&
     ' number of nodal degrees-of-freedom  . . . (ndofP ) = ',i10//5x,&
     ' number of nodal degrees-of-freedom  . . . (ndofV ) = ',i10//5x,&
     ' number of nodal degrees-of-freedom  . . . (ndofD ) = ',i10//5x,&
     ' number of load vectors  . . . . . . . . . (nlvectP) = ',i10//5x,&
     ' number of load vectors  . . . . . . . . . (nlvectV) = ',i10//5x,&
     ' number of load vectors  . . . . . . . . . (nlvectD) = ',i10//5x)
 5000 FORMAT(5X,  &
     &' MESH DATA FOR RESERVOIR AND OVERBUDEN DOMAINS:     '//5X,     &
     &'  ELEMENTS IN X-DIRECTION GLOBAL DOMAIN. . (NELX ) = ',I10//5X, &
     &'  ELEMENTS IN Y-DIRECTION GLOBAL DOMAIN. . (NELY ) = ',I10//5X, &
     &'  ELEMENTS IN Z-DIRECTION GLOBAL DOMAIN. . (NELZ ) = ',I10//5X, &
     &'  NODAL POINTS FOR GLOBAL DOMAIN . . . . . (NUMNP) = ',I10//5X, &
     &'  ELEMENTS IN X-DIRECTION RESERVOIR. . . . (NELX1) = ',I10//5X, &
     &'  ELEMENTS IN Y-DIRECTION RESERVOIR. . . . (NELY1) = ',I10//5X, &
     &'  ELEMENTS IN Z-DIRECTION RESERVOIR. . . . (NELZ1) = ',I10//5X, &
     &'  NODAL POINTS FOR RESERVOIR . . . . . . . (NUMNP1)= ',I10//5X) 
!
 6000 FORMAT(5X,  &
     &' MESH STRUCTURE AND READ INPUT DATA OPTIONS:           '//5X, &
     &' Reservoir Fine Mesh and External Gross  (IrregMesh) = ',I2/5X,&
     &'    eq. 0, Not Read External Files                     ',  /5x,&
     &'    eq. 1, Read Files with Coordinates and Conectivies ',  //5x,&
     &' Read Geomechanical Dirichlet Conditions (Dirichlet) = ',I2/5X,&
     &'    eq. 0, Not Read External Files                     ',  /5x,&
     &'    eq. 1, Read External File                          ',  //5x,&
     &' Read Geomechanical Neumann   Conditions (Neumann  ) = ',I2/5X,&
     &'    eq. 0, Not Read External Files                     ',  /5x,&
     &'    eq. 1, Read External File                          ',  //5x,&
     &' Multiply Sea Load on Neumann Conditions (I4SeaLoad) = ',I2/5X,&
     &'    eq. 0, Not Multiply by External Load               ',   /5x,&
     &'    eq. 1, Multiply by External Load                   ',  //)
 9001  FORMAT("Solver para a solucao do sistema de equacoes: Hidro:",A7, ", Geo:", A7)
!5x,&
!     &' Integer to Read Curve Dome Profile      (IDome    ) = ',I2/5X,&
!     &'    eq. 0, Not Read Sismic Profile for Dome            ',   /5x,&
!     &'    eq. 1, Read Sismic Profile for Dome                '//)
!
      end subroutine preprocessador_DS
!
!**** new *******************************************************************
!
      SUBROUTINE STRESS_INIT()
!
      use mGlobaisEscalares
      use mGeomecanica, only: ndofD, nlvectD 
      use mHidrodinamicaRT, only: ndofV, nlvectV, nlvectP, ndofP 
      use mGeomecanica,      only : ALHSD, idDesloc
      use mLeituraEscritaSimHidroGeoMec,   only : PRINT_DXMESH, PRINT_DXINFO, IFEDX
      use mLeituraEscrita,   only : iflag_tipoPrint
      use mMalha,            only : nsd, numel, nen, numnp
      use mMalha,            only : numelReserv, numnpReserv 
      use mMalha,            only : x, xc, conecNodaisElem
      use mPropGeoFisica,    only : lerPropriedadesFisicas, nelXReserv 
      use mPropGeoFisica,    only : YOUNG, PERMKX, MASCN
      use mPropGeoFisica,    only : GEOFORM, PORE, PWELL
      use mTransporte,       only : satElem
      use mGeomecanica,      only : NROWB,NED2,NINTD,NROWB2,IOPT
      use mGeomecanica,      only : SIGMAT, SIGMA0
      use mGeomecanica,      only : GEOSETUP, DIS, DIS0, AVSTRS
      use mGeomecanica,      only : POS4ITER, COMPTRACE, PRINT_DX
      use mGeomecanica,      only : DIVU, STRSS0, GEOPRSR, CYLINDER
      use mHidrodinamicaRT,  only : pressaoElem
!
      implicit none
!
      INTEGER  :: K, NODE, NEL
!
      NNP = 0
!
!.... READ STOCHASTIC FIELDS 
!    
      CALL lerPropriedadesFisicas()     
!
      IF (SALTCREEP) THEN
         CALL GEOSETUP(YOUNG,GEOFORM,NUMEL,NROWB,NINTD,NROWB2,IOPT)
      ENDIF
!
!.... SETUP INITIAL HIDROSTATIC PORE-PRESSURE 
!
      CALL GEOMECHANIC('INIT_GEOMECH_ARRAY')
!       
      IF (CYLINDER) THEN
         CALL CREEP_EXAMPLE()
         RETURN
      ENDIF
!
!.... MOUNT STIFFNESS MATRIX OF GEOMECHANIC 
!
      CALL GEOMECHANIC('ELASTIC_BBAR_MATRX')
      !
!.... MOUNTAGEM RIGHT HAND VECTOR FORCE ARRAY AND SOLVE
!
      IF (INITS3) CALL GEOMECHANIC('RIGHT_SIDE_2_SOLVE')
!
!.... COMPUTE INITIAL STRESS AND VOLUMETRIC DEFORMATION
!
      IF (INITS3) CALL POS4ITER(X, conecNodaisElem, STRSS0, DIVU)
!
!.... MOVE COMPUTED DISPLACEMENTS (DIS) TO INITIAL DISPLACEMENTS (DIS0)
!
      DIS0 = DIS
!
!.... COMPUTE TRACE OF TOTAL STRESS'S (SIGMAT)
!
      CALL COMPTRACE(pressaoElem,STRSS0,SIGMAT,NROWB, &
     &               NUMEL,numelReserv)
     
!
!.... CLEAR GEOMECHANICAL VECTOR FORCE ARRAY
!
      IF (TypeProcess.NE.'MANDEL') THEN
         CALL GEOMECHANIC('RESETS_FORCE_VECTR')
      ENDIF
!
!.... PRINT GEOMECHANICAL INITIAL CONDITIONS 
! 
      if(iflag_tipoPrint.eq.3)then
         CALL PRINT_DXMESH(X,DIS0,GEOPRSR,STRSS0,conecNodaisElem, & 
     &     YOUNG,PERMKX,CYLINDER,NUMDX)
      end if
!
      DIS  = 0.0D0
      DIVU = 0.0D0
!
      WRITE(*,*) "  "
      WRITE(*,*) "*** ************* ************* ************* ***"
      WRITE(*,*) "***                                           ***"
      WRITE(*,*) "*** INITIAL DISPLACEMENTS AND STRESS COMPUTED ***"
      WRITE(*,*) "***                                           ***"
      WRITE(*,*) "*** ************* ************* ************* ***"
      WRITE(*,*) "  "
!
      RETURN
!
 4500 FORMAT(I8,X,40(1PE15.8,2X))
!
      END SUBROUTINE
!
!**** new *******************************************************************
!
      SUBROUTINE CREEP_EXAMPLE()
!
      use mGlobaisEscalares
      use mLeituraEscrita,   only : iflag_tipoPrint
      use mGeomecanica, only: ndofD, nlvectD 
      use mHidrodinamicaRT, only: ndofV, nlvectV, nlvectP, ndofP 
      use mLeituraEscritaSimHidroGeoMec,   only : PRINT_DXMESH

!
      use mMalha,            only : x, conecNodaisElem
!
      use mPropGeoFisica,    only : YOUNG, PERMKX, MASCN
      use mPropGeoFisica,    only : TOLCREEP, PORE, PWELL
!
      use mGeomecanica,      only : DIS, DTRL, DIS0, CYLINDER
      use mGeomecanica,      only : DIVU, STRSS0, GEOPRSR, POS4STRS
      use mGeomecanica,      only : NWTNITER, RESIDUAL, RESMAX
!
      implicit none
!
      LOGICAL  :: LJUMP
      REAL(8)  :: ERRSIZE
!
      LJUMP = .FALSE.
!
!.... MOUNT STIFFNESS MATRIX OF GEOMECHANIC 
!
 50   CONTINUE
!....
        NWTNITER = NWTNITER + 1
        WRITE(*,2500) NWTNITER
!....
        IF (NWTNITER.GT.MAXITERC) LJUMP = .TRUE.
!
        CALL GEOMECHANIC('GEOMECHANICS_CREEP')
        ERRSIZE = RESIDUAL/RESMAX 
        write(*,5000) ERRSIZE 
! !
!.... ..TEST RESIDUAL NORM WITH TOLERANCE CRITERIA 4 CREEP
!
        IF ((ERRSIZE.GT.TOLCREEP).AND.(.NOT.LJUMP)) GOTO 50
       
!.... 
 70   CONTINUE
!
!.... UPDATE DISPLACEMENT 
!
      DIS = DTRL
!
!.... POST-PROCESS STRESS FIELD 
!
      CALL POS4STRS(X, conecNodaisElem, STRSS0, DIVU)
!
!.... MOVE COMPUTED DISPLACEMENTS (DIS) TO INITIAL DISPLACEMENTS (DIS0)
!
      DIS0 = DIS
!
!.... PRINT GEOMECHANICAL INITIAL CONDITIONS 
! 
      if(iflag_tipoPrint.eq.3)then
         CALL PRINT_DXMESH(X,DIS0,GEOPRSR,STRSS0,conecNodaisElem, & 
     &     YOUNG,PERMKX,CYLINDER,NUMDX)
      end if
!
      DIS  = 0.0D0
      DIVU = 0.0D0
!
      WRITE(*,*) "  "
      WRITE(*,*) "*** ************* ************* ************* ***"
      WRITE(*,*) "***                                           ***"
      WRITE(*,*) "***    END OF NON-LINEAR CREEP EXAMPLE        ***"
      WRITE(*,*) "***  SEE FILE 'dxcreep01.ht01/cilindr.dat'    ***"
      WRITE(*,*) "***     FOR RADIAL STRESS OUTPUT              ***"
      WRITE(*,*) "***                                           ***"
      WRITE(*,*) "*** ************* ************* ************* ***"
      WRITE(*,*) "  "
!
      RETURN
!
 2500 FORMAT('    NEWTON ITERATION COUNTER =',I5)
 4500 FORMAT(I8,X,40(1PE15.8,2X))
 5000 FORMAT('    RESIDUAL/RMAX = ',1PE15.8)
!
      END SUBROUTINE
!
!**** new *******************************************************************
!
      subroutine processamento_elast()
!
      use mGlobaisEscalares
      use mLeituraEscrita, only: iflag_tipoPrint
      use mGeomecanica, only: ndofD, nlvectD 
      use mHidrodinamicaRT, only: ndofV, nlvectV, nlvectP, ndofP 
      use mGeomecanica,      only : idDesloc, ALHSD
      use mHidroDinamicaRT,  only : nAlhsV, alhsV, brhsV,neqV
      use mGlobaisArranjos,  only : mat, c
!
      use mLeituraEscritaSimHidroGeoMec,   only : imprimirCondicoesIniciais
      use mLeituraEscritaSimHidroGeoMec,   only : imprimirSolucaoNoTempo
      use mLeituraEscritaSimHidroGeoMec,   only : isat,escreverArqParaviewIntermed
      use mLeituraEscritaSimHidroGeoMec,   only : iflag_sat,iflag_mass
      use mLeituraEscritaSimHidroGeoMec,   only : PRINT_DXINFO, IFEDX
!
      use mMalha,            only : nsd, numel,numelReserv, nen
      use mMalha,            only : numnp, numnpReserv
      use mMalha,            only : numLadosReserv, numLadosElem
      use mMalha,            only : x, xc
      use mMalha,            only : conecNodaisElem, conecLadaisElem
!
      use mPropGeoFisica,    only : phi, hx,hy,hz,nelx, nely, nelz
      use mPropGeoFisica,    only : perm, permkx, permky, permkz, phi0
      use mPropGeoFisica,    only : iflag_prod, tprt_prod
      use mPropGeoFisica,    only : np_rand_prod, dtprt_prod
      use mPropGeoFisica,    only : lerPropriedadesFisicas
      use mPropGeoFisica,    only : PORE, PORE0, PHIEULER
      use mPropGeoFisica,    only : YOUNG, GEOFORM, MASCN0, MASCN
      use mPropGeoFisica,    only : RHOW, RHOO
!
      use mTransporte,       only : transport, calc_prod      
      use mTransporte,       only : satElem, satElemAnt
      use mTransporte,       only : satElemL, satElemL0, satElem0
!
      use mGeomecanica,      only : SIGMAT, SIGMA0, GEOTIME, VDP, NED2
      use mGeomecanica,      only : NROWB, IOPT, NINTD, NROWB2
      use mGeomecanica,      only : DIS, DIS0, AVSTRS, PRINT_DX
      use mGeomecanica,      only : DIVU0, DIVU, STRSS0, TMANDEL
!
      use mHidrodinamicaRT,  only : hidroGeomecanicaRT,  pressaoElem
      use mHidrodinamicaRT,  only : velocLadal, pressaoElemAnt, vc
      use mHidroDinamicaRT,  only : NCONDP, NPRESPRODWELL, PCONDP, PRESMEDIAINICIAL, PRESPROD
      USE mPropGeoFisica,    only : GEOFORM,HX,HY
      
      use mMCMC,             only : imprimirCondicoesIniciaisMCMC,imprimirSolucaoNoTempoMCMC
!
      IMPLICIT NONE
!
!.... solution driver program 
!
      LOGICAL       :: JUMPINDXK, JUMPINDXL 
      INTEGER       :: I, INDEXL, INDEXK
      REAL*8        :: T1, T2
      character(21) :: labelTransp
!
      REAL(8)       :: ERRSZSIG, ERRSZVEL, ERRSIG0, ERRVEL0
      REAL(8)       :: SIGNORM, VELNORM
      REAL*8        :: TIMEINJ, TIMELOAD, AUX
      REAL(8), external :: DESPRESSURIZAR
      REAL(8), external :: DESPRESSURIZAR_INIT
      REAL(8), DIMENSION(numelReserv)      :: SIGMAK
      REAL(8), DIMENSION(1,numLadosReserv) :: VELOCITYL
      
!
!.... imprime condicoes iniciais
! 
      call imprimirCondicoesIniciais(pressaoElem, velocLadal, phi, &
           permkx, satElem, YOUNG, DIS, PORE, PRESPROD,ndofV, ndofP, ndofD)
           
      call imprimirCondicoesIniciaisMCMC(phi, satElem)
!
      IF(iflag_tipoPrint.EQ.3)THEN
         IF (NUMDX.GT.0) THEN
            IF (MOD(NNP,NUMDX).EQ.0) THEN 
               CALL PRINT_DX(DIS,PORE,pressaoElem,satElem,MASCN, &
     &                       STRSS0,VC,AVSTRS,NDOFD,NUMEL,NROWB,NUMNP )
            ENDIF
         ENDIF
      END IF
!
!.... tempo de simulacao para cada bloco do transporte
!
      dtBlocoTransp=tt/DFLOAT(nvel)
!
!.... INICIALIZACAO DO TEMPO DA GEOMECANICA
!
      GEOTIME = TT/DFLOAT(NVEL)
!
      IF (TypeProcess.EQ.'MANDEL') TMANDEL = 0.0D0
!
!-----------------------------------------------------------------------
!
!.... DIMINUICAO GRADUAL DA PRESSAO NO POCO 
!
!-----------------------------------------------------------------------
      AUX = DESPRESSURIZAR_INIT(PRESPROD,dtBlocoTransp)
!
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!
!     SIMULACAO
!
!-----------------------------------------------------------------------
!
      DO NNP=1,NVEL

! DESPRESURIZACAO LENTA DO POCO DE PRODUCAO
         PRESPROD = DESPRESSURIZAR(NPRESPRODWELL,NNP, &
     &                             PRESMEDIAINICIAL,AUX,PRESPROD)
!
         write(*,1000) tTransporte,nnp,nvel
!...  ..SETUP SATURATION FOR ITERATIVE HIDRODINAMICS AND TRANSPORT
         satElem0       = satElem
         satElemL       = satElem
!...  ..UPDATE FOR MACROTIME EVOLUTION
         pressaoElemAnt = pressaoElem(1,:)
         SIGMA0         = SIGMAT
         PORE0          = PORE
         DIVU0          = DIVU
         MASCN0         = MASCN 
!
         IF (TypeProcess.EQ.'MANDEL') THEN 
            TMANDEL = TMANDEL + GEOTIME
            CALL GEOMECHANIC('MANDEL_DATA_EXAMPL')
         ENDIF
!
!...  ..COMPUTE INJECTION GEOMECHANICAL TIME
         TIMEINJ   = TIMELOAD(GEOTIME*DFLOAT(NNP))
!...  ..BEGIN LOOP FOR ITERATIVE HIDRODINAMICS AND TRANSPORT (INDEX L)
         VELOCITYL = velocLadal
         INDEXL    = 0
         JUMPINDXL = .FALSE.
100      CONTINUE
            INDEXL =  INDEXL+1
!            DO 500 INDEXL=1,NITHIDRO
            IF (INDEXL.EQ.NITHIDRO) JUMPINDXL = .TRUE.
!...  .. UPDATE ARRAY SATURATIONS FOR ITERATIVE HIDRO-TRANSPORT
            satElemL0 = satElemL 
            satElemL  = satElem
            satElem   = satElem0
!
            IF (INDEXL.EQ.1) satElemL = satElemL0
!...  .. BEGIN LOOP FOR ITERATIVE HIDRODINAMICS AND GEOMECHANIC (INDEX K)
            SIGMAK  = SIGMAT
            INDEXK  = 0
            JUMPINDXK = .FALSE.
300         CONTINUE
               INDEXK = INDEXK+1
!               DO 400 INDEXK=1,NITGEO
               IF (INDEXK.EQ.NITGEO) JUMPINDXK = .TRUE.
               IF (INDEXK.EQ.1) SIGMAT = SIGMA0
               WRITE(*,2000) INDEXL, INDEXK, NITGEO 
               PHI  = PORE
               PHI0 = PHIEULER
               call hidroGeomecanicaRT(satElemL,satElemL0,SIGMAT, & 
     &              SIGMA0,TIMEINJ)
    
!...  .. ...
               CALL GEOMECHANIC('RIGHT_SIDE_2_SOLVE')
!
               CALL UPDTINIT(X,DIS,DIS0,NDOFD,NUMNP,INITS3)
!
               CALL GEOMECHANIC('ELAST_STRSS_SIGMAT')
!
               CALL GEOMECHANIC('RESETS_FORCE_VECTR')               
!.... .. ...
               ERRSZSIG = SIGNORM(SIGMAT,SIGMAK,numelReserv)
               IF (INDEXK.EQ.1) THEN
                     ERRSIG0 = ERRSZSIG
                  ELSE
                     ERRSIG0  = DMAX1(ERRSIG0,ERRSZSIG)
               ENDIF
               ERRSZSIG = ERRSZSIG/ERRSIG0
               WRITE(*,3000) 'SIGMA',INDEXK, ERRSZSIG
               SIGMAK = SIGMAT
               IF ((ERRSZSIG.GT.TOLSIGMA).AND.(.NOT.JUMPINDXK)) GOTO 300
! 400        CONTINUE ! END LOOP FOR INDEX_K: ITERATION HIDRODINAMICS AND GEOMECHANIC
!.... .. ...
         CALL GEOMECHANIC('UPDAT_MASS_CONTENT')
!...
         call transport(velocLadal)
     
         ERRSZVEL = VELNORM(velocLadal,VELOCITYL,numLadosReserv)
         IF (INDEXL.EQ.1) THEN 
              ERRVEL0 = ERRSZVEL
            ELSE
              ERRVEL0  = DMAX1(ERRVEL0,ERRSZVEL)
         ENDIF
         ERRSZVEL = ERRSZVEL/ERRVEL0
         WRITE(*,3000) 'VELOC',indexl,ERRSZVEL
         VELOCITYL = velocLadal 
         GEOMECH = 1
!...
         IF ((ERRSZVEL.GT.TOLVELOC).AND.(.NOT.JUMPINDXL)) GOTO 100
!500   CONTINUE   ! END LOOP FOR INDEX L: ITERATION HIDRODINAMICS AND TRANSPORT 
!...
      tTransporte=tTransporte+dtBlocoTransp
!
      call calc_prod(ndofV,numLadosElem,numnpReserv,numelReserv,& 
     &     nsd,nen, conecNodaisElem, conecLadaisElem,velocLadal, & 
     &     satElem,x,tTransporte)
!
      if (iflag_sat==1) then
         if (iflag_tipoPrint==1) then
            call gerarLabel(labelTransp,tTransporte)
            call escreverArqParaviewIntermed(isat, satElem, ndofV, &
                 numel, trim(labelTransp), len(trim(labelTransp)))
         endif
      endif
!
!.... imprime solucao intermediaria no tempo
! 
      call imprimirSolucaoNoTempo(satElem,DIS,PORE,YOUNG,pressaoElem, &
     &                    velocLadal, NCONDP, PCONDP, PRESPROD, tTransporte, ndofV, ndofP, ndofD)
      call imprimirSolucaoNoTempoMCMC(satElem,DIS,PORE,pressaoElem, &
                            velocLadal,NCONDP,PCONDP,tTransporte)
!
      IF(iflag_tipoPrint.EQ.3)THEN
         IF (NUMDX.GT.0) THEN
            IF (MOD(NNP,NUMDX).EQ.0) THEN 
               CALL PRINT_DX(DIS,PORE,pressaoElem,satElem,MASCN, &
     &                       STRSS0,VC,AVSTRS,NDOFD,NUMEL,NROWB,NUMNP )
            ENDIF
         ENDIF
      END IF
!
      END DO ! nnp=1,nvel
!
!.....CLOSE SERIES OF DATA AT NODAL POINTS
!
      IF(iflag_tipoPrint.EQ.3)THEN
         IF(NUMDX.GT.0) CALL PRINT_DXINFO('CLOSE_FEDX_FILE',IFEDX, & 
     &                                     NUMNP,NUMNP)
      END IF
      
      tempoTotalVelocidade=tempoMontagemVel+tempoSolverVel
      tempoTotalGeomecanica=tempoMontagemGeo+tempoSolverGeo
!
      write(*,*) " "
      write(*,*) "**********************************************"
      write(*,*) "Tempo total da velocidade=", tempoTotalVelocidade
      write(*,*) "      Montagem velocidade=", tempoMontagemVel
      write(*,*) "      Solver velocidade=", tempoSolverVel
      write(*,*) "**********************************************"
      write(*,*) "Tempo total da pressao   =", tempoTotalPressao
      write(*,*) "**********************************************"
      write(*,*) "Tempo total do transporte =", tempoTotalTransporte
      write(*,*) "**********************************************"
      write(*,*) "Tempo total da geomecanica=", tempoTotalGeomecanica
      write(*,*) "      Montagem geomecanica=", tempoMontagemGeo
      write(*,*) "      Solver geomecanica=", tempoSolverGeo
      write(*,*) "**********************************************"
      write(*,*) "SOMATORIO DOS TEMPOS =", tempoTotalVelocidade+ &
                                           tempoTotalPressao+ &
                                           tempoTotalTransporte+ &
                                           tempoTotalGeomecanica,  " segundos "

      return

 1000 FORMAT(/,'###########################################',/, &
               '###########################################',/, &
               'tempo inicial: ',1PE15.8,5x,'Passo: ',i5, 2x,'de ',i5,/, &
               '###########################################',/)
 1010 format('###########',/,'Fim da realizacao:', &
     & i5,/,'###########',/)
 2000 FORMAT('ITERATIVE L-INDEX = ',I2,2X,'ITERATIVE K-INDEX =',I2,' OF ',I2)
 2500 FORMAT('ITERATIVE WAY_S COUNTER = ',I2,' OF ',I2)
 3000 FORMAT(/'PROPORTIONAL ERROR OF ',A5,' ON ITERATION ',I2,X,'IS ',1PE15.8/)
!
      END SUBROUTINE
!
!**** new ******************************************************************
!
      SUBROUTINE UPDTINIT(X,DIS,DIS0,NDOFD,NUMNP,LFLAG)
!
      IMPLICIT NONE
!
      LOGICAL :: LFLAG
      INTEGER :: NODE, NDOFD, NUMNP,nnp
      REAL(8), DIMENSION(NDOFD,NUMNP) :: X, DIS, DIS0
!
      IF (.NOT.LFLAG) RETURN

      DO 100 NODE=1, NUMNP
          DIS(1,NODE) = DIS(1,NODE)-DIS0(1,NODE)
          DIS(2,NODE) = DIS(2,NODE)-DIS0(2,NODE)
!          write(2020+nnp,4500) NODE, DIS(1,NODE), DIS(2,NODE)
100   CONTINUE
!
      RETURN
 4500 FORMAT(I8,X,40(1PE15.8,2X))
!
      END SUBROUTINE 
!
!**** new *******************************************************************
!
      subroutine processamento_creep()
      use mLeituraEscrita, only: iflag_tipoPrint
      
      use mHidroDinamicaRT,  only : neqV,nAlhsV, alhsV, brhsV
      use mGeomecanica,      only : ALHSD
      use mGlobaisArranjos,  only : mat, c
      use mGlobaisEscalares
      use mGeomecanica, only: ndofD, nlvectD 
      use mHidrodinamicaRT, only: ndofV, nlvectV, nlvectP, ndofP 
!
      use mLeituraEscritaSimHidroGeoMec,   only : imprimirCondicoesIniciais,imprimirSolucaoNoTempo
      use mLeituraEscritaSimHidroGeoMec,   only : escreverArqParaviewIntermed
      use mLeituraEscritaSimHidroGeoMec,   only : PRINT_DXINFO, IFEDX
      use mLeituraEscritaSimHidroGeoMec,   only : isat, iflag_sat
!
      use mMalha,            only : nsd, numel, numelReserv, nen
      use mMalha,            only : numnp, numnpReserv
      use mMalha,            only : numLadosReserv, numLadosElem
      use mMalha,            only : x, xc
      use mMalha,            only : conecNodaisElem, conecLadaisElem
!
      use mPropGeoFisica,    only : phi, phi0, hx, hy, hz
      use mPropGeoFisica,    only : nelx, nely, nelz
      use mPropGeoFisica,    only : perm, permkx,permky,permkz
      use mPropGeoFisica,    only : iflag_prod, tprt_prod
      use mPropGeoFisica,    only : np_rand_prod, dtprt_prod
      use mPropGeoFisica,    only : lerPropriedadesFisicas
      use mPropGeoFisica,    only : DTCREEP, TOLCREEP
      use mPropGeoFisica,    only : PORE, PORE0, PHIEULER
      use mPropGeoFisica,    only : YOUNG, GEOFORM, MASCN0, MASCN
      use mPropGeoFisica,    only : RHOW, RHOO
!
      use mTransporte,       only : transport, satElem, satElemAnt, calc_prod
      use mTransporte,       only : satElemL, satElemL0, satElem0
      use mGeomecanica,      only : SIGMAT, SIGMA0, GEOTIME, RESMAX
      use mGeomecanica,      only : NWTNITER, RESIDUAL, VDP, NED2
      use mGeomecanica,      only : HMTTG, NROWB, IOPT, NINTD, NROWB2
      use mGeomecanica,      only : GEOSETUP, DIS, DIS0, AVSTRS, PRINT_DX
      use mGeomecanica,      only : DIVU0, DIVU, DTRL, STRSS0
!
      use mHidrodinamicaRT,  only : hidroGeomecanicaRT, vc, velocLadal
      use mHidrodinamicaRT,  only : pressaoElem, pressaoElemAnt
      use mHidroDinamicaRT,  only : NPRESPRODWELL,PRESMEDIAINICIAL,PRESPROD,PCONDP,NCONDP
!
      USE mPropGeoFisica,    only : GEOFORM,HX,HY
!
      implicit none
!
!.... solution driver program 
!
      LOGICAL       :: JUMPINDXK, JUMPINDXL, LJUMP 
      INTEGER       :: I, INDEXL, INDEXK
      REAL*8        :: T1, T2
      character(21) :: labelTransp
!
      REAL(8)       :: ERRSZSIG, ERRSZVEL, ERRSIG0, ERRVEL0
      REAL(8)       :: ERRSIZE
      REAL(8)       :: SIGNORM, VELNORM
      REAL*8        :: TIMEINJ, TIMELOAD, AUX
      REAL(8), external :: DESPRESSURIZAR
      REAL(8), external :: DESPRESSURIZAR_INIT
      REAL(8), DIMENSION(numelReserv)      :: SIGMAK
      REAL(8), DIMENSION(1,numLadosReserv) :: VELOCITYL
!     
      LJUMP = .FALSE.
!

!.... imprime condicoes iniciais
! 
      call imprimirCondicoesIniciais(pressaoElem, velocLadal, phi, & 
     &            permkx, satElem, YOUNG, DIS, PORE,PRESPROD, NDOFV, NDOFP, NDOFD)
!
      IF(iflag_tipoPrint.EQ.3)THEN
         IF (NUMDX.GT.0) THEN
            IF (MOD(NNP,NUMDX).EQ.0) THEN 
               CALL PRINT_DX(DIS,PORE,pressaoElem,satElem,MASCN, &
     &                       STRSS0,VC,AVSTRS,NDOFD,NUMEL,NROWB,NUMNP )
            ENDIF
         ENDIF
      END IF
!
!.... tempo de simulacao para cada bloco do transporte
!
      dtBlocoTransp=tt/nvel
!
!.... INICIALIZACAO DAS VARIAVEIS DA GEOMECANICA
!
      GEOTIME = TT/DFLOAT(NVEL)
      DTCREEP = DFLOAT(NCREEP)*GEOTIME
!
!-----------------------------------------------------------------------
!
!.... DIMINUICAO GRADUAL DA PRESSAO NO POCO 
!   
      AUX = DESPRESSURIZAR_INIT(PRESPROD,dtBlocoTransp)
!
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!
!     SIMULACAO
!
!-----------------------------------------------------------------------
!
      DO NNP=1,NVEL
!
! DESPRESURIZACAO LENTA DO POCO DE PRODUCAO
         PRESPROD = DESPRESSURIZAR(NPRESPRODWELL,NNP, &
    &                              PRESMEDIAINICIAL,AUX,PRESPROD)
         write(*,1000) tTransporte,nnp,nvel
!...  ..SETUP SATURATION FOR ITERATIVE HIDRODINAMICS AND TRANSPORT
         satElem0       = satElem
         satElemL       = satElem
!...  ..UPDATE FOR MACROTIME EVOLUTION
         pressaoElemAnt = pressaoElem(1,:)
         SIGMA0         = SIGMAT
         PORE0          = PORE
         DIVU0          = DIVU
         MASCN0         = MASCN 
!...  ..COMPUTE INJECTION GEOMECHANICAL TIME
         TIMEINJ        = TIMELOAD(GEOTIME*DFLOAT(NNP))
!...  ..BEGIN LOOP FOR ITERATIVE HIDRODINAMICS AND TRANSPORT (INDEX L)
         VELOCITYL = velocLadal
         INDEXL    = 0
         JUMPINDXL = .FALSE.
100      CONTINUE
            INDEXL =  INDEXL+1
!            DO 500 INDEXL=1,NITHIDRO
            IF (INDEXL.EQ.NITHIDRO) JUMPINDXL = .TRUE.
!...  .. UPDATE ARRAY SATURATIONS FOR ITERATIVE HIDRO-TRANSPORT
            satElemL0 = satElemL 
            satElemL  = satElem
            satElem   = satElem0
!
            IF (INDEXL.EQ.1) satElemL = satElemL0
!...  .. BEGIN LOOP FOR ITERATIVE HIDRODINAMICS AND GEOMECHANIC (INDEX K)
            SIGMAK  = SIGMAT
            INDEXK  = 0
            JUMPINDXK = .FALSE.
300         CONTINUE
               INDEXK = INDEXK+1
!               DO 400 INDEXK=1,NITGEO
               IF (INDEXK.EQ.NITGEO) JUMPINDXK = .TRUE.
               IF (INDEXK.EQ.1) SIGMAT = SIGMA0
               RESMAX   = 1.0D0
               NWTNITER = 0
               WRITE(*,2000) INDEXL, INDEXK, NITGEO 
               PHI  = PORE
               PHI0 = PHIEULER
               call hidroGeomecanicaRT(satElemL,satElemL0,SIGMAT, & 
     &              SIGMA0,TIMEINJ)
!...  .. ...
350            CONTINUE
!...  .. ...
               NWTNITER = NWTNITER + 1
               WRITE(*,2500) INDEXK,NITGEO,NWTNITER
!...  .. ...
               IF (NWTNITER.GT.MAXITERC) LJUMP = .TRUE.
!
               CALL GEOMECHANIC('GEOMECHANICS_CREEP')
               ERRSIZE = residual/resmax
               write(*,5000) ERRSIZE
               !bidustop 'depois de GEOMECHANICS_CREEP'
!
!.... .. ... NEXT TEST RESIDUAL NORM WITH TOLERANCE CRITERIA 4 CREEP
!
               IF ((ERRSIZE.GT.TOLCREEP).AND.(.NOT.LJUMP)) GOTO 350
!.... .. ...
!.... .. ...  UPDATE NON-LINEAR DISPLACEMENT 
!.... .. ...
               DIS = DTRL
!
!.... .. ...  CORRECTION WITH INITIAL DISPLACEMENTS
! 
               CALL UPDTINIT(X,DIS,DIS0,NDOFD,NUMNP,INITS3)
!.... .. ...
               call GEOMECHANIC('CREEP_STRSS_SIGMAT')
!.... .. ...
               ERRSZSIG = SIGNORM(SIGMAT,SIGMAK,numelReserv)
               IF (INDEXK.EQ.1) THEN
                     ERRSIG0 = ERRSZSIG
                  ELSE
                     ERRSIG0  = DMAX1(ERRSIG0,ERRSZSIG)
               ENDIF
               ERRSZSIG = ERRSZSIG/ERRSIG0
               WRITE(*,3000) 'SIGMA',INDEXK, ERRSZSIG
               SIGMAK = SIGMAT
               IF ((ERRSZSIG.GT.TOLSIGMA).AND.(.NOT.JUMPINDXK)) GOTO 300
! 400       CONTINUE ! END LOOP FOR INDEX_K: ITERATION HIDRODINAMICS AND GEOMECHANIC
!.... .. ...
         call GEOMECHANIC('UPDAT_MASS_CONTENT')
!...
         call transport(velocLadal)

         !bidustop 'aqui'
         ERRSZVEL = VELNORM(velocLadal,VELOCITYL,numLadosReserv)
         IF (INDEXL.EQ.1) THEN 
              ERRVEL0 = ERRSZVEL
            ELSE
              ERRVEL0  = DMAX1(ERRVEL0,ERRSZVEL)
         ENDIF
         ERRSZVEL = ERRSZVEL/ERRVEL0
         WRITE(*,3000) 'VELOC',indexl,ERRSZVEL
         VELOCITYL = velocLadal 
         GEOMECH = 1
!...
         IF ((ERRSZVEL.GT.TOLVELOC).AND.(.NOT.JUMPINDXL)) GOTO 100
!...
!500   CONTINUE   ! END LOOP FOR INDEX L: ITERATION HIDRODINAMICS AND TRANSPORT 
!...
      tTransporte=tTransporte+dtBlocoTransp
!
      call calc_prod(ndofV,numLadosElem,numnpReserv,numelReserv, & 
    &      nsd,nen, conecNodaisElem, conecLadaisElem,velocLadal, &
    &      satElem, x,tTransporte)
!
      if (iflag_sat==1) then
         if (iflag_tipoPrint==1) then
            call gerarLabel(labelTransp,tTransporte)
            call escreverArqParaviewIntermed(isat, satElem, ndofV, &
                 numel, trim(labelTransp), len(trim(labelTransp)))
         endif
      endif

!
!.... imprime solucao intermediaria no tempo
! 
      call imprimirSolucaoNoTempo(satElem,DIS,PORE,YOUNG,pressaoElem, &
    &       velocLadal, NCONDP, PCONDP, PRESPROD, tTransporte, NDOFV,NDOFP, NDOFD)
!
      IF(iflag_tipoPrint.EQ.3)THEN
         IF (NUMDX.GT.0) THEN
            IF (MOD(NNP,NUMDX).EQ.0) THEN
               CALL PRINT_DX(DIS,PORE,pressaoElem,satElem,MASCN, &
     &                       STRSS0,VC,AVSTRS,NDOFD,NUMEL,NROWB,NUMNP )
            ENDIF
         ENDIF
      END IF
!
      end do ! nnp=1,nvel
!
!.....CLOSE SERIES OF DATA AT NODAL POINTS
!
      IF(iflag_tipoPrint.EQ.3)THEN
         IF(NUMDX.GT.0) CALL PRINT_DXINFO('CLOSE_FEDX_FILE',IFEDX, & 
     &                                    NUMNP,NUMNP)
      END IF

      
      tempoTotalVelocidade=tempoMontagemVel+tempoSolverVel
      tempoTotalGeomecanica=tempoMontagemGeo+tempoSolverGeo
!
      write(*,*) " "
      write(*,*) "**********************************************"
      write(*,*) "Tempo total da velocidade=", tempoTotalVelocidade
      write(*,*) "      Montagem velocidade=", tempoMontagemVel
      write(*,*) "      Solver velocidade=", tempoSolverVel
      write(*,*) "**********************************************"
      write(*,*) "Tempo total da pressao   =", tempoTotalPressao
      write(*,*) "**********************************************"
      write(*,*) "Tempo total do transporte=", tempoTotalTransporte
      write(*,*) "**********************************************"
      write(*,*) "Tempo total da geomecanica=", tempoTotalGeomecanica
      write(*,*) "      Montagem geomecanica=", tempoMontagemGeo
      write(*,*) "      Solver geomecanica=", tempoSolverGeo
      write(*,*) "**********************************************"
      write(*,*) "TEMPO TOTAL DE EXECUCAO=", &
          tempoTotalVelocidade+tempoTotalPressao+tempoTotalTransporte+tempoTotalGeomecanica

!
1000  FORMAT(/,'###########################################',/, &
               '###########################################',/, &
               'tempo inicial: ',1PE15.8,5x,'Passo: ',i5, 2x,'de ',i5,/, &
               '###########################################',/)
1010  format('###########',/,'Fim da realizacao:', &
     & i5,/,'###########',/)
2000  FORMAT('ITERATIVE L-INDEX = ',I2,2X,'ITERATIVE K-INDEX =',I2,' OF ',I2)
2500  FORMAT('ITERATIVE WAY_S COUNTER = ',I2,' OF ',I2,2X,&
     &        '    NEWTON ITERATION COUNTER =',I5)
3000  FORMAT(/'PROPORTIONAL ERROR OF ',A5,' ON ITERATION ',I2,X,'IS ',1PE15.8/)
5000  FORMAT('     RESIDUAL/RMAX = ',1PE15.8//)


      return
    end subroutine processamento_creep
!
!**** NEW ******************************************************************
!
      FUNCTION TIMELOAD(REALTIME)
!
      use mGlobaisEscalares, only: YEARINJ
!
      REAL*8  :: X, TIMELOAD, REALTIME, TOL
!
      TOL = 1e-6
      IF (REALTIME-YEARINJ.LE.TOL) THEN 
                 X = 0.0D0
         ELSE
                 X = 1.0D0
      ENDIF
!
      TIMELOAD = X
!
      END FUNCTION
!
!**** NEW ******************************************************************
!
      FUNCTION SIGNORM(A,B,N) 
!
!.... FUNCTION TO COMPUTE SUPREMUM NORM OF ARRAY DIFFERENCE 
! 
      INTEGER :: I, N
      REAL(8) :: SIGNORM, XMAXIMO, XDIFF
      REAL(8), DIMENSION(N) :: A, B
!
      XMAXIMO = 1.0D0
!
      DO 100 I=1,N
         XDIFF = DABS(A(I)-B(I))
         XMAXIMO = DMAX1(XDIFF,XMAXIMO)
100   CONTINUE
!
      SIGNORM = XMAXIMO
!
      END FUNCTION
!
!**** NEW ******************************************************************
!
      FUNCTION VELNORM(A,B,N) 
!
!.... FUNCTION TO COMPUTE EUCLIDEAN NORM OF VELOCITIES DIFFERENCE 
! 
      INTEGER :: I, N
      REAL(8) :: VELNORM, XMAXIMO, XDIFF
      REAL(8), DIMENSION(1,N) :: A, B
!
      XMAXIMO = 0.0D0
!
      DO 100 I=1,N
         XDIFF   = DABS(A(1,I)-B(1,I))
         XMAXIMO = DMAX1(XDIFF,XMAXIMO)
100   CONTINUE
!
      VELNORM = XMAXIMO 
!
 4000 FORMAT(i5,2x,2(1PE15.8,2x))
      END FUNCTION
!
!**** new *******************************************************************
!
      subroutine GEOMECHANIC(fase)
      use mHidrodinamicaRT,  only: pressaoElemAnt, pressaoElem
      use mGeomecanica
      use mGeomecanica,      only: ALHSD, BRHSD, NALHSD, NEQD, idiagD, idDesloc
      use mGeomecanica,      only: optSolverD, ptD, iparmD, dparmD, ApGeo, AiGeo
      use mPropGeoFisica,    only: PORE0, PORE, PHI, PHI0, PHIEULER
      use mPropGeoFisica,    only: TOLCREEP, GEOFORM, YOUNG
      use mPropGeoFisica,    only: PERMKX,  MASCN0, MASCN
      use mPropGeoFisica,    only: RHOW, RHOO, XTERLOAD
      use mMalha,            only: x, xc, nsd, nen, numel, numnp
      use mMalha,            only: numelReserv, numnpReserv
      use mMalha,            only: conecNodaisElem, conecLadaisElem
      use mMalha,            only: numLadosElem, numLadosReserv
      use mGeomecanica, only: ndofD, nlvectD 
      use mHidrodinamicaRT, only: ndofV, nlvectV, nlvectP, ndofP 
      use mGlobaisEscalares, only: nnp, nrowsh
      use mGlobaisEscalares, only: tempoMontagemGeo, tempoSolverGeo
      use mTransporte,       only: satElemAnt, satElem
      use mSolverPardiso,    only: solverPardisoEsparso, escreverSistemaAlgCSRemMTX
      use mSolverGaussSkyline, only: solverGaussSkyline
      use mSolverHypre
      use mSolverHypre, only: destruirMatriz_HYPRE, criarMatriz_HYPRE  
!
      implicit none
!
      CHARACTER(LEN=18) :: FASE
!
      real*8           :: t1,t2
      integer          :: i, j, k, ierr
      real*8  ::  final_res_norm, tol
      integer :: num_iterations
!
      CHARACTER*18, DIMENSION(9) :: REFTASK
      real*8, external :: coldot
!
      DATA REFTASK(1)     ,     REFTASK(2)     ,     REFTASK(3)     /&
     &'INIT_GEOMECH_ARRAY','ELASTIC_BBAR_MATRX','RIGHT_SIDE_2_SOLVE'/&
     &     REFTASK(4)     ,     REFTASK(5)     ,     REFTASK(6)     /&
     &'GEOMECHANICS_CREEP','ELAST_STRSS_SIGMAT','CREEP_STRSS_SIGMAT'/&
     &     REFTASK(7)     ,     REFTASK(8)     ,     REFTASK(9)     /&
     &'RESETS_FORCE_VECTR','UPDAT_MASS_CONTENT','MANDEL_DATA_EXAMPL'/
!
      WRITE(*,*) ' ==> GEOMECHANIC TASK == ', FASE
!

      solver_id_G  = 1; precond_id_G = 2; tol = 1.0e-06
      solver_id_G  = 1; precond_id_G = 1; tol = 1.0e-06

      DO 10 I=1,9
         IF (FASE.EQ.REFTASK(I)) J=I
10    CONTINUE
!
      GOTO(100,200,300,400,500,600,700,800,900), J
!
100   CONTINUE     !.....  fase=='INIT_GEOMECH_ARRAY'
!
      CALL InicializacaoGEO()
!
      RETURN
!
200   CONTINUE     !....   fase=='ELASTIC_BBAR_MATRX'
!
      call timing(t1)
      call montarSistEqAlgGEO('bbarmatrix_elast',satElem)
      call timing(t2)
#ifdef mostrarTempos
      write(*,9002) t2-t1 
#endif
     tempoMontagemGeo=tempoMontagemGeo+(t2-t1)

      write(*,'(a)', ADVANCE='NO') '2, solucao do sistema de eq, GEOMECHANICS, '

      call timing(t1)

      if (optSolverD=='skyline') then
         write(*,'(2a)') ' direto ', optSolverD
         call solverGaussSkyline(alhsD,brhsD,idiagD,nalhsD,neqD, 'fact')
         call solverGaussSkyline(alhsD,brhsD,idiagD,nalhsD,neqD, 'back')
      end if
!
      if (optSolverD=='pardiso') then
         write(*,'(2a)') ' direto ', optSolverD
        if(NNP==0) &
        call solverPardisoEsparso(alhsD, brhsD, ApGeo, AiGeo, ptD, iparmD, dparmD, neqD, nalhsD, simetriaGeo, 'geo', 'reor')
        
        call solverPardisoEsparso(alhsD, brhsD, ApGeo, AiGeo, ptD, iparmD, dparmD, neqD, nalhsD, simetriaGeo, 'geo', 'fact')
        call solverPardisoEsparso(alhsD, brhsD, ApGeo, AiGeo, ptD, iparmD, dparmD, neqD, nalhsD, simetriaGeo, 'geo', 'back')
      endif
! 
      if(optSolverD=='hypre') then
         write(*,'(2a)') ' iterativo ', optSolverD

         call fecharMatriz_HYPRE    (A_HYPRE_G, parcsr_A_G )
         call fecharVetor_HYPRE     (b_HYPRE_G, par_b_G )
         call fecharVetor_HYPRE     (u_HYPRE_G, par_u_G )

         if(.not.allocated(initialGuess_G)) then
            allocate(initialGuess_G(neqD)); initialGuess_G=0.0
         endif
       
         call resolverSistemaAlgHYPRE (A_HYPRE_G, parcsr_A_G, b_HYPRE_G, par_b_G, u_HYPRE_G, par_u_G, &
                                   solver_G, solver_id_G, precond_id_G, tol,    &
                                   num_iterations, final_res_norm, initialGuess_G, brhsD, rows_G, neqD, myid, mpi_comm)

         call extrairValoresVetor_HYPRE(u_HYPRE_G, 1, neqD, rows_G,BRHSD)
         initialGuess_G=brhsD


         call destruirVetor_HYPRE(b_HYPRE_G)
         call destruirVetor_HYPRE(u_HYPRE_G)
!         call destruirMatriz_HYPRE(A_HYPRE_G)
!         call criarMatriz_HYPRE  (A_HYPRE_G, Clower_G, Cupper_G, mpi_comm )
         call criarVetor_HYPRE   (b_HYPRE_G, Clower_G, Cupper_G, mpi_comm )
         call criarVetor_HYPRE   (u_HYPRE_G, Clower_G, Cupper_G, mpi_comm )

      endif
!
      call timing(t2)
#ifdef mostrarTempos
      write(*,9003) t2-t1 
#endif
     tempoSolverGeo=tempoSolverGeo+(t2-t1)
      write(*,*) " 200 continue, valores nos extremos do vetor solucao geo,  "
      write(*,'(6e16.8)') brhsD(1    :6)
      write(*,'(6e16.8)') brhsD(neqD-5: neqD)



      RETURN      
!
300   CONTINUE     !....  fase=='RIGHT_SIDE_2_SOLVE'
!
      call timing(t1)
      CALL montarSistEqAlgGEO('right_hand_elast',satElem)
      call timing(t2)
#ifdef mostrarTempos
      write(*,9002) t2-t1 
#endif
     tempoMontagemGeo=tempoMontagemGeo+(t2-t1)

      write(*,'(a)', ADVANCE='NO') '3, solucao do sistema de eq, GEOMECHANICS, '

      call timing(t1)
!
      if (optSolverD=='skyline') then
         write(*,'(2a)') ' direto ', optSolverD
         call solverGaussSkyline(alhsD,brhsD,idiagD,nalhsD,neqD, 'back')
      end if
!
      IF (optSolverD=='pardiso') then
         write(*,'(2a)') ' direto ', optSolverD
!        call escreverSistemaAlgCSRemMTX(alhsD, brhsD, ApGeo, AiGeo,  nalhsD, neqD, "sistemaPardisoGeo300.mtx")
         call solverPardisoEsparso(alhsD, brhsD, ApGeo, AiGeo, ptD, iparmD, dparmD, neqD, nalhsD, simetriaGeo, 'geo', 'back')
      endif

      if(optSolverD=='hypre') then
         write(*,'(2a)') ' iterativo ', optSolverD

         call fecharMatriz_HYPRE    (A_HYPRE_G, parcsr_A_G )
         call fecharVetor_HYPRE     (b_HYPRE_G, par_b_G )
         call fecharVetor_HYPRE     (u_HYPRE_G, par_u_G )

         if(.not.allocated(initialGuess_G)) then
            allocate(initialGuess_G(neqD)); initialGuess_G=0.0
         endif

       !  write(*,*) "call HYPRE_IJMatrixPrint"
       !  call HYPRE_IJMatrixPrint( A_HYPRE_G, "Alhs300.out.", ierr)

         call resolverSistemaAlgHYPRE (A_HYPRE_G, parcsr_A_G, b_HYPRE_G, par_b_G, u_HYPRE_G, par_u_G, &
                                   solver_G, solver_id_G, precond_id_G, tol,    &
                                   num_iterations, final_res_norm, initialGuess_G, brhsD, rows_G, neqD, myid, mpi_comm)

         call extrairValoresVetor_HYPRE(u_HYPRE_G, 1, neqD, rows_G,BRHSD)
         initialGuess_G=brhsD


         call destruirVetor_HYPRE(b_HYPRE_G)
         call destruirVetor_HYPRE(u_HYPRE_G)
         call destruirMatriz_HYPRE(A_HYPRE_G)
         call criarMatriz_HYPRE  (A_HYPRE_G, Clower_G, Cupper_G, mpi_comm )
         call criarVetor_HYPRE   (b_HYPRE_G, Clower_G, Cupper_G, mpi_comm )
         call criarVetor_HYPRE   (u_HYPRE_G, Clower_G, Cupper_G, mpi_comm )

      endif


!
!.... UPDATE DISPLACEMENT 
! 
      CALL BTOD(idDesloc,DIS,BRHSD,NDOFD,NUMNP)
      
      call timing(t2)
#ifdef mostrarTempos
      write(*,9003) t2-t1 
#endif
     tempoSolverGeo=tempoSolverGeo+(t2-t1)

      write(*,*) " 300 continue, valores nos extremos do vetor solucao geo,  "
      write(*,'(6e16.8)') brhsD(1    :6)
      write(*,'(6e16.8)') brhsD(neqD-5: neqD)
!
      !stop
      RETURN
!
400   CONTINUE     !.... fase=='GEOMECHANICS_CREEP'


!      call destruirMatriz_HYPRE(A_HYPRE_G)
!      call criarMatriz_HYPRE (A_HYPRE_G, Clower_G, Cupper_G, mpi_comm )
!
      call timing(t1)
      call montarSistEqAlgGEO('bbarmatrix_creep',satElem)
      call timing(t2)


#ifdef mostrarTempos
      write(*,9002) t2-t1 
#endif
     tempoMontagemGeo=tempoMontagemGeo+(t2-t1)

      write(*,*) " 400 continue, valores nos extremos do vetor BRHS geo,  "
      write(*,'(6e16.8)') brhsD(1    :6)
      write(*,'(6e16.8)') brhsD(neqD-5: neqD)

!
!.... COMPUTE RESIDUAL EUCLIDEAN NORM
!
      RESIDUAL = dsqrt(COLDOT(BRHSD,BRHSD,NEQD))
!      IF ((NWTNITER.EQ.1)) RESMAX = MAX(RESIDUAL,RESMAX)
      RESMAX = DMAX1(RESIDUAL,RESMAX)
!     
      WRITE(*,4000) RESIDUAL,resmax

      write(*,'(a)', ADVANCE='NO') '4, solucao do sistema de eq, GEOMECHANICS, '
      call timing(t1)

      if (optSolverD=='skyline') then
         write(*,'(2a)') ' direto ', optSolverD
!          call solverDiretoSkyLine(alhsD, brhsD, idiagD, nalhsD, neqD, 'geo') 
         call solverGaussSkyline(alhsD,brhsD,idiagD,nalhsD,neqD, 'full')
      end if

      if (optSolverD=='pardiso') then
         write(*,'(2a)') ' direto ', optSolverD
!        call escreverSistemaAlgCSRemMTX(alhsD, brhsD, ApGeo, AiGeo,  nalhsD, neqD, "sistemaPardisoGeo400.mtx")
         call solverPardisoEsparso(alhsD, brhsD, ApGeo, AiGeo, ptD, iparmD, dparmD, neqD, nalhsD, simetriaGeo, 'geo', 'full')
      endif
      
      if(optSolverD=='hypre') then
         write(*,'(2a)') ' iterativo ', optSolverD

         call fecharMatriz_HYPRE    (A_HYPRE_G, parcsr_A_G )
         call fecharVetor_HYPRE     (b_HYPRE_G, par_b_G )
         call fecharVetor_HYPRE     (u_HYPRE_G, par_u_G )

         if(.not.allocated(initialGuess_G)) then
            allocate(initialGuess_G(neqD)); initialGuess_G=0.0
         endif

         !write(*,*) "call HYPRE_IJMatrixPrint"
         !call HYPRE_IJMatrixPrint( b_HYPRE_G, "Blhs400.out.", ierr)
         !call HYPRE_IJMatrixPrint( A_HYPRE_G, "Alhs400.out.", ierr)

       
         call resolverSistemaAlgHYPRE (A_HYPRE_G, parcsr_A_G, b_HYPRE_G, par_b_G, u_HYPRE_G, par_u_G, &
                                   solver_G, solver_id_G, precond_id_G, tol,    &
                                   num_iterations, final_res_norm, initialGuess_G, brhsD, rows_G, neqD, myid, mpi_comm)

         call extrairValoresVetor_HYPRE(u_HYPRE_G, 1, neqD, rows_G,BRHSD)
         initialGuess_G=brhsD


         call destruirVetor_HYPRE(b_HYPRE_G)
         call destruirVetor_HYPRE(u_HYPRE_G)
         call destruirMatriz_HYPRE(A_HYPRE_G)
         call criarMatriz_HYPRE  (A_HYPRE_G, Clower_G, Cupper_G, mpi_comm )
         call criarVetor_HYPRE   (b_HYPRE_G, Clower_G, Cupper_G, mpi_comm )
         call criarVetor_HYPRE   (u_HYPRE_G, Clower_G, Cupper_G, mpi_comm )
      endif

      write(*,*) " 400 continue, valores nos extremos do vetor solucao geo,  "
      write(*,'(6e16.8)') brhsD(1    :6)
      write(*,'(6e16.8)') brhsD(neqD-5: neqD)

     ! write(911,*) BRHSD


      call timing(t2)
#ifdef mostrarTempos
      write(*,9003) t2-t1 
#endif
     tempoSolverGeo=tempoSolverGeo+(t2-t1)
     
!
!.... UPDATE TRIAL DISPLACEMENT WITH INCREMENT 
!
      CALL UPDATEINCR(idDesloc,BRHSD,NDOFD,NUMNP)
     
!
!.... NEXT LINE VISCO-ELASTIC EVOLUTION AND TANGENT MATRIX UPDATE
!
      CALL POS4CREEP(x, conecNodaisElem)
!
      !stop
      RETURN
!
500   CONTINUE     !.... fase='ELAST_STRSS_SIGMAT'
!
!.... POST-PROCESS STRESS FIELD 
!
      CALL POS4ITER(X, conecNodaisElem, AVSTRS, DIVU) 
!
!... COMPUTE TRACE OF TOTAL STRESS'S (SIGMAT)
!
      CALL COMPTRACE(pressaoElem,AVSTRS,SIGMAT,NROWB,&
     &               NUMEL,numelReserv)
!
      RETURN
!
600   CONTINUE     !.... fase=='CREEP_STRSS_SIGMAT'
!
      VDP   = 0.0D0
!
!.... POST-PROCESS STRESS FIELD 
!
      CALL POS4STRS(x, conecNodaisElem, AVSTRS, DIVU)
!
!.... COMPUTE TRACE OF TOTAL STRESS'S (SIGMAT)
!
      CALL COMPTRACE(pressaoElem,AVSTRS,SIGMAT,NROWB,NUMEL,numelReserv)
!
      RETURN
!
700   CONTINUE     !....  fase=='RESETS_FORCE_VECTR'
!
      call timing(t1)
      CALL montarSistEqAlgGEO('right_hand_reset',satElem)
      call timing(t2)
#ifdef mostrarTempos
      write(*,9002) t2-t1 
#endif
     tempoMontagemGeo=tempoMontagemGeo+(t2-t1)
!
      RETURN
!
800   CONTINUE     !.... fase=='UPDAT_MASS_CONTENT'
!
      CALL POS4_MASSCNT(DIVU,DIVU0, pressaoElem, pressaoElemAnt, &
     &         PORE,PORE0,YOUNG,MASCN,MASCN0,PHIEULER,NUMEL,NUMELRESERV)
!
      PHI  = MASCN
      PHI0 = MASCN0
!
      RETURN
!
900   CONTINUE   !.... fase=='MANDEL_DATA_EXAMPL'
!
      CALL FTODMANDL(FDesloc,NDOFD,NUMNP,NLVECTD,XTERLOAD,TMANDEL)
!
      RETURN
!
 4000 FORMAT(3X,'RESIDUAL = ',1PE15.8,2X,'RMAX =',1PE15.8)
 9002 FORMAT( "Tempo de montagem da matriz e/ou vetor força = ",f12.5)
 9003 FORMAT( "Tempo do solver da geomecanica = ",f12.5)
!
      END SUBROUTINE
!
!**** new **********************************************************************
!
      subroutine InicializacaoGEO()
!     
      use mGeomecanica
      use mGlobaisEscalares, only:  NUMDX, NNP
      use mGeomecanica, only: ndofD, nlvectD 
      use mGlobaisEscalares, only: TypeProcess
!
      use mHidrodinamicaRT,  only: pressaoElem, vc
      use mTransporte,       only: satElem, satElemAnt
      use mPropGeoFisica,    only: phi, phi0, PERMKX, YOUNG
      use mPropGeoFisica,    only: PORE, PORE0, MASCN, MASCN0
      use mPropGeoFisica,    only: RHOW, RHOO, XTERLOAD, PHIEULER
!
      use mMalha,            only: numel, nsd, nen, numnp
      use mMalha,            only: conecNodaisElem
      use mMalha,            only: x, xc, numelReserv
!
      use mGeomecanica,     only: idDesloc, IDIAGD, NALHSD
      use mGeomecanica,     only: ALHSD, BRHSD, NEQD,LMD
!
!..NEW LINES 4 COMPRESSIBLE MODEL MASS CONTENT ARRAYS
!
      IMPLICIT NONE
!
      integer :: i
!.... SETUP SATURATION AT TIME ZERO
!
      satElemAnt = satElem
!
!.... SETUP EULERIAN AND 'CLASSICAL' POROSITIES AT TIME ZERO
!
      PHI0     = PHI
      PORE     = PHI
      PORE0    = PHI0
      MASCN    = PORE
      MASCN0   = PORE0
      PHIEULER = PORE
!
      if(.not.allocated(ALHSD)) allocate(ALHSD(NALHSD))
      if(.not.allocated(BRHSD)) allocate(BRHSD(NEQD))
!
!.... CLEAR LEFT AND RIGHT HAND SIDE FOR THE ELASTIC  PROBLEM 
! 
      ALHSD = 0.0D0
      BRHSD = 0.0D0
!
      IF (TypeProcess.EQ.'MANDEL') THEN 
         CALL FTODMANDL(FDesloc,NDOFD,NUMNP,NLVECTD,XTERLOAD,0.0D0)
      ENDIF
!
!.... ACCOUNT THE NODAL FORCES IN THE R.H.S. 
!
      IF (NLVECTD.GT.0) &
         CALL LOAD(idDesloc,fDesloc,brhsd,NDOFD,NUMNP,NLVECTD)
!
      IF (NLVECTD.GT.0) &
         CALL FTOD(idDesloc,DIS,fDesloc,NDOFD,NUMNP,NLVECTD)
         
!
!.... SETUP INITIAL HIDROSTATIC PRESSURE FIELD
!
      CALL PRSRINIT(GEOPRSR,pressaoElem,XC,NUMEL,numelReserv) 
!
      RETURN
!
      END SUBROUTINE
!
!**** NEW **********************************************************************
!
      subroutine alocarMemoria()
      use mGlobaisArranjos,  only: uTempoN, mat, grav, beta
      use mGeomecanica,      only: ndofD, nlvectD 
      use mHidrodinamicaRT,  only: ndofV, nlvectV, nlvectP, ndofP 
      use mGlobaisEscalares, only: geomech, novaMalha
!
      use mMalha,            only: nsd, numel, numelReserv, numnp
      use mMalha,            only: numLadosReserv, nen, numLadosElem
      use mMalha,            only: x, xc
!
      use mMalha,            only: listaDosElemsPorNo
      use mMalha,            only: listaDosElemsPorFace
      use mMalha,            only: conecNodaisElem, conecLadaisElem
      use mHidroDinamicaRT,  only: idVeloc, lmV
      use mGeomecanica,      only: idDesloc, lmD
      use mTransporte,       only: satElemAnt, satElem
      use mTransporte,       only: satElemL, satElemL0, satElem0
      use mGeomecanica
      use mPropGeoFisica,    only: GEOFORM, MASCN, MASCN0, PHIEULER  
      use mHidrodinamicaRT,  only: pressaoElem, pressaoElemAnt,PRESPRODWELL,NCONDP
      use mHidrodinamicaRT,  only: vc, ve, velocLadal, fVeloc
      use mSolverPardiso, only: listaDosElemsPorNoCRS
!
      implicit none
!
!campos
!
      allocate(pressaoElem(ndofP,numelReserv))
      pressaoElem = 0.0D0
      allocate(pressaoElemAnt(numelReserv))
      pressaoElemAnt = 0.0D0
      allocate(velocLadal (ndofV,numLadosReserv))
      velocLadal = 0.0D0
      ALLOCATE(PRESPRODWELL(NCONDP))
      PRESPRODWELL = 0.0D0
      allocate(ve (nsd,nen,numelReserv))
      ve = 0.0D0
      allocate(vc (nsd,numelReserv))
      vc = 0.0D0
      allocate(satElem(numelReserv))
      satElem = 0.0D0
      allocate(satElemAnt(numelReserv))
      satElemAnt  = 0.0D0
!
! NEXT 3 LINES NEW SATURATION ARRAY FOR ITERATIVE HIDRODINAMICS AND TRANSPORT 
!
      ALLOCATE(satElem0(numelReserv))
      satElem0    = 0.0D0
      ALLOCATE(satElemL(numelReserv))
      satElemL    = 0.0D0
      ALLOCATE(satElemL0(numelReserv))
      satElemL0   = 0.0D0
!
!malha
!
      allocate(x(nsd,numnp))
      x  = 0.0d0
      allocate(xc(nsd,numel))
      xc = 0.0d0
      allocate(conecNodaisElem(nen,numel))
      conecNodaisElem = 0
      allocate(conecLadaisElem(numLadosElem,numelReserv))
      conecLadaisElem = 0
      allocate(listaDosElemsPorNo(nen,numnp))
      listaDosElemsPorNo = 0
      allocate(listaDosElemsPorFace(numLadosElem,numLadosReserv))
      listaDosElemsPorFace = 0

#ifdef withcrs
      if(novaMalha.eqv..true.) then
         allocate(listaDosElemsPorNoCRS(7,numnp))
         listaDosElemsPorNoCRS = 0
      else 
         allocate(listaDosElemsPorNoCRS(nen,numnp))
         listaDosElemsPorNoCRS = 0
      endif
#endif
!
!contorno
!
      allocate(idVeloc(ndofV,numLadosReserv))
      idVeloc = 0
      allocate(lmV(ndofV,numLadosElem,numelReserv))
      lmV = 0
!
      if (nlvectV.ne.0)  then
         allocate(fVeloc(ndofV,numLadosReserv))
         fVeloc = 0.0d0
      endif
!
!material
      allocate(mat(numel)); mat=0.d0
!
!gravidade
      allocate(grav(3));   grav=0.d0
!
!tempo
      allocate(uTempoN(numel))
!
!.... GEOMECANICS 2-D MODEL ARRAYS 
!
      IF (NSD.EQ.2) THEN
            NROWB  = 4  
            NINTD  = 4
      ENDIF
      NROWB2 = NROWB*NROWB
      ALLOCATE(idDesloc(ndofD,NUMNP));    idDesloc  = 0
      ALLOCATE(LMD(ndofD,NEN,NUMEL));     LMD       = 0
!
      IF (nlvectD.ne.0) THEN
          ALLOCATE(fDesloc(ndofD,numnp));   fDesloc = 0.0d0
      ENDIF
      ALLOCATE(IDDIS(NDOFD, NUMNP));        IDDIS   = 0
      ALLOCATE(DIS(NDOFD, NUMNP));          DIS     = 0.0D0
      ALLOCATE(DIS0(NDOFD, NUMNP));         DIS0    = 0.0D0
      ALLOCATE(VDP(NDOFD,NUMNP));           VDP     = 0.0D0
      ALLOCATE(DTRL(NDOFD,NUMNP));          DTRL    = 0.0D0
      ALLOCATE(DIVU  (NUMEL));              DIVU    = 0.0D0
      ALLOCATE(DIVU0 (NUMEL));              DIVU0   = 0.0D0
!
      ALLOCATE(STRSS (NUMEL,NINTD,NROWB));  STRSS   = 0.0D0
      ALLOCATE(HMTTG (NUMEL,NINTD,NROWB2)); HMTTG   = 0.0D0
      ALLOCATE(ECREEP(NUMEL,NINTD,NROWB));  ECREEP  = 0.0D0
!
      ALLOCATE(GEOPRSR(NUMEL));             GEOPRSR = 0.0D0
      ALLOCATE(STRSS0(NROWB,NUMEL));        STRSS0  = 0.0D0
      ALLOCATE(AVSTRS(NROWB,NUMEL));        AVSTRS  = 0.0D0
      ALLOCATE(AVCREP(NUMEL,NROWB));        AVCREP  = 0.0D0
!
      ALLOCATE(SIGMAT(numelReserv));        SIGMAT  = 0.0D0
      ALLOCATE(SIGMA0(numelReserv));        SIGMA0  = 0.0D0
      ALLOCATE(GEOFORM(NUMEL));             GEOFORM = 'NONE'
!
!... ALLOCATE MEMORY 4 MASS CONTENT
!
      allocate(MASCN (numelReserv));        MASCN   = 0.0D0
      allocate(MASCN0(numelReserv));        MASCN0  = 0.0D0
!
!... END GEOMECANICA
!
      allocate(beta(numel))
!
      end subroutine
      
!
!**** NEW **** MODIFIED 4 HIERARCH MESH  ************************************
!
      SUBroutine TOPologiaMALhaSistEQUAcoesDS(NALHSV,NEQV,NALHSD,NEQD)
      use mGlobaisEscalares
      use mGeomecanica, only: ndofD, nlvectD, optSolverD, LMStencilGeo
      use mHidrodinamicaRT, only: ndofV, nlvectV, nlvectP, ndofP , optSolverV, LMStencilVel
      use mGlobaisArranjos
      use mLeituraEscrita, only: iin, iecho, nprint, prntel
      use mLeituraEscritaSimHidroGeoMec, only: GEOREGION_DS 
!
      use mHidroDinamicaRT, only: lmV, idVeloc, idiagV, nedV, ApVel, AiVel
      use mHidrodinamicaRT, only: fVeloc, simetriaVel
      use mHidroDinamicaRT, only: numCoefPorLinhaVel,  meanbwV
      use mGeomecanica,     only: LMD, IDIAGD,ALHSD, BRHSD, CLHSD, IDDESLOC, ApGeo, AiGeo
      use mGeomecanica,     only: AUXM,IOPT,NED,NED2,NEE,NEE2,NEESQ,NEESQ2,NESD,NSTR,simetriaGeo
      use mGeomecanica,     only: numCoefPorLinhaGeo
!
      use mMalha,          only: IrregMesh
      use mMalha,          only: numel, numnp, nsd, nen
      use mMalha,          only: numLadosElem,numLadosReserv
      use mMalha,          only: numelReserv,numnpReserv
      use mMalha,          only: x, xc, local
      use mMalha,          only: conecNodaisElem, conecLadaisElem
      use mMalha,          only: listaDosElemsPorNo
      use mMalha,          only: listaDosElemsPorFace
      use mMalha,          only: criarListaVizinhos
      use mMalha,          only: genel, genelFaces
      use mMalha,          only: formlm, renumerarMalha
!
      use mPropGeoFisica,  only: hx,hy,hz,nelx,nely,nelz, calcdim
      use mPropGeoFisica,  only: nelxReserv, nelyReserv, nelzReserv
      use mPropGeoFisica,  only: GEOFORM
!
      use mSolverPardiso, only: criarPonteirosMatEsparsa_CRS
      use mSolverPardiso, only: listaDosElemsPorNoCRS
      use mSolverPardiso, only: criarListaVizinhosCRS
      
      use mInputReader,      only: readNodeElementsDS, readMaterialPropertiesDS, readConstantBodyForcesDS
      use mInputReader,      only: genelFacesDS, leituraGeracaoConectividadesDS

      use mSolverHypre
!

!
      implicit none
!
!.... program to set arrays storage
!
      integer :: NALHSV, NEQV, NALHSD, NEQD
      integer :: numConexoesPorElem
      integer :: i, m, n, nel
      real*8, allocatable :: xl(:,:) 
      real*8 :: t1, t2
!
!... NEW FOR IRREGULAR  MESH
!
      INTEGER, DIMENSION(NEN,NUMEL)  :: OLDIEN
      INTEGER :: INIEN, NODE
      CHARACTER*30 NAMEIN
!... END FOR IRREGULAR MESH

      character(len=50) keyword_name
      integer :: ierr
      
      integer*4 :: localSize_V, localSize_G
!
!.... calculate hx, hy, hz
!
       call calcdim(nsd,numnpReserv,x)
       WRITE(*,1001) hx, hy, hz
!
!.... set element parameters
!
      allocate(npar(numParElem))
      call readNodeElementsDS
!
      if(nsd==2)nrowsh = 3
      if(nsd==3)nrowsh = 4
      NROWSH = NSD + 1
!
      nedV   = ndofV
      nprint = 0
!
      NED    = 1
      NEE    = NEN*NED
      NEESQ  = NEE*NEE
      NESD   = 2
      NNP    = 0
!..Bbar.. BEGIN
!..Bbar..  SET ELEMENT PARAMETERS VISCOELASTIC INCOMPRESSIVEL MODEL
!..Bbar..
      NED2   = 2
      NEE2   = NEN*NED2 
      NEESQ2 = NEE2*NEE2 
      NSTR   = 3 
      IOPT   = 1
!
      allocate(AUXM(NUMEL,NEE2,NEE2)); AUXM = 0.0D0
!
!..Bbar.. END
!
      if (nicode.eq.0) nicode=nen
      npint   = nicode 
!
!....... set memory pointers
! 
      allocate(c(6,numat)); c=0.d0
!
      write(iecho,1000) ntype,numel,numat,nen,npint
!
!      read material properties
!
      keyword_name = "material_properties"
      call readMaterialPropertiesDS(keyword_name, iecho, ierr)
!
!     constant body forces
!
      keyword_name = "constant_body_forces"
      call readConstantBodyForcesDS(keyword_name, iecho, ierr)
!
!    generation of conectivities
!
      keyword_name = "conectividades_nodais"
      call leituraGeracaoConectividadesDS(keyword_name, conecNodaisElem,mat,nen, ierr)
!
      keyword_name = "conectividades_ladais"
      if (nsd==2) then
         call leituraGeracaoConectividadesDS(keyword_name, conecLadaisElem, mat, numLadosELem, ierr)
      else 
         call genelFacesDS(keyword_name, conecLadaisElem, numLadosElem, nelxReserv, nelyReserv, nelzReserv, ierr)
      endif
      
      listaDosElemsPorNo=0
      call criarListaVizinhos(nen,numnp,numel,conecNodaisElem, & 
     &     listaDosElemsPorNo  )
      listaDosElemsPorFace=0
      call criarListaVizinhos(numladosElem,numLadosReserv, & 
     &     numelReserv,conecLadaisElem,listaDosElemsPorFace)


#ifdef withcrs
       if(novaMalha.eqv..true.)then
          listaDosElemsPorNoCRS=0
          call criarListaVizinhosCRS(nen, numnp, numel,7,& 
     &         conecNodaisElem,listaDosElemsPorNoCRS  )
       endif
#endif

#ifdef debug
      if(iprtin.eq.0) then
         call prntel(mat,conecNodaisElem,nen,numel,1)
         call prntel(mat,conecLadaisElem,numLadosELem,numelReserv,2)
      end if
#endif
!
!.... generation of velocity lm array
!
      call formlm(idVeloc,conecLadaisElem,lmV,nedV,nedV,& 
     &     numLadosElem,numelReserv)
!
!.... modification of idiag array
!
      call colht(idiagV,lmV,nedV,numLadosElem,numelReserv,neqV)
!
      if(optSolverV=='skyline') then
         call diag(idiagV,neqV,nalhsV)
      endif
!
      meanbwV = nalhsV/neqV
      write(iecho,9000) neqV, nalhsV, meanbwV, 8.0*nalhsV/1000/1000

!     calculo de xc (centros dos elementos)
      allocate(xl(nsd,nen)); xl=0.0
      do i=1,numel
         call local(conecNodaisElem(1,i),x,xl,nen,nsd,nsd)
         xc(1,i) = sum(xl(1,1:nen))/nen
         xc(2,i) = sum(xl(2,1:nen))/nen
         if(nsd==3)xc(3,i) = sum(xl(3,1:nen))/nen
      end do
!    
!..Bbar.. BEGIN
!-----------------------------------------------------------------------
!     B-BAR FORMULATION FOR DISPLACEMENTS 
!-----------------------------------------------------------------------
! 
!.... CLEAR IDIAG ARRAY 
! 
      IDIAGD=0
! 
!.... GENERATION OF LMD ARRAY 
!
      CALL FORMLM(idDesloc,conecNodaisElem,LMD,NDOFD,NED2,NEN,NUMEL) 
!
!.... MODIFICATION OF IDIAGD ARRAY 
!
      CALL COLHT(IDIAGD,LMD,NED2,NEN,NUMEL,NEQD)
!
!.... DETERMINE ADDRESSES OF DIAGONALS IN LEFT-HAND-SIDE MATRIX 
!
      if(optSolverD=='skyline') then
         CALL DIAG(IDIAGD,NEQD,NALHSD)
      endif
!
!.... SETUP GEO-MECHANICAL REGIONS   
!
      CALL GEOREGION_DS(XC,NUMEL)
!

      if(optSolverV=='pardiso'.or.optSolverV=='hypre') then
         call timing(t1)
         if(nsd==2) numCoefPorLinhaVel=7
         if(nsd==3) numCoefPorLinhaVel=11
         numConexoesPorElem=numLadosElem
         call criarPonteirosMatEsparsa_CRS(nsd, ndofV, neqV, & 
     &        numCoefPorLinhaVel, conecLadaisElem, & 
     &        listaDosElemsPorFace, idVeloc, LMStencilVel, ApVel, AiVel, numLadosReserv, &
     &        numLadosElem, numConexoesPorElem, nalhsV, & 
     &        simetriaVel)
!
         call timing(t2)
#ifdef mostrarTempos
         write(*,9002) t2-t1 
#endif
     endif !if(optSolverV=='pardiso'.or.optSolverV=='hypre') then
!     

     if(optSolverD=='pardiso'.or.optSolverD=='hypre') then
         call timing(t1)
         if(geomech==1) then
            if(novaMalha.eqv..true.) then
               numCoefPorLinhaGeo=26
               numConexoesPorElem=7
               call criarPonteirosMatEsparsa_CRS(nsd, ndofD, neqD, &
     &              numCoefPorLinhaGeo, conecNodaisElem, & 
     &              listaDosElemsPorNoCRS, idDesloc, LMStencilGeo, ApGeo, AiGeo, numnp, nen, & 
     &              numConexoesPorElem, nalhsD, simetriaGeo)
            else
               numCoefPorLinhaGeo=18
               numConexoesPorElem=4
               call criarPonteirosMatEsparsa_CRS(nsd, ndofD, neqD, & 
     &              numCoefPorLinhaGeo, conecNodaisElem, & 
     &              listaDosElemsPorNo, idDesloc, LMStencilGeo, ApGeo, AiGeo, numnp, nen, & 
     &              numConexoesPorElem, nalhsD, simetriaGeo)

            endif

            call timing(t2)
#ifdef mostrarTempos
            write(*,9003) t2-t1
#endif
         endif
     endif !if(optSolverD=='pardiso'.or.optSolverD=='hypre') then

!
      print*, "NALHSV=", NALHSV, "NALHSD=", NALHSD
      
      if(optSolverV=="hypre") then
      
         Clower_V = 1 - 1
         Cupper_V = neqV-1
         localSize_V=CUpper_V-Clower_V+1
         if(.not.allocated(rows_V)) allocate(rows_V(localSize_V))
         do i = 1, localSize_V
            rows_V(i) = i - 1
         end do
         call criarMatriz_HYPRE    (A_HYPRE_V, Clower_V, Cupper_V, mpi_comm )
         call criarVetor_HYPRE     (b_HYPRE_V, Clower_V, Cupper_V, mpi_comm )
         call criarVetor_HYPRE     (u_HYPRE_V, Clower_V, Cupper_V, mpi_comm )

      endif

      if(optSolverD=="hypre") then
      
         Clower_G = 1 - 1
         Cupper_G = neqD-1
         localSize_G=CUpper_G-Clower_G+1
         if(.not.allocated(rows_G)) allocate(rows_G(localSize_G))
         do i = 1, localSize_G
            rows_G(i) = i - 1
         end do
         call criarMatriz_HYPRE    (A_HYPRE_G, Clower_G, Cupper_G, mpi_comm )
         call criarVetor_HYPRE     (b_HYPRE_G, Clower_G, Cupper_G, mpi_comm )
         call criarVetor_HYPRE     (u_HYPRE_G, Clower_G, Cupper_G, mpi_comm )

     endif 
     
!
      return
!
 1000 format(//,&
      ' two/three-n o d e    e l e m e n t s ',//,5x,&
      ' element type number . . . . . . . . . . (ntype ) = ',i10,//5x,&
      ' number of elements  . . . . . . . . . . (numel ) = ',i10,//5x,&
      ' number of element material sets . . . . (numat ) = ',i10,//5x,&
      ' number of element nodes . . . . . . . . (nen   ) = ',i10,//5x,&
      ' number of integration points. . . . . . (npint  ) = ',i10)
 1001  FORMAT(' hx= ',F12.5,2X,'hy= ',F12.5,2X,' hz= ',F12.5)
 1500  FORMAT(I10)
 4000  format(///,&
     ' m a t e r i a l   s e t   d a t a             ',  //5x,&
     ' number of material sets . . . . . . .(numat ) = ',i10,//,2x,&
     & 'set',4x,'Kx ',4x,'Ky',4x,'Kz')
 4500  FORMAT(27I8)
 5000  format(i10,5x,5f10.0)
 6000  format(2x,i3,1x,5(1x,1pe11.4))
 7000 format(8f10.0)
 8000 format(///,&
      ' g r a v i t y   v e c t o r   c o m p o n e n t s     ',//5x,&
      ' exemplo 1. . . . . . . . . . . . . .  = ',      1pe15.8,//5x,&
      ' exemplo 2 . . . . . . . . . . . . . . = ',      1pe15.8,//5x,&
      ' exemplo 3............................ = ',      1pe15.8,//)
 9000 format('Calculo das velocidades nodais' //&
      ' e q u a t i o n    s y s t e m    d a t a              ',  //5x,&
      ' number of equations . . . . . . . . . . . . (neq    ) = ',i8//5x,&
      ' number of terms in left-hand-side matrix  . (nalhs  ) = ',i12//5x,&
      ' mean half bandwidth . . . . . . . . . . . . (meanbw ) = ',i8//5x,&
      ' memoria necessaria para a matriz do sistema  (Mbytes)  = ',e10.2)
9002  FORMAT( "Tempo de pre-processamento da Velocidade para solver externo",f12.5)
9003  FORMAT( "Tempo de pre-processamento da Geomecanic para solver externo",f12.5)
9500  FORMAT(5X,  &
     &' ** ** ** ** ** ** ** ** *** ** ** ** ** ** ** ** ** **',/5X,&
     &' ** ** ** ** ** ** ** ** *** ** ** ** ** ** ** ** ** **',/5X,&
     &' **                                                  **',/5X,&
     &' **                                                  **',/5X,&
     &' **  CONECTIVITIES STRUCTURE:                        **',/5X,&
     &' **  OBTAINED FROM ', A30 ,                     '    **',/5X,&
     &' **                                                  **',/5X,&
     &' **                                                  **',/5X,&
     &' **                                                  **',/5X,&
     &' ** ** ** ** ** ** ** ** *** ** ** ** ** ** ** ** ** **',/5X,&
     &' ** ** ** ** ** ** ** ** *** ** ** ** ** ** ** ** ** **',/5X)
!
      END  SUBROUTINE 
!
!**************** ********* ************ *********** ********************
!
    FUNCTION DESPRESSURIZAR(NPWELL,NNP,PREINIC,AUX,BWP)
!
      IMPLICIT NONE
      INTEGER :: NPWELL,NNP
      REAL(8) :: AUX,PREINIC,BWP,DESPRESSURIZAR
      IF(NPWELL.GE.NNP)THEN
         BWP = PREINIC-REAL(NNP)*AUX
         write(*,*)'%%%%%%%%%%%%%%%%%%%%%%%%%%%%'
         write(*,*)'PRESSAO NO FUNDO DO POCO.....:',BWP
         WRITE(*,*)'PRESSAO MINIMA INICIAL.......:',PREINIC
         write(*,*)'%%%%%%%%%%%%%%%%%%%%%%%%%%%%'
      END IF
      DESPRESSURIZAR = BWP
    END FUNCTION DESPRESSURIZAR
!
!===========================================================================
!
    FUNCTION DESPRESSURIZAR_INIT(BWP,DT)
!
      use mMalha,            only : numel
      USE mPropGeoFisica,    only : GEOFORM
      use mHidroDinamicaRT,  only : pressaoElem
      use mHidroDinamicaRT,  ONLY : TPRESPRODWELL,NPRESPRODWELL,PRESMEDIAINICIAL
!
      IMPLICIT NONE
!
      REAL(8) :: BWP,DT,AUX
      REAL(8) :: DESPRESSURIZAR_INIT
      INTEGER :: NEL
!
      PRESMEDIAINICIAL = 1E30
      DO NEL=1,NUMEL
         IF(GEOFORM(NEL).EQ.'RESERVATORIO')THEN
            IF(pressaoElem(1,NEL).LE.PRESMEDIAINICIAL)THEN
               PRESMEDIAINICIAL=pressaoElem(1,NEL)
            end IF
         END IF
      END DO
      WRITE(*,*)'#############################################################'
      WRITE(*,*)'PRESSAO INICIAL NO FUNDO DO POCO DE PRODUCAO: ###############'
      WRITE(*,*)PRESMEDIAINICIAL
      WRITE(*,*)'#############################################################'
!
      NPRESPRODWELL = NINT(TPRESPRODWELL/DT)
      IF(NPRESPRODWELL.EQ.0)THEN
         AUX = 0.0
      ELSE
         AUX = (PRESMEDIAINICIAL-BWP)/REAL(NPRESPRODWELL)
      END IF
!
      WRITE(*,*)'#############################################################'
      WRITE(*,*)'TEMPO TOTAL PARA CHEGAR A PRESSAO MINIMA NO POCO DE PRODUCAO:'
      WRITE(*,*)REAL(NPRESPRODWELL)*DT
      WRITE(*,*)'#############################################################'
!
      DESPRESSURIZAR_INIT = AUX
!
    END FUNCTION DESPRESSURIZAR_INIT
