!==============================================================================
!         programa de elementos finitos  
!         baseado em: The Finite Element Method, Hughes, T. J. R., (2003)
! 
!         + implementacoes de Abimael Loula
!
!         + novo projeto: modular e fortran 90, por
!         Eduardo Garcia,        bidu@lncc.br 
!         Tuane Lopes,           tuane@lncc.br
!
!         LNCC/MCT
!         Petropolis, 07.2013
!==============================================================================
      module mLeituraEscrita
!
      implicit none
!
      integer :: iin,iecho,icoords,iconects,iconectsL
      integer :: ignuplot,iparaviewS,iparaviewP,iparaviewV
      integer :: isatTransiente
      integer :: nprint
      integer :: qtdImpSat

      integer :: ipres  = 20
      integer :: isat   = 21
      integer :: iphi   = 22
      integer :: ivel   = 23
      integer :: imass  = 24
      integer :: ifiles = 25
      integer :: imasl  = 27
      integer :: iyou   = 28
      integer :: idis   = 29
!
      integer :: ifdata = 150
      integer :: ifrand = 151
      integer :: ipress = 152
      integer :: iperm  = 153
      integer :: iporo  = 154
      integer :: iveloc = 155
      integer :: isaturacao = 156
!
!... NEW FOR GEOFORMATIONS SETUP FILES COUNTER GT 500
!
      INTEGER :: IFEDX = 501
      INTEGER :: IFNOD = 502
      INTEGER :: IFMSH = 503
!
      INTEGER :: ISTDX = 504
!
      INTEGER :: ISTYN = 505
      INTEGER :: ISTPS = 506
      INTEGER :: ISTPE = 507
      INTEGER :: ISTPH = 508
      INTEGER :: INLAY = 509
      INTEGER :: ISTSV = 510
      INTEGER :: IS3XX = 511
      INTEGER :: IS3YY = 512
      INTEGER :: IS3XY = 513
      INTEGER :: IS3ZZ = 514
!
      integer :: iflag_pres, iflag_sat, iflag_phi, iflag_vel
      integer :: iflag_disp
      integer :: iflag_mass, iflag_masl
      integer :: iflag_tipoPrint
!
      character(len=128) :: ifpres_out
      character(len=128) :: ifsat_out
      character(len=128) :: ifphi_out
      character(len=128) :: ifperm_out
      character(len=128) :: ifvel_out
      character(len=128) :: ifdis_out
      character(len=128) :: ifmass_out
      character(len=128) :: ifmasl_out
      character(len=128) :: ifyou_out
!
      real(8) :: tprt_pres, dtprt_pres
      real(8) :: tprt_sat , dtprt_sat
      real(8) :: tprt_phi , dtprt_phi
      real(8) :: tprt_vel , dtprt_vel
      real(8) :: tprt_dis , dtprt_dis
      real(8) :: tprt_mass, dtprt_mass
      real(8) :: tprt_masl, dtprt_masl
      real(8) :: TPRT_MASG, DTPRT_MASG
      real(8) :: TPRT_SIGT, DTPRT_SIGT
      real(8) :: TPRT_YOU, DTPRT_YOU
      integer :: nppres, npsat, npphi
      integer :: npvel, npmass, npmasl, NPMASG, NPSIGT, npdis
      integer :: npcontpres,npcontsat,npcontphi,npcontvel,npcontdis
!
      logical :: apenasReservatorio=.false.
!
      contains
!
!**** NEW *******************************************************************
!
      subroutine lerDataIn
!
      use mPropGeoFisica
      use mGlobaisEscalares
!
      implicit none
      character(len=128) :: flag
!
      open(unit=ifdata, file= 'data.in'  )
!
      flag="# linear(1) nao-linear(2)"
      call ireadstat(flag,ifdata,iflag_linear)
!
      flag="# viscosidade da agua"
      call readstat(flag,ifdata,xmiw)
!
      flag="# viscosidade do oleo"
      call readstat(flag,ifdata,xmio)
!
      flag="# massa especifica da agua na pressao de referencia"
      call readstat(flag,ifdata,rhow)
!
      flag="# massa especifica do oleo na pressao de referencia"
      call readstat(flag,ifdata,rhoo)
!
      FLAG="# VOLUME FORMATION"
      CALL READSTAT(FLAG,IFDATA,FORMVOL)
!
      FLAG="# WATER BULK MODULUS"
      CALL READSTAT(FLAG,IFDATA,BULKWATER)
!
      FLAG="# OIL BULK MODULUS"
      CALL READSTAT(FLAG,IFDATA,BULKOIL)
!
      flag="# saturacao residual da agua"
      call readstat(flag,ifdata,srw)
!
      flag="# saturacao residual do oleo"
      call readstat(flag,ifdata,sro)
!
      flag="# saturacao inicial da agua"
      call readstat(flag,ifdata,sinicial)
!
      flag="# saturacao na injecao"
      call readstat(flag,ifdata,sinj)
!
      flag="# saturacao do bloco"
      call readstat(flag,ifdata,sbloco)
!
      flag="# x central do bloco de sat"
      call readstat(flag,ifdata,xcbloco)
!
      flag="# y central do bloco de sat"
      call readstat(flag,ifdata,ycbloco)
!
      flag="# z central do bloco de sat"
      call readstat(flag,ifdata,zcbloco)
!
      flag="# lx do bloco de sat"
      call readstat(flag,ifdata,xlbloco)
!
      flag="# ly do bloco de sat"
      call readstat(flag,ifdata,ylbloco)
!
      flag="# lz do bloco de sat"
      call readstat(flag,ifdata,zlbloco)
!
      flag="# PERMEABILITY RESERVOIR BOTTOM-REGION"
      call readstat(flag,ifdata,perminicial)
!
      flag="# permeabilidade do bloco"
      call readstat(flag,ifdata,permbloco)
!
      flag="# x central do bloco de perm"
      call readstat(flag,ifdata,xcbloco_perm)
!
      flag="# y central do bloco de perm"
      call readstat(flag,ifdata,ycbloco_perm)
!
      flag="# z central do bloco de perm"
      call readstat(flag,ifdata,zcbloco_perm)
!
      flag="# lx do bloco de perm"
      call readstat(flag,ifdata,xlbloco_perm)
!
      flag="# ly do bloco de perm"
      call readstat(flag,ifdata,ylbloco_perm)
!
      flag="# lz do bloco de perm"
      call readstat(flag,ifdata,zlbloco_perm)
!
      flag="# porosidade inicial"
      call readstat(flag,ifdata,phiinicial)
!
      flag="# porosidade do bloco"
      call readstat(flag,ifdata,phibloco)
!
      flag="# x central do bloco de phi"
      call readstat(flag,ifdata,xcbloco_phi)
!
      flag="# y central do bloco de phi"
      call readstat(flag,ifdata,ycbloco_phi)
!
      flag="# z central do bloco de phi"
      call readstat(flag,ifdata,zcbloco_phi)
!
      flag="# lx do bloco de phi"
      call readstat(flag,ifdata,xlbloco_phi)
!
      flag="# ly do bloco de phi"
      call readstat(flag,ifdata,ylbloco_phi)
!
      flag="# lz do bloco de phi"
      call readstat(flag,ifdata,zlbloco_phi)
!
      FLAG="# PASSO PARA IMPRESSAO OPENDX-FILES"
      CALL IREADSTAT(flag,ifdata,NUMDX)
!
      flag="# impressao para matlab(0), paraview(1), paraview mrb(2) ou DX(3)"
      call ireadstat(flag,ifdata,iflag_tipoPrint) 
!
      flag="# impressao da saturacao"
      call ireadstat(flag,ifdata,iflag_sat) 
      read(ifdata,"(a)") ifsat_out
      read(ifdata,    *) npsat
! 
      flag="# impressao da pressao"
      call ireadstat(flag,ifdata,iflag_pres)
      read(ifdata,"(a)") ifpres_out
      read(ifdata,    *) nppres
! 
      flag="# impressao da velocidade"
      call ireadstat(flag,ifdata,iflag_vel)
      read(ifdata,"(a)") ifvel_out
      read(ifdata,    *) npvel
!
      flag="# impressao da porosidade"
      call ireadstat(flag,ifdata,iflag_phi)
      read(ifdata,"(a)") ifphi_out
      read(ifdata,    *) npphi
!
      flag="# impressao dos deslocamentos"
      call ireadstat(flag,ifdata,iflag_disp)
      read(ifdata,"(a)") ifdis_out
      read(ifdata,    *) npdis
!
      flag="# calculo da massa de agua"
      call ireadstat(flag,ifdata,iflag_mass)
      read(ifdata,"(a)") ifmass_out
      read(ifdata,    *) npmass
!
      flag="# balanco local da massa de agua"
      call ireadstat(flag,ifdata,iflag_masl)
      read(ifdata,"(a)") ifmasl_out
      read(ifdata,    *) npmasl
!
      !!! provisorio
      rmi=xmiw/xmio
      nsw=1
      nso=1

      if(geomech==1) then
         nsw=2
         nso=2
      endif
!
1000  FORMAT(A6)
!
      end subroutine lerDataIn
!
!===========================================================================
!
      subroutine lerRandfilesIn
!
      use mPropGeoFisica
      use mMalha, only: nsd
      use mGlobaisEscalares, only: geomech, iflag_beta
!
      implicit none

      character(len=128) :: flag
!
      open(unit=ifrand, file= 'randfiles.in'  )
!
!.... Numero de realizacoes
!
      flag="# realizacoes"
      call ireadstat(flag,ifrand,nrand)
!
!.... Leitura da permeabilidade
!
      flag="# permeabilidade"
      call ireadstat(flag,ifrand,iflag_read_perm) 
      read(ifrand,"(a)") perm_inx
      read(ifrand,"(a)") perm_iny
      if(nsd==3) read(ifrand,"(a)") perm_inz
      perm_inx=trim(perm_inx)    
      perm_iny=trim(perm_iny) 
      if(nsd==3)  perm_inz=trim(perm_inz)
      read(ifrand,*) kg, rho
!
!.... Leitura da porosidade
!
      flag="# porosidade"
      call ireadstat(flag,ifrand,iflag_read_phi) 
      read(ifrand,"(a)") phi_in
      phi_in=trim(phi_in)
      read(ifrand,*) kgphi, rhophi
!
!.... Leitura da porosidade
!
      flag="# beta"
      call ireadstat(flag,ifrand,iflag_beta) 
      read(ifrand,"(a)") beta_in
      beta_in=trim(beta_in)
      read(ifrand,*) kgbeta, rhobeta
!
!.... Leitura de Modulo de Young do Reservatorio
!
      if(geomech==1) then
         flag="# Modulo de Young Reservatorio"
         CALL IREADSTAT(FLAG,IFRAND,IFLAG_READ_YNG) 
         READ(IFRAND,"(A)") YNG_IN
         YNG_IN=trim(YNG_IN)
         READ(IFRAND,*) KGYNG, RHOYNG
      endif
!
!.... Leitura de Modulo de Young Dominio
!
      if(geomech==1) then
         flag="# Modulo de Young Dominio (Meio)"
         CALL IREADSTAT(FLAG,IFRAND,IFLAG_READ_YNG2) 
         READ(IFRAND,"(A)") YNG2_IN
         YNG2_IN=trim(YNG2_IN)
         READ(IFRAND,*) KGYNG2, RHOYNG2
         READ(IFRAND,*) YNGX1,YNGX2
         flag="# Modulo de Young Dominio (Top)"
         CALL IREADSTAT(FLAG,IFRAND,IFLAG_READ_YNG3) 
         READ(IFRAND,"(A)") YNG3_IN
         YNG3_IN=trim(YNG3_IN)
         READ(IFRAND,*) KGYNG3, RHOYNG3
         READ(IFRAND,*) YNG3X1,YNG3X2
      endif
!
!.... Saida do mixing
!
      flag="# mixing"
      call ireadstat(flag,ifrand,iflag_mix)
      read(ifrand,"(a)") mixing_out
      mixing_out=trim(mixing_out)
      read(ifrand,*)npmix
!
!.... Saida da producao
!
      flag="# producao"
      call ireadstat(flag,ifrand,iflag_prod)
      read(ifrand,"(a)") prod_out
      prod_out=trim(prod_out)
      read(ifrand,*)npprod
!
    end subroutine lerRandfilesIn
!
!===========================================================================
!
      subroutine lerGeoMechParam 
!
!      use mPropGeoFisica
!
      IMPLICIT NONE
!
      INTEGER I, NREGION, INGEOMECH
! 
      CHARACTER(LEN=128) :: FLAG, NAMEFILE, TEXT
      CHARACTER(LEN=128), DIMENSION(6) :: GEONAME
      INTEGER, DIMENSION(6) :: GEOINT
!
      GEOINT = 0 
      INGEOMECH = 520
!                  123456789+123456789+123456789+
      NAMEFILE  = 'geoMechFiles.in'
      OPEN(UNIT=INGEOMECH, FILE= NAMEFILE,STATUS='OLD')
!
!.... HEADER OF FILE 
!
      FLAG="# Files with Geological Formations"
      READ(INGEOMECH,"(A)") TEXT
      IF (TRIM(FLAG).NE.TRIM(TEXT)) THEN
           CALL CODERROR(6,NAMEFILE//'   ')
      ENDIF
!
      DO 100 I=1,6
         READ(INGEOMECH,"(A)") GEONAME(I)
         GEOINT(I) = 520+I
         write(*,*) geoint(i), geoname(i)
100   CONTINUE
!
      DO 200 I=1,6
         OPEN(UNIT=GEOINT(I), FILE=GEONAME(I), STATUS='OLD')
         CALL GEOREADER(GEOINT(I))
         CLOSE(GEOINT(I))
200   CONTINUE
!
      CLOSE(INGEOMECH)

      RETURN
!
      END SUBROUTINE
!
!**** *******************************************************************
!
      SUBROUTINE GEOREADER(IFILE)
!
      use mGlobaisEscalares, only: S3DIM
      use mPropGeoFisica,    only: YUNGVECT, POISVECT, RHODVECT
      use mPropGeoFisica,    only: GRAINBLK, PORELAGR, BULK
!
      IMPLICIT NONE
!
      INTEGER :: IFILE, IREGION
! 
      CHARACTER(LEN=128) :: FLAG, TEXT
!
      CHARACTER*18, DIMENSION(6) :: REGION 
      real(8) :: YOUNG, POISSON, GRBULK
!
      DATA   REGION(1)    ,     REGION(2)     ,      REGION(3)       /&
     & '#RESERVOIR--REGION','#RIFT-UNDER-REGION','#RIGHT-SIDE-BURDEN'/&
      &      REGION(4)    ,     REGION(5)     ,      REGION(6)       /&
     & '#SALINE-CAP-REGION','#LEFT-SIDE--BURDEN','#POST-SALT--REGION'/ 
!
      IREGION = IFILE-520
!
      READ(IFILE,1000) TEXT
!
      IF (TRIM(TEXT).NE.REGION(IREGION)) THEN 
         CALL CODERROR(6, REGION(IREGION))
      ENDIF
!
      FLAG="# YOUNG MODULUS"
      CALL READSTAT(flag,IFILE,YUNGVECT(IREGION))
      YOUNG = YUNGVECT(IREGION)
!
      FLAG="# POISSON RATIO"
      CALL READSTAT(FLAG,IFILE,POISVECT(IREGION)) 
      POISSON = POISVECT(IREGION)
!
      FLAG="# ROCK DENSITY"
      CALL READSTAT(FLAG,IFILE,RHODVECT(IREGION)) 
!
      FLAG="# BULK SOLID GRAIN"
      CALL READSTAT(FLAG,IFILE,GRAINBLK(IREGION))
      GRBULK = GRAINBLK(IREGION)
!
      FLAG="# POROSITY"
      CALL READSTAT(FLAG,IFILE,PORELAGR(IREGION))
!
      IF (GRBULK.LE.BULK(YOUNG,POISSON,S3DIM)) THEN 
         CALL CODERROR(1,REGION(IREGION))
      ENDIF

      RETURN

1000  FORMAT(A18)
!
      END SUBROUTINE
!
!===========================================================================
!
      subroutine lerNumericParam 
!
      use mPropGeoFisica,     only : TOLCREEP
      use mGlobaisEscalares,  only : NITGEO, NITHIDRO, SPLITT, IBBAR
!
      IMPLICIT NONE
!
      INTEGER I, NUMPARAM
! 
      CHARACTER(LEN=128) :: FLAG, NAMEFILE, TEXT
!
      NUMPARAM = 520
!                  123456789+123456789+123456789+
      NAMEFILE  = 'NumericParam.in'
      OPEN(UNIT=NUMPARAM, FILE= NAMEFILE, STATUS='OLD')
!
!.... HEADER OF FILE 
!
      FLAG="# Numerical Parameters"
      READ(NUMPARAM,"(A)") TEXT
      IF (TRIM(FLAG).NE.TRIM(TEXT)) THEN
           CALL CODERROR(6,NAMEFILE//'   ')
      ENDIF
!
      FLAG="# ITERATIONS FOR TRANSPORT<-->GEOMEC"
      CALL IREADSTAT(FLAG,NUMPARAM,NITGEO)
!
      FLAG="# ITERATIONS FOR HIDRODI <--> TRANSP"
      CALL IREADSTAT(FLAG,NUMPARAM,NITHIDRO)
!
      FLAG="# SPLITTING PROCEDURE"
      READ(NUMPARAM,"(A)") TEXT
      IF (TRIM(FLAG).NE.TRIM(TEXT)) THEN
           CALL CODERROR(6,NAMEFILE//'   ')
      ENDIF
      READ(NUMPARAM,1000) SPLITT
!
      FLAG="# TOLERANCE FOR ALL ITERATIVE" 
      CALL READSTAT(FLAG,NUMPARAM,TOLCREEP)
!
      FLAG="# BBAR METHOD (NOT=0/YES=1)"
      CALL IREADSTAT(FLAG,NUMPARAM,IBBAR)
!
      CLOSE(NUMPARAM)
!
      RETURN
!
1000  FORMAT(A6)
!
      END SUBROUTINE
!
!===========================================================================
!
      subroutine lerCreepParam
!
      use mPropGeoFisica,     only: DTCREEP, CREEPZERO, POWERN, SIGMAREF
      use mGlobaisEscalares,  only: NCREEP
!
      IMPLICIT NONE
!
      INTEGER I, ICREEP
! 
      CHARACTER(LEN=128) :: FLAG, NAMEFILE, TEXT
!
      ICREEP    = 520
!                  
      NAMEFILE  = 'CreepParameters.in'
      OPEN(UNIT=ICREEP, FILE= NAMEFILE, STATUS='OLD')
!
!.... HEADER OF FILE 
!           123456789+123456789+123456789+12
      FLAG="# Creep Data Parameters and Info"
      READ(ICREEP,"(A)") TEXT
      IF (TRIM(FLAG).NE.TRIM(TEXT)) THEN
           CALL CODERROR(6,NAMEFILE//'   ')
      ENDIF
!
      FLAG="# CREEP DEFORMATION REFERENCE" 
      CALL READSTAT(FLAG,ICREEP,CREEPZERO)
!
      FLAG="# TIME INCREMENT FOR CREEP" 
      CALL READSTAT(FLAG,ICREEP,DTCREEP)
!
      FLAG="# STRESS REFERENCE FOR CREEP" 
      CALL READSTAT(FLAG,ICREEP,SIGMAREF)
!
      FLAG="# POWER OF CREEP LAW" 
      CALL READSTAT(FLAG,ICREEP,POWERN)
!
      FLAG="# CREEP TIME AS MULTIPLE OF GEOTIME"
      CALL IREADSTAT(flag,ICREEP,NCREEP)
!
      CLOSE(ICREEP)
!
      RETURN
!
      END SUBROUTINE
!
!===========================================================================
!
      subroutine lerSimulatorParam 
!
      use mPropGeoFisica
      use mGlobaisEscalares
      use mMCMC,             only: PRESPROD,TPRESPRODWELL
      use mMCMC,             only: ZFUNDOPOCCO
!
      IMPLICIT NONE
!
      INTEGER I, ISimulator
! 
      CHARACTER(LEN=128) :: FLAG, NAMEFILE, TEXT, FLAG1
!
      ISimulator  = 520
!                    123456789+123456789+12
      NAMEFILE    = 'SimulatorParam.in'
      OPEN(UNIT=ISimulator, FILE= NAMEFILE, STATUS='OLD')
!
!.... HEADER OF FILE 
!           123456789+123456789+123456789+12
      FLAG="# Simulator Parameters"
      READ(ISimulator,"(A)") TEXT
      IF (TRIM(FLAG).NE.TRIM(TEXT)) THEN
           CALL CODERROR(6,NAMEFILE//'   ')
      ENDIF
!
      flag="# SEA DEPTH IN METERS"
      call readstat(flag,ISimulator,SEADEPTH)
!
      SEALOAD = SEADEPTH*RHOW*10.0D0
!
      flag="# instante inicial"
      call readstat(flag,ISimulator,tzero)
      tTransporte = tzero
!
      flag="# tempo total de simulacao"
      call readstat(flag,ISimulator,tt)
!
      flag="# numero de calculos da velocidade"
      call ireadstat(flag,ISimulator,nvel)
!
      FLAG="# TEMPO DE INICIO DE INJECAO"
      CALL readstat(flag,ISimulator,YEARINJ)
!
      WRITE(*,*)'########################################'
      FLAG='# PRESSAO NO FUNCO DO POCCO DE PRODUCAO'
      READ(ISimulator,"(a)")FLAG1
      if(trim(flag).eq.trim(flag1))THEN
         READ(ISimulator,*)PRESPROD
         WRITE(*,*)'PRESSAO NO POCO DE PRODUCAO:',PRESPROD
      ELSE
         WRITE(*,*) 'ERRO NA LEITURA DA PRESSAO DO POCCO'
         STOP
      END IF
      WRITE(*,*)'########################################'
!
      FLAG='# COTA DO FUNDO DO POCCO'
      READ(ISimulator,"(a)")FLAG1
      if(trim(flag).eq.trim(flag1))THEN
         READ(ISimulator,*)ZFUNDOPOCCO
         WRITE(*,*)'COTA DO FUNDO DO POCCO:',ZFUNDOPOCCO
      ELSE
         WRITE(*,*) 'ERRO NA LEITURA DA COTA DO FUNDO DO POCCO'
         STOP
      END IF
      WRITE(*,*)'########################################'
!
      FLAG='# TEMPO PARA A DESCOMPRESSAO DO POCCO'
      READ(ISimulator,"(a)")FLAG1
      if(trim(flag).eq.trim(flag1))THEN
         READ(ISimulator,*)TPRESPRODWELL
         WRITE(*,*)'TEMPO PARA DESCOMPRESSAO DO POCCO:',TPRESPRODWELL
      ELSE
         WRITE(*,*) 'ERRO NA LEITURA DO TEMPO PARA ADESCOMPRESSAO'
         STOP
      END IF
      WRITE(*,*)'########################################'
!
      FLAG='# Pressao de referencia'
      READ(ISimulator,"(a)")FLAG1
      if(trim(flag).eq.trim(flag1))THEN
         READ(ISimulator,*)PRESSAOREF
         WRITE(*,*)'PRESSAO DE REFERENCIA:',PRESSAOREF
      ELSE
         WRITE(*,*) 'ERRO NA LEITURA DA PRESSAO DE REFERENCIA'
         STOP
      END IF
      WRITE(*,*)'########################################'
!
      FLAG='# Cota de referencia'
      READ(ISimulator,"(a)")FLAG1
      if(trim(flag).eq.trim(flag1))THEN
         READ(ISimulator,*)COTAREF
         WRITE(*,*)'COTA DE REFERENCIA:',COTAREF
      ELSE
         WRITE(*,*) 'ERRO NA LEITURA DA COTA DE REFERENCIA'
         STOP
      END IF
      WRITE(*,*)'########################################'
!
      CLOSE(ISimulator)
!
      RETURN
!
      END SUBROUTINE
!
!**** new ***************************************************************
!
      subroutine readstat(flag,ifile,x)
!
      implicit none
      integer :: ifile
      real(8) :: x
      character(len=128) :: flag,flag1
!
      read(ifile,"(a)") flag1
      if(trim(flag).eq.trim(flag1)) then
      read(ifile,*) x
      else
      write(*,*) "Erro na leitura de ", flag
      stop
      end if
!
      end subroutine
!
!=========================================================================
!
      subroutine ireadstat(flag,ifile,n)
!
      implicit none

      integer, intent(in)  :: ifile
      integer, intent(out) :: n
      character(len=128), intent(in) :: flag
      character(len=128) :: flag1
!
      read(ifile,"(a)") flag1
      if (trim(flag).eq.trim(flag1)) then
         read(ifile,*) n
      else
         write(*,*) "Erro na leitura de ", flag
      stop
      end if
!
      end subroutine
!
!==========================================================================
!
      subroutine abrirArquivosInformacoesMalha()
!  
!        iin    = input unit number
!        iecho  = output unit of input data
!        iouter  = output unit of error norms
!
      character(len=2) :: nRank
      character(len=20) :: nomeIecho
!
      iin           = 11
      iecho         = 12
      icoords       = 13
      iconects      = 14
      iconectsL     = 15
      isatTransiente= 16
!
      open(unit=iin   , file= 'input.in')
!
      nomeIecho='echo.dat'
      open(unit=iecho , file= nomeIecho)
!
#ifdef debug
      open(unit=icoords    , file= 'coordenadas.dat')
      open(unit=iconects   , file= 'conectsNodais.dat')
      open(unit=iconectsL  , file= 'conectsLadais.dat')
#endif

   end subroutine 
!
!========================================================================
!
      subroutine abrirArquivosResultados
      use mPropGeoFisica
      use mGlobaisEscalares
!
      implicit none
!
      character(128) :: aux

      if(iflag_sat==1)then
         if(iflag_tipoPrint==0) then
            ifsat_out=trim(ifsat_out)//'sat_amostra.res'
            open(unit=isat,file=ifsat_out,status='unknown')
         end if
         if(iflag_tipoPrint==1) then         
            ifsat_out=trim(ifsat_out)//'resultadoSat.vtk'
            open(unit=isat,file=ifsat_out,status='unknown')
         end if
         if(iflag_tipoPrint==2) then         
            ifsat_out=trim(ifsat_out)//'resultadoSat-'
         end if
      end if

      if(iflag_pres==1)then
         if(iflag_tipoPrint==0) then
            ifpres_out=trim(ifpres_out)//'pres_amostra.res'
            open(unit=ipres,file=ifpres_out,status='unknown')
         end if
         if(iflag_tipoPrint==1) then
            ifpres_out=trim(ifpres_out)//'resultadoPressao.vtk'
            open(unit=ipres,file=ifpres_out,status='unknown')
         end if
         if(iflag_tipoPrint==2) then
            ifpres_out=trim(ifpres_out)//'resultadoPressao-'
         end if
      end if

      if(iflag_vel==1)then
         if(iflag_tipoPrint==0) then
            ifvel_out=trim(ifvel_out)//'vel_amostra.res'
            open(unit=ivel,file=ifvel_out,status='unknown')
         end if
         if(iflag_tipoPrint==1) then
            ifvel_out=trim(ifvel_out)//'resultadoVel.vtk'
            open(unit=ivel,file=ifvel_out,status='unknown')
         end if
         if(iflag_tipoPrint==2) then
            ifvel_out=trim(ifvel_out)//'resultadoVel-'
         end if
      end if

      if(iflag_disp==1)then
         if(iflag_tipoPrint==0) then
            ifdis_out=trim(ifdis_out)//'dis_amostra.res'
            open(unit=idis,file=ifdis_out,status='unknown')
         end if
         if(iflag_tipoPrint==1) then
            ifdis_out=trim(ifdis_out)//'resultadoDis.vtk'
            open(unit=idis,file=ifdis_out,status='unknown')
         end if
         if(iflag_tipoPrint==2) then
            ifdis_out=trim(ifdis_out)//'resultadoDis-'
         end if
      end if

      if(iflag_phi==1)then
         if(iflag_tipoPrint==0) then
            ifphi_out=trim(ifphi_out)//'phi_amostra.res'
            ifperm_out=trim(ifphi_out)//'perm_amostra.res'
            open(unit=iphi, file=ifphi_out, status='unknown')
            open(unit=iperm,file=ifperm_out,status='unknown')
         end if
         if(iflag_tipoPrint==1) then
            ifphi_out= trim(ifphi_out)//'resultadoPhi.vtk'
            ifperm_out=trim(ifphi_out)//'resultadoPerm.vtk'
            open(unit=iphi, file=ifphi_out, status='unknown')
            open(unit=iperm,file=ifperm_out,status='unknown')
         end if
         if(iflag_tipoPrint==2) then
            aux = ifphi_out
            ifphi_out= trim(ifphi_out)//'resultadoPhi-'
            ifyou_out= trim(aux      )//'resultadoYoung-'
            ifperm_out=trim(aux      )//'resultadoPerm-'
         end if
      end if

     if(iflag_mass==1)then
         ifmass_out=trim(ifmass_out)//'mass_amostra'
         open(unit=imass,file=ifmass_out,status='unknown')
      end if

      if(iflag_masl==1)then
         ifmasl_out=trim(ifmasl_out)//'bal_mass_amostra.res'
         open(unit=imasl,file=ifmasl_out,status='unknown')
      end if
!
      end subroutine
!
!**** new ********************************************************************
!
      subroutine fecharArquivos()
!
      close(iin   )
      close(iecho )
      close(icoords )
      close(ignuplot)
      close(iconectsL)
      close(iconects )
      close(isat)
      close(ipres)
      close(ivel)
      close(iphi)
      close(iperm)
      close(imass)
      close(imasl)
! 
      end subroutine fecharArquivos
!      
! =======================================================================
!   
       subroutine inittime
! 
       use mPropGeoFisica
       use mGlobaisEscalares, only: tt, tzero, NVEL

       implicit none
!
       real(8) :: tempo,TOL=1e-8
! 
! .... ajustes para impressao estocastico
! 
!        t0 = tzero
! 
       npcontpres = 0
       npcontsat  = 0
       npcontphi  = 0
       npcontvel  = 0
       npcontdis  = 0
!
       np_rand_mix  = 1
       np_rand_prod = 1
       np_rand_conc = 1
       np_rand_prodF= 1
       ninit_prodF  = 1
! 
! .... tamanhos dos intervalos de impressoes
! 
       dtprt_pres = (tt-tzero)/nppres
       dtprt_sat  = (tt-tzero)/npsat
       dtprt_phi  = (tt-tzero)/npphi
       dtprt_vel  = (tt-tzero)/npvel
       dtprt_dis  = (tt-tzero)/npdis
       dtprt_mass = (tt-tzero)/npmass
       dtprt_masl = (tt-tzero)/npmasl
       dtprt_prod = (tt-tzero)/npprod
       dtprt_mix  = (tt-tzero)/npmix
       dtprt_prodF= (tt-tzero)/npprodF
       dtprt_conc = (tt-tzero)/npconc
! 
! .... inicializa os contadores de impressao
! 
       if(npsat.ge.nvel)then
          npsat = nvel
	  dtprt_sat  = (tt-tzero)/npsat
       else
          tempo=1.0
          do while(tempo.gt.TOL)
             npsat=npsat+1
	     dtprt_sat  = (tt-tzero)/npsat
             tempo = dmod(dtprt_sat,(tt-tzero)/real(nvel))
          end do
      endif
      write(*,*)'###################################'
      write(*,*)'NUMERO de IMPRESSOES SAT ',npsat
      write(*,*)'Delta t IMPRESSAO        ',dtprt_sat
 ! 
       if(npphi.ge.nvel)then
          npphi = nvel
          dtprt_phi  = (tt-tzero)/npphi
       else
          tempo=1.0
          do while(tempo.gt.TOL)
             npphi=npphi+1
	     dtprt_phi  = (tt-tzero)/npphi
             tempo = dmod(dtprt_phi,(tt-tzero)/real(nvel))
          end do
       endif
       write(*,*)'###################################'
       write(*,*)'NUMERO de IMPRESSOES PHI ',npphi
       write(*,*)'Delta t IMPRESSAO        ',dtprt_phi
! 
       if(nppres.ge.nvel)then
          nppres = nvel
	  dtprt_pres  = (tt-tzero)/nppres
       else
          tempo=1.0
          do while(tempo.gt.TOL)
             nppres=nppres+1
	     dtprt_pres  = (tt-tzero)/nppres
             tempo = dmod(dtprt_pres,(tt-tzero)/real(nvel))
          end do
       endif
       write(*,*)'###################################'
       write(*,*)'NUMERO de IMPRESSOES PRES',nppres
       write(*,*)'Delta t IMPRESSAO        ',dtprt_pres
! 
       if(npvel.ge.nvel)then
          npvel = nvel
          dtprt_vel  = (tt-tzero)/npvel
       else
          tempo=1.0
          do while(tempo.gt.TOL)
             npvel=npvel+1
	     dtprt_vel  = (tt-tzero)/npvel
             tempo = dmod(dtprt_vel,(tt-tzero)/real(nvel))
          end do
       endif
       write(*,*)'###################################'
       write(*,*)'NUMERO de IMPRESSOES VEL ',npvel
       write(*,*)'Delta t IMPRESSAO        ',dtprt_vel
! 
       if(npdis.ge.nvel)then
          npdis = nvel
          dtprt_dis  = (tt-tzero)/npdis
       else
          tempo=1.0
          do while(tempo.gt.TOL)
             npdis=npdis+1
	     dtprt_dis  = (tt-tzero)/npdis
             tempo = dmod(dtprt_dis,(tt-tzero)/real(nvel))
          end do
       endif
       write(*,*)'###################################'
       write(*,*)'NUMERO de IMPRESSOES DIS ',npdis
       write(*,*)'Delta t IMPRESSAO        ',dtprt_dis
! 
       if(npmass.ge.nvel)then
          npmass = nvel
          dtprt_mass  = (tt-tzero)/npmass
       else
          tempo=1.0
          do while(tempo.gt.TOL)
             npmass=npmass+1
	     dtprt_mass  = (tt-tzero)/npmass
             tempo = dmod(dtprt_mass,(tt-tzero)/real(nvel))
          end do
       endif
       write(*,*)'###################################'
       write(*,*)'NUMERO de IMPRESSOES MASSA ',npmass
       write(*,*)'Delta t IMPRESSAO          ',dtprt_mass
       write(*,*)'###################################'
!! 
       tprt_pres= dtprt_pres
       tprt_sat = dtprt_sat
       tprt_phi = dtprt_phi
       tprt_vel = dtprt_vel
       tprt_dis = dtprt_dis
       tprt_mass= dtprt_mass
       tprt_masl= dtprt_masl
       tprt_prod= dtprt_prod
       tprt_mix = dtprt_mix
       tprt_prodF= dtprt_prodF
       tprt_conc= dtprt_conc
! 
       end subroutine
!
!*****************************************************
!
      subroutine imprimirCondicoesIniciais(pressaoElem, velocLadal, phi, perm, satElem, YOUNG, DIS, PORE)
      use mGlobaisEscalares, only: ndofP, ndofV, ndofD, tTransporte, numdx
      use mMalha,            only: x, conecNodaisElem, conecLadaisElem, nen, nsd
      use mMalha,            only: numel, numnp, numLadosElem
      use mMalha,            only: numelReserv, numLadosReserv
!
      implicit none
!
      real*8, intent(in) :: pressaoElem(ndofP, numelReserv)
      real*8, intent(in) :: velocLadal(ndofV,numLadosReserv)
      real*8, intent(in) :: phi(numelReserv), perm(numelReserv), satElem(numelReserv)
      real*8, intent(in) :: YOUNG(numel), DIS(ndofd,NUMNP), PORE(numelReserv)
      integer :: i, k, DZERO = 0
      LOGICAL :: SIM,NAO
      real*8  :: ZERO=0.D0
!
      SIM=.TRUE.
      NAO=.FALSE.
!
!.... imprime a condicao inicial da massa de agua
!
       if(iflag_mass==1)then
          call MONITORMASSATOTAL(satElem,phi,zero)
       end if
!
!.... imprime a condicao inicial: pressao
!
       if(iflag_pres==1) then
          if(iflag_tipoPrint==0) then
             call prt(nsd,numelReserv,tTransporte,pressaoElem,ipres)
          end if
!
          if(iflag_tipoPrint==1) then
             if(apenasReservatorio.eqv..true.) then
                call escreverArqParaviewReservatorio(ipres, pressaoElem, ndofP, &
                     numelReserv, nen, conecNodaisElem, 1, 't=0.0', len('t=0.0')) 
             else
                call escreverArqParaview(ipres, pressaoElem, ndofP, numel, nen, &
                     conecNodaisElem, 1, 't=0.0', len('t=0.0')) 
             endif
          endif
!
          if(iflag_tipoPrint==2) then
             call escreverArqParaview_escalar(ipres,pressaoElem,ZERO,ifpres_out,nen, &
                  NAO,DZERO,'PRES')
          end if
       endif
!  
!.... imprime a condicao inicial: velocidade
!
      if(iflag_vel==1) then
         if(iflag_tipoPrint==0) then
            call prtvB(nsd,numelReserv,tTransporte,velocLadal,ndofV, conecLadaisElem,&
                                           numLadosElem,ivel) 
         end if
         if(iflag_tipoPrint==1) then
            write(ivel,*) "impressao da velocidade ladal nao implementada para o paraview"
         endif
         if(iflag_tipoPrint==2) then
            write(ivel,*) "impressao da velocidade ladal nao implementada para o paraview2"
         endif         
      end if
!  
!.... imprime a condicao inicial: deslocamentos
!
      if(iflag_disp==1) then
         if(iflag_tipoPrint==0) then
            write(idis,*) "impressao dos deslocamentos nao implementada para o matlab"
         end if
         if(iflag_tipoPrint==1) then
            write(idis,*) "impressao dos deslocamentos implementada para o paraview"
         endif
         if(iflag_tipoPrint==2) then
            call escreverArqParaview_vetor(idis,DIS,ZERO,ifdis_out,nen,SIM,DZERO,'DISP',ndofD)
         endif         
      end if
!     
!.... imprime a condicao inicial: saturacao
!
      if(iflag_sat==1) then    
         if(iflag_tipoPrint==0) then 
            call prt (nsd,numelReserv,tTransporte,satElem,isat)
         end if
         if(iflag_tipoPrint==1) then 
            if(apenasReservatorio.eqv..true.) then
               call escreverArqParaviewReservatorio(isat, satElem, ndofV,&
                     numelReserv, nen, conecNodaisElem, 1, 't=0.0', len('t=0.0')) 
            else
               call escreverArqParaview(isat, satElem, ndofV, numel, nen, &
                                  conecNodaisElem, 1, 't=0.0', len('t=0.0')) 
            endif
         endif
         if(iflag_tipoPrint==2) then 
             call escreverArqParaview_escalar(isat,satElem,ZERO,ifsat_out,nen, &
                  NAO,DZERO,' SAT')
          end if
      endif
!
!.... imprime a condicao inicial: porosidade
!
       if(iflag_phi==1) then
          if(iflag_tipoPrint==0) then
             call prt(nsd,numelReserv,tTransporte,phi,iphi)
             call prt(nsd,numelReserv,tTransporte,perm,iperm)
          end if
          if(iflag_tipoPrint==1) then
             if(apenasReservatorio.eqv..true.) then
                call escreverArqParaviewReservatorio(iphi, phi, ndofV, numelReserv, nen, & 
                    conecNodaisElem, 1, 'porosidade', len('porosidade')) 
             else
                call escreverArqParaview(iphi, phi, ndofV, numel, nen, & 
                    conecNodaisElem, 1, 'porosidade', len('porosidade')) 
             endif
             if(apenasReservatorio.eqv..true.) then
                call escreverArqParaviewReservatorio(iperm, perm, ndofV, numelReserv, nen, &
                     conecNodaisElem, 1, 'permeabilidade', len('permeabilidade'))
             else
                call escreverArqParaview(iperm, perm, ndofV, numel, nen, &
                     conecNodaisElem, 1, 'permeabilidade', len('permeabilidade'))
             endif
          endif
          if(iflag_tipoPrint==2) then
             call escreverArqParaview_escalar(iphi,phi,ZERO,ifphi_out,nen, &
                  NAO,DZERO,'PORE')
             call escreverArqParaview_escalar(iperm,perm,ZERO,ifperm_out,nen, &
                  NAO,DZERO,'PERM')
             call escreverArqParaview_escalar(iyou,YOUNG,ZERO,ifyou_out,nen, &
                  SIM,DZERO,'YOUN')
          end if
       endif
!
       RETURN
!
     END SUBROUTINE imprimirCondicoesIniciais
!
!*****************************************************
!
      subroutine imprimirSolucaoNoTempo(sat,DIS,PORE,YOUNG,pressaoElem, velocLadal, tempo)
       use mGlobaisEscalares, only : ndofP, ndofV, ndofD
       use mMalha,            only : nsd, numel, numelReserv, numLadosReserv, numLadosElem, conecLadaisElem
       use mMalha,            only : nen, numnp
       use MMCMC
!
      implicit none
!
      real*8, intent(in) :: pressaoElem(ndofP,numelReserv), velocLadal(ndofV,numLadosReserv)
      real*8, intent(in) :: sat(numelReserv),DIS(ndofD,numnp),PORE(numelReserv),YOUNG(numel)
      real*8, intent(in) :: tempo
!
      character(21) :: labelTransp
      character(21) :: label
      character(21) :: num
      real*8        :: TOL=1e-6
      integer       :: UM=1
      LOGICAL       :: SIM,NAO
! MCMC !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
      SIM=.TRUE.
      NAO=.FALSE.
! 
!
!.... imprime a condicao inicial da massa de agua
!
       IF(IFLAG_MASS==1)THEN
          IF(ABS(tempo-TPRT_MASS).LE.TOL)THEN
             TPRT_MASS = TPRT_MASS+DTPRT_MASS
             call MONITORMASSATOTAL(sat,PORE,tempo)
          END IF
       end if
!
!
      IF(IFLAG_CONC)THEN
         IF(ABS(tempo-TPRT_CONC).LE.TOL)THEN
	   TPRT_CONC = TPRT_CONC+DTPRT_CONC
	   CALL MONITOR(CONC_OUT,PCONDC,NCONDC,ELEM_CONDC,sat,tempo)
         ENDIF
      ENDIF
      IF(IFLAG_PRODF)THEN
         IF(ABS(tempo-TPRT_PRODF).LE.TOL)THEN
	   TPRT_PRODF = TPRT_PRODF+DTPRT_PRODF
	   CALL MONITORPROD(PRODF_OUT,PCONDP,NCONDP,ELEM_CONDP,velocLadal,sat,tempo)
         ENDIF
      ENDIF
      IF(IFLAG_PRESF)THEN
         IF(ABS(tempo-TPRT_PRESF).LE.TOL)THEN
	   TPRT_PRESF = TPRT_PRESF+DTPRT_PRESF
	   CALL MONITOR(PRESF_OUT,PCONDPR,NCONDPR,ELEM_CONDPR,pressaoElem,tempo)
         ENDIF
      ENDIF
      IF(IFLAG_DISF)THEN
         IF(ABS(tempo-TPRT_DISF).LE.TOL)THEN
	   TPRT_DISF = TPRT_DISF+DTPRT_DISF
	   CALL MONITORV(DISF_OUT,PCONDD,NCONDD,ELEM_CONDD,DIS,tempo)
         ENDIF
      ENDIF
! MCMC !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!.... imprime a velocidade no centro
!
      if(iflag_vel==1) then
         if(iflag_tipoPrint==0) then
            call prtvB(nsd,numelReserv,tempo,velocLadal,ndofV, conecLadaisElem, numLadosElem,ivel) 
         endif
         if(iflag_tipoPrint==1) then
            write(ivel,*) "impressao da velocidade ladal nao implementada para o paraview"
         endif
         if(iflag_tipoPrint==2) then
            if(abs(tempo-tprt_vel).le.TOL)then
               tprt_vel=tprt_vel+dtprt_vel
	       npcontvel=npcontvel+1
               write(ivel,*) "impressao da velocidade ladal nao implementada para o paraview2"
            endif
         endif  
      end if
!
!.... imprime os deslocamentos
!
      if(iflag_disp==1) then
         if(iflag_tipoPrint==0) then
            write(idis,*) "impressao dos deslocamentos nao implementada para o matlab"
         endif
         if(iflag_tipoPrint==1) then
            write(idis,*) "impressao da velocidade ladal nao implementada para o paraview"
         endif
         if(iflag_tipoPrint==2) then
            if(abs(tempo-tprt_dis).le.TOL)then
               tprt_dis=tprt_dis+dtprt_dis
	       npcontdis=npcontdis+1
               call escreverArqParaview_vetor(idis,DIS,tempo,ifdis_out,nen,SIM,npcontdis,'DISP',ndofD)
            endif
         endif  
      end if
!
!.... imprime a saturacao no tempo
      if (iflag_sat==1) then
         if(iflag_tipoPrint==0) then
            call prt(nsd,numelReserv,tempo,sat,isat)
         end if
         if (iflag_tipoPrint==1) then
            call gerarLabel(labelTransp,tempo)
            call escreverArqParaviewIntermed(isat, sat, ndofV, &
                 numel, trim(labelTransp), len(trim(labelTransp)))
         endif
         if(iflag_tipoPrint==2) then
            if(abs(tempo-tprt_sat).le.TOL)then
               tprt_sat=tprt_sat+dtprt_sat
	       npcontsat=npcontsat+1
               call escreverArqParaview_escalar(isat,sat,tempo,ifsat_out,nen,NAO,npcontsat,' SAT')
            end if
         endif
      endif
!
!.... imprime a pressao no tempo 
!
      if(iflag_pres==1) then
         if(iflag_tipoPrint==0) then
            call prt(nsd,numelReserv,tempo,pressaoElem,ipres)
         end if
         if(iflag_tipoPrint==1) then
            write(num,'(f12.5)') tempo
            label="t="//ADJUSTL(num)
            call escreverArqParaviewIntermed(ipres, pressaoElem, ndofP, numelReserv, trim(label), len(trim(label)))
         end if
         if(iflag_tipoPrint==2) then
            if(abs(tempo-tprt_pres).le.TOL)then
               tprt_pres=tprt_pres+dtprt_pres
	       npcontpres=npcontpres+1
               call escreverArqParaview_escalar(ipres,pressaoElem,tempo,ifpres_out,nen,NAO,npcontpres,'PRES')
            end if
         endif
      end if
!
!.... imprime a porosidade no tempo 
!
      if(iflag_phi==1) then
         if(iflag_tipoPrint==0) then
            call prt(nsd,numelReserv,tempo,PORE,iphi)
         end if
         if(iflag_tipoPrint==1) then
            write(num,'(f12.5)') tempo
            label="t="//ADJUSTL(num)
            call escreverArqParaviewIntermed(iphi, PORE, ndofP, numelReserv, trim(label), len(trim(label)))
         end if
         if(iflag_tipoPrint==2) then
            if(abs(tempo-tprt_phi).le.TOL)then
               tprt_phi=tprt_phi+dtprt_phi
	       npcontphi=npcontphi+1
               call escreverArqParaview_escalar(iphi,PORE,tempo, &
                    ifphi_out,nen,NAO,npcontphi,'PORE')
!               call escreverArqParaview_escalar(iyou,YOUNG,tempo, &
!                    ifyou_out,nen,SIM,npcontphi,'YOUN')
            end if
         endif
      end if
!
      RETURN
!
      END SUBROUTINE 
!
!**** new ******************************************************************
!
      subroutine echo
!
      implicit none
!
!.... program to echo input data
!
      character*4 ia(20)
      integer :: iech, i
!
!     cabeçalho
      write(iecho,500)

      read(iin,1000) iech
      if (iech.eq.0) return
!
      write(iecho,2000) iech
      backspace iin
!
      do 100 i=1,100000
         read(iin,3000,end=200) ia
         if (mod(i,50).eq.1) write(iecho,4000)
         write(iecho,5000) ia
  100 continue
!
  200 continue
      rewind iin
      read(iin,1000) iech
!
      return
!
 500  format('programa de elementos finitos em fortran 90 baseado em:',// &
      'The Finite Element Method, Hughes, T. J. R., (2003)'//)
 1000 format(16i10)
 2000 format('1',' i n p u t   d a t a   f i l e               ',  //5x,&
     ' echo print code . . . . . . . . . . . . . . (iecho ) = ',i10//5x,&
     '    eq. 0, no echo of input data                        ',   /5x,&
     '    eq. 1, echo input data                              ',   ///)
 3000 format(20a4)
 4000 format(' ',8('123456789*'),//)
 5000 format(' ',20a4)

      end subroutine echo
!
!**** new *******************************************************************
!
      subroutine leituraGeracaoCoordenadas_OLD(x, nsd, numnp, iin, icoords, iprtin)
      use mMalha, only: genfl
!
!.... program to read, generate and write coordinate data
!
      implicit none
!
!.... remove above card for single-precision operation
!
      real*8, intent(inout) ::  x(nsd,*)
      integer, intent(in)   :: nsd, numnp, iin, icoords, iprtin
!
      integer :: i, n
!      
      call genfl(x,nsd,iin)
!
      if (iprtin.eq.1) return
!
#ifdef debug
         write(icoords,*) "# Coordenadas ", nsd
         do n=1,numnp
            write(icoords,2000) n,(x(i,n),i=1,nsd) 
         end do
#endif
!
      return
!
!  1000 format('1',' n o d a l   c o o r d i n a t e   d a t a '///5x,&
!      ' node no.',3(13x,' x',i1,' ',:)//)
 2000 format(6x,i12,10x,3(1pe15.8,2x))
      end subroutine
!
!**** NEW ** MODIFIED FOR IRREGULAR MESH ***********************************
!
      SUBRoutine LEITuraGERAcaoCOORdenadas(x, nsd, numnp, iin, icoords, iprtin)
      use mMalha, only: IrregMesh, genfl
      use mMalha, only: POSTLINE, SALTLINE, RTOPLINE, RBTTLINE, RIFTLINE
      use mMalha, only: LEFTLINE, RGHTLINE, IGEOFORM, IDome
!
!.... program to read, generate and write coordinate data
!
      implicit none
!
      INTEGER, INTENT(IN)    :: NSD, NUMNP, IIN, ICOORDS, IPRTIN
      REAL(8), INTENT(INOUT) ::  X(NSD,NUMNP)
!
      INTEGER      :: INXCORD
      CHARACTER*30 :: NAMEIN
      CHARACTER(LEN=128) :: FLAG, ROTULO
!
      INTEGER :: I, N, NODE
      REAL(8) :: XLEFT,XRIGHT,YBOTTOM,YTOP
!
      call genfl(x,nsd,iin)
!
      IF (IrregMesh.EQ.0) THEN
            WRITE(*,5000) 'REGULAR MESH  : input.in'
         ELSE
            INXCORD  = 530
            NAMEIN   = 'nodes_rsrv_outsburden.in'
            WRITE(*,5000) NAMEIN 
            OPEN(UNIT=INXCORD,FILE=NAMEIN,STATUS='OLD')
            DO 100 NODE=1,NUMNP
               READ(INXCORD,1500) X(1,NODE), X(2,NODE)
!              WRITE(*,1500) X(1,NODE), X(2,NODE)
100         CONTINUE
            CLOSE(INXCORD)
      ENDIF
!
      XLEFT   = 0.0D0
      YTOP    = 0.0D0
      XRIGHT  = 0.0D0
      YBOTTOM = 0.0D0
!
      DO 200 NODE=1,NUMNP
         XLEFT   = DMIN1(X(1,NODE),XLEFT)
         XRIGHT  = DMAX1(X(1,NODE),XRIGHT)
         YTOP    = DMAX1(X(2,NODE),YTOP)
         YBOTTOM = DMIN1(X(2,NODE),YBOTTOM)
200   CONTINUE
!
      INXCORD  = 531
!                 123456789+123456789+123456789+
      NAMEIN   = 'layer_reference.in'
      OPEN(UNIT=INXCORD,FILE=NAMEIN,STATUS='OLD')
!
!.... HEADER OF FILE 
!
      FLAG = '#REFERENCE LINES OF GEO-FORMATIONS'
      READ(INXCORD,'(A)') ROTULO
      IF (TRIM(FLAG).NE.TRIM(ROTULO)) THEN
           CALL CODERROR(6,NAMEIN)
      ENDIF
!
!            123456789+123456789+123456789+123456789+123456789+123456789+12345
      FLAG ='#READ GEOFORMATIONS FROM elmnt_geoformation.in FILE (NOT=0/YES=1)'
      CALL IREADSTAT(FLAG,INXCORD,IGEOFORM)
!
      FLAG = '#TOP POST-SALT HORIZONTAL Y-AXIS'
      CALL READSTAT(FLAG,INXCORD,POSTLINE)
!      READ(INXCORD,3700) POSTLINE 
      FLAG = '#SALT-DOME TOP HORIZONTAL Y-AXIS'
      CALL READSTAT(FLAG,INXCORD,SALTLINE)
!      READ(INXCORD,3700) SALTLINE 
      FLAG = '#TOP RESERVOIR HORIZONTAL Y-AXIS'
      CALL READSTAT(FLAG,INXCORD,RTOPLINE)
!      READ(INXCORD,3700) RTOPLINE 
      FLAG = '#BOTTOM RESERVOIR HORIZONTAL Y-AXIS'
      CALL READSTAT(FLAG,INXCORD,RBTTLINE)
!     READ(INXCORD,3700) RBTTLINE 
      FLAG ='#BOTTOM RIFT HORIZONTAL Y-AXIS'
      CALL READSTAT(FLAG,INXCORD,RIFTLINE)
!      READ(INXCORD,3700) RIFTLINE 
      FLAG ='#LEFT RESERVOIR VERTICAL X-AXIS'
      CALL READSTAT(FLAG,INXCORD,LEFTLINE)
!      READ(INXCORD,3700) LEFTLINE
      FLAG ='#RIGHT RESERVOIR VERTICAL X-AXIS'
      CALL READSTAT(FLAG,INXCORD,RGHTLINE)
!      READ(INXCORD,3700) RGHTLINE 
      FLAG ='#REFERENCE FOR SISMIC DOME (NOT=0/YES=1)'
      CALL IREADSTAT(FLAG,INXCORD,IDOME)
!
      CLOSE(INXCORD)
!
      IF (YTOP.NE.POSTLINE)    POSTLINE=YTOP
      IF (YBOTTOM.NE.RIFTLINE) RIFTLINE=YBOTTOM
!                                            123456789+12345678
      IF (XLEFT.GT.LEFTLINE) CALL CODERROR(6,NAMEIN)
!                                              123456789+12345678
      IF (XRIGHT.LT.RGHTLINE) CALL CODERROR(6,NAMEIN)
!
      IF (RIFTLINE.GT.RBTTLINE) CALL CODERROR(6,NAMEIN)
      IF (RBTTLINE.GT.RTOPLINE) CALL CODERROR(6,NAMEIN)
      IF (RTOPLINE.GT.SALTLINE) CALL CODERROR(6,NAMEIN)
      IF (SALTLINE.GT.POSTLINE) CALL CODERROR(6,NAMEIN)
!
      IF (IPRTIN.EQ.1) RETURN
!
#ifdef debug
         write(icoords,*) "# Coordenadas ", nsd
         do n=1,numnp
            write(icoords,2000) n,(x(i,n),i=1,nsd) 
         end do
#endif
!
      RETURN
!
 1000 FORMAT(I10)
 1500 FORMAT(2X,40(1PE15.8,2X))
 2000 FORMAT(6x,i12,10x,3(1pe15.8,2x))
 3700 FORMAT(35X,1PE15.8)
 5000 FORMAT(5X,  &
     &' ** ** ** ** ** ** ** ** *** ** ** ** ** ** ** ** ** **',/5X,&
     &' ** ** ** ** ** ** ** ** *** ** ** ** ** ** ** ** ** **',/5X,&
     &' **                                                  **',/5X,&
     &' **                                                  **',/5X,&
     &' **  COORDINATES MESH STRUCTURE:                     **',/5X,&
     &' **  OBTAINED FROM ', A30 ,                     '    **',/5X,&
     &' **                                                  **',/5X,&
     &' **                                                  **',/5X,&
     &' **                                                  **',/5X,&
     &' ** ** ** ** ** ** ** ** *** ** ** ** ** ** ** ** ** **',/5X,&
     &' ** ** ** ** ** ** ** ** *** ** ** ** ** ** ** ** ** **',/5X)
!
      END SUBROUTINE
!
!**** NEW **** MODIFIED FOR IRREGULAR MESH  ***************** 
!
      SUBROUTINE GEOREGION(XC,NUMEL)
!
      use mMalha,         only: IrregMesh, IGEOFORM
      use mPropGeoFisica, only: GEOINDIC, GEOYLOC
      use mPropGeoFisica, only: GEOFORM, MECLAW
!
!.... PROGRAM TO SETUP GEOMECHANICAL REGIONS: RESERVOIR-PRE-SAL, DOMO, POS-SAL
! 
      IMPLICIT NONE 
!  
      INTEGER :: NEL, NUMEL, INGEOFOR
!
      CHARACTER*30 :: NINGEOFOR
!
      INTEGER,      DIMENSION(NUMEL)   :: IREG, IMEC
      REAL(8),      DIMENSION(2,NUMEL) :: XC
!
      IF (IGEOFORM.EQ.0) THEN 
            DO 200 NEL=1,NUMEL
               GEOFORM(NEL) = GEOYLOC(XC(1,NEL),XC(2,NEL))
               IF (GEOFORM(NEL).EQ.'SALT_CAPROCK') THEN
                    MECLAW(NEL) = 2
                  ELSE
                    MECLAW(NEL) = 1
               ENDIF
 200        CONTINUE
         ELSE
            INGEOFOR  = 530
            NINGEOFOR = 'elmnt_geoformation.in'
            OPEN(UNIT=INGEOFOR,FILE=NINGEOFOR,STATUS='OLD')
            DO 300 NEL=1,NUMEL
               READ(INGEOFOR,1010) GEOFORM(NEL)
 300        CONTINUE
            CLOSE(INGEOFOR)
      ENDIF
!
      RETURN
!
1000  FORMAT(I10)
1010  FORMAT(A12)
!
      END SUBROUTINE
!
!**** new **********************************************************************
!
      subroutine leituraCodigosCondContorno(id, ndof, numnp, neq, iin, iecho,iprtin)
!
!.... program to read, generate and write boundary condition data
!        and establish equation numbers
!
      use mMalha, only: igen

      integer, intent(in) :: ndof, numnp, iin, iecho, iprtin
      integer:: neq
      integer, intent(inout) :: id(ndof,numnp)
!
      integer :: nn, n, i
      logical pflag
!
      id = 0
      call igen(id,ndof, iin)
!
      if (iprtin.eq.0) then
         nn=0
         do 200 n=1,numnp
         pflag = .false.
!
         do 100 i=1,ndof
         if (id(i,n).ne.0) pflag = .true.
  100    continue
!
         if (pflag) then      
            nn = nn + 1
            if (mod(nn,50).eq.1) write(iecho,1000) (i,i=1,ndof)
            write(iecho,2000) n,(id(i,n),i=1,ndof)
         endif
  200    continue
      endif
!
!.... establish equation numbers
!
      neq = 0
!
      do 400 n=1,numnp
!
      do 300 i=1,ndof
      if (id(i,n).eq.0) then
         neq = neq + 1
         id(i,n) = neq
      else
         id(i,n) = 1 - id(i,n)
      endif
!
  300 continue
!
  400 continue
!
      return
!
 1000 format('1',' n o d a l   b o u n d a r y   c o n d i t i o n & 
     &         c o  d e s'/// &
      5x,' node no.',3x,6(6x,'dof',i1:)//)
 2000 format(6x,i10,5x,6(5x,i10))
!
      end subroutine
!
!**** NEW *** MODIFIED leituraCodigosCondContorno SUBROUTINE *****************
!
      SUBROUTINE ReadDirichlet4GeoMechanic(id, ndof, numnp, neq, IIN, iecho,iprtin)
!
!.... PROGRAM TO READ/IMPORT DATA, GENERATE AND WRITE BOUNDARY CONDITION DATA
!        AND ESTABLISH EQUATION NUMBERS
!
      use mMalha, only: igen
      use mMalha, only: Dirichlet
!
      INTEGER, INTENT(IN) :: NDOF, NUMNP, IIN, IECHO, IPRTIN
      INTEGER, INTENT(INOUT) :: ID(NDOF,NUMNP)
      INTEGER  :: NEQ
!
      INTEGER :: NN, N, I
      LOGICAL PFLAG
!.... NEXT TWO LINES NEW VARS FOR IMPORT DIRICHLET BOUNDARY CONDITIONS
      INTEGER  :: INPUTFILE
      CHARACTER*30 :: ROTULO, NAMEIN
!
      ID = 0
!
      INPUTFILE = IIN 
      CALL IGEN(ID,NDOF,INPUTFILE)
!
      IF (Dirichlet.EQ.0) THEN 
            WRITE(*,5000) 'DEFAULT FILE  : input.in'
         ELSE
            INPUTFILE = 531
            NAMEIN   = 'dirichlet_for_geomech.in'
            WRITE(*,5000) NAMEIN
            OPEN(UNIT=INPUTFILE,FILE=NAMEIN,STATUS='OLD')
            READ(INPUTFILE,*) ROTULO
            CALL IGEN(ID,NDOF,INPUTFILE)
            CLOSE(INPUTFILE)
      ENDIF
!
      IF (IPRTIN.EQ.0) THEN
         NN=0
         DO 200 N=1,NUMNP
            PFLAG = .FALSE.
            DO 100 I=1,NDOF
               IF (ID(I,N).NE.0) PFLAG = .TRUE.
  100       CONTINUE
!
            IF (PFLAG) THEN      
               NN = NN + 1
               IF (MOD(NN,50).EQ.1) WRITE(IECHO,1000) (I,I=1,NDOF)
               WRITE(IECHO,2000) N,(ID(I,N),I=1,NDOF)
            ENDIF
  200    CONTINUE
      ENDIF
!
!.... ESTABLISH EQUATION NUMBERS
!
      NEQ = 0
!
      DO 400 N=1,NUMNP
         DO 300 I=1,NDOF
            IF (ID(I,N).EQ.0) THEN
                 NEQ = NEQ + 1
                 ID(I,N) = NEQ
              ELSE
                 ID(I,N) = 1 - ID(I,N)
            ENDIF
  300    CONTINUE
  400 CONTINUE
!
      RETURN
!
 1000 FORMAT('1',' N O D A L   B O U N D A R Y   C O N D I T I O N & 
     &         C O  D E S'/// &
      5X,' NODE NO.',3X,6(6X,'DOF',I1:)//)
 2000 FORMAT(6X,I10,5X,6(5X,I10))
!
 5000 FORMAT(5X,  &
     &' ** ** ** ** ** ** ** ** *** ** ** ** ** ** ** ** ** **',/5X,&
     &' ** ** ** ** ** ** ** ** *** ** ** ** ** ** ** ** ** **',/5X,&
     &' **                                                  **',/5X,&
     &' **                                                  **',/5X,&
     &' **  DIRICHLET GEOMECHANICAL CONDITIONS:             **',/5X,&
     &' **  OBTAINED FROM ', A30 ,                     '    **',/5X,&
     &' **                                                  **',/5X,&
     &' **                                                  **',/5X,&
     &' **                                                  **',/5X,&
     &' ** ** ** ** ** ** ** ** *** ** ** ** ** ** ** ** ** **',/5X,&
     &' ** ** ** ** ** ** ** ** *** ** ** ** ** ** ** ** ** **',/5X)
!
      END SUBROUTINE
!
!**** new **********************************************************************
!
      subroutine leituraValoresCondContorno(f,ndof,numnp,j,nlvect,iprtin)
!
!.... program to read, generate and write nodal input data
!
!        f(ndof,numnp,nlvect) = prescribed forces/kinematic data (j=0)
!                             = nodal body forces(j=1)
!
      use mMalha, only : genfl
      implicit none
!
!.... remove above card for single-precision operation
!
      integer :: ndof, numnp, j, nlvect, iprtin
      real*8 f(ndof,numnp,nlvect)

      logical lzero
      integer nlv
      character(len=35) :: rotulo
!
!     call clear(f,nlvect*numnp*ndof)
!erro bidu      f(1:nlvect,1:numnp,1:ndof)=0.0
!
      F(1:NDOF,1:NUMNP,1:NLVECT)=0.0D0
!
      do 100 nlv=1,nlvect
      call genfl(f(1,1,nlv),ndof,iin)
!       
      call ztest(f(1,1,nlv),ndof*numnp,lzero)
! 
      if (iprtin.eq.0) then
!
         if (lzero) then
            if (j.eq.0) write(iecho,1000) nlv
            if (j.eq.1) write(iecho,2000)
         else
            if (j.eq.0) call printf(f,ndof,numnp,nlv)
!
            if (j.eq.1) then
               rotulo=" n o d a l  b o d y  f o r c e s"
               call printd (rotulo, f,ndof,numnp,iecho)
            end if
!
         endif
      endif
!
  100 continue
!
      return
 1000 format('1'//,' there are no nonzero prescribed forces and ',&
         'kinematic boundary conditions for load vector number ',i10)
 2000 format('1'//,' there are no nonzero nodal body forces')
      end subroutine
!
!**** NEW *** MODIFIED leituraValoresCondContorno SUBROUTINE ******************
!
      SUBROUTINE ReadNeumann4GeoMechanic(f,ndof,numnp,j,nlvect,IIN,iprtin)
!
!.... PROGRAM TO READ/IMPORT INFO TO GENERATE AND WRITE NODAL INPUT DATA
!....  MODIFIED ONLY FOR ONE AND ONLY ONE FORCE VECTOR LOAD: NLVECT = 1
!
!....    F(NDOF,NUMNP,NLVECT) = PRESCRIBED FORCES/KINEMATIC DATA (J=0)
!                             = NODAL BODY FORCES(J=1)
!
      use mMalha, only : NEUMANN, genfl
      IMPLICIT NONE
!
!.... remove above card for single-precision operation
!
      INTEGER :: NDOF, NUMNP, J, NLVECT, IIN, IPRTIN, INPUTFILE
      REAL*8 F(NDOF,NUMNP,NLVECT)
!
      LOGICAL LZERO
      INTEGER NLV
      CHARACTER(LEN=35) :: ROTULO, NAMEIN
!
!     call clear(f,nlvect*numnp*ndof)
!erro bidu      f(1:nlvect,1:numnp,1:ndof)=0.0
!
      F(1:NDOF,1:NUMNP,1:NLVECT)=0.0D0
!
      DO 100 NLV=1,NLVECT
         INPUTFILE = IIN
         CALL GENFL(F(1,1,NLV),NDOF,INPUTFILE)
         IF ((NEUMANN.EQ.1).AND.(NLVECT.EQ.1)) THEN
               INPUTFILE = 531
               NAMEIN = 'neumann_for_geomech.in'
               WRITE(*,5000) NAMEIN
               OPEN(UNIT=INPUTFILE,FILE=NAMEIN,STATUS='OLD')
               READ(INPUTFILE,*) ROTULO
               CALL GENFL(F(1,1,NLV),NDOF,INPUTFILE)
               CLOSE(INPUTFILE)
            ELSE
               WRITE(*,5000) 'DEFAULT FILE  : input.in'
         ENDIF
!
         CALL ZTEST(F(1,1,NLV),NDOF*NUMNP,LZERO)
! 
         IF (IPRTIN.EQ.0) THEN
            IF (LZERO) then
                 IF (J.EQ.0) WRITE(IECHO,1000) NLV
                 IF (J.EQ.1) WRITE(IECHO,2000)
               ELSE
                 IF (J.EQ.0) CALL PRINTF(F,NDOF,NUMNP,NLV)
                 IF (J.EQ.1) THEN
                    ROTULO=" N O D A L  B O D Y  F O R C E S"
                    CALL PRINTD (ROTULO, F,NDOF,NUMNP,IECHO)
                 END IF
            ENDIF
         ENDIF
  100 CONTINUE
!
      RETURN
 1000 FORMAT('1'//,' THERE ARE NO NONZERO PRESCRIBED FORCES AND ',&
         'KINEMATIC BOUNDARY CONDITIONS FOR LOAD VECTOR NUMBER ',I10)
 2000 FORMAT('1'//,' THERE ARE NO NONZERO NODAL BODY FORCES')
!
 5000 FORMAT(5X,  &
     &' ** ** ** ** ** ** ** ** *** ** ** ** ** ** ** ** ** **',/5X,&
     &' ** ** ** ** ** ** ** ** *** ** ** ** ** ** ** ** ** **',/5X,&
     &' **                                                  **',/5X,&
     &' **                                                  **',/5X,&
     &' **  NEUMANN GEOMECHANICAL CONDITIONS:               **',/5X,&
     &' **  OBTAINED FROM ', A30 ,                     '    **',/5X,&
     &' **                                                  **',/5X,&
     &' **                                                  **',/5X,&
     &' **                                                  **',/5X,&
     &' ** ** ** ** ** ** ** ** *** ** ** ** ** ** ** ** ** **',/5X,&
     &' ** ** ** ** ** ** ** ** *** ** ** ** ** ** ** ** ** **',/5X)
!
      END SUBROUTINE
!
!**** new ********************************************************************
!
      subroutine printf(f,ndof,numnp,nlv)
!
!.... program to print prescribed force and boundary condition data
!
      implicit none
!
!.... remove above card for single precision operation
!
      integer ndof, numnp, nlv
      real*8 :: f(ndof,numnp,*)
!
      logical lzero
      integer :: nn, n, i
!
      nn = 0
!
      do 100 n=1,numnp
      call ztest(f(1,n,nlv),ndof,lzero)
      if (.not.lzero) then
         nn = nn + 1
         if (mod(nn,50).eq.1) write(iecho,1000) nlv,(i,i=1,ndof)
         write(iecho,2000) n,(f(i,n,nlv),i=1,ndof)
      endif
  100 continue
!
      return
!
 1000 format('1',&
     ' p r e s c r i b e d   f o r c e s   a n d   k i n e m a t i c ',&
     '  b o u n d a r y   c o n d i t i o n s'//5x,&
     ' load vector number = ',i10///5x,&
     ' node no.',6(13x,'dof',i1,:)/)
 2000 format(6x,i10,10x,6(1pe15.8,2x))
      end subroutine
!      
!**** new ********************************************************************
!
      subroutine printd(name,dva,ndof,numnp,icode)
!
!.... program to print kinematic data
!
      implicit none
!
!.... remove above card for single precision operation
!
      integer :: ndof,numnp, icode
      character (LEN=*) ::  name
      real*8 dva(ndof,*)
!
      logical lzero
      integer nn, n, i
!
      nn = 0
!
      do 100 n=1,numnp
!       call ztest(dva(1,n),ndof,lzero)
!       if (.not.lzero) then
!          nn = nn + 1
!          if (mod(nn,50).eq.1) &
!            write(icode,1000) name,(i,i=1,ndof)
         write(icode,2000) n,(dva(i,n),i=1,ndof)
!       endif
  100 continue
!
      return
!
 1000 format('1',11a4//1x,'node',6(11x,'dof',i1)/)
 2000 format(1x,i10,2x,6(1pe30.10,2x))
      end subroutine
!
!**** new *******************************************************************
!
      subroutine printp(a,idiag,neq,nsq,*)
      use mGlobaisEscalares
!
!.... program to print array d after Crout factorization 
!        a = u(transpose) * d * u
!
      implicit none
!
!.... remove above card for single precision operation
!
      integer :: neq, nsq
      real*8 :: a(*)
      integer :: idiag(*)
!
      integer :: n, i
!
      do 100 n=1,neq
         if (mod(n,50).eq.1) write(iecho,1000) nsq
         write(iecho,1000)
         i = idiag(n)
         write(iecho,2000) n,a(i)
  100 continue
!
      return 1
!
 1000 format('1',' array d of factorization',/&
     ' a = u(transpose) * d * u ',                                //5x,&
     ' time sequence number   . . . . . . . . . . . . (nsq) = ',i10//5x)
 2000 format(1x,i10,4x,1pe20.8)
      end subroutine
!
!**** new *******************************************************************
!
      subroutine prntel(mat,conectElem,nen,numel,tipo)
      implicit none
!
!.... program to print data for element with "nen" nodes
!
!        note: presently the label formats are limited to
!              elements with one to nine nodes
!
      integer :: nen, numel
      integer :: mat(*),conectElem(nen,*)
      integer :: tipo
!
      integer n, i
!
      if(tipo==1) then
      write(iconects,*) "# Conectividades nodais"
      do n=1,numel
        write(iconects,2000) n,mat(n),(conectElem(i,n),i=1,nen)
      end do
      end if

      if(tipo==2) then
      write(iconectsL,*) "# Conectividades ladais"
      do  n=1,numel
        write(iconectsL,3000) n,mat(n),(conectElem(i,n),i=1,nen)
      end do
      end if
!
      return
!
 2000 format(1x,i10,9(2x,i10))
 3000 format(1x,i10,7(2x,i10))
      end subroutine
!
!**** new *******************************************************************
!
      subroutine printResultado(dva, ndof, numnp, inicio, fim, icode)
!
!.... program to print kinematic data
!
      implicit none
!
!.... remove above card for single precision operation
!
      integer :: ndof, numnp, inicio, fim, icode
      real*8  :: dva(ndof,numnp)
!
      integer :: n, i
!
      write(icode,*) "# Solucao"
      do 100 n=inicio,fim
         write(icode,2000) n,(dva(i,n),i=1,ndof)
         !write(*,*) n,(dva(i,n),i=1,ndof)
  100 continue
!
      return
 2000 format(1x,i10,2x,6(1pe13.6,2x))
      end subroutine
!
!**** new *******************************************************************
!
      subroutine prtgnup(name,x,dva,nsd,ndof,numnp,icode)
!
!.... program to print kinematic data
!
      implicit none
!
!.... remove above card for single precision operation
!
      integer :: nsd, ndof, numnp, icode
      character*4 name(11)
      real*8 :: x(nsd,*),dva(ndof,*)
!
      integer :: n, j, i
!
      write(icode,*) name
      do 100 n=1,numnp
         write(icode,2000) (x(j,n),j=1,nsd), (dva(i,n),i=1,ndof)
  100 continue
!
      return
!
 2000 format(6(1pe13.6,2x))
      end subroutine
!
!**** new *******************************************************************
!
    subroutine escreverArqParaviewReservatorio(arquivo, campo, dim1, dim2, nen, conectElem, tipo, rotulo, tamRot)
    use mMalha, only: x, nsd, numelReserv, numnpReserv

    implicit none
    integer, intent(in) :: arquivo,dim1, dim2
    double precision, intent(in) :: campo(dim1, dim2)
    integer :: nen
    integer :: conectElem(nen,numelReserv)
    integer :: tipo  !tipo=1 para elemento, e tipo=2 para no
    integer :: tamRot

    character(len=tamRot) :: rotulo

  
    write(arquivo,'(a)')'# vtk DataFile Version 3.0'
    write(arquivo,'(a)')'vtk output'
    write(arquivo,'(a)')'ASCII'
    write(arquivo,'(a)')'DATASET UNSTRUCTURED_GRID'
    write(arquivo,'(a,i10,a)')'POINTS', numnpReserv,' float '

    call escreverPontosNodais  (arquivo,x, numnpReserv, nsd)
! 
    write(arquivo,'(a,i10,i10)')'CELLS', numelReserv , (nen+1) * numelReserv
    call escreverConectividades(arquivo,conectElem, numelReserv, nen, nsd) !apenas reservatório
! 
    write(arquivo,'(a,i10)')'CELL_TYPES ', numelReserv
    call escreverTiposElementos(arquivo,numelReserv,nsd)
! 

    if(tipo==1) write(arquivo,'(a,i10)')'CELL_DATA ', numelReserv

    if(tipo==2) write(arquivo,'(a,i10)')'POINT_DATA',  dim1*dim2

    write(arquivo,'(3a)')'SCALARS ', trim(rotulo), ' float '
    write(arquivo,'(a)')'LOOKUP_TABLE default'

    call escreverEscalaresPorElemento(arquivo,campo, dim2,tamRot)

    end subroutine escreverArqParaviewReservatorio
!
!**** new *******************************************************************
!
    subroutine escreverArqParaview(arquivo, campo, dim1, dim2, nen, conectElem, tipo, rotulo, tamRot)
    use mMalha, only: x, nsd, numel,numelReserv, numnp, numnpReserv
!
    implicit none
    integer, intent(in) :: arquivo,dim1, dim2
    double precision, intent(in) :: campo(dim1, dim2)
    integer :: nen
    integer :: conectElem(nen,numel)
    integer :: tipo  !tipo=1 para elemento, e tipo=2 para no
    integer :: tamRot

    character(len=tamRot) :: rotulo

  
    write(arquivo,'(a)')'# vtk DataFile Version 3.0'
    write(arquivo,'(a)')'vtk output'
    write(arquivo,'(a)')'ASCII'
    write(arquivo,'(a)')'DATASET UNSTRUCTURED_GRID'
    write(arquivo,'(a,i10,a)')'POINTS', numnp,' float '

    call escreverPontosNodais  (arquivo,x, numnp, nsd)
! 
    write(arquivo,'(a,i10,i10)')'CELLS', numel , (nen+1) * numel
    call escreverConectividades(arquivo,conectElem, numel, nen, nsd) !todo o domínio
! 
    write(arquivo,'(a,i10)')'CELL_TYPES ', numel
    call escreverTiposElementos(arquivo,numel,nsd)
! 

    if(tipo==1) write(arquivo,'(a,i10)')'CELL_DATA ', numel

    if(tipo==2) write(arquivo,'(a,i10)')'POINT_DATA',  dim1*dim2

    write(arquivo,'(3a)')'SCALARS ', trim(rotulo), ' float '
    write(arquivo,'(a)')'LOOKUP_TABLE default'

    call escreverEscalaresPorElemento(arquivo,campo, dim2,tamRot)
!
    end subroutine escreverArqParaview

!**** new *******************************************************************
      subroutine escreverPontosNodais  (arquivo,coords, numnp, nsd)
      implicit none
      integer, intent(in) :: arquivo,numnp, nsd
      real*8,  intent(in) :: coords(nsd,numnp)
!
      real*8  :: coordZ = 0.0 
      integer :: d, i
!
      if(nsd==2) then 
        do i=1,numnp
            write(arquivo,'(3(1x, 1pe15.8))') (coords(d,i),d=1,nsd), coordZ 
        end do
      end if

      if(nsd==3) then
        do i=1,numnp
            write(arquivo,'(3(1x, 1pe15.8))') (coords(d,i),d=1,nsd)
        end do
      end if
      end subroutine escreverPontosNodais


!**** new *******************************************************************
      subroutine escreverConectividades(arquivo,conectElem, numel, nen, nsd)
      implicit none
      integer, intent(in)  :: arquivo,numel, nen, nsd
      integer, intent(in)  :: conectElem(nen,numel)
!
      integer n, i
!
      if(nsd==2) then
      do  n=1,numel
        write(arquivo,'(i10,9(2x,i10))') nen, (conectElem(i,n)-1, i = 1, nen) 
      end do
      end if

      if(nsd==3) then
      do  n=1,numel
        write(arquivo,'(i10,18(2x,i10))') nen, (conectElem(i,n)-1, i = 1, nen) 
      end do
      end if

 end subroutine escreverConectividades

!**** new *******************************************************************
      subroutine escreverTiposElementos(arquivo,numel, nsd)
      implicit none
      integer, intent(in)   :: arquivo, numel, nsd
!
      integer :: i
!
      if(nsd==2) then
      do  i =1,numel
        write(arquivo,'(a)') '9'!trim(adjustl(tipo))
      end do
      end if 
!
      if(nsd==3) then
      do  i =1,numel
        write(arquivo,'(a)') '12'!trim(adjustl(tipo))
      end do
      end if 

      end subroutine escreverTiposElementos
!
!**** new *******************************************************************
!
      subroutine escreverEscalaresNodais(arquivo,v, tam1, tam2, rotulo, tamRot)
!
      use mMalha, only: numelReserv
!
      implicit none
      integer, intent(in)  :: arquivo,tam1,tam2
      real*8, intent(in)   :: v(tam1,tam2)
      integer :: tamRot
      character(len=tamRot) :: rotulo
!
      character(len=tamRot+5) ::  rotuloN
      integer :: i,j
      character(len=5):: eixo
      real*8 :: limite,zero
!
      limite=1.e-15
      zero=0.0d0
      do i=1,tam1

        if(i>1) then
           write(eixo,'(i0)') i
           rotuloN=trim(rotulo)//trim(eixo)
           write(arquivo,'(3a)')'SCALARS ', trim(rotuloN), ' float '
           write(arquivo,'(a)')'LOOKUP_TABLE default'
        endif
!
      if(apenasReservatorio.eqv..true.) then
        do j=1, tam2
             if(v(i,j).lt.limite) then
                write(arquivo,*) zero
             else 
                write(arquivo,*) v(i,j)
             end if
        end do
!
      else

       do j=1, tam2
          if(tam2<=numelReserv) then
             if(v(i,j).lt.limite) then
                write(arquivo,*) zero
             else
                write(arquivo,*) v(i,j)
             end if
          else
             write(arquivo,*) zero
          endif
      end do
!
        endif
!
      end do
!
      end subroutine escreverEscalaresNodais
!
!**** new *******************************************************************
!
      subroutine escreverEscalaresPorElemento(arquivo,v, tam, tamRot)
!
      use mMalha, only: numelReserv
!
      implicit none
      integer, intent(in)  :: arquivo,tam
      real*8, intent(in)   :: v(tam)
      integer :: tamRot
!
      integer :: j
      real*8 :: limite,zero

      limite=1.e-15
      zero=0.0d0

      if(apenasReservatorio.eqv..true.) then
      do j=1, tam
              if(v(j).lt.limite) then
                write(arquivo,*) zero
             else
                write(arquivo,*) v(j)
             end if
      end do
!
      else
!
      do j=1, tam
          if(tam<=numelReserv) then
             if(v(j).lt.limite) then
                write(arquivo,*) zero
             else
                write(arquivo,*) v(j)
             end if
          else
             write(arquivo,*) zero
          endif
      end do
      endif
!
      end subroutine 
!
!**** new *******************************************************************
!
    subroutine escreverArqParaviewIntermed(arquivo, campo, dim1, dim2, rotulo, tamRot)
    use mMalha, only: x, nsd, numel, numnp

    implicit none
    integer, intent(in) :: arquivo,dim1, dim2
    double precision, intent(in) :: campo(dim1, dim2)

    integer :: tamRot
    character(len=tamRot) :: rotulo

    write(arquivo,'(3a)')'SCALARS ', trim(rotulo), ' float '
    write(arquivo,'(a)')'LOOKUP_TABLE default'

     call escreverEscalaresNodais(arquivo,campo, dim1, dim2, rotulo, tamRot)

     end subroutine escreverArqParaviewIntermed
!
!----------------------------------------------------------------------
!
      subroutine paraview_geraCase(steps)

      implicit none

      integer :: steps
!
      integer :: numInicial, incremento
      real*8  :: incTempo
      integer :: i

      numInicial=0
      incremento=1
      incTempo =0.0

      open(unit=124,file="./out/transiente.case",status="unknown")

      write(124, "('FORMAT',/,'type:',2x,'ensight')")
      write(124, *)
      write(124, "('GEOMETRY',/,'model:',2x,'solucao.geo')")
      write(124, *)
      write(124, "('VARIABLE',/,'scalar per element:', 2x, 'Saturacao', 2x, 'solucao.***' )")
      write(124, *)
      write(124, "('TIME',/,'time set: 1')")
      write(124, "('number of steps:', i10)"), steps
      write(124, "('filename start number:', i10)"), numInicial
      write(124, "('filename increment:', i10)"), incremento
      write(124, "('time values:')")

      do i=1, steps+1
           write(124, *), incTempo
           incTempo=incTempo+1.0
      end do

      end subroutine
!
!----------------------------------------------------------------------
!
    subroutine paraview_geometria(numel,numnp,nsd,x,conecNodaisElem)

    implicit none

      integer :: numel,numnp,nsd
      real(8), dimension(nsd,*) :: x   
      integer :: conecNodaisElem(8,numel)
!
      integer :: i
      open(unit=125,file="./out/solucao.geo",status="unknown")

      write(125,'(a)')'Title1'
      write(125,'(a)')'Title2'
      write(125,'(a)')'node id given'
      write(125,'(a)')'element id given'
      write(125,'(a)')'coordinates'
      write(125,'(i8)')  numnp

      do i = 1, numnp
      WRITE (125,'(I8,3E12.5)') I,x(1,i),x(2,i),x(3,i)
      enddo

      WRITE (125,'(A,/,A,/,A,/,I8)')                     &
                                'part 1'           ,    &
                                'malha'            ,    &
                                'hexa8'            ,    &
                                 numel

      WRITE (125,'(9I8)')  (I,conecNodaisElem(1,i),conecNodaisElem(2,i),conecNodaisElem(3,i), &
                                        conecNodaisElem(4,i),conecNodaisElem(5,i),conecNodaisElem(6,i), & 
                                        conecNodaisElem(7,i),conecNodaisElem(8,i),i=1, numel ) 

     end subroutine paraview_geometria
!
!*****************************************************
!
      subroutine imprimirCaseParaview(x, conecNodaisElem, pressaoElem, satElem, phi, perm)

      use mMalha,               only: numnp,nsd,numel,nen,numelReserv

      implicit none

      real*8 :: x(nsd, numnp)
      integer :: conecNodaisElem(nen,*)
      real*8 ::  pressaoElem(*), satElem(*), phi(*), perm(*)
!
      if(iflag_sat==1) then    
         if(iflag_tipoPrint==1) then 
            open(unit=ipress    , file= './out/solucao.P0001')
            open(unit=iperm     , file= './out/solucao.K0001')
            open(unit=iporo     , file= './out/solucao.PHI0001')
!             open(unit=iveloc      , file= 'solucao.V0001')
            open(unit=isaturacao, file= './out/solucao.S0001')

            call paraview_geraCase(qtdImpSat)
            call paraview_geometria(numel,numnp,nsd,x, conecNodaisElem)
            call paraview_escalarPorElemento(numel, pressaoElem,ipress)
            call paraview_escalarPorElemento(numel, perm, iperm)
            call paraview_escalarPorElemento(numel, phi, iporo)
            call paraview_escalarPorElemento(numelReserv, satElem,isaturacao)
         endif
      endif
!
      end subroutine
!
!----------------------------------------------------------------------
!
      subroutine paraview_vetorPorElemento(numel,campo,iarq)

! ainda precisa implementar
      implicit none
!
      integer :: numel
      real(8), dimension(*) :: campo
      integer :: iarq

      integer :: i

!
      write(iarq,"('Ensight Scalar passo     1')")
      write(iarq,"('part 1')")
      write(iarq,"('hexa8')")
!
      write(iarq,"(6e12.5)") (campo(i),i=1,numel)
!      
      close(iarq)
!
      end subroutine    
!
!----------------------------------------------------------------------
!
      subroutine paraview_escalarPorElemento(numel,campo,iarq)
      implicit none
!
      integer :: numel
      real(8), dimension(*) :: campo
      integer :: iarq

      integer :: i

      print*, "gerando", iarq
!
      write(iarq,"('Ensight Scalar passo     1')")
      write(iarq,"('part 1')")
      write(iarq,"('hexa8')")
!
      write(iarq,"(6e12.5)") (campo(i),i=1,numel)
!      
      close(iarq)
!
      end subroutine    
!
!----------------------------------------------------------------------
!
      subroutine paraview_escalarPorElementoTransiente(numel,campo,passo,iarq)
      implicit none
!
      integer :: numel,passo,iarq
      real(8), dimension(*) :: campo   
      character(len=128) :: name,sol
      character(len=8)   :: c
      integer :: i
      real(4) :: x     
!      
      x=0.001

      sol="solucao"

      write(c,"(f7.3)") x*passo
      c=adjustl(c)
      name='./out/'//trim(sol)//c 
!      
      open(unit=iarq,file=name,status="unknown")
!
      write(iarq,"('Ensight Scalar passo ',i5)") passo
      write(iarq,"('part 1')")
      write(iarq,"('hexa8')")
!      
!     imprime as coordenadas
!
       write(iarq,"(6(e12.5))") (real(campo(i)),i=1,numel)
!      
      close(iarq)
!
      passo=passo+1
!
      end subroutine          
!
!=======================================================================
!     
      subroutine prt(nsd,numel,t0,u,iunit)
!      
      use mMalha, only: xc
!      
      implicit none
!     
!     imprime campos escalares para o gnuplot ou para o matlab
!
      integer                   :: nsd,numel
      real(8), dimension(*)     :: u
      real(8)                   :: t0
      integer :: iunit
!     
      integer :: nel
!     
      write(iunit,"('#TIMESTEP PRINT OUT = ',f15.8)") t0
      write(iunit,*)
!     
      do nel=1,numel
         write(iunit,"(5(f25.15,2x))") xc(1:nsd,nel),u(nel)
      end do
!     
      write(iunit,*)
!     
      end subroutine
!     
!=======================================================================
!
      subroutine prtvB(nsd,numel,t0,velocLadal,ndofV,  conecLadaisElem, numLadosElem, iunit)
!
      use mMalha, only: xc
!
      implicit none
!
!     imprime campos vetoriais para o gnuplot ou para o matlab
!
      integer                   :: numel,nsd, ndofV ,numLadosElem
      real(8), dimension(ndofV,*) :: velocLadal
      real(8)                   :: t0
      integer                   :: conecLadaisElem(numLadosElem,numel)
!
      integer :: nel
      real(8) :: vc(nsd)
      real*8 :: mediaCentro
!
      integer :: iunit

      vc=0.0
!
      write(iunit,"('#TIMESTEP PRINT OUT = ',f15.8)") t0
!
      write(iunit,*)
!
      do nel=1,numel
!
        vc(1) = (velocLadal(1,conecLadaisElem(2,nel))+velocLadal(1,conecLadaisElem(4,nel)))/2.0
        vc(2) = (velocLadal(1,conecLadaisElem(1,nel))+velocLadal(1,conecLadaisElem(3,nel)))/2.0
        if(nsd==3) vc(3) = (velocLadal(1,conecLadaisElem(5,nel))+velocLadal(1,conecLadaisElem(6,nel)))/2.0
        mediaCentro=sum(vc)/nsd
        write(iunit,"(6(f25.15,2x))")xc(1:nsd,nel),mediaCentro
!
      end do
!
      write(iunit,*)
!
      end subroutine
!     
!=======================================================================
!     
      subroutine prtv(nsd,numel,ndofV,numLadosReserv,t0,u,iunit)
!      
      use mMalha, only: xc
!
      implicit none
!     
!     imprime campos vetoriais para o gnuplot ou para o matlab
!
      integer                   :: numel,nsd, ndofV, numLadosReserv
      real(8), dimension(ndofV,numLadosReserv) :: u
      real(8)                   :: t0
!
      integer :: nel
!     
      integer :: iunit
!     
      write(iunit,"('#TIMESTEP PRINT OUT = ',f15.8)") t0
!
      write(iunit,*)
!     
      do nel=1,numel
         write(iunit,"(5(f25.15,2x))") xc(1:nsd,nel),u(1,nel),u(2,nel)
      end do
!     
      write(iunit,*)
!     
      end subroutine
!
!*** NEW ***************************************************************
!
      SUBROUTINE SETUPDX()
!
      use mGlobaisEscalares, only: NUMDX, NITGEO, NITHIDRO, PATHDX
!..
!...  PROGRAM TO SETUP MANAGER FILES FOR GRAPHICAL INTERFACE OPEN-DX
!
      CHARACTER*30 NIFEDX,NIFNOD,NIFMSH,NIS3XX,NIS3YY,NIS3XY,NIS3ZZ
      CHARACTER*30 NISTDX,NISTYN,NISTPS,NISTPE,NISTSV,NISTPH
!
      CHARACTER*2 ASTEP1
      CHARACTER*2 ASTEP2
      CHARACTER*20 MKDIR
!
      IF (NUMDX.EQ.0) RETURN
!
      WRITE(ASTEP1,'(I2.2)') NITGEO
      WRITE(ASTEP2,'(I2.2)') NITHIDRO
!
      PATHDX = 'dxelast'//ASTEP1
      PATHDX = 'dxelast'//ASTEP1//'.ht'//ASTEP2
!      PATHDX = 'dxcreep'//ASTEP1
!      PATHDX = 'dxcreep'//ASTEP1//'.ht'//ASTEP2
      MKDIR  = 'rm -r '//PATHDX
      CALL SYSTEM(mkdir)
      MKDIR  = 'mkdir '//PATHDX
      CALL SYSTEM(mkdir)
!
      NIFEDX = PATHDX//'/nodestoc.dx'
      NIFNOD = PATHDX//'/fnodes.stoc'
      NIFMSH = PATHDX//'/femesh.stoc'
      NISTDX = PATHDX//'/parametr.dx'
      NISTYN = PATHDX//'/youngmd.dat'
      NISTPS = PATHDX//'/poisson.dat'
      NISTPE = PATHDX//'/permeab.dat'
      NISTPH = PATHDX//'/displac.dat'
      NISTSV = PATHDX//'/geoprsr.dat'
      NIS3XX = PATHDX//'/strs0xx.dat'
      NIS3YY = PATHDX//'/strs0yy.dat'
      NIS3XY = PATHDX//'/strs0xy.dat'
      NIS3ZZ = PATHDX//'/strs0zz.dat'
!
      OPEN(UNIT=IFEDX, FILE= NIFEDX)
      OPEN(UNIT=IFNOD, FILE= NIFNOD)
      OPEN(UNIT=IFMSH, FILE= NIFMSH)
      OPEN(UNIT=ISTDX, FILE= NISTDX)
      OPEN(UNIT=ISTYN, FILE= NISTYN)
      OPEN(UNIT=ISTPS, FILE= NISTPS)
      OPEN(UNIT=ISTPE, FILE= NISTPE)
      OPEN(UNIT=ISTPH, FILE= NISTPH)
      OPEN(UNIT=ISTSV, FILE= NISTSV)
      OPEN(UNIT=IS3XX, FILE= NIS3XX)
      OPEN(UNIT=IS3YY, FILE= NIS3YY)
      OPEN(UNIT=IS3XY, FILE= NIS3XY)
      OPEN(UNIT=IS3ZZ, FILE= NIS3ZZ)
!
      END SUBROUTINE
!
!**** NEW **** FOR DATA EXPLORER OUT PUT ************************************* 
!
      SUBROUTINE PRINT_DXMESH(X,DIS,GEOPRSR,STRSS0,conecNodaisElem,YOUNG,PERM)
!
      use mGlobaisEscalares, only: NUMDX
      use mMalha,            only: NEN, NSD, numel, numnp, numelReserv
      use mPropGeoFisica,    only: GEOFORM, GEOINDIC
!
      INTEGER :: K, NODE, NEL
      INTEGER, DIMENSION(NEN,NUMEL)   :: conecNodaisElem
      REAL(8), DIMENSION(NSD,NUMNP)   :: X, DIS
      REAL(8), DIMENSION(NUMEL)       :: YOUNG, GEOPRSR
      REAL(8), DIMENSION(numelReserv) :: PERM
      REAL(8), DIMENSION(4,NUMEL)     :: STRSS0
      REAL(8)                         :: XPRINT
! 
!.... PRINT DISPLACEMENTS DATA FOR OPEN-DX FILE
! 
      DO 10 NODE=1,NUMNP
 	 WRITE(IFNOD,1900) X(1,NODE), X(2,NODE)
 	 WRITE(ISTPH,1900) DIS(1,NODE), DIS(2,NODE)
 10   CONTINUE 
!
      CLOSE(IFNOD)
      CLOSE(ISTPH)
! 
!.... PRINT CONECTIVITIES DATA FOR OPEN-DX FILE
! 
      DO 20 NEL=1,NUMEL 
         WRITE(IFMSH,2000) conecNodaisElem(1,NEL)-1, & 
     &                     conecNodaisElem(2,NEL)-1, &
     &                     conecNodaisElem(4,NEL)-1, &
     &                     conecNodaisElem(3,NEL)-1 
!
         XPRINT = 1.0D-16
         IF (GEOFORM(NEL).EQ.'RESERVATORIO') XPRINT = PERM(NEL)
         WRITE(ISTPE,1900) XPRINT
         WRITE(ISTYN,1900) YOUNG(NEL)
         WRITE(ISTPS,1900) GEOINDIC('POISSON',GEOFORM(NEL))
         WRITE(ISTSV,1900) -GEOPRSR(NEL)
         WRITE(IS3XX,1900) STRSS0(1,NEL)
         WRITE(IS3YY,1900) STRSS0(2,NEL)
         WRITE(IS3XY,1900) STRSS0(3,NEL)
         WRITE(IS3ZZ,1900) STRSS0(4,NEL)
 20   CONTINUE
!
! WRITE(IS3ZZ,1900)GEOINDIC('POISSON',GEOFORM(NEL))*(STRSS0(1,NEL)+STRSS0(2,NEL)),STRSS0(4,NEL)
      CLOSE(IFMSH)
      CLOSE(ISTPE)
      CLOSE(ISTYN)
      CLOSE(ISTPS)
      CLOSE(ISTSV)
      CLOSE(IS3XX)
      CLOSE(IS3YY)
      CLOSE(IS3XY)
      CLOSE(IS3ZZ)
!                     
!.... OPEN-DX DATA FILES FOR INITIAL FIELDS AND PARAMETERS
!
      CALL PRINT_DXINFO('OPEN_FEMDX_FILE',ISTDX,NUMNP,NUMEL)
      CALL PRINT_DXINFO('WRITE_STDX_FILE',ISTDX,NUMNP,NUMEL)
!
      IF (NUMDX.EQ.0) RETURN
!
!.... OPEN-DX DATA FILES FOR TRANSIENT FIELDS
!
      CALL PRINT_DXINFO('OPEN_FEMDX_FILE',IFEDX,NUMNP,NUMEL)
!
      RETURN
!
1900  FORMAT(2X,6(1PE15.8,2X)) 
1901  FORMAT(a12,x,6(1PE15.8,2X)) 
2000  FORMAT(27I6) 
! 
      END SUBROUTINE
!
!**** NEW **** FOR DATA EXPLORER OUT PUT ************************************* 
!
      SUBROUTINE PRINT_DXINFO(TASK,IFILE,NUMNP,NUMEL)
!
      use mGlobaisEscalares, only: NUMDX, NNP, NVEL
!
!..... PROGRAM TO SET-UP AND WRITE DATA ON OPEN-DX FORMAT
!
      IMPLICIT REAL*8 (A-H,O-Z)
!
      CHARACTER*15 TASK
!
      INTEGER JJ, IFILE, NUMNP, NUMEL, NINDX, NNSTEP, NTINDX
!
      IF (NUMDX.EQ.0) RETURN

      IF(TASK=='OPEN_FEMDX_FILE') THEN
! 
!.... PRINT NODAL AND MESH INFORMATION FOR DATA EXPLORER FILE
! 
	WRITE(IFILE,1000) '## OpenDX format File' 
	WRITE(IFILE,1000) '## OutPut Data  at Nodal Points in the' 
	WRITE(IFILE,1000) '## sense of Finite Element Method'
        WRITE(IFILE,1000) '##=========================================='
        WRITE(IFILE,1000) '##=========================================='
! 
!.... Nodes of finite element mesh 
!
	WRITE(IFILE,1000) '# ' 
	WRITE(IFILE,1000) '## Nodes locations'
	WRITE(IFILE,1000) '# '  
	WRITE(IFILE,1500)  &
     & 'object 1 class array type float rank 1 shape 2 items ', &
     & NUMNP,' data file "fnodes.stoc"' 
!
!..... Conectivity of finite element mesh
!
	WRITE(IFILE,1000) '# ' 
	WRITE(IFILE,1000) '## Connectivity' 
	WRITE(IFILE,1000) '# ' 
	WRITE(IFILE,1500) &
     & 'object 2 class array type int rank 1 shape 4 items ', &
     & NUMEL,' data file "femesh.stoc"' 
	WRITE(IFILE,1000) 'attribute "element type" string "quads"' 
	WRITE(IFILE,1000) 'attribute "ref" string "positions"' 
	WRITE(IFILE,1000) '#  '
! 
      ENDIF
!
      IF(TASK=='WRITE_STDX_FILE') THEN
! 
	WRITE(IFILE,1000)'# Scalar field : Young Modulus'
	WRITE(IFILE,1800) 3, NUMEL,'youngmd'
!
	WRITE(IFILE,1000)'# Scalar field : Poisson Ratio'
	WRITE(IFILE,1800) 4, NUMEL, 'poisson'
!
	WRITE(IFILE,1000)'# Scalar field : Permeability'
	WRITE(IFILE,1800) 5, NUMEL, 'permeab'
!
	WRITE(IFILE,1000)'# Scalar field : Geo Prssure'
	WRITE(IFILE,1800) 6, NUMEL, 'geoprsr'
! 
	WRITE(IFILE,1000)'# Vector field : Displacements' 
	WRITE(IFILE,1900) 7, NUMNP, 'displac'
!
	WRITE(IFILE,1000)'# Scalar field : Init. Stress XX'
	WRITE(IFILE,1800) 8, NUMEL, 'strs0xx'
! 
	WRITE(IFILE,1000)'# Scalar field : Init. Stress YY'
	WRITE(IFILE,1800) 9, NUMEL, 'strs0yy'
! 
	WRITE(IFILE,1000)'# Scalar field : Init. Stress XY'
	WRITE(IFILE,1800) 10, NUMEL, 'strs0xy'
!
	WRITE(IFILE,1000)'# Scalar field : Init. Stress ZZ'
	WRITE(IFILE,1800) 11, NUMEL, 'strs0zz'
! 
!...... YOUNG FIELD INFORMATION 
! 
	WRITE(IFILE,1000)'# Next object is a member of the: '
	WRITE(IFILE,1000)'#   Scalar YOUNG series'
	WRITE(IFILE,3000) 12, 3
! 
!...... POISSON FIELD INFORMATION 
! 
	WRITE(IFILE,1000)'# Next object is a member of the: '
	WRITE(IFILE,1000)'#   Scalar POISSON series'
	WRITE(IFILE,3000) 13, 4
! 
!...... PERMEABILITY FIELD INFORMATION 
! 
	WRITE(IFILE,1000)'# Next object is a member of the: '
	WRITE(IFILE,1000)'#   Scalar PERMEABILITY series'
	WRITE(IFILE,3000) 14, 5
!
!...... GEO PRESSURE FIELD INFORMATION 
! 
	WRITE(IFILE,1000)'# Next object is a member of the: '
	WRITE(IFILE,1000)'#   Scalar GEO PRESSURE series'
	WRITE(IFILE,3000) 15, 6
!
!...... DISPLACEMENTS FIELD INFORMATION 
! 
	WRITE(IFILE,1000)'# Next object is a member of the: '
	WRITE(IFILE,1000)'#   Scalar DISPLACEMENTS series'
	WRITE(IFILE,3000) 16, 7
!
!...... INITIAL STRESS XX FIELD INFORMATION 
! 
	WRITE(IFILE,1000)'# Next object is a member of the: '
	WRITE(IFILE,1000)'#   Scalar INIT STRESS XX series'
	WRITE(IFILE,3000) 17, 8
!
!...... INITIAL STRESS YY FIELD INFORMATION 
! 
	WRITE(IFILE,1000)'# Next object is a member of the: '
	WRITE(IFILE,1000)'#   Scalar INIT STRESS YY series'
	WRITE(IFILE,3000) 18, 9
!
!...... INITIAL STRESS XX FIELD INFORMATION 
! 
	WRITE(IFILE,1000)'# Next object is a member of the: '
	WRITE(IFILE,1000)'#   Scalar INIT STRESS XY series'
	WRITE(IFILE,3000) 19, 10
!
!...... INITIAL STRESS XX FIELD INFORMATION 
! 
	WRITE(IFILE,1000)'# Next object is a member of the: '
	WRITE(IFILE,1000)'#   Scalar INIT STRESS ZZ series'
	WRITE(IFILE,3000) 20, 11
!
        WRITE(IFILE,1000)'#  ' 
        WRITE(IFILE,1000)'# Here we create the YOUNG MODULUS serie object'
        WRITE(IFILE,1000)'object "young" class series'
        WRITE(IFILE,7000) 0,12,0
!
        WRITE(IFILE,1000)'#  ' 
        WRITE(IFILE,1000)'# Here we create the POISSON RATIO serie object'
        WRITE(IFILE,1000)'object "poisson" class series'
        WRITE(IFILE,7000) 0,13,0
!
        WRITE(IFILE,1000)'#  ' 
        WRITE(IFILE,1000)'# Here we create the PERMEABILITY serie object'
        WRITE(IFILE,1000)'object "permeability" class series'
        WRITE(IFILE,7000) 0,14,0
!
        WRITE(IFILE,1000)'#  ' 
        WRITE(IFILE,1000)'# Here we create the GEO PRESSURE serie object'
        WRITE(IFILE,1000)'object "geoprsr" class series'
        WRITE(IFILE,7000) 0,15,0
!
        WRITE(IFILE,1000)'#  ' 
        WRITE(IFILE,1000)'# Here we create the DISPLACEMENTS serie object'
        WRITE(IFILE,1000)'object "displacements" class series'
        WRITE(IFILE,7000) 0,16,0
!
        WRITE(IFILE,1000)'#  ' 
        WRITE(IFILE,1000)'# Here we create the INIT STRESS XX serie object'
        WRITE(IFILE,1000)'object "strs0xx" class series'
        WRITE(IFILE,7000) 0,17,0
!
        WRITE(IFILE,1000)'#  ' 
        WRITE(IFILE,1000)'# Here we create the INIT STRESS YY serie object'
        WRITE(IFILE,1000)'object "strs0yy" class series'
        WRITE(IFILE,7000) 0,18,0
!
        WRITE(IFILE,1000)'#  ' 
        WRITE(IFILE,1000)'# Here we create the INIT STRESS XY serie object'
        WRITE(IFILE,1000)'object "strs0xy" class series'
        WRITE(IFILE,7000) 0,19,0
!
        WRITE(IFILE,1000)'#  ' 
        WRITE(IFILE,1000)'# Here we create the INIT STRESS ZZ serie object'
        WRITE(IFILE,1000)'object "strs0zz" class series'
        WRITE(IFILE,7000) 0,20,0
!
        WRITE(IFILE,1000)'#  ' 
        WRITE(IFILE,1000)'# Structure of VARIAVEL OF DATA FILE'
        WRITE(IFILE,1000)'object "campos" class group'
        WRITE(IFILE,1000)'member "young" value "young"' 
        WRITE(IFILE,1000)'member "poisson" value "poisson"' 
        WRITE(IFILE,1000)'member "permeability" value "permeability"' 
        WRITE(IFILE,1000)'member "geoprsr" value "geoprsr"' 
        WRITE(IFILE,1000)'member "displacements" value "displacements"'
        WRITE(IFILE,1000)'member "strs0xx" value "strs0xx"' 
        WRITE(IFILE,1000)'member "strs0yy" value "strs0yy"' 
        WRITE(IFILE,1000)'member "strs0xy" value "strs0xy"' 
        WRITE(IFILE,1000)'member "strs0zz" value "strs0zz"' 
        WRITE(IFILE,1000)'#  ' 
        WRITE(IFILE,1000)'end'
        WRITE(IFILE,1000)'#  ' 
!
      ENDIF
!
      IF(TASK=='WRITE_FEDX_DATA') THEN
!
!..... PRINT NODAL DATA FOR OPEN-DX FILE
! 
        NINDX = NNP/NUMDX
!
        NNSTEP=24*(NINDX)+3
! 
!...... DISPLACEMENTS: VECTOR FIELD
! 
	WRITE(IFILE,1000)'# Vector field : DISPLACEMENTS' 
	WRITE(IFILE,2000) NNSTEP, NUMNP, NINDX 
! 
!...... DISPLACEMENT SERIES INFORMATION 
! 
	WRITE(IFILE,1000)'# Next object is a member of the: '
	WRITE(IFILE,1000)'#   DISPLACEMENT series'
	WRITE(IFILE,3000) NNSTEP+1, NNSTEP
! 
!.....  PRESSURE: SCALAR FIELD
! 
	WRITE(IFILE,1000)'# Scalar field : PRESSURE' 
	WRITE(IFILE,4000) NNSTEP+2, NUMEL, NINDX 
! 
!...... PRESSURE SERIES INFORMATION 
! 
	WRITE(IFILE,1000)'# Next object is a member of the: '
	WRITE(IFILE,1000)'#   Scalar PRESSURE series'
	WRITE(IFILE,3000) NNSTEP+3, NNSTEP+2
! 
!.....  GEOMECHANIC POROSITY: PORE SCALAR FIELD
! 
	WRITE(IFILE,1000)'# Scalar field : PORE' 
	WRITE(IFILE,4010) NNSTEP+4, NUMEL, NINDX 
! 
!...... GEOMECHANIC POROSITY: PORE SERIES INFORMATION 
! 
	WRITE(IFILE,1000)'# Next object is a member of the: '
	WRITE(IFILE,1000)'#   Scalar PORE series'
	WRITE(IFILE,3000) NNSTEP+5, NNSTEP+4
! 
!.....  CREEP SCALAR FIELD
! 
	WRITE(IFILE,1000)'# Scalar field : ' 
	WRITE(IFILE,4020) NNSTEP+6, NUMEL, NINDX 
! 
!...... CREEP SERIES INFORMATION 
! 
	WRITE(IFILE,1000)'# Next object is a member of the: '
	WRITE(IFILE,1000)'#   Scalar CREEP series'
	WRITE(IFILE,3000) NNSTEP+7, NNSTEP+6
! 
!.....  SATURATION: SCALAR FIELD
! 
	WRITE(IFILE,1000)'# Scalar field : SATURATION' 
	WRITE(IFILE,4030) NNSTEP+8, NUMEL, NINDX 
! 
!...... SATURATION SERIES INFORMATION 
! 
	WRITE(IFILE,1000)'# Next object is a member of the: '
	WRITE(IFILE,1000)'#   Scalar SATURATION series'
	WRITE(IFILE,3000) NNSTEP+9, NNSTEP+8
! 
!.....  VELOCITY VECTOR FIELD
! 
	WRITE(IFILE,1000)'# Vector field: VELOCITY' 
	WRITE(IFILE,4040) NNSTEP+10, NUMEL, NINDX 
! 
!...... VELOCITY SERIES INFORMATION 
! 
	WRITE(IFILE,1000)'# Next object is a member of the: '
	WRITE(IFILE,1000)'#   Vector Field VELOCITY series'
	WRITE(IFILE,3000) NNSTEP+11, NNSTEP+10
! 
!.....  PERMEABILITY: SCALAR FIELD
! 
	WRITE(IFILE,1000)'# Scalar field : MASS CONTENT' 
	WRITE(IFILE,4050) NNSTEP+12, NUMEL, NINDX 
! 
!...... PERMEABILITY SERIES INFORMATION 
! 
	WRITE(IFILE,1000)'# Next object is a member of the: '
	WRITE(IFILE,1000)'#   Scalar MASS CONTENT series'
	WRITE(IFILE,3000) NNSTEP+13, NNSTEP+12
! 
!.....  STRESS ON X DIRECTION
! 
	WRITE(IFILE,1000)'# Scalar field : Stress_X' 
	WRITE(IFILE,4060) NNSTEP+14, NUMEL, NINDX 
! 
!...... STRESS_X SERIES INFORMATION 
! 
	WRITE(IFILE,1000)'# Next object is a member of the: '
	WRITE(IFILE,1000)'#   Scalar Stress_X series'
	WRITE(IFILE,3000) NNSTEP+15, NNSTEP+14
! 
!.....  STRESS ON Y DIRECTION
! 
	WRITE(IFILE,1000)'# Scalar field : Stress_Y' 
	WRITE(IFILE,4070) NNSTEP+16, NUMEL, NINDX 
! 
!...... STRESS_Y SERIES INFORMATION 
! 
	WRITE(IFILE,1000)'# Next object is a member of the: '
	WRITE(IFILE,1000)'#   Scalar Stress_Y series'
	WRITE(IFILE,3000) NNSTEP+17, NNSTEP+16
! 
!.....  SHEAR STRESS XY
! 
	WRITE(IFILE,1000)'# Scalar field : Stress_XY'
	WRITE(IFILE,4080) NNSTEP+18, NUMEL, NINDX 
! 
!...... SHEAR STRESS_XY SERIES INFORMATION 
! 
	WRITE(IFILE,1000)'# Next object is a member of the: '
	WRITE(IFILE,1000)'#   Scalar Stress_XY series'
	WRITE(IFILE,3000) NNSTEP+19, NNSTEP+18
! 
!.....  STRESS Z 
! 
	WRITE(IFILE,1000)'# Scalar field : Stress_Z'
	WRITE(IFILE,4090) NNSTEP+20, NUMEL, NINDX 
! 
!...... SHEAR STRESS_XY SERIES INFORMATION 
! 
	WRITE(IFILE,1000)'# Next object is a member of the: '
	WRITE(IFILE,1000)'#   Scalar Stress_Z series'
	WRITE(IFILE,3000) NNSTEP+21, NNSTEP+20
! 
!.....  PRINCIPAL STRESS RATIO S2/S1 
! 
	WRITE(IFILE,1000)'# Scalar field : Princ. Stress'
	WRITE(IFILE,4100) NNSTEP+22, NUMEL, NINDX 
! 
!...... PRINCIPAL STRESS SERIES INFORMATION 
! 
	WRITE(IFILE,1000)'# Next object is a member of the: '
	WRITE(IFILE,1000)'#   Scalar Princ. Stress series'
	WRITE(IFILE,3000) NNSTEP+23, NNSTEP+22
!
      ENDIF
!
      IF(TASK=='CLOSE_FEDX_FILE') THEN
!
!.... PRINT SERIES LINKS INFORMATION FOR DATA EXPLORER FILE
!
      NTINDX=(NVEL/NUMDX)+1
!
      WRITE(IFILE,1000)'#  ' 
      WRITE(IFILE,1000)'# Here we create the DISPLACEMENT series object'
      WRITE(IFILE,1000)'object "displacement" class series'
      DO 301 JJ=1,NTINDX
         WRITE(IFILE,7000) JJ-1,24*JJ-20,JJ-1
 301  CONTINUE
!
      WRITE(IFILE,1000)'#  ' 
      WRITE(IFILE,1000)'# Here we create the PRESSURE series object'
      WRITE(IFILE,1000)'object "pressure" class series'
      DO 302 JJ=1,NTINDX
         WRITE(IFILE,7000) JJ-1,24*JJ-18,JJ-1
 302  CONTINUE
!
      WRITE(IFILE,1000)'#  ' 
      WRITE(IFILE,1000)'# Here we create the PORE series object'
      WRITE(IFILE,1000)'object "pore" class series'
!
      DO 303 JJ=1,NTINDX
         WRITE(IFILE,7000) JJ-1,24*JJ-16,JJ-1
 303  CONTINUE
      WRITE(IFILE,1000)'#  ' 
      WRITE(IFILE,1000)'# Here we create the CREEP series object'
      WRITE(IFILE,1000)'object "creep" class series'
      DO 304 JJ=1,NTINDX
         WRITE(IFILE,7000) JJ-1,24*JJ-14,JJ-1
 304  CONTINUE
!
      WRITE(IFILE,1000)'#  ' 
      WRITE(IFILE,1000)'# Here we create the SATURATION series object'
      WRITE(IFILE,1000)'object "saturation" class series'
      DO 305 JJ=1,NTINDX
         WRITE(IFILE,7000) JJ-1,24*JJ-12,JJ-1
 305  CONTINUE
!
      WRITE(IFILE,1000)'#  ' 
      WRITE(IFILE,1000)'# Here we create the VELOCITY series object'
      WRITE(IFILE,1000)'object "velocity" class series'
      DO 306 JJ=1,NTINDX
         WRITE(IFILE,7000) JJ-1,24*JJ-10,JJ-1
 306  CONTINUE
!
      WRITE(IFILE,1000)'#  ' 
      WRITE(IFILE,1000)'# Here we create the MASS CONTENT series object'
      WRITE(IFILE,1000)'object "mass_content" class series'
      DO 307 JJ=1,NTINDX
         WRITE(IFILE,7000) JJ-1,24*JJ-8,JJ-1
 307  CONTINUE
!
      WRITE(IFILE,1000)'#  ' 
      WRITE(IFILE,1000)'# Here we create the STRESS X series object'
      WRITE(IFILE,1000)'object "stress_x" class series'
      DO 308 JJ=1,NTINDX
         WRITE(IFILE,7000) JJ-1,24*JJ-6,JJ-1
 308  CONTINUE
!
      WRITE(IFILE,1000)'#  ' 
      WRITE(IFILE,1000)'# Here we create the STRESS Y series object'
      WRITE(IFILE,1000)'object "stress_y" class series'
      DO 309 JJ=1,NTINDX
        WRITE(IFILE,7000) JJ-1,24*JJ-4,JJ-1
 309  CONTINUE
!
      WRITE(IFILE,1000)'#  ' 
      WRITE(IFILE,1000)'# Here we create the SHEAR STRESS series object'
      WRITE(IFILE,1000)'object "stress_xy" class series'
      DO 310 JJ=1,NTINDX
         WRITE(IFILE,7000) JJ-1,24*JJ-2,JJ-1
 310  CONTINUE
!
      WRITE(IFILE,1000)'#  ' 
      WRITE(IFILE,1000)'# Here we create the STRESS Z series object'
      WRITE(IFILE,1000)'object "stress_z" class series'
      DO 311 JJ=1,NTINDX
         WRITE(IFILE,7000) JJ-1,24*JJ,JJ-1
 311  CONTINUE
!
      WRITE(IFILE,1000)'#  ' 
      WRITE(IFILE,1000)'# Here we create the PRINC. STRESS serie object'
      WRITE(IFILE,1000)'object "s2s1" class series'
      DO 312 JJ=1,NTINDX
         WRITE(IFILE,7000) JJ-1,24*JJ+2,JJ-1
 312  CONTINUE
!
      WRITE(IFILE,1000)'#  ' 
      WRITE(IFILE,1000)'# Structure of VARIAVEL OF DATA FILE'
      WRITE(IFILE,1000)'object "campos" class group'
      WRITE(IFILE,1000)'member "displacement" value "displacement"'
      WRITE(IFILE,1000)'member "pressure" value "pressure"'
      WRITE(IFILE,1000)'member "pore" value "pore"'
      WRITE(IFILE,1000)'member "creep" value "creep"'
      WRITE(IFILE,1000)'member "saturation" value "saturation"'
      WRITE(IFILE,1000)'member "velocity" value "velocity"' 
      WRITE(IFILE,1000)'member "mass_content" value "mass_content"' 
      WRITE(IFILE,1000)'member "stress_x" value "stress_x"' 
      WRITE(IFILE,1000)'member "stress_y" value "stress_y"' 
      WRITE(IFILE,1000)'member "stress_xy" value "stress_xy"' 
      WRITE(IFILE,1000)'member "stress_z" value "stress_z"' 
      WRITE(IFILE,1000)'member "s2s1" value "s2s1"' 
      WRITE(IFILE,1000)'#  ' 
      WRITE(IFILE,1000)'end'
      WRITE(IFILE,1000)'#  ' 
!
      ENDIF
!
      IF(TASK=='OTHERS__ANOTHER') THEN
!
      ENDIF
!
!.... FORMATOS DE SAIDA  OPEN-DX
!
 1000 FORMAT(A) 
 1500 FORMAT(A,I7,A)
 1800 FORMAT('object ',I3,' class array type float rank 0 items ',I8, &
     &' data file "',A7,'.dat"'/                    &
     &' attribute "dep" string "connections"'/'#  ') 
 1900 FORMAT('object ',I3,                                 &
     & ' class array type float rank 1 shape 2 items', I8, &
     &' data file "',A7,'.dat"'/                    &
     &'attribute "dep" string "positions"'/'#  ') 
!
 2000 FORMAT('object ',I5,                                 &
     & ' class array type float rank 1 shape 2 items', I8, &
     &' data file "disp',I3.3,'.stoc"'/                    &
     &'attribute "dep" string "positions"'/'#  ') 
!
 3000 FORMAT('object ',I5,' class field'/    &
     &'component "positions" value 1'/       &
     &'component "connections" value 2'/     &
     &'component "data" value ',I5/'#  ')
!
 4000 FORMAT('object ',I5,                     &
     &' class array type float rank 0 items ', &
     & I8,' data file "prsr',I3.3,'.stoc"'/    &
     &'attribute "dep" string "connections"'/'#  ') 
!
 4010 FORMAT('object ',I5,                       &
     &' class array type float rank 0 items ',   &
     & I8,' data file "pore',I3.3,'.stoc"'/      &
     &'attribute "dep" string "connections"'/'#  ') 
!
 4020 FORMAT('object ',I5,                       &
     &' class array type float rank 0 items ',   &
     & I8,' data file "crep',I3.3,'.stoc"'/      &
     &'attribute "dep" string "connections"'/'#  ') 
!
 4030 FORMAT('object ',I5,                       &
     &' class array type float rank 0 items ',   &
     & I8,' data file "satr',I3.3,'.stoc"'/      &
     &'attribute "dep" string "connections"'/'#  ') 
!
 4040 FORMAT('object ',I5,                                 &
     & ' class array type float rank 1 shape 2 items', I8, &
     &' data file "velt',I3.3,'.stoc"'/                    &
     &'attribute "dep" string "connections"'/'#  ') 
!
 4050 FORMAT('object ',I5,                       &
     &' class array type float rank 0 items ',   &
     & I8,' data file "masc',I3.3,'.stoc"'/      &
     &'attribute "dep" string "connections"'/'#  ') 
!
 4060 FORMAT('object ',I5,                       &
     &' class array type float rank 0 items ',   &
     & I8,' data file "sigx',I3.3,'.stoc"'/      &
     &'attribute "dep" string "connections"'/'#  ') 
!
 4070 FORMAT('object ',I5,                       &
     &' class array type float rank 0 items ',   &
     & I8,' data file "sigy',I3.3,'.stoc"'/      &
     &'attribute "dep" string "connections"'/'#  ') 
!
4080  FORMAT('object ',I5,                         &
     & ' class array type float rank 0 items', I8, &
     &' data file "sigt',I3.3,'.stoc"'/            &
     &'attribute "dep" string "connections"'/'#  ') 
!
4090  FORMAT('object ',I5,                         &
     & ' class array type float rank 0 items', I8, &
     &' data file "sigz',I3.3,'.stoc"'/            &
     &'attribute "dep" string "connections"'/'#  ') 
!
4100  FORMAT('object ',I5,                          &
     & ' class array type float rank 0 items', I8,  &
     &' data file "s2s1',I3.3,'.stoc"'/             &
     &'attribute "dep" string "connections"'/'#  ') 
7000  FORMAT('member ',I5,' value ',I5,' position ',I5)
!
      END SUBROUTINE
!
!**** NEW **************************************************************** 
!
      SUBROUTINE CODERROR(IERRO,TASK)
!
!.... PROGRAM TO WRITE OUT ERROR
!
      IMPLICIT NONE
!
      INTEGER I, IERRO
!
      CHARACTER(LEN=18) :: TASK
      CHARACTER*18, DIMENSION(6 ) :: REFTASK
!
      DATA   REFTASK(1)     ,     REFTASK(2)     ,     REFTASK(3)     / &
      & 'BULK_MODULUS_GRAIN','YOUNG_MODULUS_LESS','POROSITY__NEGATIVE'/ &
      &     REFTASK(4)      ,     REFTASK(5)     ,     REFTASK(6)     / &
      & 'MASS_CONTENT_NEGAT','GEOMECH_PARAMETERS','INCONSISTENT__DATA'/
!
       WRITE(*,*) '---------------'
       WRITE(*,*) 'LOG ERROR FILE '
       WRITE(*,*) '---------------'
!
!..... TEST STEP
!
       GOTO (100,200,300,400,500,600), IERRO
!
 100  CONTINUE
       WRITE(*,*) '   '
       WRITE(*,*) 'PROBLEM: '//REFTASK(IERRO)//' IN '//TASK
       WRITE(*,*) '   '
!
       STOP
!
 200  CONTINUE
       WRITE(*,*) '   '
       WRITE(*,*) 'PROBLEM: '//REFTASK(IERRO)//' IN '//TASK
       WRITE(*,*) '   '
!
       STOP
!
 300  CONTINUE
!
       WRITE(*,*) '   '
       WRITE(*,*) 'LINEAR APPROXIMATION OF POROSITY FURNISH'
       WRITE(*,*) 'NEGATIVE VALUES, ORIGIN: EXTERNAL LOAD TO GREAT'
       WRITE(*,*) '   '
!
       STOP
!
400   CONTINUE
!
       WRITE(*,*) '   '
       WRITE(*,*) 'LINEAR APPROXIMATION OF MASS CONTENT OR POROSITY'
       WRITE(*,*) 'FURNISH NEGATIVE VALUES, PROBLEM ORIGIN: '
       WRITE(*,*) 'EXTERNAL LOAD GREATER THAN YOUNG MODULUS;'
       WRITE(*,*) 'INJECTION RATE TO GREATER'
       WRITE(*,*) '   '
!
       STOP
!
500   CONTINUE
!
       WRITE(*,*) '   '
       WRITE(*,*) 'PROBLEM READING FILES THAT CONTAINING'
       WRITE(*,*) 'GEOMECHANICAL PARAMETERS FOR EACH FORMATIONS'
       WRITE(*,*) '   '
!
       STOP
!
600   CONTINUE
!
       WRITE(*,*) '   '
       WRITE(*,*) 'PROBLEM: '//REFTASK(IERRO)//' IN '//TASK
       WRITE(*,*) '   '
!
       STOP
!
       RETURN
!
      END SUBROUTINE
!
      END MODULE
!
!******** ************ ************ ************ *********** **********
!
