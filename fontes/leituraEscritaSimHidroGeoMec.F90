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
      module mLeituraEscritaSimHidroGeoMec
!
      implicit none
!
      integer :: iparaviewS,iparaviewP,iparaviewV
      integer :: isatTransiente
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
      
      logical :: apenasReservatorio=.false.
!
      contains
      
      subroutine fecharArquivosSimHidroGeoMec()
!     
      close(isat)
      close(ipres)
      close(ivel)
      close(iphi)
      close(iperm)
      close(imass)
      close(imasl)
! 
      end subroutine fecharArquivosSimHidroGeoMec
!
!**** NEW *******************************************************************
!
      subroutine lerDataIn_DS
!
      use mPropGeoFisica
      use mGlobaisEscalares
      use mInputReader,    only: readIntegerKeywordValue,readRealKeywordValue,readOutFlagKeyword
      use mLeituraEscrita, only: iflag_tipoPrint
!
      implicit none
      character(len=50) :: keyword_name
      integer :: ierr
!
        keyword_name = "iflag_linear"
        call readIntegerKeywordValue(keyword_name, iflag_linear, iflag_linear, ierr)

        keyword_name = "viscosidade_da_agua"
        call readRealKeywordValue(keyword_name, xmiw, xmiw, ierr)

        keyword_name = "viscosidade_do_oleo"
        call readRealKeywordValue(keyword_name, xmio, xmio, ierr)

!         keyword_name = "massa_especifica_da_agua"
!         call readRealKeywordValue(keyword_name, rhow, rhow, ierr)
! !
!         keyword_name = "massa_especifica_do_oleo"
!         call readRealKeywordValue(keyword_name, rhoo, rhoo, ierr)

        keyword_name = "saturacao_residual_da_agua"
        call readRealKeywordValue(keyword_name, srw, srw, ierr)
!
        keyword_name = "saturacao_residual_do_oleo"
        call readRealKeywordValue(keyword_name, sro, sro, ierr)
!
        keyword_name = "saturacao_inicial_da_agua"
        call readRealKeywordValue(keyword_name, sinicial, sinicial, ierr)
!
        keyword_name = "saturacao_na_injecao"
        call readRealKeywordValue(keyword_name, sinj, sinj, ierr)
!
        keyword_name = "saturacao_do_bloco"
        call readRealKeywordValue(keyword_name, sbloco, sbloco, ierr)
!
        keyword_name = "x_central_do_bloco_de_sat"
        call readRealKeywordValue(keyword_name, xcbloco, xcbloco, ierr)
!
        keyword_name = "y_central_do_bloco_de_sat"
        call readRealKeywordValue(keyword_name, ycbloco, ycbloco, ierr)
!
        keyword_name = "z_central_do_bloco_de_sat"
        call readRealKeywordValue(keyword_name, zcbloco, zcbloco, ierr)
!
        keyword_name = "lx_do_bloco_de_sat"
        call readRealKeywordValue(keyword_name, xlbloco, xlbloco, ierr)
!
        keyword_name = "ly_do_bloco_de_sat"
        call readRealKeywordValue(keyword_name, ylbloco, ylbloco, ierr)
!
        keyword_name = "lz_do_bloco_de_sat"
        call readRealKeywordValue(keyword_name, zlbloco, zlbloco, ierr)
!
        keyword_name = "permeability_reservoir_bottom_region"
        call readRealKeywordValue(keyword_name, perminicial, perminicial, ierr)
!
        keyword_name = "permeabilidade_do_bloco"
        call readRealKeywordValue(keyword_name, permbloco, permbloco, ierr)
!
        keyword_name = "x_central_do_bloco_de_perm"
        call readRealKeywordValue(keyword_name, xcbloco_perm, xcbloco_perm, ierr)
!
        keyword_name = "y_central_do_bloco_de_perm"
        call readRealKeywordValue(keyword_name, ycbloco_perm, ycbloco_perm, ierr)
!
        keyword_name = "z_central_do_bloco_de_perm"
        call readRealKeywordValue(keyword_name, zcbloco_perm, zcbloco_perm, ierr)
!
        keyword_name = "lx_do_bloco_de_perm"
        call readRealKeywordValue(keyword_name, xlbloco_perm, xlbloco_perm, ierr)
!
        keyword_name = "ly_do_bloco_de_perm"
        call readRealKeywordValue(keyword_name, ylbloco_perm, ylbloco_perm, ierr)
!
        keyword_name = "lz_do_bloco_de_perm"
        call readRealKeywordValue(keyword_name, zlbloco_perm, zlbloco_perm, ierr)
!
        keyword_name = "porosidade_inicial"
        call readRealKeywordValue(keyword_name, phiinicial, phiinicial, ierr)
!
        keyword_name = "porosidade_do_bloco"
        call readRealKeywordValue(keyword_name, phibloco, phibloco, ierr)
!
        keyword_name = "x_central_do_bloco_de_phi"
        call readRealKeywordValue(keyword_name, xcbloco_phi, xcbloco_phi, ierr)
!
        keyword_name = "y_central_do_bloco_de_phi"
        call readRealKeywordValue(keyword_name, ycbloco_phi, ycbloco_phi, ierr)
!
        keyword_name = "z_central_do_bloco_de_phi"
        call readRealKeywordValue(keyword_name, zcbloco_phi, zcbloco_phi, ierr)
!
        keyword_name = "lx_do_bloco_de_phi"
        call readRealKeywordValue(keyword_name, xlbloco_phi, xlbloco_phi, ierr)
!
        keyword_name = "ly_do_bloco_de_phi"
        call readRealKeywordValue(keyword_name, ylbloco_phi, ylbloco_phi, ierr)
!
        keyword_name = "lz_do_bloco_de_phi"
        call readRealKeywordValue(keyword_name, zlbloco_phi, zlbloco_phi, ierr)
!
        keyword_name = "passo_para_impressao_opendx_files"
        call readIntegerKeywordValue(keyword_name, NUMDX, NUMDX, ierr)
!
!       flag="# impressao para matlab(0), paraview(1), paraview mrb(2) ou DX(3)"
        keyword_name = "impressao_para_matlab_paraview_paraviewMrb_DX"
        call readIntegerKeywordValue(keyword_name, iflag_tipoPrint, iflag_tipoPrint, ierr)
!
        keyword_name = "impressao_da_saturacao"
        call readOutFlagKeyword(keyword_name, iflag_sat, ifsat_out, npsat, ierr)
!        
        keyword_name = "impressao_da_pressao"
        call readOutFlagKeyword(keyword_name, iflag_pres, ifpres_out, nppres, ierr)
!
        keyword_name = "impressao_da_velocidade"
        call readOutFlagKeyword(keyword_name, iflag_vel, ifvel_out, npvel, ierr)
!
        keyword_name = "impressao_da_porosidade"
        call readOutFlagKeyword(keyword_name, iflag_phi, ifphi_out, npphi, ierr)
!
        keyword_name = "impressao_dos_deslocamentos"
        call readOutFlagKeyword(keyword_name, iflag_disp, ifdis_out, npdis, ierr)
!
        !Read npmass
        keyword_name = "calculo_da_massa_de_agua"
        call readOutFlagKeyword(keyword_name, iflag_mass, ifmass_out, npmass, ierr)

        !Read npmasl
        keyword_name = "balanco_local_da_massa_de_agua"
        call readOutFlagKeyword(keyword_name, iflag_masl, ifmasl_out, npmasl, ierr)
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
      end subroutine lerDataIn_DS
!
!===========================================================================
!
      subroutine lerRandfilesIn_DS
!
      use mPropGeoFisica
      use mMalha, only: nsd
      use mGlobaisEscalares, only: geomech, iflag_beta
      use mInputReader, only: readIntegerKeywordValue,readRealKeywordValue
      use mInputReader, only: readOutFlagKeyword
!
      implicit none

      character(len=50) :: keyword_name
      integer :: ierr
!
!.... Numero de realizacoes
!
      keyword_name="realizacoes"
      call readIntegerKeywordValue(keyword_name, nrand, nrand, ierr)
!
!.... Kozeny-Carman relation
      call readRandKozenyCarman
!
!.... Leitura da permeabilidade
!
      call readRandPermeability
!
!.... Leitura da porosidade
!
      call readRandPorosity
!
!.... Leitura da porosidade
!
      call readRandBeta
!
!.... Leitura de Modulo de Young do Reservatorio
!
      if(geomech==1) then
         keyword_name="modulo_de_young_reservatorio"
         call readRandYoungModule(keyword_name, IFLAG_READ_YNG,YNG_IN,KGYNG, RHOYNG)
      endif
!
!.... Leitura de Modulo de Young Dominio
!
      if(geomech==1) then
         keyword_name="modulo_de_young_dominio_meio"
         call readRandYoungModuleCAP(keyword_name, IFLAG_READ_YNG2,YNG2_IN,KGYNG2,RHOYNG2,YNGX1,YNGX2)
         
         keyword_name="modulo_de_young_dominio_top"
         call readRandYoungModuleCAP(keyword_name, IFLAG_READ_YNG3,YNG3_IN,KGYNG3,RHOYNG3,YNG3X1,YNG3X2)
      endif
!
!.... Saida do mixing
!
      call readRandMixing
!
!.... Saida da producao
!
      call readRandProduction
!
    contains

            !> Efetua a leitura de dados de permeabilidade
    subroutine readRandKozenyCarman
!
        use mInputReader,    only: file_lines, findKeyword      
        
        implicit none
        
        integer keyword_line
        character(len=50) keyword_name
        keyword_name = "Kozeny-Carman"
        keyword_line = findKeyword(keyword_name)
        if (keyword_line.eq.0) then
            return
        end if
        read(file_lines(keyword_line:), *) nkozenycarman
        keyword_line = keyword_line + 1
        read(file_lines(keyword_line:),*) sup_KC,const_KC
    end subroutine readRandKozenyCarman 
    
        !> Efetua a leitura de dados de permeabilidade
    subroutine readRandPermeability
!
        use mInputReader,    only: file_lines, findKeyword      
        
        implicit none
        
        integer keyword_line
        character(len=50) keyword_name
        keyword_name = "permeabilidade"
        keyword_line = findKeyword(keyword_name)
        if (keyword_line.eq.0) then
            return
        end if
        read(file_lines(keyword_line:), *) iflag_read_perm
        keyword_line = keyword_line + 1
        read(file_lines(keyword_line:),"(a)") perm_inx
        keyword_line = keyword_line + 1
        read(file_lines(keyword_line:),"(a)") perm_iny
        keyword_line = keyword_line + 1
        if (nsd==3) then
            read(file_lines(keyword_line:),"(a)") perm_inz
            keyword_line = keyword_line + 1
        endif
        perm_inx=trim(perm_inx)
        perm_iny=trim(perm_iny)
        if(nsd==3)  perm_inz=trim(perm_inz)
        read(file_lines(keyword_line:),*) kg, rho
    end subroutine readRandPermeability !*****************************************************************************

        !> Efetua a leitura de dados de porosidade
    subroutine readRandPorosity
        use mInputReader,    only: file_lines, findKeyword
        integer keyword_line
        character(len=50) keyword_name
        keyword_name = "porosidade"
        keyword_line = findKeyword(keyword_name)
        if (keyword_line.eq.0) then
            return
        end if
        read(file_lines(keyword_line:), *) iflag_read_phi
        keyword_line = keyword_line + 1
        read(file_lines(keyword_line:),"(a)") phi_in
        keyword_line = keyword_line + 1
        phi_in=trim(phi_in)
        read(file_lines(keyword_line:),*) kgphi, rhophi
        keyword_line = keyword_line + 1
        read(file_lines(keyword_line:),*) normalphi
    end subroutine readRandPorosity !*****************************************************************************
    
        !> Efetua a leitura de dados de beta
    subroutine readRandBeta
        use mInputReader,    only: file_lines, findKeyword
        integer keyword_line
        character(len=50) keyword_name
        keyword_name = "beta"
        keyword_line = findKeyword(keyword_name)
        if (keyword_line.eq.0) then
            return
        end if
        read(file_lines(keyword_line:), *) iflag_beta
        keyword_line = keyword_line + 1
        read(file_lines(keyword_line:),"(a)") beta_in
        keyword_line = keyword_line + 1
        beta_in=trim(beta_in)
        read(file_lines(keyword_line:),*) kgbeta, rhobeta
    end subroutine readRandBeta !*****************************************************************************

    !> Efetua a leitura de dados de Modulo de Young
    subroutine readRandYoungModule(keyword_name, IFLAG_READ_YNG,yng_in,KGYNG, RHOYNG)
        use mInputReader,    only: file_lines, findKeyword
        integer keyword_line
        character(len=50) keyword_name
        integer :: IFLAG_READ_YNG
        real*8 :: KGYNG, RHOYNG
        CHARACTER(len=128) :: YNG_IN
        
        keyword_line = findKeyword(keyword_name)
        if (keyword_line.eq.0) then
            return
        end if
        read(file_lines(keyword_line:), *) IFLAG_READ_YNG
        keyword_line = keyword_line + 1
        read(file_lines(keyword_line:),"(a)") yng_in
        keyword_line = keyword_line + 1
        yng_in=trim(yng_in)
        read(file_lines(keyword_line:),*) KGYNG, RHOYNG
    end subroutine readRandYoungModule !*****************************************************************************

    
    !> Efetua a leitura de dados de Modulo de Young
    subroutine readRandYoungModuleCAP(keyword_name, IFLAG_READ_YNG,yng_in,KGYNG, RHOYNG,YNGX1,YNGX2)
        use mInputReader,    only: file_lines, findKeyword
        integer keyword_line
        character(len=50) keyword_name
        integer :: IFLAG_READ_YNG
        CHARACTER(len=128) :: yng_in
        real*8 :: KGYNG, RHOYNG,YNGX1,YNGX2
        
        keyword_line = findKeyword(keyword_name)
        if (keyword_line.eq.0) then
            return
        end if
        read(file_lines(keyword_line:), *) IFLAG_READ_YNG
        keyword_line = keyword_line + 1
        read(file_lines(keyword_line:),"(a)") yng_in
        keyword_line = keyword_line + 1
        yng_in=trim(yng_in)
        read(file_lines(keyword_line:),*) KGYNG, RHOYNG
        keyword_line = keyword_line + 1
        read(file_lines(keyword_line:),*) YNGX1,YNGX2
    end subroutine readRandYoungModuleCAP !*****************************************************************************
    
    !> Efetua a leitura de saida mixing
    subroutine readRandMixing
        use mInputReader,    only: file_lines, findKeyword
        integer keyword_line
        character(len=50) keyword_name
        keyword_name = "mixing"
        keyword_line = findKeyword(keyword_name)
        if (keyword_line.eq.0) then
            return
        end if
        read(file_lines(keyword_line:), *) iflag_mix
        keyword_line = keyword_line + 1
        read(file_lines(keyword_line:),"(a)") mixing_out
        keyword_line = keyword_line + 1
        mixing_out=trim(mixing_out)
        read(file_lines(keyword_line:),*) npmix
    end subroutine readRandMixing !************************************************************************************

    !> Efetua a leitura de saida mixing
    subroutine readRandProduction
        use mInputReader,    only: file_lines, findKeyword       
        integer keyword_line
        character(len=50) keyword_name
        keyword_name = "producao"
        keyword_line = findKeyword(keyword_name)
        if (keyword_line.eq.0) then
            return
        end if
        read(file_lines(keyword_line:), *) iflag_prod
        keyword_line = keyword_line + 1
        read(file_lines(keyword_line:),"(a)") prod_out
        keyword_line = keyword_line + 1
        prod_out=trim(prod_out)
        read(file_lines(keyword_line:),*) npprod
    end subroutine readRandProduction !********************************************************************************


    end subroutine lerRandfilesIn_DS
!
!**** *******************************************************************
!
      SUBROUTINE lerGeoMechParam_DS
!
      use mGlobaisEscalares, only: S3DIM, SALTCREEP
      use mPropGeoFisica,    only: YUNGVECT, POISVECT, RHODVECT
      use mPropGeoFisica,    only: GRAINBLK, PORELAGR, BULK
      use mPropGeoFisica,    only: MEANSATR, GEOMECLAW
      use mPropgeoFisica,    only: MEANRHOW, MEANRHOO, MEANDENS
      use mPropgeoFisica,    only: MEANBLKW, MEANBLKO, MEANBULK
      use mPropGeoFisica,    only: RHOW, RHOO, BULKWATER, BULKOIL
      use mInputReader, only: readIntegerKeywordValue,readRealKeywordValue
      use mInputReader, only: readStringKeywordValue,readOutFlagKeyword

!
      IMPLICIT NONE
!
      INTEGER :: IFILE, IREGION
! 
      CHARACTER*6 :: TEXT
      CHARACTER*18, DIMENSION(6) :: REGION 
!
      real(8) :: YOUNG, POISSON, GRBULK, MEANWATER, MEANOIL, SATURA
      
      character(len=50) :: keyword_name
      integer :: ierr
!
      REGION(1)="reservoir_region"
      REGION(2)="rift_under_region"
      REGION(3)="right_side_burden"
      REGION(4)="saline_cap_region"
      REGION(5)="left_side_burden"
      REGION(6)="post_salt_region"

      DO IREGION=1, 6

      keyword_name = "young_modulus_"//trim(REGION(IREGION))
      call readRealKeywordValue(keyword_name, YUNGVECT(IREGION), YUNGVECT(IREGION), ierr)
      YOUNG = YUNGVECT(IREGION)
!
      keyword_name = "poisson_ratio_"//trim(REGION(IREGION))
      call readRealKeywordValue(keyword_name, POISVECT(IREGION), POISVECT(IREGION), ierr)    
      POISSON = POISVECT(IREGION)
!
      keyword_name = "rock_density_"//trim(REGION(IREGION))
      call readRealKeywordValue(keyword_name, RHODVECT(IREGION), RHODVECT(IREGION), ierr)    
!
      keyword_name = "bulk_solid_grain_"//trim(REGION(IREGION))
      call readRealKeywordValue(keyword_name, GRAINBLK(IREGION), GRAINBLK(IREGION), ierr)    
      GRBULK = GRAINBLK(IREGION)
!
!..TEST BULK MODULUS ORDER OF SCHELETON AND GRAIN
!
      IF (GRBULK.LE.BULK(YOUNG,POISSON,S3DIM)) THEN 
         CALL CODERROR(1,REGION(IREGION))
      ENDIF
!
      keyword_name = "mean_water_density_"//trim(REGION(IREGION))
      call readRealKeywordValue(keyword_name, MEANRHOW(IREGION), MEANRHOW(IREGION), ierr)       
      IF (IREGION.EQ.1) RHOW=MEANRHOW(IREGION)
      MEANWATER = MEANRHOW(IREGION)
!
      keyword_name = "mean_oil_density_"//trim(REGION(IREGION))
      call readRealKeywordValue(keyword_name, MEANRHOO(IREGION), MEANRHOO(IREGION), ierr)   
      IF (IREGION.EQ.1) RHOO=MEANRHOO(IREGION)
      MEANOIL = MEANRHOO(IREGION)
!
      keyword_name = "mean_water_saturation_"//trim(REGION(IREGION))
      call readRealKeywordValue(keyword_name, MEANSATR(IREGION), MEANSATR(IREGION), ierr)
      SATURA = MEANSATR(IREGION)
!
      MEANDENS(IREGION) = SATURA*MEANWATER+(1.0D0-SATURA)*MEANOIL
      
!
      keyword_name = "mean_bulk_water_"//trim(REGION(IREGION))
      call readRealKeywordValue(keyword_name, MEANBLKW(IREGION), MEANBLKW(IREGION), ierr)
      MEANWATER = MEANBLKW(IREGION)
      IF (IREGION.EQ.1) BULKWATER = MEANBLKW(IREGION)
!
      keyword_name = "mean_bulk_oil_"//trim(REGION(IREGION))
      call readRealKeywordValue(keyword_name, MEANBLKO(IREGION), MEANBLKO(IREGION), ierr)
      MEANOIL = MEANBLKO(IREGION)
      IF (IREGION.EQ.1) BULKOIL = MEANBLKW(IREGION)
!
      MEANBULK(IREGION) = SATURA/MEANWATER+(1.0D0-SATURA)/MEANOIL
!
      keyword_name = "porosity_"//trim(REGION(IREGION))
      call readRealKeywordValue(keyword_name, PORELAGR(IREGION), PORELAGR(IREGION), ierr)
!..TEST BIOT COEFICIENT GREATER THAN POROSITY 
      GRBULK=GRBULK*(1.0D0-PORELAGR(IREGION))
!
      IF (GRBULK.LE.BULK(YOUNG,POISSON,S3DIM)) THEN 
         WRITE(*,*) '..TEST BIOT COEFICIENT GREATER THAN POROSITY '
         CALL CODERROR(1,REGION(IREGION))
      ENDIF
!
      keyword_name = "stress-strain_relation_"//trim(REGION(IREGION))
      call readStringKeywordValue(keyword_name, TEXT, ierr)

      IF (TRIM(TEXT(1:5)).EQ.'CREEP') THEN 
            SALTCREEP = .TRUE.
            GEOMECLAW(IREGION) = TEXT(1:5)
         ELSE
           IF (TRIM(TEXT(1:5)).NE.'ELAST') THEN
              WRITE(*,*) 'UNKNOWN MODEL: ',TEXT
              CALL CODERROR(5,REGION(IREGION))
           ENDIF
           GEOMECLAW(IREGION) = TEXT(1:5)
      ENDIF
      
      END DO
!
      RETURN
!
1000  FORMAT(A18)
!
    END SUBROUTINE lerGeoMechParam_DS
!
!===========================================================================
!
      subroutine lerNumericParam_DS
!
      use mPropGeoFisica,     only : TOLCREEP
      use mGlobaisEscalares,  only : NITGEO, NITHIDRO, SPLITT, IBBAR
      use mGlobaisEscalares,  only : MAXITERC, TOLSIGMA, TOLVELOC
      use mInputReader, only: readIntegerKeywordValue,readRealKeywordValue
      use mInputReader, only: readStringKeywordValue
!
      IMPLICIT NONE
!
      INTEGER I, NUMPARAM
! 
      CHARACTER(LEN=128) :: NAMEFILE, TEXT
      
      character(len=50) :: keyword_name
      integer :: ierr

!
      keyword_name="max_number_for_hidro-geomechanic_iterations"
      call readIntegerKeywordValue(keyword_name, NITGEO, NITGEO, ierr)
!
      keyword_name="tolerance_for_hidro-geomechacnic_convergenc" 
      call readRealKeywordValue(keyword_name, TOLSIGMA, TOLSIGMA, ierr)
!
      keyword_name="max_number_for_hidro-transport_iterations"
      call readIntegerKeywordValue(keyword_name, NITHIDRO, NITHIDRO, ierr)
!
      keyword_name="tolerance_for_hidro-transport_convergenc" 
      call readRealKeywordValue(keyword_name, TOLVELOC, TOLVELOC, ierr)
!
      keyword_name="splitting_procedure"
      call readStringKeywordValue(keyword_name, SPLITT, ierr)

      keyword_name="maximum_number_of_creep_iterations"
      call readIntegerKeywordValue(keyword_name, MAXITERC, MAXITERC, ierr)
!
      keyword_name="tolerance_for_creep_convergence" 
      call readRealKeywordValue(keyword_name, TOLCREEP, TOLCREEP, ierr)
!
      keyword_name="bbar_method"
      call readIntegerKeywordValue(keyword_name, IBBAR, IBBAR, ierr)
!
      RETURN
!
1000  FORMAT(A6)
!
    END SUBROUTINE lerNumericParam_DS
!
!===========================================================================
!
      subroutine lerCreepParam_DS
!
      use mPropGeoFisica,     only: DTCREEP, CREEPZERO, POWERN, SIGMAREF
      use mGlobaisEscalares,  only: NCREEP
      use mInputReader,       only: readIntegerKeywordValue,readRealKeywordValue
!
      IMPLICIT NONE
!
      INTEGER I, ICREEP
! 
      CHARACTER(LEN=128) :: NAMEFILE, TEXT
      character(len=50) :: keyword_name
      integer :: ierr
!
      keyword_name="creep_deformation_reference" 
      call readRealKeywordValue(keyword_name, CREEPZERO, CREEPZERO, ierr)
!
      keyword_name="time_increment_for_creep" 
      call readRealKeywordValue(keyword_name, DTCREEP, DTCREEP, ierr)
!
      keyword_name="stress_reference_for_creep" 
      call readRealKeywordValue(keyword_name, SIGMAREF, SIGMAREF, ierr)
!
      keyword_name="power_of_creep_law" 
      call readRealKeywordValue(keyword_name, POWERN, POWERN, ierr)
!
      keyword_name="creep_time_as_multiple_of_geotime"
      call readIntegerKeywordValue(keyword_name, NCREEP, NCREEP, ierr)
!
      RETURN
!
      END SUBROUTINE
!
!===========================================================================
!
      subroutine lerSimulatorParam_DS
!
      use mPropGeoFisica
      use mGlobaisEscalares
      use mInputReader, only: readIntegerKeywordValue,readRealKeywordValue
!
      IMPLICIT NONE
!
      INTEGER I, ISimulator
! 
      CHARACTER(LEN=128) :: NAMEFILE, TEXT, FLAG1
      character(len=50)  :: keyword_name
      integer :: ierr
!
      keyword_name="top_external_load_em_pascal"
      call readRealKeywordValue(keyword_name, XTERLOAD, XTERLOAD, ierr)
!
      keyword_name="sea_depth_in_meters"
      call readRealKeywordValue(keyword_name, SEADEPTH, SEADEPTH, ierr)
!
      keyword_name="instante_inicial"
      call readRealKeywordValue(keyword_name, tzero, tzero, ierr)
      tTransporte = tzero
!
      keyword_name="tempo_total_de_simulacao"
      call readRealKeywordValue(keyword_name, tt, tt, ierr)
!
      keyword_name="numero_de_calculos_da_velocidade"
      call readIntegerKeywordValue(keyword_name, nvel, nvel, ierr)
!
      keyword_name="tempo_de_inicio_de_injecao"
      call readRealKeywordValue(keyword_name, YEARINJ, YEARINJ, ierr)
!
!       keyword_name='pressao_de_referencia'
!       call readRealKeywordValue(keyword_name, PRESSAOREF, PRESSAOREF, ierr)
!       WRITE(*,*)'PRESSAO DE REFERENCIA:',PRESSAOREF
!       IF(ierr==1) THEN  
!          WRITE(*,*) 'ERRO NA LEITURA DA PRESSAO DE REFERENCIA'
!          STOP
!       END IF
!       WRITE(*,*)'########################################'
! !
!       keyword_name='cota_de_referencia'
!       call readRealKeywordValue(keyword_name, COTAREF, COTAREF, ierr)
!       IF(ierr.ne.1) THEN  
!          WRITE(*,*)'COTA DE REFERENCIA:',COTAREF
!       ELSE
!          WRITE(*,*) 'ERRO NA LEITURA DA COTA DE REFERENCIA'
!          STOP
!       END IF
! !
!       WRITE(*,*)'########################################'
! !
!       keyword_name='pressao_no_fundo_do_poco_de_producao'
!       call readRealKeywordValue(keyword_name, PRESPROD, PRESPROD, ierr)
!       IF(ierr.ne.1) THEN  
!          WRITE(*,*)'PRESSAO NO POCO DE PRODUCAO:',PRESPROD     
!       ELSE
!          WRITE(*,*) 'ERRO NA LEITURA DA PRESSAO DO POCCO'
!          STOP
!       END IF
!       WRITE(*,*)'########################################'
! !
!       keyword_name='cota_do_fundo_do_poco'
!       call readRealKeywordValue(keyword_name, ZFUNDOPOCCO, ZFUNDOPOCCO, ierr)
!       IF(ierr.ne.1) THEN  
!          WRITE(*,*)'COTA DO FUNDO DO POCCO:',ZFUNDOPOCCO
!       ELSE
!          WRITE(*,*) 'ERRO NA LEITURA DA COTA DO FUNDO DO POCCO'
!          STOP
!       END IF
!       WRITE(*,*)'########################################'
! !
!       keyword_name='tempo_para_a_descompressao_do_poco'
!       call readRealKeywordValue(keyword_name, TPRESPRODWELL, TPRESPRODWELL, ierr)
!       IF(ierr.ne.1) THEN  
!          WRITE(*,*)'TEMPO PARA DESCOMPRESSAO DO POCCO:',TPRESPRODWELL
!       ELSE
!          WRITE(*,*) 'ERRO NA LEITURA DO TEMPO PARA ADESCOMPRESSAO'
!          STOP
!       END IF
!       WRITE(*,*)'########################################'
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
!========================================================================
!
      subroutine abrirArquivosResultados
      use mPropGeoFisica
      use mGlobaisEscalares
      use mLeituraEscrita, only: iflag_tipoPrint
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
      subroutine imprimirCondicoesIniciais(pressaoElem, velocLadal, phi, perm, satElem, YOUNG, DIS, &
     &                                              PORE, PRESPROD, ndofV, ndofP, ndofD)
      use mGlobaisEscalares, only: tTransporte, numdx
      use mMalha,            only: x, conecNodaisElem, conecLadaisElem, nen, nsd
      use mMalha,            only: numel, numnp, numLadosElem
      use mMalha,            only: numelReserv, numLadosReserv
      USE mPropGeoFisica,    only: SRW,PHIINICIAL,PERMINICIAL
      use mLeituraEscrita,   only: iflag_tipoPrint
      use mLeituraEscrita,   only: prt
      
!
      implicit none
!
      real*8, intent(in) :: pressaoElem(ndofP, numelReserv)
      real*8, intent(in) :: velocLadal(ndofV,numLadosReserv)
      real*8, intent(in) :: phi(numelReserv), perm(numelReserv), satElem(numelReserv)
      real*8, intent(in) :: YOUNG(numel), DIS(ndofd,NUMNP), PORE(numelReserv), PRESPROD

      integer, intent(in)  ::  ndofV, ndofP, ndofD
      integer :: i, k, DZERO = 0
      LOGICAL :: SIM,NAO
      real*8  :: ZERO=0.D0
      real*8  :: YOUNGINICIAL=5e19
!
      SIM=.TRUE.
      NAO=.FALSE.
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
                  NAO,DZERO,'PRES',ndofP,PRESPROD)
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
                  NAO,DZERO,' SAT',ndofP,SRW)
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
                  NAO,DZERO,'PORE',ndofP,PHIINICIAL)
             call escreverArqParaview_escalar(iperm,perm,ZERO,ifperm_out,nen, &
                  NAO,DZERO,'PERM',ndofP,PERMINICIAL)
             call escreverArqParaview_escalar(iyou,YOUNG,ZERO,ifyou_out,nen, &
                  SIM,DZERO,'YOUN',ndofP,YOUNGINICIAL)
          end if
       endif
!
       RETURN
!
     END SUBROUTINE imprimirCondicoesIniciais
!
!*****************************************************
!
      subroutine imprimirSolucaoNoTempo(sat,DIS,PORE,YOUNG,pressaoElem,velocLadal,NCONDP, &
                                        PCONDP,PRESPROD,tempo, ndofV, ndofP, ndofD)
!
       use mMalha,            only: nsd, numel, numelReserv, numLadosReserv, numLadosElem, conecLadaisElem
       use mMalha,            only: nen, numnp
       USE mPropGeoFisica,    only: SRW,PHIINICIAL
       use mLeituraEscrita,   only: iflag_tipoPrint,prt
! 

      implicit none
!
      real*8, intent(in) :: pressaoElem(ndofP,numelReserv), velocLadal(ndofV,numLadosReserv)
      real*8, intent(in) :: sat(numelReserv),DIS(ndofD,numnp),PORE(numelReserv),YOUNG(numel)
      REAL(8), DIMENSION(2,100) :: PCONDP
      real*8, intent(in) :: PRESPROD, Tempo
      integer, intent(in)  ::  ndofV, ndofP, ndofD
      integer :: NCONDP
!
      character(21) :: labelTransp
      character(21) :: label
      character(21) :: num
      real*8        :: TOL=1e-6
      integer       :: UM=1
      LOGICAL       :: SIM,NAO
!
      SIM=.TRUE.
      NAO=.FALSE.
!
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
               call escreverArqParaview_escalar(isat,sat,tempo,ifsat_out,nen,NAO,npcontsat,' SAT',ndofP, SRW)
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
               call escreverArqParaview_escalar(ipres,pressaoElem,tempo,ifpres_out,nen,NAO,npcontpres,'PRES',ndofP, PRESPROD)
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
                    ifphi_out,nen,NAO,npcontphi,'PORE',ndofP,PHIINICIAL)
            end if
         endif
      end if
!
      RETURN
!
    END SUBROUTINE imprimirSolucaoNoTempo
!
!**** NEW ** MODIFIED FOR IRREGULAR MESH ***********************************
!
      SUBRoutine leituraGeoformations_DS(x, nsd, numnp)
      use mMalha, only: POSTLINE, SALTLINE, RTOPLINE, RBTTLINE, RIFTLINE
      use mMalha, only: LEFTLINE, RGHTLINE, IGEOFORM, IDome
      use mInputReader, only: readIntegerKeywordValue,readRealKeywordValue
!
!.... program to read, generate and write coordinate data
!
      implicit none
!
      INTEGER, INTENT(IN)    :: NSD, NUMNP
      REAL(8), INTENT(INOUT) ::  X(NSD,NUMNP)
!
      INTEGER      :: INXCORD
      CHARACTER*30 :: NAMEIN
      CHARACTER(LEN=128) :: FLAG, ROTULO
!
      INTEGER :: I, N, NODE
      REAL(8) :: XLEFT,XRIGHT,YBOTTOM,YTOP
      
      integer*4:: keyword_line
      character(len=50) keyword_name
      integer :: ierr
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
      keyword_name='read_geoformations_from_elmnt_geoformation'
      call readIntegerKeywordValue(keyword_name, IGEOFORM, IGEOFORM, ierr)      
!
      keyword_name = 'top_post_salt_horizontal_y_axis'
      call readRealKeywordValue(keyword_name, POSTLINE, POSTLINE, ierr)

      keyword_name = 'salt_dome_top_horizontal_y_axis'
      call readRealKeywordValue(keyword_name, SALTLINE, SALTLINE, ierr)

      keyword_name = 'top_reservoir_horizontal_y_axis'
      call readRealKeywordValue(keyword_name, RTOPLINE, RTOPLINE, ierr)

      keyword_name = 'bottom_reservoir_horizontal_y_axis'
      call readRealKeywordValue(keyword_name, RBTTLINE, RBTTLINE, ierr)

      keyword_name ='bottom_rift_horizontal_y_axis'
      call readRealKeywordValue(keyword_name, RIFTLINE, RIFTLINE, ierr)

      keyword_name ='left_reservoir_vertical_x_axis'
      call readRealKeywordValue(keyword_name, LEFTLINE, LEFTLINE, ierr)

      keyword_name ='right_reservoir_vertical_x_axis'
      call readRealKeywordValue(keyword_name, RGHTLINE, RGHTLINE, ierr)

      keyword_name ='reference_for_sismic_dome'
      call readIntegerKeywordValue(keyword_name, IDOME, IDOME, ierr)
!
      IF (YTOP.NE.POSTLINE)    POSTLINE=YTOP
      IF (YBOTTOM.NE.RIFTLINE) RIFTLINE=YBOTTOM
!                                            
      IF (XLEFT.GT.LEFTLINE) CALL CODERROR(6,NAMEIN)
!                                            
      IF (XRIGHT.LT.RGHTLINE) CALL CODERROR(6,NAMEIN)
!
      IF (RIFTLINE.GT.RBTTLINE) CALL CODERROR(6,NAMEIN)
      IF (RBTTLINE.GT.RTOPLINE) CALL CODERROR(6,NAMEIN)
      IF (RTOPLINE.GT.SALTLINE) CALL CODERROR(6,NAMEIN)
      IF (SALTLINE.GT.POSTLINE) CALL CODERROR(6,NAMEIN)
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
      SUBROUTINE GEOREGION_DS(XC,NUMEL)
!
      use mMalha,         only: IrregMesh, IGEOFORM
      use mPropGeoFisica, only: GEOFORM, GEOINDIC, GEOYLOC
      use mInputReader,   only: leituraRegiaoDs
!
!.... PROGRAM TO SETUP GEOMECHANICAL REGIONS: RESERVOIR-PRE-SAL, DOMO, POS-SAL
! 
      IMPLICIT NONE 
!  
      INTEGER :: NEL, NUMEL, INGEOFOR
!
      CHARACTER*30 :: NINGEOFOR
!
      INTEGER,  DIMENSION(NUMEL)   :: IREG
      REAL(8),  DIMENSION(2,NUMEL) :: XC
      
      character(len=50) keyword_name      
      integer :: ierr
!
      IF (IGEOFORM.EQ.0) THEN 
            DO 200 NEL=1,NUMEL
               GEOFORM(NEL) = GEOYLOC(XC(1,NEL),XC(2,NEL))
 200        CONTINUE
         ELSE
            keyword_name = "elmnt_geoformation"
            call leituraRegiaoDs(keyword_name, GEOFORM, NUMEL, ierr)
      ENDIF
!
      RETURN
!
1000  FORMAT(I10)
1010  FORMAT(A12)

!
      END SUBROUTINE
!
!**** new *******************************************************************
!     
     SUBROUTINE escreverArqParaview_vetor(arquivo,campo,passo,fname,nen, &
                              NRESERV,nprint,LABEL,dim)
!
    USE mMalha,            only: x,nsd,numel,numnp
    USE mPropGeoFisica,    only: GEOFORM
    use mMalha,            only: conecNodaisElem, conecLadaisElem
!
    IMPLICIT NONE
!
    INTEGER :: ARQUIVO,ISTAT,d,i,n,nen,dim,j
    REAL(8), DIMENSION(dim,*) :: CAMPO
    CHARACTER(LEN=128)    :: fname,NAME
    CHARACTER(LEN=10)     :: TEMP
    REAL(8)               :: COORDZ = 0.0,PASSO
    integer               :: tipo=1
    CHARACTER(len=4)      :: EXT,LABEL
    CHARACTER(len=5)      :: C
    INTEGER               :: VARIOSARQ,ZERO,NPRINT
    LOGICAL               :: NRESERV
!
    VARIOSARQ = 1
    ZERO = 0
!
    IF(VARIOSARQ.EQ.1)THEN
       WRITE(C,300)nprint
       C=ADJUSTL(C)
       EXT='.vtk'
       NAME=ADJUSTL(TRIM(fname))//TRIM(C)//TRIM(EXT)
       OPEN(UNIT=ARQUIVO,FILE=name,STATUS='UNKNOWN',&
            ACTION='READWRITE',IOSTAT=ISTAT)
    ELSE
       WRITE(C,300)ZERO
       C=ADJUSTL(C)
       EXT='.vtk'
       NAME=ADJUSTL(TRIM(fname))//TRIM(C)//TRIM(EXT)
       OPEN(UNIT=ARQUIVO,FILE=name,STATUS='UNKNOWN',&
            ACTION='READWRITE',IOSTAT=ISTAT,POSITION='APPEND')
    END IF
    WRITE(*,*)'ARQUIVO DE IMPRESSAO:',NAME
!
    IF(ISTAT.NE.0)THEN
       WRITE(*,*)'ERROR ON OPENING INPUT FILE: ',fname
       STOP
    END IF
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    IF(VARIOSARQ.EQ.1)THEN
       write(arquivo,'(a)')'# vtk DataFile Version 3.0'
       write(arquivo,'(a)')'vtk output'
       write(arquivo,'(a)')'ASCII'
       write(arquivo,'(a)')'DATASET UNSTRUCTURED_GRID'
       write(arquivo,'(a,i10,a)')'POINTS', numnp,' float '
!    write(arquivo,'(a,i10,a)')'POINTS', numnpReserv,' float '
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!  Escreve as coordenadas nodais
!
       if(nsd==2) then 
          do i=1,numnp
             write(arquivo,'(3(1x, 1pe15.8))') (x(d,i),d=1,nsd), coordZ 
          end do
       end if
!
       if(nsd==3) then
          do i=1,numnp
             write(arquivo,'(3(1x, 1pe15.8))') (x(d,i),d=1,nsd)
          end do
       end if
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Escreve as conectividades
       write(arquivo,'(a,i10,i10)')'CELLS', numel , (nen+1) * numel
       if(nsd==2) then
          do  n=1,numel
             write(arquivo,'(i10,9(2x,i10))') nen, (conecNodaisElem(i,n)-1, i = 1, nen) 
          end do
       end if
!
       if(nsd==3) then
          do  n=1,numel
             write(arquivo,'(i10,18(2x,i10))') nen, (conecNodaisElem(i,n)-1, i = 1, nen) 
          end do
       end if
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Escreve o tipo de celula
       write(arquivo,'(a,i10)')'CELL_TYPES ', numel
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
!
       if(tipo==1) write(arquivo,'(a,i10)')'POINT_DATA ', numnp

!    if(tipo==2) write(arquivo,'(a,i10)')'POINT_DATA',  numel*ndofP
    END IF
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    WRITE(TEMP,100)PASSO
    TEMP=ADJUSTL(TRIM(TEMP))
    WRITE(*,*)'TEMPO da IMPRESSAO:',TEMP
    IF(VARIOSARQ.EQ.1)THEN
       write(arquivo,'(3a,i10)')'SCALARS ', 'DISP' , ' float ',dim
    ELSE
       write(arquivo,'(3a,i10)')'SCALARS ', 't='//TRIM(TEMP)//'' , ' float ',dim
    END IF
    write(arquivo,'(2a)')'LOOKUP_TABLE ','default'
!
    DO I=1,NUMNP
       WRITE(ARQUIVO,*)(CAMPO(J,I),J=1,DIM)
    ENDDO
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    CLOSE(ARQUIVO) 
100 FORMAT(F10.5)
200 FORMAT(E15.7)
300 FORMAT(I5)
!
  END SUBROUTINE escreverArqParaview_vetor
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
SUBROUTINE escreverArqParaview_escalar(arquivo,campo,passo,fname,nen, &
     NRESERV,nprint,LABEL,dim,MINIMO)
!
    USE mMalha,            only: x,nsd,numel,numnp,numelReserv,numnpReserv
    USE mPropGeoFisica,    only: GEOFORM
    use mMalha,            only: conecNodaisElem, conecLadaisElem

!
    IMPLICIT NONE
!
    INTEGER :: ARQUIVO,ISTAT,d,i,n,nen,nprint,NPOSN,dim
    REAL(8), DIMENSION(dim,*) :: CAMPO
    REAL(8) :: MINIMO
    CHARACTER(LEN=128)    :: fname,NAME
    CHARACTER(LEN=21)     :: TEMPO
    REAL(8)               :: COORDZ = 0.0,PASSO,PROP
    integer               :: tipo=1
    INTEGER               :: NUMNPLOCAL,NUMELLOCAL
    LOGICAL               :: NRESERV
    CHARACTER(len=4)      :: EXT,LABEL
    CHARACTER(len=5)      :: C
    INTEGER               :: VARIOSARQ,ZERO
    INTEGER, DIMENSION(numel) :: ELEMRESERV
!
!     MINIMO = 1E30
    NUMNPLOCAL = numnp
    NUMELLOCAL = NUMEL
    VARIOSARQ = 1
    IF(VARIOSARQ.EQ.1)THEN
       WRITE(C,300)nprint
       C=ADJUSTL(C)
       EXT='.vtk'
       NAME=ADJUSTL(TRIM(fname))//TRIM(C)//TRIM(EXT)
    ELSE
       WRITE(C,300)ZERO
       C=ADJUSTL(C)
       EXT='.vtk'
       NAME=ADJUSTL(TRIM(fname))//TRIM(C)//TRIM(EXT)
    END IF
    WRITE(*,*)'ARQUIVO DE IMPRESSAO:',NAME
!
    IF(VARIOSARQ.EQ.1)THEN
       OPEN(UNIT=ARQUIVO,FILE=name,STATUS='UNKNOWN',&
            ACTION='WRITE',IOSTAT=ISTAT)
    ELSE
       OPEN(UNIT=ARQUIVO,FILE=name,STATUS='UNKNOWN',&
            ACTION='READWRITE',IOSTAT=ISTAT,POSITION='APPEND')
    END IF
    IF(ISTAT.NE.0)THEN
       WRITE(*,*)'ERROR ON OPENING INPUT FILE: ',fname
       STOP
    END IF
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
    IF(VARIOSARQ.EQ.1)THEN
       write(arquivo,'(a)')'# vtk DataFile Version 3.0'
       write(arquivo,'(a)')'vtk output'
       write(arquivo,'(a)')'ASCII'
       write(arquivo,'(a)')'DATASET UNSTRUCTURED_GRID'
       write(arquivo,'(a,i10,a)')'POINTS', NUMNPLOCAL,' float '
       !    write(arquivo,'(a,i10,a)')'POINTS', numnpReserv,' float '
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!  Escreve as coordenadas nodais
!
       if(nsd==2) then 
          do i=1,NUMNPLOCAL
             write(arquivo,'(3(1x, 1pe15.8))') (x(d,i),d=1,nsd), coordZ 
          end do
       end if
!
       if(nsd==3) then
          do i=1,NUMNPLOCAL
             write(arquivo,'(3(1x, 1pe15.8))') (x(d,i),d=1,nsd)
          end do
       end if
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Escreve as conectividades
       write(arquivo,'(a,i10,i10)')'CELLS', NUMELLOCAL , (nen+1) * NUMELLOCAL
       if(nsd==2) then
          do  n=1,NUMELLOCAL
             write(arquivo,'(i10,9(2x,i10))') nen, (conecNodaisElem(i,n)-1, i = 1, nen) 
          end do
       end if
!
       if(nsd==3) then
          do  n=1,NUMELLOCAL
             write(arquivo,'(i10,18(2x,i10))') nen, (conecNodaisElem(i,n)-1, i = 1, nen) 
          end do
       end if
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Escreve o tipo de celula
       write(arquivo,'(a,i10)')'CELL_TYPES ', NUMELLOCAL
!
       if(nsd==2) then
          do  i =1,NUMELLOCAL
             write(arquivo,'(a)') '9'!trim(adjustl(tipo))
          end do
       end if
!
       if(nsd==3) then
          do  i =1,NUMELLOCAL
             write(arquivo,'(a)') '12'!trim(adjustl(tipo))
          end do
       end if
!
       if(tipo==1) write(arquivo,'(a,i10)')'CELL_DATA ', NUMELLOCAL
    ENDIF
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    WRITE(TEMPO,'(F20.4)')PASSO
    TEMPO=ADJUSTL(TRIM(TEMPO))
    WRITE(*,*)'TEMPO da IMPRESSAO:',TEMPO
    IF(VARIOSARQ.EQ.1)THEN
       write(arquivo,'(3a)')'SCALARS ', ' '//TRIM(LABEL)//'' , ' float '
    ELSE
       write(arquivo,'(3a)')'SCALARS ', 't='//TRIM(TEMPO)//'' , ' float '
    ENDIF
    write(arquivo,'(2a)')'LOOKUP_TABLE ','default'
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!     if(LABEL.EQ.' SAT')MINIMO = SRW
!     if(LABEL.EQ.'PRES')MINIMO = PRESPROD!PRESMEDIAINICIAL
!     IF(LABEL.EQ.'PORE')MINIMO = PHIINICIAL
!     IF(LABEL.EQ.'PERM')MINIMO = PERMINICIAL
!     IF(LABEL.EQ.'YOUN')MINIMO = 5e19
!
    DO I=1,NUMELLOCAL
       IF(NRESERV)THEN
          PROP = CAMPO(dim,I)
       ELSE
          IF(GEOFORM(I).EQ.'RESERVATORIO')THEN
             PROP=CAMPO(dim,I)
          ELSE
             PROP=MINIMO
          END IF
       END IF
       WRITE(ARQUIVO,*)PROP
    ENDDO
!
    CLOSE(ARQUIVO) 
100 FORMAT(F10.5)
200 FORMAT(E15.7)
300 FORMAT(I5)
!
END SUBROUTINE escreverArqParaview_escalar

!
!*****************************************************
!
      subroutine imprimirCaseParaview(x, conecNodaisElem, pressaoElem, satElem, phi, perm)

      use mMalha,          only: numnp,nsd,numel,nen,numelReserv
      use mLeituraEscrita, only: iflag_tipoPrint

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
      use mGlobaisEscalares, only: NITGEO, NITHIDRO
      use mGlobaisEscalares, only: NUMDX, PATHDX, SALTCREEP
!..
!...  PROGRAM TO SETUP MANAGER FILES FOR GRAPHICAL INTERFACE OPEN-DX
!
      CHARACTER*30 NIFEDX,NIFNOD,NIFMSH,NIS3XX,NIS3YY,NIS3XY,NIS3ZZ
      CHARACTER*30 NISTDX,NISTYN,NISTPS,NISTPE,NISTSV,NISTPH
!
      CHARACTER*2 ASTEP1
      CHARACTER*2 ASTEP2
      CHARACTER*5 MECLAW
      CHARACTER*20 MKDIR
!
      IF (NUMDX.EQ.0) RETURN
!
      WRITE(ASTEP1,'(I2.2)') NITGEO
      WRITE(ASTEP2,'(I2.2)') NITHIDRO
!
      MECLAW = 'elast'
      IF (SALTCREEP) MECLAW = 'creep'
!
      PATHDX = 'dx'//MECLAW//ASTEP1//'.ht'//ASTEP2
!
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
      use mGlobaisEscalares, only: CYLINDER, NUMDX, PATHDX
      use mMalha,            only: XC, NEN, NSD, numel, numnp, numelReserv
      use mPropGeoFisica,    only: GEOFORM, GEOINDIC
!
      INTEGER :: K, NODE, NEL, ICENTER
      INTEGER, DIMENSION(NEN,NUMEL)   :: conecNodaisElem
      REAL(8), DIMENSION(NSD,NUMNP)   :: X, DIS
      REAL(8), DIMENSION(NUMEL)       :: YOUNG, GEOPRSR
      REAL(8), DIMENSION(numelReserv) :: PERM
      REAL(8), DIMENSION(4,NUMEL)     :: STRSS0
      REAL(8) :: XPRINT, RADIAL, ANGLE, PI
      REAL(8) :: SGXXCOS2, SGYYSIN2, SGXYSIN2, SIGMAR
      CHARACTER*30                    :: CENTER
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
      IF (CYLINDER) THEN
         PI = 4.0D0*DATAN(1.0D0)
         ICENTER = 516
         CENTER = PATHDX//'/cilindr.dat'
         OPEN(UNIT=ICENTER,FILE=CENTER)    
      ENDIF
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
         WRITE(ISTSV,1900) GEOPRSR(NEL)
         WRITE(IS3XX,1900) STRSS0(1,NEL) 
         WRITE(IS3YY,1900) STRSS0(2,NEL)
         WRITE(IS3XY,1900) STRSS0(3,NEL)
         WRITE(IS3ZZ,1900) STRSS0(4,NEL)
         IF (CYLINDER) THEN
            radial=DSQRT(XC(1,NEL)**2+XC(2,NEL)**2)
            angle=DATAN(XC(2,NEL)/XC(1,NEL))
! 
            sgxxcos2 = (dcos(angle))**2*STRSS0(1,NEL)
            sgyysin2 = (dsin(angle))**2*STRSS0(2,NEL)
            sgxysin2 = (dsin(2.0d0*angle))*STRSS0(3,NEL)
!
            sigmar = sgxxcos2+sgyysin2+sgxysin2 
!  
            WRITE(ICENTER,1900) Radial,ANGLE*180.0D0/PI,SIGMAR

         ENDIF
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
      IF (CYLINDER) CLOSE(ICENTER)
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
      
    subroutine readSetupPhaseDS(nlvectV, nlvectP, nlvectD, optSolverV, optSolverD)
        use mInputReader,      only:readStringKeywordValue, readIntegerKeywordValue
        use mGlobaisArranjos,  only: title
        use mGlobaisEscalares, only: exec, iprtin, optCC, TypeProcess
        use mMalha,            only: nsd, numel, numnp
        use mMalha,            only: numnpReserv, numelReserv, numelReserv
        use mPropGeoFisica,    only: nelx, nely, nelz
        use mPropGeoFisica,    only: nelxReserv, nelyReserv, nelzReserv
        use mGlobaisEscalares, only: geomech
        use mMalha,            only: IrregMesh, Dirichlet, Neumann, I4SeaLoad
        
        implicit none
        integer :: nlvectV, nlvectP, nlvectD
        character(len=*) :: optSolverV, optSolverD
        
        character(len=50) keyword_name
        integer :: ierr

        keyword_name = "type_process"
        call readStringKeywordValue(keyword_name, TypeProcess, ierr)
!
        keyword_name = "title"
        call readStringKeywordValue(keyword_name, title, ierr)
!
        keyword_name = "exec"
        call readIntegerKeywordValue(keyword_name, exec, 0_4, ierr)
!
        keyword_name = "iprtinc"
        call readIntegerKeywordValue(keyword_name, iprtin, 0_4, ierr)
!
        keyword_name = "nsd"
        call readIntegerKeywordValue(keyword_name, nsd, 0_4, ierr)
!
        keyword_name = "geomech"
        call readIntegerKeywordValue(keyword_name, geomech, 0_4, ierr)
        
        keyword_name = "irregMesh"
        call readIntegerKeywordValue(keyword_name, IrregMesh, 0_4, ierr)        

        keyword_name = "Dirichlet"
        call readIntegerKeywordValue(keyword_name, Dirichlet, 0_4, ierr) 

        keyword_name = "Neumann"
        call readIntegerKeywordValue(keyword_name, Neumann, 0_4, ierr)  

        keyword_name = "I4SeaLoad"
        call readIntegerKeywordValue(keyword_name, I4SeaLoad, 0_4, ierr)  
        
        keyword_name = "numnp"
        call readIntegerKeywordValue(keyword_name, numnp, 0_4, ierr)
        !Reads numel
        keyword_name = "numel"
        call readIntegerKeywordValue(keyword_name, numel, 0_4, ierr)
        !Reads nelx
        keyword_name = "nelx"
        call readIntegerKeywordValue(keyword_name, nelx, 0_4, ierr)
        !Reads nely
        keyword_name = "nely"
        call readIntegerKeywordValue(keyword_name, nely, 0_4, ierr)
        !Reads nelz
        keyword_name = "nelz"
        call readIntegerKeywordValue(keyword_name, nelz, 0_4, ierr)
        !Reads numnpReserv
        keyword_name = "numnpReserv"
        call readIntegerKeywordValue(keyword_name, numnpReserv, 0_4, ierr)
        !Reads numelReserv
        keyword_name = "numelReserv"
        call readIntegerKeywordValue(keyword_name, numelReserv, 0_4, ierr)
        !Reads nelxReserv
        keyword_name = "nelxReserv"
        call readIntegerKeywordValue(keyword_name, nelxReserv, 0_4, ierr)
        !Reads nelyReserv
        keyword_name = "nelyReserv"
        call readIntegerKeywordValue(keyword_name, nelyReserv, 0_4, ierr)
        !Reads nelzReserv
        keyword_name = "nelzReserv"
        call readIntegerKeywordValue(keyword_name, nelzReserv, 0_4, ierr)
        !Reads nlvectP
        keyword_name = "nlvectP"
        call readIntegerKeywordValue(keyword_name, nlvectP, 0_4, ierr)
        !Reads nlvectV
        keyword_name = "nlvectV"
        call readIntegerKeywordValue(keyword_name, nlvectV, 0_4, ierr)
        !Reads nlvectD
        keyword_name = "nlvectD"
        call readIntegerKeywordValue(keyword_name, nlvectD, 0_4, ierr)

        keyword_name = "solver_hidrodinamica"
        call readStringKeywordValue(keyword_name, optSolverV, ierr)
        
        keyword_name = "solver_geomecanica"
        call readStringKeywordValue(keyword_name, optSolverD, ierr)
       
        
        return
    end subroutine readSetupPhaseDS 
!
!=======================================================================
!
      subroutine prtvB3D(x,nen,nsd,numel,conecNodaisElem,t0,velocLadal,ndofV, numLados, conecLadaisElem, numLadosElem, iunit)
!
      use mMalha, only: local 
!
      implicit none
!
!     imprime campos vetoriais para o gnuplot ou para o matlab
!
      real*8 :: x(nsd,*)
      integer                   :: nen,numel,nsd, ndofV ,numLados,numLadosElem
      real(8), dimension(ndofV,numLados) :: velocLadal
      real(8)                   :: t0
      integer                   :: conecNodaisElem(nen,numel),conecLadaisElem(numLadosElem,numel)
!
      integer  :: nel
      real*8   :: xl(nsd,nen)
      real(8)  :: xg,yg,zg,vcx,vcy,vcz
      integer :: lado1,lado2,lado3,lado4,lado5,lado6
!
      integer :: iunit
!
      write(iunit,"('#TIMESTEP PRINT OUT = ',f15.8)") t0
!
      write(iunit,*)
!
      do nel=1,numel
!
        call local(conecNodaisElem(1,nel),x,xl,nen,nsd,nsd)
        xg = sum(xl(1,1:nen))/nen
        yg = sum(xl(2,1:nen))/nen
        zg = sum(xl(3,1:nen))/nen

        lado1 = conecLadaisElem(1,nel);  lado2 = conecLadaisElem(2,nel);
        lado3 = conecLadaisElem(3,nel);  lado4 = conecLadaisElem(4,nel);
        lado5 = conecLadaisElem(5,nel);  lado6 = conecLadaisElem(6,nel);
        vcx = (velocLadal(1,lado2)+velocLadal(1,lado4))/2.0
        vcy = (velocLadal(1,lado1)+velocLadal(1,lado3))/2.0
        vcz = (velocLadal(1,lado5)+velocLadal(1,lado6))/2.0
        write(iunit,"(6(f25.15,2x))") xg,yg,zg,vcx,vcy,vcz

      end do
!
      write(iunit,*)
!
      end subroutine
      
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
      END MODULE
!
!******** ************ ************ ************ *********** **********
!
